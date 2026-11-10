"""Correlated-k (CKD) coefficients of SOS-ABS as NetCDF / xarray.

The ASCII files ``fic/COEFF_CKD/<res>cmm1/coef_<GAS>_<numax>_<numin>_<res>cmm1``
(one file per gas and 500 cm-1 chunk) are gathered in one NetCDF file per
spectral resolution. The ASCII files stay the input of the Fortran code;
the NetCDF file is a faster, self-describing copy for Python analyses.

Content of the dataset (``n_term`` <= 5 terms of the CKD expansion)::

    a(gas, wavenumber, term)                         weights a_i
    k(gas, wavenumber, term, pressure, temperature)  k_i(P, T), all gases but H2O
    k_h2o(wavenumber, term, h2o_fraction, pressure, temperature)
    n_term(gas, wavenumber)                          number of terms

The transmission of a gas over an interval for a column amount ``u``
(molecules cm-2) is ``sum_i a_i exp(-k_i(P, T) u)``. ``n_term = 0`` means
no absorption: the Fortran code then uses a single term with ``a = 1`` and
``k = 0``; in the dataset, ``a`` and ``k`` are NaN for the unused terms.

Example::

    import sosabs
    sosabs.ckd.to_netcdf(10)                 # fic/COEFF_CKD/ckd_10cmm1.nc
    ds = sosabs.ckd.open_ckd(10)
    ds.k.sel(gas="O2").sel(wavenumber=13100, method="nearest")
"""

import glob
import os
import re
import time

import numpy as np
import xarray as xr

from .runner import sos_root

#: gases of the CKD files, in the order of the Fortran arrays (H2O CO2 O3 N2O
#: CO CH4 O2 NO2); N2 files are provided but not used by SOS-ABS
GASES = ("H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "NO2", "N2")
RESOLUTIONS = (1, 5, 10)
N_TERM_MAX = 5  # CTE_CKD_NAI_MAX in inc/SOS.h


def ckd_dir(resolution, root=None):
    """Directory of the ASCII CKD files for ``resolution`` (1, 5 or 10 cm-1)."""
    if resolution not in RESOLUTIONS:
        raise ValueError("resolution must be one of %s" % (RESOLUTIONS,))
    return os.path.join(sos_root(root), "fic", "COEFF_CKD",
                        "%dcmm1" % resolution)


def default_path(resolution, root=None):
    """Default NetCDF file: ``fic/COEFF_CKD/ckd_<res>cmm1.nc``."""
    return os.path.join(sos_root(root), "fic", "COEFF_CKD",
                        "ckd_%dcmm1.nc" % resolution)


# ---------------------------------------------------------------------------
# Reading of the ASCII files
# ---------------------------------------------------------------------------

_NAME = re.compile(r"coef_([A-Za-z0-9]+)_(\d+)_(\d+)_(\d+)cmm1$")


def read_ascii(path):
    """Read one ASCII CKD file (same logic as ``READ_CKD_COEFF``).

    Returns a dict with the grids (``temperature``, ``pressure``,
    ``h2o_fraction`` or None), the intervals (``numin``, ``numax``,
    ``wl_min``, ``wl_mean``, ``wl_max``, ``n_term``) and the arrays
    ``a`` (interval, term) and ``k`` (interval, term, [conc,] pressure,
    temperature), NaN-padded to ``N_TERM_MAX`` terms.
    """
    with open(path) as f:
        lines = f.read().splitlines()
    # header: text between the first two lines of dashes
    seps = [i for i, l in enumerate(lines) if l.startswith("-----")]
    if len(seps) < 2:
        raise ValueError("%s: header not found" % path)
    tok = np.array(" ".join(lines[seps[1] + 1:]).split(), dtype=float)

    pos = 0

    def take(n):
        nonlocal pos
        out = tok[pos:pos + n]
        if out.size != n:
            raise ValueError("%s: unexpected end of file" % path)
        pos += n
        return out

    numax, numin, step = take(3)
    nt = int(take(1)[0])
    temp = take(nt)
    npres = int(take(1)[0])
    pres = take(npres)
    is_h2o = os.path.basename(path).startswith("coef_H2O_")
    if is_h2o:
        nc = int(take(1)[0])
        conc = take(nc)
        shape = (nc, npres, nt)
    else:
        conc = None
        shape = (npres, nt)

    n_int = int(round((numax - numin) / step))
    info = np.empty((n_int, 6))
    a = np.full((n_int, N_TERM_MAX), np.nan)
    k = np.full((n_int, N_TERM_MAX) + shape, np.nan)
    nrow = int(np.prod(shape[:-1]))          # lines per term
    ncol = len(shape) + nt                   # indices + nt values
    for i in range(n_int):
        info[i] = take(6)
        nexp = int(info[i, 5])
        if nexp == 0:
            continue
        if nexp > N_TERM_MAX:
            raise ValueError("%s: %d terms > %d" % (path, nexp, N_TERM_MAX))
        a[i, :nexp] = take(nexp)
        blk = take(nexp * nrow * ncol).reshape(nexp, nrow, ncol)
        # check the indices written on each line (term, [conc,] pressure)
        idx = blk[:, :, :ncol - nt].astype(int)
        expect = np.stack(np.meshgrid(*[np.arange(1, s + 1) for s in
                                        (nexp,) + shape[:-1]],
                                      indexing="ij"), -1)
        if not np.array_equal(idx, expect.reshape(nexp, nrow, -1)):
            raise ValueError("%s: unexpected line order in interval %d"
                             % (path, i))
        k[i, :nexp] = blk[:, :, ncol - nt:].reshape((nexp,) + shape)
    if pos != tok.size:
        raise ValueError("%s: %d values not read" % (path, tok.size - pos))
    return dict(temperature=temp, pressure=pres, h2o_fraction=conc,
                wl_min=info[:, 0], wl_mean=info[:, 1], wl_max=info[:, 2],
                # the files give the larger wavenumber first
                numin=np.minimum(info[:, 3], info[:, 4]),
                numax=np.maximum(info[:, 3], info[:, 4]),
                n_term=info[:, 5].astype(int), a=a, k=k, step=step)


def _files(resolution, gas, root):
    files = []
    for p in glob.glob(os.path.join(ckd_dir(resolution, root),
                                    "coef_%s_*" % gas)):
        m = _NAME.search(os.path.basename(p))
        if m and m.group(1) == gas:
            files.append((int(m.group(3)), p))
    return [p for _, p in sorted(files)]


def read_gas(resolution, gas, root=None):
    """Read all the files of one gas, sorted by increasing wavenumber."""
    parts = [read_ascii(p) for p in _files(resolution, gas, root)]
    if not parts:
        raise FileNotFoundError("no CKD file for %s at %d cm-1"
                                % (gas, resolution))
    for key in ("temperature", "pressure", "h2o_fraction"):
        ref = parts[0][key]
        for p in parts[1:]:
            if (ref is None) != (p[key] is None) or (
                    ref is not None and not np.array_equal(ref, p[key])):
                raise ValueError("%s: %s grid differs between files"
                                 % (gas, key))
    out = {key: parts[0][key] for key in
           ("temperature", "pressure", "h2o_fraction")}
    for key in ("wl_min", "wl_mean", "wl_max", "numin", "numax", "n_term",
                "a", "k"):
        arr = np.concatenate([p[key] for p in parts])
        out[key] = arr
    order = np.argsort(out["numin"])          # files list intervals downward
    for key in ("wl_min", "wl_mean", "wl_max", "numin", "numax", "n_term",
                "a", "k"):
        out[key] = out[key][order]
    return out


# ---------------------------------------------------------------------------
# Dataset and NetCDF
# ---------------------------------------------------------------------------


def to_dataset(resolution, gases=GASES, root=None):
    """Read the ASCII CKD files of ``resolution`` into an xarray Dataset."""
    data = {g: read_gas(resolution, g, root) for g in gases}
    ref = data[gases[0]]
    for g in gases:
        for key in ("numin", "numax"):
            if not np.array_equal(data[g][key], ref[key]):
                raise ValueError("%s: spectral intervals differ" % g)
    others = [g for g in gases if g != "H2O"]
    o = data[others[0]] if others else None
    for g in others:
        for key in ("temperature", "pressure"):
            if not np.array_equal(data[g][key], o[key]):
                raise ValueError("%s: %s grid differs" % (g, key))

    nu = 0.5 * (ref["numin"] + ref["numax"])
    coords = {
        "gas": ("gas", np.array(gases, dtype=object)),
        "wavenumber": ("wavenumber", nu, {
            "long_name": "central wavenumber of the interval",
            "units": "cm-1"}),
        "term": ("term", np.arange(1, N_TERM_MAX + 1), {
            "long_name": "index i of the CKD expansion term"}),
        "wavenumber_min": ("wavenumber", ref["numin"], {"units": "cm-1"}),
        "wavenumber_max": ("wavenumber", ref["numax"], {"units": "cm-1"}),
        "wavelength": ("wavenumber", ref["wl_mean"], {
            "long_name": "mean wavelength of the interval",
            "units": "micrometer"}),
        "wavelength_min": ("wavenumber", ref["wl_min"],
                           {"units": "micrometer"}),
        "wavelength_max": ("wavenumber", ref["wl_max"],
                           {"units": "micrometer"}),
    }
    if others:
        coords["pressure"] = ("pressure", o["pressure"],
                              {"units": "hPa", "long_name": "pressure"})
        coords["temperature"] = ("temperature", o["temperature"],
                                 {"units": "K", "long_name": "temperature"})
    ds = xr.Dataset(coords=coords)
    ds["n_term"] = (("gas", "wavenumber"),
                    np.stack([data[g]["n_term"] for g in gases]).astype("i1"),
                    {"long_name": "number of terms of the CKD expansion "
                                  "(0: no absorption)"})
    ds["a"] = (("gas", "wavenumber", "term"),
               np.stack([data[g]["a"] for g in gases]),
               {"long_name": "weights a_i of the CKD expansion"})
    if others:
        ds["k"] = (("gas", "wavenumber", "term", "pressure", "temperature"),
                   np.stack([data[g]["k"] if g != "H2O" else
                             np.full_like(o["k"], np.nan) for g in gases]),
                   {"long_name": "absorption coefficients k_i(P, T) "
                                 "(NaN for H2O, see k_h2o)",
                    "units": "cm2 molecule-1"})
    if "H2O" in gases:
        h = data["H2O"]
        ds = ds.assign_coords(
            h2o_pressure=("h2o_pressure", h["pressure"], {"units": "hPa"}),
            h2o_temperature=("h2o_temperature", h["temperature"],
                             {"units": "K"}),
            h2o_fraction=("h2o_fraction", h["h2o_fraction"], {
                "long_name": "H2O mole fraction", "units": "1"}))
        ds["k_h2o"] = (("wavenumber", "term", "h2o_fraction", "h2o_pressure",
                        "h2o_temperature"), h["k"],
                       {"long_name": "H2O absorption coefficients "
                                     "k_i(C, P, T)",
                        "units": "cm2 molecule-1"})
    ds.attrs = {
        "title": "SOS-ABS correlated-k distribution coefficients",
        "spectral_resolution": "%d cm-1" % resolution,
        "source": "fic/COEFF_CKD/%dcmm1 (ASCII files read by the SOS-ABS "
                  "routine READ_CKD_COEFF)" % resolution,
        "references": "CKD coefficients computed by HYGEOS for SOS-ABS "
                      "(CNES); Lacis and Oinas, J. Geophys. Res. 96, "
                      "9027-9063, 1991",
        "Conventions": "CF-1.8",
        "comment": "transmission = sum_i a_i exp(-k_i u); n_term = 0 means "
                   "no absorption (one term a = 1, k = 0 in SOS-ABS); "
                   "N2 is provided but not used by SOS-ABS",
        "history": "%s: created by sosabs.ckd.to_netcdf"
                   % time.strftime("%Y-%m-%d"),
    }
    return ds


def to_netcdf(resolution, path=None, gases=GASES, root=None,
              dtype="float32", complevel=4):
    """Convert the ASCII CKD files of ``resolution`` to one NetCDF file.

    ``dtype="float32"`` (default) is lossless in practice: the ASCII values
    are themselves single-precision numbers (smallest value 1.4e-45, the
    float32 denormal minimum); ``"float64"`` doubles the size. Returns the
    path.
    """
    ds = to_dataset(resolution, gases, root)
    path = path or default_path(resolution, root)
    enc = {}
    for v in ds.data_vars:
        e = {"zlib": True, "complevel": complevel, "shuffle": True}
        if ds[v].dtype.kind == "f":
            e["dtype"] = dtype
            e["_FillValue"] = np.nan
        enc[v] = e
    ds.to_netcdf(path, encoding=enc)
    return path


def open_ckd(resolution=10, path=None, root=None, create=False):
    """Open the NetCDF CKD file of ``resolution`` (lazy xarray Dataset).

    With ``create=True`` the file is built from the ASCII files if missing.
    """
    path = path or default_path(resolution, root)
    if not os.path.exists(path):
        if not create:
            raise FileNotFoundError(
                "%s not found: run sosabs.ckd.to_netcdf(%d)"
                % (path, resolution))
        to_netcdf(resolution, path, root=root)
    return xr.open_dataset(path)


def transmittance(ds, gas, amount, pressure, temperature, h2o_fraction=None):
    """Band transmission sum_i a_i exp(-k_i u) of one gas and one layer.

    ``amount`` in molecules cm-2; ``pressure`` (hPa), ``temperature`` (K)
    and, for H2O, ``h2o_fraction`` are interpolated linearly in the tables
    (the Fortran code uses a spline in temperature). Returns a DataArray
    versus wavenumber.
    """
    a = ds.a.sel(gas=gas)
    if gas == "H2O":
        if h2o_fraction is None:
            raise ValueError("h2o_fraction is required for H2O")
        k = ds.k_h2o.interp(h2o_fraction=h2o_fraction, h2o_pressure=pressure,
                            h2o_temperature=temperature)
    else:
        k = ds.k.sel(gas=gas).interp(pressure=pressure,
                                     temperature=temperature)
    t = (a * np.exp(-k * amount)).sum("term", min_count=1)
    return t.where(ds.n_term.sel(gas=gas) > 0, 1.0).rename("transmittance")
