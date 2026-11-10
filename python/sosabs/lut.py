"""Fast generation of SOS-ABS look-up tables (LUT) in NetCDF.

A LUT is the set of SOS-ABS simulations for all the combinations of the
values of a few parameters (the *grid*), the other parameters being fixed.
The simulations are run in parallel and gathered in one
:class:`xarray.Dataset`:

* one dimension per grid parameter (``wl``, ``h2o``, ``sza``...),
* ``vza`` and ``azi`` for the Stokes parameters,
* the fixed parameters, the code version and the conventions stored as
  NetCDF attributes.

Speed comes from:

* the single-pass azimuth recombination of the Fourier series in the code
  (``SOS_TRPHI_MULTI``, all azimuths from one reading of the OS file),
* the in-memory f2py backend (no ASCII result file),
* the parallel runs, each worker keeping its own cache of Mie and surface
  matrix files (no concurrent writing of the same file).

Example
-------
>>> from sosabs import build_lut
>>> lut = build_lut(
...     grid={"-SOS_Main.Wa": [0.90, 0.94, 0.98], "-AP.H2O": [0.5, 2.0]},
...     fixed={"-ANG.Thetas": 40., "-ANG.Rad.NbGauss": 24, ...},
...     outfile="lut.nc", n_workers=8)
"""

import datetime
import itertools
import multiprocessing as mp
import os
import queue
import shutil
import subprocess
import tempfile
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import xarray as xr

from . import params as P
from . import runner

# keys of the parameters set by the LUT builder itself
_RESERVED = ("resroot", "dir_mie", "dir_surf", "itrphi")


# ---------------------------------------------------------------------------
# Workers
# ---------------------------------------------------------------------------

_SLOTS = None  # queue of the cache slots (set in each worker process)


def _init_worker(slots):
    global _SLOTS
    _SLOTS = slots


def _run_job(job):
    """Run one simulation in a cache slot; return (index, dataset|error)."""
    idx, par, backend, root, cache_dir, slots = job
    slots = slots if slots is not None else _SLOTS
    slot = slots.get()
    try:
        sdir = os.path.join(cache_dir, "slot_%03d" % slot)
        par = dict(par)
        par["dir_mie"] = os.path.join(sdir, "MIE")
        par["dir_surf"] = os.path.join(sdir, "SURF")
        os.makedirs(sdir, exist_ok=True)
        work = tempfile.mkdtemp(prefix="run_", dir=sdir)
        try:
            ds = runner.run(par, backend=backend, root=root, workdir=work,
                            isolate=False)
            return idx, ds.load()
        finally:
            shutil.rmtree(work, ignore_errors=True)
    except Exception as exc:  # reported by build_lut
        return idx, exc
    finally:
        slots.put(slot)


# ---------------------------------------------------------------------------
# LUT builder
# ---------------------------------------------------------------------------


def build_lut(grid, fixed=None, outfile=None, n_workers=None,
              backend="auto", root=None, cache_dir=None, common_vza=True,
              dtype="float32", attrs=None, verbose=True):
    """Run SOS-ABS on a parameter grid and gather the results.

    Parameters
    ----------
    grid : dict
        Parameters of the LUT and their values, e.g.
        ``{"-SOS_Main.Wa": [0.90, 0.91], "-AP.H2O": [0.5, 1, 2]}``
        (keywords or Python names). The order gives the order of the
        dimensions.
    fixed : dict, optional
        Other SOS-ABS parameters (same for all the runs). A parameter of
        ``grid`` given in ``fixed`` too takes the grid values.
    outfile : str, optional
        NetCDF file to write.
    n_workers : int, optional
        Number of parallel runs (default: number of CPUs, max. 32).
    backend : {"auto", "f2py", "exe"}
        See :func:`sosabs.run`.
    root : str, optional
        SOS-ABS root directory (default: ``$SOS_ABS_ROOT``).
    cache_dir : str, optional
        Directory of the Mie / surface files cache, one sub-directory per
        worker. Reusing it between LUTs avoids recomputing the Mie files.
        Default: a temporary directory removed at the end.
    common_vza : bool
        Keep only the viewing angles common to all the runs (Gauss and
        user angles): the solar zenith angle is added by SOS-ABS to the
        output angles, so it differs between runs when ``sza`` is a LUT
        dimension. Ignored if ``-ANG.Thetas`` is not in ``grid``.
    dtype : str
        Storage type of the Stokes parameters in the NetCDF file.
    attrs : dict, optional
        Extra global attributes.

    Returns
    -------
    xarray.Dataset
    """
    t0 = time.time()
    root = runner.sos_root(root)
    os.environ["SOS_ABS_ROOT"] = root
    if backend == "auto":
        backend = runner.available_backend(root)

    grid = P.normalize(grid)
    fixed = P.normalize(fixed or {})
    for name in grid:
        if name in _RESERVED:
            raise ValueError("%s cannot be a LUT dimension" % name)
        fixed.pop(name, None)       # the grid values override fixed
    names = list(grid)
    values = [np.atleast_1d(np.asarray(grid[n])) for n in names]
    shape = tuple(len(v) for v in values)
    if n_workers is None:
        n_workers = min(os.cpu_count() or 1, 32)
    n_jobs = int(np.prod(shape))
    n_workers = max(1, min(n_workers, n_jobs))

    tmp_cache = cache_dir is None
    cache_dir = os.path.abspath(cache_dir or tempfile.mkdtemp(
        prefix="sosabs_cache_"))
    os.makedirs(cache_dir, exist_ok=True)

    jobs = []
    for idx in itertools.product(*[range(n) for n in shape]):
        par = dict(fixed)
        for name, v, i in zip(names, values, idx):
            par[name] = v[i].item()
        jobs.append([idx, par, backend, root, cache_dir])

    if verbose:
        print("SOS-ABS LUT: %d runs (%s), %d workers, backend %s" % (
            n_jobs, " x ".join("%s:%d" % (P.dim_name(n), s)
                               for n, s in zip(names, shape)),
            n_workers, backend))

    results = {}
    try:
        if backend == "f2py":
            # one process per run (fresh Fortran state), slots shared
            # through a manager queue
            ctx = mp.get_context("fork")
            with ctx.Manager() as manager:
                slots = manager.Queue()
                for k in range(n_workers):
                    slots.put(k)
                with ctx.Pool(n_workers, maxtasksperchild=1,
                              initializer=_init_worker,
                              initargs=(slots,)) as pool:
                    it = pool.imap_unordered(
                        _run_job, [tuple(j) + (None,) for j in jobs])
                    _collect(it, results, n_jobs, verbose, t0)
        else:
            slots = queue.Queue()
            for k in range(n_workers):
                slots.put(k)
            with ThreadPoolExecutor(n_workers) as ex:
                it = ex.map(_run_job, [tuple(j) + (slots,) for j in jobs])
                _collect(it, results, n_jobs, verbose, t0)
    finally:
        if tmp_cache:
            shutil.rmtree(cache_dir, ignore_errors=True)

    lut = _assemble(names, values, shape, results, common_vza, fixed,
                    backend, root, time.time() - t0, attrs)
    if outfile:
        write_netcdf(lut, outfile, dtype=dtype)
        if verbose:
            print("LUT written to %s" % outfile)
    return lut


def _collect(iterator, results, n_jobs, verbose, t0):
    lock = threading.Lock()
    for k, (idx, res) in enumerate(iterator, 1):
        with lock:
            results[idx] = res
        if isinstance(res, Exception):
            warnings.warn("run %s failed: %s" % (idx, res))
        if verbose and (k == n_jobs or k % max(1, n_jobs // 20) == 0):
            print("  %d/%d runs done (%.0f s)" % (k, n_jobs,
                                                 time.time() - t0),
                  flush=True)


def _assemble(names, values, shape, results, common_vza, fixed, backend,
              root, elapsed, extra_attrs):
    ok = [ds for ds in results.values() if isinstance(ds, xr.Dataset)]
    if not ok:
        raise RuntimeError("all the SOS-ABS runs failed:\n%s" % "\n".join(
            str(e) for e in results.values()))
    failed = [idx for idx, ds in results.items()
              if not isinstance(ds, xr.Dataset)]

    def select(ds):
        if common_vza and "tetas" in names and "gauss_weight" in ds.coords:
            keep = (ds.gauss_weight > 0) | (ds.user_angle == 1)
            ds = ds.isel(vza=np.where(keep.values)[0])
        return ds

    ref = select(ok[0])
    vza, azi = ref.vza.values, ref.azi.values

    dims = [P.dim_name(n) for n in names]
    coords = {}
    for name, dim, v in zip(names, dims, values):
        a = dict(P.DIM_ATTRS.get(name, {}))
        a["sos_keyword"] = P.KEYWORD[name]
        coords[dim] = (dim, v, a)
    coords["vza"] = ("vza", vza, ref.vza.attrs)
    coords["azi"] = ("azi", azi, ref.azi.attrs)
    for c in ("gauss_weight", "user_angle"):
        if c in ref.coords:
            coords[c] = ("vza", ref[c].values, ref[c].attrs)

    data = {}
    for var in ref.data_vars:
        vshape = shape + ref[var].shape
        data[var] = np.full(vshape, np.nan)
    for idx, ds in results.items():
        if not isinstance(ds, xr.Dataset):
            continue
        ds = select(ds)
        if ds.sizes["vza"] != len(vza) or not np.allclose(ds.vza, vza,
                                                         atol=1e-3):
            raise RuntimeError(
                "viewing angles differ between runs: use common_vza=True "
                "or a fixed -ANG.Thetas")
        for var in ref.data_vars:
            data[var][idx] = ds[var].values

    dvars = {var: (tuple(dims) + ref[var].dims, data[var], ref[var].attrs)
             for var in ref.data_vars}
    lut = xr.Dataset(dvars, coords=coords)

    # global attributes: those of the runs without the grid parameters
    attrs = {k: v for k, v in ref.attrs.items()
             if not any(k == "sos" + P.KEYWORD[n].replace("-", "_")
                        .replace(".", "_") for n in names)}
    now = datetime.datetime.now(datetime.timezone.utc)
    attrs.update(
        title="SOS-ABS look-up table",
        Conventions="CF-1.8",
        date_created=now.isoformat(timespec="seconds"),
        history="%s: created by sosabs.build_lut (%d runs, %.0f s)" % (
            now.strftime("%Y-%m-%d %H:%M:%S UTC"), int(np.prod(shape)),
            elapsed),
        lut_dimensions=" ".join(dims),
        lut_keywords=" ".join(P.KEYWORD[n] for n in names),
        sos_code_version=_code_version(root),
        sos_backend=backend,
        n_failed_runs=len(failed),
    )
    if failed:
        attrs["failed_runs"] = str(failed)
    attrs.update(extra_attrs or {})
    lut.attrs = attrs
    return lut


def _code_version(root):
    """'SOS-ABS version x.y' and git commit of the code, if available."""
    version = ""
    try:
        with open(os.path.join(root, "info_Version.txt")) as f:
            version = f.readline().strip()
    except OSError:
        pass
    try:
        commit = subprocess.run(
            ["git", "-C", root, "describe", "--always", "--dirty"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
            timeout=10).stdout.strip()
        if commit:
            version += " (git %s)" % commit
    except (OSError, subprocess.SubprocessError):
        pass
    return version


def write_netcdf(lut, outfile, dtype="float32", complevel=4):
    """Write a LUT in a compressed NetCDF4 file."""
    enc = {}
    for var in lut.data_vars:
        e = dict(zlib=True, complevel=complevel)
        if lut[var].ndim > len(lut.attrs.get("lut_dimensions", "")
                               .split()):
            e["dtype"] = dtype
        enc[var] = e
    # NetCDF attributes cannot be booleans or None
    lut = lut.copy()
    for obj in [lut] + [lut[v] for v in lut.variables]:
        for k, v in list(obj.attrs.items()):
            if isinstance(v, bool) or v is None:
                obj.attrs[k] = str(v)
    lut.to_netcdf(outfile, encoding=enc, engine="netcdf4")
    return outfile
