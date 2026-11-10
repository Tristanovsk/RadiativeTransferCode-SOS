"""Run one SOS-ABS simulation and return it as an :class:`xarray.Dataset`.

Two backends are available:

``"f2py"``
    In-memory call of ``SOS_PROC`` through the f2py module ``sos``
    (``binding/compilation_f2py``). The Stokes parameters are returned in
    double precision, without writing nor parsing the ASCII result files.
``"exe"``
    Call of the executable ``exe/SOS_ABS_MAIN.exe`` and reading of its
    ASCII result files (6 significant digits).

``"auto"`` (default) uses f2py when the module can be imported.

The radiation field is always recombined as a polar diagram
(``-SOS.View 2``) on the azimuth grid 0, ``pas_phi``, ..., 360 degrees.
"""

import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import xarray as xr

from . import params as P

# unbuffered Fortran console output (read when libgfortran is loaded),
# needed to redirect it to the log file of each run
os.environ.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")

LOGFILE = "sos_console.log"

# ---------------------------------------------------------------------------
# Location of the code
# ---------------------------------------------------------------------------


def sos_root(root=None):
    """Root directory of SOS-ABS (``root`` or ``$SOS_ABS_ROOT``)."""
    root = root or os.environ.get("SOS_ABS_ROOT")
    if not root:
        # python/sosabs/runner.py -> repository root
        root = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            "..", ".."))
    if not os.path.isdir(os.path.join(root, "fic")):
        raise RuntimeError(
            "SOS-ABS root directory not found (%s): set SOS_ABS_ROOT "
            "to the directory containing fic/, exe/ and binding/" % root)
    return os.path.abspath(root)


def _import_f2py(root):
    """Import the f2py module ``sos`` from ``root/binding``."""
    path = os.path.join(root, "binding")
    if path not in sys.path:
        sys.path.insert(0, path)
    import sos  # noqa: E402  (built by binding/compilation_f2py)
    return sos


def available_backend(root=None):
    """'f2py' if the f2py module can be imported, else 'exe'."""
    try:
        _import_f2py(sos_root(root))
        return "f2py"
    except ImportError:
        return "exe"


# ---------------------------------------------------------------------------
# Single run
# ---------------------------------------------------------------------------


def run(params, backend="auto", root=None, workdir=None, keep=False,
        verbose=False, isolate=True):
    """Run SOS-ABS for one set of parameters.

    Parameters
    ----------
    params : dict
        SOS-ABS parameters, indexed by keyword (``"-AP.H2O"``) or by
        Python name (``"h2o"``). Missing parameters take the default
        values of :data:`sosabs.params.PARAMETERS`. ``-AER.DirMie`` and
        ``-SURF.Dir`` can point to shared cache directories (the Mie and
        surface files are reused if they exist).
    backend : {"auto", "f2py", "exe"}
    root : str, optional
        SOS-ABS root directory (default: ``$SOS_ABS_ROOT``).
    workdir : str, optional
        Directory of the run (``-SOS_Main.ResRoot``). A temporary
        directory is used by default.
    keep : bool
        Keep the temporary working directory.
    verbose : bool
        Print the output of the code (exe backend).
    isolate : bool
        f2py backend: call the code in a child process (fresh Fortran
        state at each call, the calling process, e.g. a Jupyter kernel,
        is protected from a Fortran STOP).

    Returns
    -------
    xarray.Dataset
        Stokes parameters of the upward (``*_up``) and downward
        (``*_down``) fields versus ``vza`` and ``azi``, fluxes and the
        input parameters as attributes.
    """
    root = sos_root(root)
    os.environ["SOS_ABS_ROOT"] = root
    if backend == "auto":
        backend = available_backend(root)

    par = P.full_parameters(params)
    par["itrphi"] = 2                      # polar diagram
    if par["pas_phi"] in (None, P.NOT_INT):
        par["pas_phi"] = 5

    tmp = None
    if workdir is None:
        tmp = workdir = tempfile.mkdtemp(prefix="sosabs_")
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    par["resroot"] = workdir
    if not par["dir_mie"]:
        par["dir_mie"] = os.path.join(workdir, "MIE")
    if not par["dir_surf"]:
        par["dir_surf"] = os.path.join(workdir, "SURF")

    try:
        if backend == "f2py" and isolate:
            ctx = mp.get_context("fork")
            with ctx.Pool(1) as pool:
                ds = pool.apply(_run_f2py, (par, root))
        elif backend == "f2py":
            ds = _run_f2py(par, root)
        elif backend == "exe":
            ds = _run_exe(par, root, verbose)
        else:
            raise ValueError("unknown backend %r" % backend)
        angles = _read_used_angles(
            os.path.join(workdir, "SOS", par["ficangles_res_lum"]))
        ds = _add_angle_info(ds, angles)
    finally:
        if tmp is not None and not keep:
            shutil.rmtree(tmp, ignore_errors=True)

    ds.attrs.update(_global_attrs(par, backend))
    return ds


class _console_to_file:
    """Redirect the file descriptors 1 and 2 (Fortran output) to a file."""

    def __init__(self, path):
        self.path = path

    def __enter__(self):
        sys.stdout.flush()
        sys.stderr.flush()
        self.saved = (os.dup(1), os.dup(2))
        fd = os.open(self.path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        os.dup2(fd, 1)
        os.dup2(fd, 2)
        os.close(fd)
        return self

    def __exit__(self, *exc):
        os.dup2(self.saved[0], 1)
        os.dup2(self.saved[1], 2)
        for fd in self.saved:
            os.close(fd)
        return False


def _log_tail(path, n=15):
    try:
        with open(path, errors="replace") as f:
            lines = [ln.rstrip() for ln in f if ln.strip()]
        return "\n".join(lines[-n:])
    except OSError:
        return ""


def _run_f2py(par, root):
    sos = _import_f2py(root)
    args = {name: par[name] for name, _, _ in P.PARAMETERS}
    log = os.path.join(par["resroot"], LOGFILE)
    with _console_to_file(log):
        out = sos.sos_proc(**args, ier=0, trace=0)
    (nblum, ind_angout, phi, vza,
     sca_up, i_up, q_up, u_up, _pa_up, _pr_up, _lp_up,
     sca_dn, i_dn, q_dn, u_dn, _pa_dn, _pr_dn, _lp_dn,
     tdir, fdiff_down, ftot_down, fdiff_up, coef_tronca) = out
    if nblum <= 0 or not np.isfinite(i_up[0, 0]):
        raise RuntimeError("SOS_PROC failed:\n" + _log_tail(log))
    nphi = len(range(0, 361, int(par["pas_phi"])))
    sl = (slice(0, nphi), slice(0, nblum))
    data = dict(
        I_up=i_up[sl], Q_up=q_up[sl], U_up=u_up[sl],
        I_down=i_dn[sl], Q_down=q_dn[sl], U_down=u_dn[sl],
        scat_angle_up=sca_up[sl], scat_angle_down=sca_dn[sl])
    fluxes = dict(flux_dir_down=tdir, flux_diff_down=fdiff_down,
                  flux_tot_down=ftot_down, flux_diff_up=fdiff_up,
                  coef_tronca=coef_tronca)
    return _make_dataset(phi[:nphi], vza[:nblum], data, fluxes)


def _cmd_value(v):
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _run_exe(par, root, verbose):
    exe = os.path.join(root, "exe", "SOS_ABS_MAIN.exe")
    if not os.path.isfile(exe):
        raise RuntimeError("executable not found: %s (run make)" % exe)
    cmd = [exe]
    for name, kw, default in P.PARAMETERS:
        value = par[name]
        # give only the defined parameters, the executable has the
        # same defaults
        if value is None or value == default:
            if name not in ("resroot", "itrphi", "pas_phi"):
                continue
        cmd += [kw, _cmd_value(value)]
    cmd += ["-SOS.ResFileUp", "SOS_Up.txt", "-SOS.ResFileDown",
            "SOS_Down.txt", "-SOS.Flux", "SOS_Flux.txt"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    if verbose:
        print(proc.stdout)
    sosdir = os.path.join(par["resroot"], "SOS")
    up = os.path.join(sosdir, "SOS_Up.txt")
    if proc.returncode != 0 or not os.path.isfile(up):
        raise RuntimeError("SOS_ABS_MAIN failed:\n" + proc.stdout[-3000:])

    def load(fname):
        a = np.loadtxt(os.path.join(sosdir, fname), comments="#")
        phi = np.unique(a[:, 0])
        nvza = a.shape[0] // len(phi)
        a = a.reshape(len(phi), nvza, -1)
        return phi, a[0, :, 1], a

    phi, vza, up = load("SOS_Up.txt")
    _, _, dn = load("SOS_Down.txt")
    data = dict(
        I_up=up[..., 3], Q_up=up[..., 4], U_up=up[..., 5],
        I_down=dn[..., 3], Q_down=dn[..., 4], U_down=dn[..., 5],
        scat_angle_up=up[..., 2], scat_angle_down=dn[..., 2])
    fluxes = _read_flux_file(os.path.join(sosdir, "SOS_Flux.txt"))
    return _make_dataset(phi, vza, data, fluxes)


def _read_flux_file(fname):
    txt = open(fname).read()

    def get(pattern):
        m = re.search(pattern + r"\s*:?\s*([-+0-9.Ee]+)", txt)
        return float(m.group(1)) if m else np.nan

    return dict(
        flux_dir_down=get(r"Downward direct flux at BOA"),
        flux_diff_down=get(r"Downward diffuse flux at BOA"),
        flux_tot_down=get(r"Downward total flux at BOA"),
        flux_diff_up=get(r"Upward diffuse flux at TOA \(normalized by "
                         r"TOA solar flux\)"),
        coef_tronca=np.nan)


def _read_used_angles(fname):
    """Cosines, weights and user flags of the angles (SOS_UsedAngles.txt)."""
    rows = []
    with open(fname) as f:
        for line in f:
            parts = line.replace("D", "E").split()
            if len(parts) == 4 and parts[0].isdigit():
                rows.append((float(parts[1]), float(parts[2]),
                             int(parts[3])))
    return np.array(rows)


# ---------------------------------------------------------------------------
# Dataset construction
# ---------------------------------------------------------------------------

VAR_ATTRS = {
    "I": ("Stokes parameter I", "normalised radiance pi*L/E0"),
    "Q": ("Stokes parameter Q", "normalised radiance pi*Q/E0"),
    "U": ("Stokes parameter U", "normalised radiance pi*U/E0"),
}

NORMALISATION = (
    "Stokes parameters normalised to the extraterrestrial solar "
    "irradiance E0 (on a surface normal to the solar beam): X = pi*L/E0. "
    "The reflectance is rho = X / cos(sza).")

AZIMUTH_CONVENTION = (
    "Relative azimuth between the viewing direction (line of sight) and "
    "the solar azimuth: 0 deg <-> sensor and Sun in opposite half-planes "
    "(forward scattering, sun glint for the upward field, solar aureole "
    "for the downward field); 180 deg <-> sensor and Sun in the same "
    "half-plane (backscattering).")


def _make_dataset(phi, vza, data, fluxes):
    coords = dict(
        vza=("vza", np.asarray(vza, float),
             dict(units="degree", long_name="viewing zenith angle",
                  description="zenith angle of the line of sight: from "
                  "nadir for the upward field, from zenith for the "
                  "downward field")),
        azi=("azi", np.asarray(phi, float),
             dict(units="degree", long_name="relative azimuth angle",
                  description=AZIMUTH_CONVENTION)))
    dvars = {}
    for name, arr in data.items():
        arr = np.asarray(arr, float).T          # (vza, azi)
        if name.startswith("scat_angle"):
            attrs = dict(units="degree", long_name="scattering angle (%s "
                         "field)" % name.split("_")[-1])
        else:
            stokes, field = name.split("_")
            ln, desc = VAR_ATTRS[stokes]
            attrs = dict(units="1",
                         long_name="%s, %sward field" % (ln, field),
                         description=desc + "; " + NORMALISATION)
        dvars[name] = (("vza", "azi"), arr, attrs)
    flux_attrs = dict(
        flux_dir_down="direct downward flux at the surface (direct "
                      "transmittance)",
        flux_diff_down="diffuse downward flux at the surface",
        flux_tot_down="total downward flux at the surface",
        flux_diff_up="diffuse upward flux at TOA",
        coef_tronca="truncation coefficient of the aerosol phase "
                    "function")
    for name, value in fluxes.items():
        attrs = dict(long_name=flux_attrs[name], units="1")
        if name.startswith("flux"):
            attrs["description"] = ("normalised to the solar flux at TOA "
                                    "(E0*cos(sza))")
        dvars[name] = ((), float(value), attrs)
    return xr.Dataset(dvars, coords=coords)


def _add_angle_info(ds, angles):
    """Add the Gauss weights and user flags along vza."""
    if angles.size == 0 or angles.shape[0] < ds.sizes["vza"]:
        return ds
    # SOS_UsedAngles.txt is sorted by decreasing cosine = increasing vza
    a = angles[:ds.sizes["vza"]]
    ds = ds.assign_coords(
        gauss_weight=("vza", a[:, 1],
                      dict(long_name="Gauss quadrature weight",
                           description="0 for the angles added to the "
                           "Gauss angles (solar zenith angle, user "
                           "angles)")),
        user_angle=("vza", a[:, 2].astype("int8"),
                    dict(long_name="1 for the user angles")))
    return ds


# file names of the code (not useful as metadata)
_INTERNAL_FILES = (
    "ficmain_log", "ficangles_res_lum", "ficangles_res_mie", "ficanglog",
    "ficgranu_log", "ficmie_log", "ficgranu", "ficprofil_log",
    "ficsurf_log", "ficsos_log", "ficsos_res_bin", "fictrans", "ficflux")


def _global_attrs(par, backend):
    from . import __version__
    attrs = dict(
        title="SOS-ABS radiative transfer simulation",
        source="SOS-ABS (Successive Orders of Scattering coupled with "
               "gaseous ABSorption), CNES",
        references="Lenoble et al. (2007), JQSRT 107(3), "
                   "doi:10.1016/j.jqsrt.2007.03.010",
        sosabs_python_version=__version__,
        sos_backend=backend,
        normalisation=NORMALISATION,
        azimuth_convention=AZIMUTH_CONVENTION,
    )
    for name, kw, default in P.PARAMETERS:
        if name in ("resroot", "dir_mie", "dir_surf") + _INTERNAL_FILES:
            continue
        if name == "ficsurf" and par[name] == "DEFAULT":
            continue
        value = par[name]
        if (value is None or value in (P.NOT_INT, P.NOT_DBLE)
                or (isinstance(value, str) and value.startswith("NO_"))):
            continue
        attrs["sos" + kw.replace("-", "_").replace(".", "_")] = (
            value if isinstance(value, (int, float, str)) else str(value))
    return attrs
