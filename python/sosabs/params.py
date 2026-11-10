"""SOS-ABS input parameters.

Every parameter of the code has a command-line keyword (``-AP.H2O``) and a
Python name (``h2o``, the argument name of the f2py routine ``sos_proc``).
The table :data:`PARAMETERS` gives both names and the default value used
when the parameter is not given (the same defaults as ``SOS_ABS_MAIN``).

In :mod:`sosabs`, parameters can be given with either name::

    {"-SOS_Main.Wa": 0.94, "-AP.H2O": 2.0}
    {"wa_simu": 0.94, "h2o": 2.0}
"""

NOT_INT = -999      # CTE_NOT_DEFINED_VALUE_INT
NOT_DBLE = -999.0   # CTE_NOT_DEFINED_VALUE_DBLE

# (python / f2py name, command-line keyword, default)
# The order is the argument order of the f2py routine sos.sos_proc.
PARAMETERS = [
    ("resroot", "-SOS_Main.ResRoot", None),
    ("ficmain_log", "-SOS_Main.Log", "NO_LOG_FILE"),
    ("wa_simu", "-SOS_Main.Wa", NOT_DBLE),
    ("nbmu_gauss_lum", "-ANG.Rad.NbGauss", NOT_INT),
    ("ficangles_user_lum", "-ANG.Rad.UserAngFile", "NO_USER_ANGLES"),
    ("tetas", "-ANG.Thetas", NOT_DBLE),
    ("ficangles_res_lum", "-ANG.Rad.ResFile", "SOS_UsedAngles.txt"),
    ("nbmu_gauss_mie", "-ANG.Aer.NbGauss", NOT_INT),
    ("ficangles_user_mie", "-ANG.Aer.UserAngFile", "NO_USER_ANGLES"),
    ("ficangles_res_mie", "-ANG.Aer.ResFile", "Aer_UsedAngles.txt"),
    ("ficanglog", "-ANG.Log", "NO_LOG_FILE"),
    ("waref_aot", "-AER.Waref", NOT_DBLE),
    ("aot_ref", "-AER.AOTref", NOT_DBLE),
    ("itronc_aer", "-AER.Tronca", 1),
    ("ficgranu_log", "-AER.Log", "NO_LOG_FILE"),
    ("ficmie_log", "-AER.MieLog", "NO_LOG_FILE"),
    ("dir_mie", "-AER.DirMie", None),
    ("ficgranu", "-AER.ResFile", "Aerosols.txt"),
    ("imod_aer", "-AER.Model", NOT_INT),
    ("rn_wa", "-AER.MMD.MRwa", NOT_DBLE),
    ("in_wa", "-AER.MMD.MIwa", NOT_DBLE),
    ("rn_waref", "-AER.MMD.MRwaref", NOT_DBLE),
    ("in_waref", "-AER.MMD.MIwaref", NOT_DBLE),
    ("igranu", "-AER.MMD.SDtype", NOT_INT),
    ("lnd_radius_mmd_aer", "-AER.MMD.LNDradius", NOT_DBLE),
    ("lnd_lnvar_mmd_aer", "-AER.MMD.LNDvar", NOT_DBLE),
    ("jd_slope_mmd_aer", "-AER.MMD.JD.slope", NOT_DBLE),
    ("jd_rmin_mmd_aer", "-AER.MMD.JD.rmin", NOT_DBLE),
    ("jd_rmax_mmd_aer", "-AER.MMD.JD.rmax", 50.0),
    ("imodele_wmo", "-AER.WMO.Model", NOT_INT),
    ("c_wmo_dl", "-AER.WMO.DL", NOT_DBLE),
    ("c_wmo_ws", "-AER.WMO.WS", NOT_DBLE),
    ("c_wmo_oc", "-AER.WMO.OC", NOT_DBLE),
    ("c_wmo_so", "-AER.WMO.SO", NOT_DBLE),
    ("imodele_sf", "-AER.SF.Model", NOT_INT),
    ("rh", "-AER.SF.RH", NOT_DBLE),
    ("mode_param_bilnd", "-AER.BMD.VCdef", NOT_INT),
    ("user_cv_coarse", "-AER.BMD.CoarseVC", NOT_DBLE),
    ("user_cv_fine", "-AER.BMD.FineVC", NOT_DBLE),
    ("rtauct_waref", "-AER.BMD.RAOT", NOT_DBLE),
    ("bmd_cm_mrwa", "-AER.BMD.CM.MRwa", NOT_DBLE),
    ("bmd_cm_miwa", "-AER.BMD.CM.MIwa", NOT_DBLE),
    ("bmd_cm_mrwaref", "-AER.BMD.CM.MRwaref", NOT_DBLE),
    ("bmd_cm_miwaref", "-AER.BMD.CM.MIwaref", NOT_DBLE),
    ("bmd_cm_rmodal", "-AER.BMD.CM.SDradius", NOT_DBLE),
    ("bmd_cm_var", "-AER.BMD.CM.SDvar", NOT_DBLE),
    ("bmd_fm_mrwa", "-AER.BMD.FM.MRwa", NOT_DBLE),
    ("bmd_fm_miwa", "-AER.BMD.FM.MIwa", NOT_DBLE),
    ("bmd_fm_mrwaref", "-AER.BMD.FM.MRwaref", NOT_DBLE),
    ("bmd_fm_miwaref", "-AER.BMD.FM.MIwaref", NOT_DBLE),
    ("bmd_fm_rmodal", "-AER.BMD.FM.SDradius", NOT_DBLE),
    ("bmd_fm_var", "-AER.BMD.FM.SDvar", NOT_DBLE),
    ("ficextdata_aer", "-AER.ExtData", "NO_USER_AEROSOLS_PHAZE_FCT"),
    ("ficmixture_aer", "-AER.DefMixture", "NO_USER_AEROSOLS_MIXTURE"),
    ("ficuser_aer", "-AER.UserFile", "NO_USER_AEROSOLS"),
    ("ficprofil_log", "-AP.Log", "NO_LOG_FILE"),
    ("tr", "-AP.MOT", NOT_DBLE),
    ("hr", "-AP.HR", NOT_DBLE),
    ("ha", "-AP.AerHS.HA", NOT_DBLE),
    ("iprofil", "-AP.AerProfile.Type", NOT_INT),
    ("zmin", "-AP.AerLayer.Zmin", NOT_DBLE),
    ("zmax", "-AP.AerLayer.Zmax", NOT_DBLE),
    ("psurf", "-AP.Psurf", 1013.0),
    ("h2o", "-AP.H2O", NOT_DBLE),
    ("o3", "-AP.O3", NOT_DBLE),
    ("co2", "-AP.CO2", NOT_DBLE),
    ("ch4", "-AP.CH4", NOT_DBLE),
    ("absprofil", "-AP.AbsProfile.Type", NOT_INT),
    ("ficabsprofil", "-AP.AbsProfile.UserFile", "NO_USER_ABS_PROFILE_FILE"),
    ("nustep", "-AP.SpectralResol", NOT_DBLE),
    ("isurf", "-SURF.Type", NOT_INT),
    ("dir_surf", "-SURF.Dir", None),
    ("ficsurf_log", "-SURF.Log", "NO_LOG_FILE"),
    ("surf_ind", "-SURF.Ind", NOT_DBLE),
    ("wind", "-SURF.Glitter.Wind", NOT_DBLE),
    ("k0_roujean", "-SURF.Roujean.K0", NOT_DBLE),
    ("k1_roujean", "-SURF.Roujean.K1", NOT_DBLE),
    ("k2_roujean", "-SURF.Roujean.K2", NOT_DBLE),
    ("alpha_nadal", "-SURF.Nadal.Alpha", NOT_DBLE),
    ("beta_nadal", "-SURF.Nadal.Beta", NOT_DBLE),
    ("coef_c_maignan", "-SURF.Maignan.C", NOT_DBLE),
    ("rho", "-SURF.Alb", NOT_DBLE),
    ("ficsurf", "-SURF.File", "DEFAULT"),
    ("ficsos_log", "-SOS.Log", "NO_LOG_FILE"),
    ("ficsos_res_bin", "-SOS.ResBin", "SOS_Result.bin"),
    ("fictrans", "-SOS.Trans", "NO_OUTPUT"),
    ("ficflux", "-SOS.Flux", "NO_OUTPUT"),
    ("zout", "-SOS.OutputAlt", -1.0),
    ("igmax", "-SOS.IGmax", 100),
    ("ipolar", "-SOS.Ipolar", 1),
    ("itrphi", "-SOS.View", 2),
    ("phios", "-SOS.View.Phi", NOT_DBLE),
    ("pas_phi", "-SOS.View.Dphi", 5),
    ("imode_ckd_calcul", "-SOS.AbsModeCKD", NOT_INT),
]

KEYWORD = {name: kw for name, kw, _ in PARAMETERS}
NAME = {kw: name for name, kw, _ in PARAMETERS}
DEFAULTS = {name: default for name, _, default in PARAMETERS}

# Short names of the LUT dimensions for the usual parameters
DIM_NAMES = {
    "wa_simu": "wl",
    "tetas": "sza",
    "h2o": "h2o",
    "o3": "o3",
    "aot_ref": "aot_ref",
    "wind": "wind",
    "psurf": "pressure",
    "zout": "z_out",
    "surf_ind": "n_surf",
    "rho": "albedo",
    "imodele_wmo": "wmo_model",
    "imodele_sf": "sf_model",
    "rh": "rh",
    "imode_ckd_calcul": "ckd_mode",
    "nustep": "ckd_resolution",
}

# Units and long names of the LUT dimensions
DIM_ATTRS = {
    "wa_simu": dict(units="um", long_name="wavelength of the simulation"),
    "tetas": dict(units="degree", long_name="solar zenith angle"),
    "h2o": dict(units="g cm-2", long_name="total column water vapour"),
    "o3": dict(units="Dobson", long_name="total column ozone"),
    "aot_ref": dict(units="1",
                    long_name="aerosol optical thickness at -AER.Waref"),
    "wind": dict(units="m s-1", long_name="wind speed (Cox and Munk)"),
    "psurf": dict(units="hPa", long_name="surface pressure"),
    "zout": dict(units="km", long_name="output altitude"),
    "surf_ind": dict(units="1",
                     long_name="refractive index surface / air"),
    "rho": dict(units="1", long_name="lambertian surface albedo"),
    "imode_ckd_calcul": dict(
        units="1", long_name="CKD absorption mode",
        description="1: radiance computed for each CKD term; "
                    "2: one radiance computation for the CKD-equivalent "
                    "absorption optical thickness"),
    "nustep": dict(units="cm-1",
                   long_name="spectral resolution of the CKD coefficients"),
}


def normalize(params):
    """Return a dict indexed by Python names from keywords or names.

    Raises KeyError for an unknown parameter.
    """
    out = {}
    for key, value in (params or {}).items():
        if key in DEFAULTS:
            out[key] = value
        elif key in NAME:
            out[NAME[key]] = value
        else:
            raise KeyError("unknown SOS-ABS parameter: %r" % key)
    return out


def full_parameters(params):
    """Defaults updated by ``params`` (Python names)."""
    full = dict(DEFAULTS)
    full.update(normalize(params))
    return full


def dim_name(name):
    """Short dimension name of a parameter (Python name or keyword)."""
    name = NAME.get(name, name)
    return DIM_NAMES.get(name, name)
