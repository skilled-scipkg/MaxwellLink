import numpy as np

fs_to_au = 41.341373335  # 1 fs in atomic units


def meep_to_atomic_units_E(Emu_vec3, time_units_fs: float = 0.1):
    """
    Convert the regularized E-field integral from MEEP units to atomic units.
    Kept here for backward compatibility; MEEP backend uses it directly.
    """
    factor = 1.2929541569381223e-6 / (time_units_fs**2)
    return np.asarray(Emu_vec3, dtype=float) * factor


def atomic_to_meep_units_SourceAmp(amp_au_vec3, time_units_fs: float = 0.1):
    """Convert source amplitude (a.u.) to MEEP units. (dipole/time)"""
    factor = 0.002209799779149953
    return np.asarray(amp_au_vec3, dtype=float) * factor
