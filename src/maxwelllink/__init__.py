'''
Lazy-loading package for MaxwellLink to avoid importing heavy dependencies, such as FDTD libraries
'''

from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("maxwelllink")
except PackageNotFoundError:
    __version__ = "0.0.0"

# Declare names for linters/docs, but don't import:
__all__ = [
    "TLSMolecule", "SocketMolecule", "update_molecules",
    "update_molecules_no_mpi", "update_molecules_no_socket",
    "SocketHub", "get_available_host_port",
    "mxl_driver_main", "launch_driver", "terminate_driver",
    "mxl_lammps_main",
]

# Lazy attribute loader: import submodules *only when accessed*.
def __getattr__(name):
    if name in {"TLSMolecule", "SocketMolecule",
                "update_molecules", "update_molecules_no_mpi", "update_molecules_no_socket"}:
        from .molecule_fast import (
            TLSMolecule, SocketMolecule, update_molecules,
            update_molecules_no_mpi, update_molecules_no_socket
        )
        return locals()[name]
    if name in {"SocketHub", "get_available_host_port"}:
        from .sockets import SocketHub, get_available_host_port
        return locals()[name]
    if name in {"mxl_driver_main", "launch_driver", "terminate_driver"}:
        from .mxl_drivers.python.mxl_driver import (
            mxl_driver_main, launch_driver, terminate_driver
        )
        return locals()[name]
    if name == "mxl_lammps_main":
        from .mxl_drivers.lammps.install import mxl_lammps_main
        return mxl_lammps_main
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
