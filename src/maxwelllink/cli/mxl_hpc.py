"""Manage persistent MaxwellLink HPC profile settings.

This module provides the implementation behind ``mxl hpc set``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from . import mxl_init

_REQUIRED_HPC_KEYS = (
    "slurm_default_partition",
    "slurm_defaults",
    "slurm_resource_policy",
)


def _load_hpc_profile(path: Path) -> dict:
    """Load an HPC profile JSON file.

    Parameters
    ----------
    path : pathlib.Path
        Path to JSON file.

    Returns
    -------
    dict
        Parsed JSON object.
    """
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"HPC profile must be a JSON object: {path}")
    return data


def _validate_hpc_profile(data: dict) -> None:
    """Validate required keys and value types of an HPC profile.

    Parameters
    ----------
    data : dict
        Parsed profile object.

    Raises
    ------
    ValueError
        If required keys are missing or values are not strings.
    """
    missing = [key for key in _REQUIRED_HPC_KEYS if key not in data]
    if missing:
        raise ValueError(f"Missing required HPC profile keys: {', '.join(missing)}")

    bad = [key for key in _REQUIRED_HPC_KEYS if not isinstance(data[key], str)]
    if bad:
        raise ValueError("HPC profile keys must map to strings: " + ", ".join(bad))


def set_hpc_profile(source_file: Path, destination_file: Path | None = None) -> Path:
    """Install a persistent global HPC profile.

    Parameters
    ----------
    source_file : pathlib.Path
        User-provided JSON profile.
    destination_file : pathlib.Path or None, default=None
        Install destination. When ``None``, uses
        ``~/.maxwelllink/HPC_PROFILE.json``.

    Returns
    -------
    pathlib.Path
        Installed destination path.
    """
    if not source_file.exists():
        raise FileNotFoundError(f"Profile file not found: {source_file}")

    data = _load_hpc_profile(source_file)
    _validate_hpc_profile(data)

    destination = destination_file or mxl_init._global_hpc_profile_path()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return destination


def mxl_hpc_main(argv: list[str] | None = None) -> int:
    """Run the ``mxl hpc`` CLI command family.

    Parameters
    ----------
    argv : list of str or None, default=None
        Optional command-line arguments. When ``None``, uses ``sys.argv``.

    Returns
    -------
    int
        ``0`` on success, otherwise non-zero on failure.
    """
    parser = argparse.ArgumentParser(
        prog="mxl hpc",
        description="Manage persistent HPC profile settings for MaxwellLink.",
    )
    subparsers = parser.add_subparsers(dest="hpc_command", required=True)

    set_parser = subparsers.add_parser(
        "set",
        help="Install a user HPC profile to ~/.maxwelllink/HPC_PROFILE.json",
    )
    set_parser.add_argument(
        "file",
        type=str,
        help="Path to a JSON profile file.",
    )

    args = parser.parse_args(argv)

    if args.hpc_command == "set":
        try:
            src = Path(args.file).expanduser().resolve()
            dest = set_hpc_profile(src)
        except Exception as exc:
            print(f"[mxl-hpc] ERROR: {exc}", file=sys.stderr)
            return 2

        print("[mxl-hpc] Installed profile at", dest)
        return 0

    parser.error(f"Unknown hpc command: {args.hpc_command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(mxl_hpc_main())
