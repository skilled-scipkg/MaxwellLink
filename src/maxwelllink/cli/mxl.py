# --------------------------------------------------------------------------------------#
# Copyright (c) 2026 MaxwellLink                                                       #
# This file is part of MaxwellLink. Repository: https://github.com/TaoELi/MaxwellLink  #
# If you use this code, always credit and cite arXiv:2512.06173.                       #
# See AGENTS.md and README.md for details.                                             #
# --------------------------------------------------------------------------------------#

"""Top-level ``mxl`` command dispatcher for MaxwellLink CLI actions."""

from __future__ import annotations

import argparse

from .mxl_clean import mxl_clean_main
from .mxl_hpc import mxl_hpc_main
from .mxl_init import mxl_init_main


def main(argv: list[str] | None = None) -> int:
    """Run the ``mxl`` command dispatcher.

    Parameters
    ----------
    argv : list of str or None, default=None
        Optional command-line arguments. When ``None``, uses ``sys.argv``.

    Returns
    -------
    int
        Exit status code from the selected subcommand.
    """
    parser = argparse.ArgumentParser(
        prog="mxl",
        description="MaxwellLink convenience CLI.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser(
        "init",
        help="Initialize current folder for MaxwellLink agent workflows.",
    )
    init_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite conflicting files/symlinks in the destination directory.",
    )

    clean_parser = subparsers.add_parser(
        "clean",
        help="Clean workspace artifacts created by `mxl init` / `mxl-init`.",
    )
    clean_parser.add_argument(
        "--force",
        action="store_true",
        help="Remove conflicting managed paths even when they were modified.",
    )

    hpc_parser = subparsers.add_parser(
        "hpc",
        help="Manage persistent HPC profile settings.",
    )
    hpc_subparsers = hpc_parser.add_subparsers(dest="hpc_command", required=True)
    hpc_set_parser = hpc_subparsers.add_parser(
        "set",
        help="Install a global HPC profile from JSON.",
    )
    hpc_set_parser.add_argument(
        "file",
        type=str,
        help="Path to a JSON profile file.",
    )

    args = parser.parse_args(argv)
    if args.command == "init":
        passthrough = ["--force"] if args.force else []
        return mxl_init_main(passthrough)
    if args.command == "clean":
        passthrough = ["--force"] if args.force else []
        return mxl_clean_main(passthrough)
    if args.command == "hpc":
        if args.hpc_command == "set":
            return mxl_hpc_main(["set", args.file])
        return mxl_hpc_main([args.hpc_command])

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
