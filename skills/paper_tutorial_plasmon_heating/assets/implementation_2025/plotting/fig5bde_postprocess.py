#!/usr/bin/env python3
"""Scoped postprocess for Fig.5b/5d/5e reproduction.

Expected runtime layout:
  implementation_2025/
    meep_plasmon_empty/
      vac/flux0_a2.79.dat
      no_mol_with_dielectric/flux_a2.79_r1.11.dat
    meep_plasmon_HCN_excitation_tls_strong/nmol_256_with_dielectric/mol_*.npz
    meep_plasmon_HCN_excitation_bomd_strong/nmol_256_with_dielectric/mol_*.npz
    plotting/fig5bde_postprocess.py
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

HCN_STRETCH_CM = 3465.5930
MAP_EXTENT_UM = 1.4


def load_flux_spectrum(base_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    vac_file = base_dir / "meep_plasmon_empty" / "vac" / "flux0_a2.79.dat"
    geom_file = (
        base_dir
        / "meep_plasmon_empty"
        / "no_mol_with_dielectric"
        / "flux_a2.79_r1.11.dat"
    )
    if not vac_file.is_file():
        raise FileNotFoundError(f"Missing spectrum reference file: {vac_file}")
    if not geom_file.is_file():
        raise FileNotFoundError(f"Missing spectrum geometry file: {geom_file}")

    f0 = np.genfromtxt(vac_file, delimiter=",")
    f = np.genfromtxt(geom_file, delimiter=",")

    if f0.ndim != 2 or f.ndim != 2 or f0.shape[1] < 2 or f.shape[1] < 2:
        raise ValueError("Flux files must be two-column arrays with frequency and flux columns.")
    if f0.shape[0] != f.shape[0]:
        raise ValueError(
            f"Flux grid mismatch: vac points={f0.shape[0]}, geometry points={f.shape[0]}"
        )

    refl = -f[:, 1] / f0[:, 1]
    absorption = 1.0 - refl
    freq_cm = f0[:, 0] * 1e4
    return freq_cm, absorption


def gather_molecular_gains(
    base_dir: Path, branch: str, nmol: int
) -> Tuple[np.ndarray, List[int], Path]:
    case_dir = base_dir / branch / f"nmol_{nmol}_with_dielectric"
    side = int(round(math.sqrt(nmol)))
    if side * side != nmol:
        raise ValueError(f"nmol must be a perfect square, got nmol={nmol}")

    gains: List[float] = []
    missing: List[int] = []
    for idx in range(nmol):
        path = case_dir / f"mol_{idx}_data.npz"
        if not path.is_file():
            missing.append(idx)
            continue
        data = np.load(path)
        if "energy_au" not in data:
            raise KeyError(f"energy_au not found in {path}")
        energy = np.asarray(data["energy_au"], dtype=float)
        gains.append(float(np.mean(energy - energy[0])))

    return np.asarray(gains, dtype=float), missing, case_dir


def gains_to_map(gains: np.ndarray, nmol: int) -> np.ndarray:
    side = int(round(math.sqrt(nmol)))
    grid = gains.reshape((side, side)).transpose()
    grid = grid - np.min(grid)
    return grid


def anisotropy_ratio(energy_map: np.ndarray) -> Tuple[float, float, float]:
    side = energy_map.shape[0]
    coords = np.linspace(-MAP_EXTENT_UM, MAP_EXTENT_UM, side)
    xx, yy = np.meshgrid(coords, coords, indexing="xy")

    y_gap_mask = (np.abs(xx) < 0.2) & (np.abs(yy) > 0.8)
    x_gap_mask = (np.abs(yy) < 0.2) & (np.abs(xx) > 0.8)

    y_gap_mean = float(np.mean(energy_map[y_gap_mask]))
    x_gap_mean = float(np.mean(energy_map[x_gap_mask]))
    if x_gap_mean <= 0.0:
        return float("inf"), y_gap_mean, x_gap_mean
    return y_gap_mean / x_gap_mean, y_gap_mean, x_gap_mean


def dominant_hotspot(energy_map: np.ndarray) -> Dict[str, float | bool]:
    side = energy_map.shape[0]
    coords = np.linspace(-MAP_EXTENT_UM, MAP_EXTENT_UM, side)
    flat_idx = int(np.argmax(energy_map))
    y_idx, x_idx = np.unravel_index(flat_idx, energy_map.shape)
    x_peak = float(coords[x_idx])
    y_peak = float(coords[y_idx])
    hotspot_yedge_pass = bool(abs(y_peak) >= 0.8)
    hotspot_xcenter_pass = bool(abs(x_peak) <= 0.35)
    return {
        "hotspot_x_um": x_peak,
        "hotspot_y_um": y_peak,
        "hotspot_yedge_pass": hotspot_yedge_pass,
        "hotspot_xcenter_pass": hotspot_xcenter_pass,
        "hotspot_yedge_centered_pass": bool(hotspot_yedge_pass and hotspot_xcenter_pass),
    }


def summarize_spectrum(freq_cm: np.ndarray, absorption: np.ndarray) -> Dict[str, float]:
    peak_idx = int(np.argmax(absorption))
    return {
        "n_points": int(freq_cm.size),
        "peak_cm": float(freq_cm[peak_idx]),
        "max_absorption": float(absorption[peak_idx]),
        "absorption_at_hcn": float(np.interp(HCN_STRETCH_CM, freq_cm, absorption)),
    }


def summarize_map(
    gains: np.ndarray, energy_map: np.ndarray, threshold: float
) -> Dict[str, float | bool]:
    ratio, y_gap_mean, x_gap_mean = anisotropy_ratio(energy_map)
    summary: Dict[str, float | bool] = {
        "mean_gain": float(np.mean(gains)),
        "min_gain": float(np.min(gains)),
        "max_gain": float(np.max(gains)),
        "anisotropy_ratio_y_over_x": float(ratio),
        "anisotropy_threshold": float(threshold),
        "anisotropy_pass": bool(ratio >= threshold),
        "y_gap_mean": float(y_gap_mean),
        "x_gap_mean": float(x_gap_mean),
    }
    summary.update(dominant_hotspot(energy_map))
    return summary


def add_spectrum_checks(metrics: Dict[str, object]) -> Dict[str, bool]:
    s = metrics["spectrum"]
    checks = {
        "spectrum_grid_has_200_points": bool(s["n_points"] == 200),
        "spectrum_peak_in_3450_3535_cm": bool(3450.0 <= s["peak_cm"] <= 3535.0),
        "spectrum_max_absorption_gt_0p90": bool(s["max_absorption"] > 0.90),
        "spectrum_absorption_at_hcn_gt_0p70": bool(s["absorption_at_hcn"] > 0.70),
    }
    return checks


def make_figure(
    mode: str,
    figure_out: Path,
    spectrum_data: Tuple[np.ndarray, np.ndarray] | None,
    tls_map: np.ndarray | None,
    bomd_map: np.ndarray | None,
) -> None:
    panels: List[str] = []
    if spectrum_data is not None:
        panels.append("spectrum")
    if tls_map is not None:
        panels.append("tls")
    if bomd_map is not None:
        panels.append("bomd")

    if not panels:
        return

    fig, axes = plt.subplots(1, len(panels), figsize=(4.4 * len(panels), 3.8))
    if len(panels) == 1:
        axes = [axes]

    for ax, panel in zip(axes, panels):
        if panel == "spectrum":
            freq_cm, absorption = spectrum_data
            ax.plot(freq_cm, absorption, color="tab:red", lw=1.8)
            ax.axvline(HCN_STRETCH_CM, color="tab:green", ls="--", lw=1.2)
            ax.set_xlabel("frequency (cm^-1)")
            ax.set_ylabel("absorption")
            ax.set_xlim(2000, 5000)
            ax.set_ylim(0.0, 1.05)
            ax.set_title("Figure 5b scope")
        elif panel == "tls":
            im = ax.imshow(
                tls_map,
                origin="lower",
                extent=[-MAP_EXTENT_UM, MAP_EXTENT_UM, -MAP_EXTENT_UM, MAP_EXTENT_UM],
                cmap="viridis",
            )
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="energy gain [a.u.]")
            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)")
            ax.set_title("Figure 5d scope (TLS)")
        elif panel == "bomd":
            im = ax.imshow(
                bomd_map,
                origin="lower",
                extent=[-MAP_EXTENT_UM, MAP_EXTENT_UM, -MAP_EXTENT_UM, MAP_EXTENT_UM],
                cmap="viridis",
            )
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="energy gain [a.u.]")
            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)")
            ax.set_title("Figure 5e scope (BOMD)")

    fig.tight_layout()
    fig.savefig(figure_out, dpi=220)


def evaluate(mode: str, base_dir: Path, nmol: int) -> Tuple[Dict[str, object], Tuple[np.ndarray, np.ndarray] | None, np.ndarray | None, np.ndarray | None]:
    metrics: Dict[str, object] = {
        "mode": mode,
        "nmol": nmol,
        "base_dir": str(base_dir),
        "checks": {},
    }

    spectrum_data = None
    tls_map = None
    bomd_map = None

    include_spectrum = mode in {"all", "fig5b"}
    include_tls = mode in {"all", "fig5d"}
    include_bomd = mode in {"all", "fig5e"}

    if include_spectrum:
        freq_cm, absorption = load_flux_spectrum(base_dir)
        spectrum_data = (freq_cm, absorption)
        metrics["spectrum"] = summarize_spectrum(freq_cm, absorption)
        metrics["checks"].update(add_spectrum_checks(metrics))

    tls_mean_gain = None
    if include_tls or mode == "fig5e":
        tls_gains, tls_missing, tls_case_dir = gather_molecular_gains(
            base_dir, "meep_plasmon_HCN_excitation_tls_strong", nmol
        )
        metrics["tls"] = {
            "case_dir": str(tls_case_dir),
            "present_count": int(tls_gains.size),
            "missing_count": int(len(tls_missing)),
            "missing_indices_head": tls_missing[:10],
        }
        metrics["checks"]["tls_family_complete"] = bool(len(tls_missing) == 0)
        if len(tls_missing) == 0:
            tls_map = gains_to_map(tls_gains, nmol)
            tls_summary = summarize_map(tls_gains, tls_map, threshold=3.0)
            metrics["tls"].update(tls_summary)
            metrics["checks"]["tls_nonnegative_mean_gains"] = bool(tls_summary["min_gain"] >= -1e-12)
            metrics["checks"]["tls_anisotropy_ge_3"] = bool(tls_summary["anisotropy_pass"])
            metrics["checks"]["tls_hotspot_yedge_centered"] = bool(
                tls_summary["hotspot_yedge_centered_pass"]
            )
            tls_mean_gain = float(tls_summary["mean_gain"])

    bomd_mean_gain = None
    if include_bomd:
        bomd_gains, bomd_missing, bomd_case_dir = gather_molecular_gains(
            base_dir, "meep_plasmon_HCN_excitation_bomd_strong", nmol
        )
        metrics["bomd"] = {
            "case_dir": str(bomd_case_dir),
            "present_count": int(bomd_gains.size),
            "missing_count": int(len(bomd_missing)),
            "missing_indices_head": bomd_missing[:10],
        }
        metrics["checks"]["bomd_family_complete"] = bool(len(bomd_missing) == 0)
        if len(bomd_missing) == 0:
            bomd_map = gains_to_map(bomd_gains, nmol)
            bomd_summary = summarize_map(bomd_gains, bomd_map, threshold=4.0)
            metrics["bomd"].update(bomd_summary)
            metrics["checks"]["bomd_nonnegative_mean_gains"] = bool(
                bomd_summary["min_gain"] >= -1e-12
            )
            metrics["checks"]["bomd_anisotropy_ge_4"] = bool(bomd_summary["anisotropy_pass"])
            metrics["checks"]["bomd_hotspot_yedge_centered"] = bool(
                bomd_summary["hotspot_yedge_centered_pass"]
            )
            bomd_mean_gain = float(bomd_summary["mean_gain"])

    if include_bomd and tls_mean_gain is not None and bomd_mean_gain is not None:
        metrics["checks"]["crosspanel_bomd_mean_gt_tls_mean"] = bool(
            bomd_mean_gain > tls_mean_gain
        )

    return metrics, spectrum_data, tls_map, bomd_map


def main() -> int:
    parser = argparse.ArgumentParser(description="Postprocess Figure 5b/5d/5e scoped outputs.")
    parser.add_argument(
        "--mode",
        choices=["all", "fig5b", "fig5d", "fig5e"],
        default="all",
        help="Which scoped figure set to evaluate.",
    )
    parser.add_argument(
        "--base-dir",
        default="..",
        help="Path to implementation_2025 root relative to this script run location.",
    )
    parser.add_argument("--nmol", type=int, default=256, help="Molecule count in coupled runs.")
    parser.add_argument(
        "--figure-out",
        default="fig5bde_scope.pdf",
        help="Output figure path. Use empty string to skip figure generation.",
    )
    parser.add_argument(
        "--report-out",
        default="fig5bde_metrics.json",
        help="Output JSON metrics path.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return nonzero if any recorded check is false.",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    metrics, spectrum_data, tls_map, bomd_map = evaluate(args.mode, base_dir, args.nmol)

    report_out = Path(args.report_out)
    report_out.parent.mkdir(parents=True, exist_ok=True)
    with report_out.open("w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    if args.figure_out:
        figure_out = Path(args.figure_out)
        figure_out.parent.mkdir(parents=True, exist_ok=True)
        make_figure(args.mode, figure_out, spectrum_data, tls_map, bomd_map)

    print(json.dumps(metrics, indent=2))

    if args.strict:
        checks = metrics.get("checks", {})
        failed = [name for name, passed in checks.items() if not bool(passed)]
        if failed:
            print("Failed checks:")
            for item in failed:
                print(f"- {item}")
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
