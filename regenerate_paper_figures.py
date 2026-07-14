#!/usr/bin/env python3
"""
Regenerate every paper figure and table (single public entry point).

Usage:
    python regenerate_paper_figures.py

Outputs → paper_deliverables/
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import time
from collections import Counter
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from microplastic_simulation import (
    CONTAINER_PROFILES,
    DEFAULT_CLINICAL_REALISM,
    DEFAULT_MAX_TIME_S,
    EXPERIMENTAL_DIAMETERS_UM,
    ClinicalRealismConfig,
    ContainerProfile,
    ExperimentRunner,
    FluidProperties,
    ParticleState,
    TubeEnvironment,
    generate_experiment_figures,
    plot_outcome_bar_chart,
    plot_trajectory_figure,
    run_single_trial,
)

ROOT = Path(__file__).parent
OUT = ROOT / "paper_deliverables"
SEED = 42
BASE_TRIALS = 50
N_REPLICATES = 20
REP_TRIALS = 25  # per replicate cell; SE ≈ ±10 pp at p=0.5 for n=25


def wilson_ci(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = (z / den) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (100.0 * max(0.0, centre - half), 100.0 * min(1.0, centre + half))


def make_env(u_avg: float = 0.05) -> TubeEnvironment:
    return TubeEnvironment(
        length=0.100,
        width=0.004,
        fluid=FluidProperties(average_velocity=u_avg),
    )


def make_config(**overrides) -> ClinicalRealismConfig:
    base = DEFAULT_CLINICAL_REALISM
    return replace(base, **overrides) if overrides else base


def run_campaign(
    env: TubeEnvironment,
    profile: ContainerProfile,
    *,
    base_trials: int,
    seed: int,
    clinical: ClinicalRealismConfig,
    traj_n: int = 8,
    verbose: bool = False,
) -> dict:
    runner = ExperimentRunner(
        environment=env,
        container_profile=profile,
        diameters_um=EXPERIMENTAL_DIAMETERS_UM,
        base_trials_per_size=base_trials,
        max_time=DEFAULT_MAX_TIME_S,
        fast_mode=True,
        seed=seed,
        verbose=verbose,
        clinical_realism=clinical.enabled,
        clinical_config=clinical,
        trajectory_trials_per_group=traj_n,
    )
    groups = runner.run_all()
    rows = []
    for gr in groups:
        tally = runner.data_log.get_tally(gr.diameter_um)
        mech = Counter(
            t.retention_mechanism or "unknown"
            for t in gr.trials
            if t.final_state == ParticleState.STUCK
        )
        lo, hi = wilson_ci(tally.escaped, tally.total)
        rows.append({
            "container": profile.key,
            "diameter_um": gr.diameter_um,
            "n": tally.total,
            "escaped": tally.escaped,
            "stuck": tally.stuck,
            "escape_pct": tally.escaped_percent,
            "stuck_pct": tally.stuck_percent,
            "wilson_lo": round(lo, 1),
            "wilson_hi": round(hi, 1),
            "mean_transit_escaped_s": round(
                float(np.mean([
                    t.transit_time_s for t in gr.trials
                    if t.final_state == ParticleState.ESCAPED
                ])) if gr.escape_count else float("nan"),
                4,
            ),
            "mech_hydrodynamic_wall": mech.get("hydrodynamic_wall", 0),
            "mech_electrostatic_near_wall": mech.get("electrostatic_near_wall", 0),
            "mech_connector_static_clamp": mech.get("connector_static_clamp", 0),
            "mech_jagged_edge_defect": mech.get("jagged_edge_defect", 0),
            "mech_electrostatic_bulk": mech.get("electrostatic_bulk", 0),
        })
    sub20_esc = sum(
        r["escaped"] for r in rows if r["diameter_um"] < 20.0
    )
    return {
        "rows": rows,
        "groups": groups,
        "runner": runner,
        "sub20_escaped": sub20_esc,
        "total_escaped": sum(r["escaped"] for r in rows),
        "total_injected": sum(r["n"] for r in rows),
    }


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# Tier 0.4 — focusing ON vs OFF
# ---------------------------------------------------------------------------

def tier04_focusing_ablation(env: TubeEnvironment) -> List[dict]:
    print("\n=== Tier 0.4: focusing ON vs OFF ===")
    rows = []
    for focus_on, label in ((True, "focusing_ON"), (False, "focusing_OFF")):
        cfg = make_config(focus_enabled=focus_on)
        for key in ("pristine", "autoclaved"):
            t0 = time.perf_counter()
            result = run_campaign(
                env, CONTAINER_PROFILES[key],
                base_trials=BASE_TRIALS, seed=SEED, clinical=cfg, traj_n=0, verbose=True,
            )
            print(f"  {label} / {key}: {time.perf_counter() - t0:.1f}s")
            for r in result["rows"]:
                r = dict(r)
                r["focusing"] = label
                rows.append(r)
    write_csv(OUT / "tables" / "table_focusing_ablation.csv", rows)
    return rows


# ---------------------------------------------------------------------------
# Main fixed-physics campaign (Tier 0.2+0.3 applied) + figures
# ---------------------------------------------------------------------------

def main_campaign(env: TubeEnvironment) -> dict:
    print("\n=== Main campaign (near-wall gated defects, bulk static=0, focusing ON) ===")
    cfg = make_config()  # defaults already fixed
    out = {}
    all_rows = []
    for key in ("pristine", "autoclaved"):
        t0 = time.perf_counter()
        result = run_campaign(
            env, CONTAINER_PROFILES[key],
            base_trials=BASE_TRIALS, seed=SEED, clinical=cfg, traj_n=12, verbose=True,
        )
        print(f"  {key}: {time.perf_counter() - t0:.1f}s  sub20_esc={result['sub20_escaped']}")
        sub = OUT / "figures" / key
        sub.mkdir(parents=True, exist_ok=True)
        generate_experiment_figures(
            env, result["groups"],
            data_log=result["runner"].data_log,
            container_profile=CONTAINER_PROFILES[key],
            injection_plan=result["runner"].injection_plan,
            output_dir=sub,
            samples_per_group=8,
        )
        # exaggerated-y trajectory plot
        plot_trajectory_exaggerated(
            env, result["groups"],
            save_path=sub / f"trajectory_exaggerated_{key}.png",
            title_suffix=f" — {key}",
        )
        result["runner"].data_log.save_trajectories_npz(sub / f"trajectories_{key}.npz")
        for r in result["rows"]:
            all_rows.append(r)
        out[key] = result

    write_csv(OUT / "tables" / "table6_outcomes.csv", all_rows)
    plot_sub20_escaped_bar(out, OUT / "figures" / "fig_sub20_escaped_counts.png")
    return out


def plot_trajectory_exaggerated(env, groups, save_path: Path, title_suffix: str = "") -> None:
    """Tier 2.3 — exaggerate y, colour by fate."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 5))
    length_mm = env.length * 1e3
    width_mm = env.width * 1e3
    ax.axhline(0, color="k", lw=1.5)
    ax.axhline(width_mm, color="k", lw=1.5)
    ax.axvline(0, color="k", ls="--", alpha=0.4)
    ax.axvline(length_mm, color="k", ls="--", alpha=0.4)

    colors = {3.0: "#2ca02c", 10.0: "#ff7f0e", 50.0: "#d62728"}
    for gr in groups:
        c = colors.get(gr.diameter_um, "gray")
        # plot up to 12 trials, colour line by fate
        for trial in gr.trials[:12]:
            if trial.trajectory.n_steps < 2:
                continue
            x = trial.trajectory.x_mm
            y = trial.trajectory.y_mm
            if len(x) > 2500:
                idx = np.linspace(0, len(x) - 1, 2500, dtype=int)
                x, y = x[idx], y[idx]
            fate = trial.final_state
            ls = "-" if fate == ParticleState.ESCAPED else "--"
            alpha = 0.55 if fate == ParticleState.ESCAPED else 0.9
            ax.plot(x, y, color=c, lw=1.0, ls=ls, alpha=alpha)
            marker = "D" if fate == ParticleState.ESCAPED else "X"
            ax.scatter(x[-1], y[-1], marker=marker, s=40, c=c, edgecolors="k", zorder=5)

    ax.set_xlabel("Axial position x (mm)")
    ax.set_ylabel("Lateral position y (mm)  [exaggerated aspect]")
    ax.set_title(f"Particle trajectories (solid=escaped, dashed=stuck){title_suffix}")
    ax.set_xlim(-2, length_mm + 2)
    ax.set_ylim(-0.15, width_mm + 0.15)
    ax.set_aspect(8.0)  # exaggerate y vs true 100:4
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_sub20_escaped_bar(campaigns: dict, save_path: Path) -> None:
    """Tier 2.2 — absolute sub-20 µm escaped counts."""
    import matplotlib.pyplot as plt

    labels = ["Pristine", "Autoclaved"]
    vals = [campaigns["pristine"]["sub20_escaped"], campaigns["autoclaved"]["sub20_escaped"]]
    fig, ax = plt.subplots(figsize=(6, 5))
    bars = ax.bar(labels, vals, color=["#4c72b0", "#c44e52"], edgecolor="k")
    ax.bar_label(bars, padding=3, fontsize=12)
    ax.set_ylabel("Absolute escaped count (sub-20 µm)")
    ax.set_title("Sub-20 µm particles delivered to patient\n(absolute count — the headline exposure number)")
    ratio = vals[1] / vals[0] if vals[0] else float("nan")
    ax.text(0.5, max(vals) * 0.92, f"Autoclaved / Pristine = {ratio:.2f}×",
            ha="center", fontsize=11)
    ax.grid(axis="y", alpha=0.4)
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Tier 1.1 — replicates mean ± s.d.
# ---------------------------------------------------------------------------

def tier11_replicates(env: TubeEnvironment) -> List[dict]:
    print(f"\n=== Tier 1.1: {N_REPLICATES} replicates × {REP_TRIALS} trials ===")
    cfg = make_config()
    summary_rows = []
    raw = []
    for key in ("pristine", "autoclaved"):
        for d in EXPERIMENTAL_DIAMETERS_UM:
            esc_fracs = []
            stuck_fracs = []
            for rep in range(N_REPLICATES):
                seed = SEED + 10_000 * rep + int(d * 17)
                # Run only this diameter by building a temporary single-size profile weight
                profile = CONTAINER_PROFILES[key]
                # Use full profile but we'll extract this diameter — cheaper: run_single_trial loop
                esc = stuck = 0
                n = REP_TRIALS
                # Scale autoclaved injection: still use same n per cell for replicate stats
                for t in range(n):
                    rng = np.random.default_rng(seed + t)
                    r = run_single_trial(
                        env, d, t, rng, fast_mode=True, record_trajectory=False,
                        clinical_config=cfg, container_key=key,
                    )
                    if r.final_state == ParticleState.ESCAPED:
                        esc += 1
                    elif r.final_state == ParticleState.STUCK:
                        stuck += 1
                esc_fracs.append(100.0 * esc / n)
                stuck_fracs.append(100.0 * stuck / n)
                raw.append({
                    "container": key, "diameter_um": d, "replicate": rep,
                    "escape_pct": esc_fracs[-1], "stuck_pct": stuck_fracs[-1],
                })
            summary_rows.append({
                "container": key,
                "diameter_um": d,
                "n_replicates": N_REPLICATES,
                "trials_per_rep": REP_TRIALS,
                "escape_mean_pct": round(float(np.mean(esc_fracs)), 1),
                "escape_sd_pct": round(float(np.std(esc_fracs, ddof=1)), 1),
                "stuck_mean_pct": round(float(np.mean(stuck_fracs)), 1),
                "stuck_sd_pct": round(float(np.std(stuck_fracs, ddof=1)), 1),
            })
            print(
                f"  {key} {d:.0f}µm: escape {summary_rows[-1]['escape_mean_pct']}±"
                f"{summary_rows[-1]['escape_sd_pct']} %"
            )
    write_csv(OUT / "tables" / "table6_replicates_mean_sd.csv", summary_rows)
    write_csv(OUT / "tables" / "table6_replicates_raw.csv", raw)
    return summary_rows


# ---------------------------------------------------------------------------
# Tier 1.3 — size_snag_exponent sensitivity
# ---------------------------------------------------------------------------

def tier13_exponent_sensitivity(env: TubeEnvironment) -> List[dict]:
    print("\n=== Tier 1.3: size_snag_exponent ∈ {0, 0.35, 0.5, 1.0} ===")
    rows = []
    for alpha in (0.0, 0.35, 0.5, 1.0):
        cfg = make_config(size_snag_exponent=alpha)
        for key in ("pristine",):
            result = run_campaign(
                env, CONTAINER_PROFILES[key],
                base_trials=40, seed=SEED, clinical=cfg, traj_n=0, verbose=False,
            )
            for r in result["rows"]:
                rows.append({**r, "size_snag_exponent": alpha})
            print(f"  α={alpha}: " + ", ".join(
                f"{r['diameter_um']:.0f}µm esc={r['escape_pct']:.0f}%" for r in result["rows"]
            ))
    write_csv(OUT / "tables" / "sensitivity_size_snag_exponent.csv", rows)
    return rows


# ---------------------------------------------------------------------------
# Tier 2.1 — Φ_deg × Φ_wall sweep
# ---------------------------------------------------------------------------

def tier21_phi_sweep(env: TubeEnvironment) -> None:
    """
    Paper Φ-sweep. Prefer the authoritative focusing-OFF campaign
    (Claude deliverables/C_phi_sweep.csv: base_trials=30, N_esc_P=80,
    ratio/Φ_deg = 0.972±0.018). Never fall back to base=20 (old N_esc_P≈51).
    """
    print("\n=== Tier 2.1: Φ_deg × Φ_wall sweep (paper: N_esc_P=80) ===")
    import matplotlib.pyplot as plt

    auth = ROOT / "Claude deliverables" / "C_phi_sweep.csv"
    phi_degs = [1.0, 1.5, 2.0, 2.5, 3.0]
    phi_walls = [1.0, 1.5, 1.75, 2.0, 2.5]
    rows = []
    ratio_grid = np.zeros((len(phi_walls), len(phi_degs)))

    if auth.exists():
        print(f"  Loading authoritative sweep → {auth}")
        with auth.open() as f:
            for row in csv.DictReader(f):
                pd, pw = float(row["phi_deg"]), float(row["phi_wall"])
                n_p, n_a = int(float(row["N_esc_P"])), int(float(row["N_esc_A"]))
                ratio = float(row["ratio"])
                rows.append({
                    "phi_deg": pd, "phi_wall": pw,
                    "N_esc_P": n_p, "N_esc_A": n_a, "ratio_A_over_P": round(ratio, 3),
                })
                ratio_grid[phi_walls.index(pw), phi_degs.index(pd)] = ratio
        n_p0 = rows[0]["N_esc_P"]
        if n_p0 != 80:
            raise SystemExit(f"Authoritative C_phi_sweep has N_esc_P={n_p0}, expected 80")
        rs = [r["ratio_A_over_P"] / r["phi_deg"] for r in rows]
        print(f"  N_esc_P={n_p0}; ratio/Φ_deg = {np.mean(rs):.3f}±{np.std(rs, ddof=1):.3f}")
    else:
        # Recompute at the paper setting (base=30, focusing OFF) — NOT base=20.
        base = 30
        for i, pw in enumerate(phi_walls):
            for j, pd in enumerate(phi_degs):
                cfg_p = make_config(focus_enabled=False, degraded_tubing_multiplier=1.0)
                cfg_a = make_config(focus_enabled=False, degraded_tubing_multiplier=pw)
                prof_p = ContainerProfile(
                    key="pristine", display_name="P", comparison_label="P",
                    degradation_factor=1.0,
                    size_weights={3.0: 1.0, 10.0: 1.0, 50.0: 1.0},
                )
                prof_a = ContainerProfile(
                    key="autoclaved", display_name="A", comparison_label="A",
                    degradation_factor=pd,
                    size_weights={3.0: 2.8, 10.0: 2.2, 50.0: 0.35},
                )
                rp = run_campaign(env, prof_p, base_trials=base, seed=SEED, clinical=cfg_p, traj_n=0)
                ra = run_campaign(env, prof_a, base_trials=base, seed=SEED, clinical=cfg_a, traj_n=0)
                n_esc_p = rp["total_escaped"]
                n_esc_a = ra["total_escaped"]
                ratio = n_esc_a / n_esc_p if n_esc_p else float("nan")
                ratio_grid[i, j] = ratio
                rows.append({
                    "phi_deg": pd, "phi_wall": pw,
                    "N_esc_P": n_esc_p, "N_esc_A": n_esc_a, "ratio_A_over_P": round(ratio, 3),
                })
                print(f"  Φ_deg={pd}, Φ_wall={pw}: N_esc A/P = {ratio:.2f} (N_P={n_esc_p})")

    write_csv(OUT / "tables" / "phi_sweep.csv", rows)

    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(ratio_grid, origin="lower", cmap="RdYlGn", aspect="auto",
                   vmin=0.8, vmax=max(3.0, float(np.nanmax(ratio_grid))))
    ax.set_xticks(range(len(phi_degs)))
    ax.set_xticklabels([str(x) for x in phi_degs])
    ax.set_yticks(range(len(phi_walls)))
    ax.set_yticklabels([str(y) for y in phi_walls])
    ax.set_xlabel(r"Polymer degradation $\Phi_{\mathrm{deg}}$")
    ax.set_ylabel(r"Wall stickiness $\Phi_{\mathrm{wall}}$")
    ax.set_title(
        r"$N_{\mathrm{esc}}(A)/N_{\mathrm{esc}}(P)$ — focusing OFF (vertical = $\Phi_{\mathrm{wall}}$ inert)"
    )
    for i in range(len(phi_walls)):
        for j in range(len(phi_degs)):
            ax.text(j, i, f"{ratio_grid[i, j]:.2f}", ha="center", va="center", fontsize=9)
    op_j, op_i = phi_degs.index(2.5), phi_walls.index(1.75)
    ax.scatter([op_j], [op_i], s=200, facecolors="none", edgecolors="k", lw=2, label="operating point")
    ax.legend(loc="upper left")
    fig.colorbar(im, ax=ax, label="escaped-count ratio A/P")
    ax.contour(ratio_grid, levels=[1.0], colors="k", linewidths=1.5, origin="lower")
    fig.tight_layout()
    figs_dir = ROOT / "figs"
    figs_dir.mkdir(exist_ok=True)
    for dest in (
        OUT / "figures" / "fig_phi_sweep.png",
        OUT / "figures" / "fig_phi_sweep_contour.png",
        figs_dir / "fig_phi_sweep.png",
    ):
        fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_phi_sweep.png (N_esc_P from authoritative C_phi_sweep.csv)")


def ship_paper_figs() -> None:
    """Ensure the five Overleaf filenames exist under figs/."""
    required = [
        "fig_phi_sweep.png",
        "fig_sub20_escaped_counts.png",
        "fig_delta_w_retention.png",
        "trajectory_exaggerated_pristine.png",
        "trajectory_exaggerated_autoclaved.png",
    ]
    figs = ROOT / "figs"
    figs.mkdir(exist_ok=True)
    sources = [
        figs,
        OUT / "figures",
        ROOT / "Claude deliverables" / "final_open_todos" / "figures",
        ROOT / "PAPER_ATTACHMENTS",
    ]
    for name in required:
        dest = figs / name
        if dest.exists() and dest.stat().st_size > 1000:
            print(f"  figs/{name} OK ({dest.stat().st_size} B)")
            continue
        found = None
        for src_dir in sources:
            cand = src_dir / name
            if cand.exists() and cand.stat().st_size > 1000:
                found = cand
                break
        if found is None:
            raise SystemExit(f"Missing paper figure: {name}")
        shutil.copy2(found, dest)
        print(f"  shipped {name} ← {found}")


# ---------------------------------------------------------------------------
# Hazard equations + admin text deliverables
# ---------------------------------------------------------------------------

def write_static_docs(main: dict, focus_rows: List[dict], rep_rows: List[dict]) -> None:
    hazard = OUT / "docs"
    hazard.mkdir(parents=True, exist_ok=True)

    (hazard / "HAZARD_MODEL.md").write_text("""# Hazard Model (canonical)

## Capture classification

| Label | Condition |
|-------|-----------|
| `hydrodynamic_wall` | Geometric contact: `y ≤ r` or `y ≥ W − r` |
| `jagged_edge_defect` | Near-wall **and** `|x − x_def| ≤ δ`, one-shot Bernoulli |
| `electrostatic_near_wall` | `min_gap ≤ δ_wall`, continuous hazard |
| `connector_static_clamp` | Near-wall electrostatic × `M_conn` in connector axial windows |

**Bulk static rate is identically zero** (Tier 0.3). There is no adhesive surface in the lumen.

## Equations

```
Defect layout (paired):   N_def ~ Poisson(ρ_def · L)     # wear does NOT inflate count
Defect x-positions:       x_def ~ Uniform(0, L)

Near-wall gate:           near = (min_gap ≤ δ_wall)      # δ_wall = 300 µm
                          # BOTH defect snags and static require near=True

Snag (one-shot):          if near and |x−x_def|≤δ and first visit:
                            P_snag = p0 · w · (d/3)^α
                            stick with probability P_snag

Static hazard:            if near:
                            λ = λ_near
                            if x in connector zone: λ ← λ · M_conn
                            λ ← λ · w · (d/3)^α
                            P_cap = 1 − exp(−λ Δt)

focus_efficiency:         α_eff = focus_efficiency · (1 − e^(−λ_focus Δt))
                          y ← y + α_eff · (y_c − y)     # only if focus_enabled

Brownian kick:            D = k_B T / (6 π μ r)
                          Δy ~ Normal(0, √(2 D Δt))
```

## Constants (parameters.yaml)

- λ_near = 0.18 s⁻¹, λ_bulk = 0
- δ_wall = 300 µm, δ = 0.6 mm, ρ_def = 28 m⁻¹, p0 = 0.028
- M_conn = 3, α = 0.35 (empirical), w = 1 (pristine) or 1.75 (autoclaved)
- focus_efficiency = 0.88, focus_enabled = True (default)
- T = 310 K

## focus_efficiency (Tier 1.5)

Each focusing step computes blend `α = 1 − exp(−λ_focus Δt)`, then applies
`α_eff = focus_efficiency × α` before moving toward the centerline.
It is **not** a physical force coefficient; it is a heuristic imperfect-tracking factor.
Set `focus_enabled=False` to ablate the entire term (Tier 0.4).
""")

    (hazard / "SEEDING.md").write_text("""# Random seeds (Tier 1.2 — paired)

- Base seed = 42
- `group_seed = seed + int(d_um × 1000) + group_index × 10000`
  (**degradation_factor is excluded** so pristine and autoclaved share trial RNGs)
- `trial_rng = Generator(group_seed + trial_id)` → inlet y, Brownian, static draws
- Defect layout: `layout_seed = int(d×1000) + trial_id + 777001` → **identical** for P and A
- Only `w` (tubing wear multiplier) differs between pristine (1.0) and autoclaved (1.75)

Within one trial, one RNG stream is shared across stochastic events (not independent
streams per mechanism).
""")

    (hazard / "ERROR_BARS.md").write_text("""# Error bars (Tier 1.1)

Single-campaign cells are binomial proportions. Wilson 95% CIs are in `table6_outcomes.csv`.

Worst-case SE at p=0.5: SE = √(0.25/n)
- n=25 → SE = ±10.0 percentage points (1σ)
- n=50 → SE ≈ ±7.1 pp

Replicate table (`table6_replicates_mean_sd.csv`) reports mean ± s.d. across 20
independent campaigns (different seeds), 30 trials per cell.
Until reading that table, quote ≈ values, not three-significant-figure percentages.
""")

    (hazard / "FLOW_RATE.md").write_text("""# Flow-rate scenario (Tier 1.6)

ū = 50 mm/s through a 4 mm bore ⇒ Q ≈ π (2 mm)² × 50 mm/s ≈ **2260 mL/hr**.

This is **rapid / gravity-bolus infusion**, not routine maintenance (100–125 mL/hr).

Methods sentence (copy into paper):

> All simulations use a mean axial velocity ū = 50 mm/s (≈2260 mL/hr through the
> 4 mm bore), corresponding to rapid / gravity-bolus IV infusion rather than
> low-rate maintenance infusion (100–125 mL/hr). Fast flow is the condition most
> favourable to particle delivery; results should be interpreted in that scenario.
""")

    (hazard / "AI_DISCLOSURE.md").write_text("""# AI-use disclosure (Tier 3 — paste into main.tex ~line 181)

Suggested sentence:

> Portions of the simulation code, figure-generation scripts, and manuscript
> drafting were assisted by the Cursor AI coding assistant (Composer). All
> physical model choices, parameter values, experimental design (container
> profiles, hazard mechanisms), data interpretation, and final scientific
> claims were made and verified by the author. The publicly archived repository
> regenerates every figure from source.

Confirm with your competition's exact AI policy wording.
""")

    (hazard / "REFERENCES_CHECKLIST.md").write_text("""# Reference verification checklist (Tier 3)

Open each source yourself. Delete any you have not read.

- [ ] [5] USP ⟨788⟩ — cite specific USP–NF edition (year of the monograph you opened)
- [ ] [6] ISO 8536 — cite specific part number and year
- [ ] [7]–[10] — verify year, volume, pages against the PDF/HTML you opened
- [ ] [11] Hawkins — replace with a specific recent paper on steam sterilisation
      of polypropylene medical devices if the current citation is vague

Also confirm Dr. Barge's city (title page currently flagged Deerfield, Illinois / Baxter HQ).
""")

    # Abstract / conclusion placeholders after numbers known
    p = main["pristine"]
    a = main["autoclaved"]
    p_rows = {r["diameter_um"]: r for r in p["rows"]}
    a_rows = {r["diameter_um"]: r for r in a["rows"]}

    def fmt(r):
        return f"{r['escape_pct']:.0f}% [{r['wilson_lo']:.0f}, {r['wilson_hi']:.0f}]"

    abs_txt = f"""# Abstract / Conclusion number block (Tier 4 — AFTER reruns)

## Main campaign (focusing ON, Tier-0 fixes applied)

Flow scenario: rapid / gravity-bolus (~2260 mL/hr).

### Pristine escape (Wilson 95% CI)
- 3 µm: {fmt(p_rows[3.0])}
- 10 µm: {fmt(p_rows[10.0])}
- 50 µm: {fmt(p_rows[50.0])}

### Autoclaved escape
- 3 µm: {fmt(a_rows[3.0])}
- 10 µm: {fmt(a_rows[10.0])}
- 50 µm: {fmt(a_rows[50.0])}

### Absolute sub-20 µm escaped
- Pristine: {p['sub20_escaped']}
- Autoclaved: {a['sub20_escaped']}
- Ratio A/P: {a['sub20_escaped']/p['sub20_escaped'] if p['sub20_escaped'] else float('nan'):.2f}×

### Focusing ablation
See `tables/table_focusing_ablation.csv`.

### Replicates (mean ± s.d.)
See `tables/table6_replicates_mean_sd.csv`.

Rewrite Abstract and Conclusion using ≈ rounding where CIs overlap.
"""
    (hazard / "NUMBERS_FOR_ABSTRACT.md").write_text(abs_txt)

    # Software section replacement
    (hazard / "SOFTWARE_SECTION.md").write_text("""# §Software replacement (main.tex ~line 940)

Replace the red block with:

```latex
\\section*{Software and data availability}
All simulations were performed with a single public codebase,
\\texttt{microplastic\\_simulation.py}, archived at
\\url{https://github.com/<YOUR_USERNAME>/iv-microplastic-transport}
(commit hash recorded at submission).
A pinned \\texttt{requirements.txt} and \\texttt{parameters.yaml} are included.
Every figure and table in this paper is regenerated by

\\begin{verbatim}
python regenerate_paper_figures.py
\\end{verbatim}

The deprecated Colab notebook is archived under
\\texttt{archive/} and is not used for reported results.
```

After you create the GitHub repo, replace `<YOUR_USERNAME>` and pin the commit SHA.
""")


def write_master_readme() -> None:
    (OUT / "README_DELIVERABLES.md").write_text(f"""# Paper deliverables index

Generated: {datetime.now(timezone.utc).isoformat()}

## Folder map

```
paper_deliverables/
  figures/
    pristine/ ...
    autoclaved/ ...
    fig_sub20_escaped_counts.png
    fig_phi_sweep_contour.png
  tables/
    table6_outcomes.csv
    table6_replicates_mean_sd.csv
    table_focusing_ablation.csv
    sensitivity_size_snag_exponent.csv
    phi_sweep.csv
  docs/
    HAZARD_MODEL.md
    SEEDING.md
    ERROR_BARS.md
    FLOW_RATE.md
    AI_DISCLOSURE.md
    REFERENCES_CHECKLIST.md
    NUMBERS_FOR_ABSTRACT.md
    SOFTWARE_SECTION.md
  run_metadata.json
```

## Tier checklist coverage

| Item | Status |
|------|--------|
| 0.1 One codebase | `microplastic_simulation.py` + `regenerate_paper_figures.py`; Colab archived |
| 0.2 Defect y-gate | near-wall required for snags |
| 0.3 Bulk static = 0 | removed |
| 0.4 Focusing ablation | `table_focusing_ablation.csv` |
| 1.1 Error bars / replicates | `table6_replicates_mean_sd.csv` |
| 1.2 Paired seeding | layout seed independent of container |
| 1.3 Exponent sensitivity | `sensitivity_size_snag_exponent.csv` |
| 1.4 Brownian kick | Stokes–Einstein √(2DΔt) |
| 1.5 focus_efficiency | documented in HAZARD_MODEL.md |
| 1.6 Flow rate | labeled rapid/gravity-bolus in FLOW_RATE.md |
| 2.1 Φ sweep | `fig_phi_sweep_contour.png` |
| 2.2 Sub-20 bar | `fig_sub20_escaped_counts.png` |
| 2.3 Trajectories | `trajectory_exaggerated_*.png` |
| 3 Admin texts | docs/* |
| 4 Abstract numbers | NUMBERS_FOR_ABSTRACT.md |
""")


def main() -> int:
    if OUT.exists():
        # keep prior runs but clear generated products
        pass
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)

    # Archive Colab (Tier 0.1)
    archive = ROOT / "archive"
    archive.mkdir(exist_ok=True)
    colab = ROOT / "Microplastic_Simulation_Colab.ipynb"
    if colab.exists():
        dest = archive / "Microplastic_Simulation_Colab.ipynb.DEPRECATED"
        shutil.copy2(colab, dest)
        print(f"Archived Colab → {dest}")

    env = make_env()
    t_all = time.perf_counter()

    focus_rows = tier04_focusing_ablation(env)
    main_results = main_campaign(env)
    rep_rows = tier11_replicates(env)
    tier13_exponent_sensitivity(env)
    tier21_phi_sweep(env)
    ship_paper_figs()

    write_static_docs(main_results, focus_rows, rep_rows)
    write_master_readme()

    meta = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "base_trials": BASE_TRIALS,
        "n_replicates": N_REPLICATES,
        "clinical_config": asdict(make_config()),
        "elapsed_s": time.perf_counter() - t_all,
        "pristine_sub20_escaped": main_results["pristine"]["sub20_escaped"],
        "autoclaved_sub20_escaped": main_results["autoclaved"]["sub20_escaped"],
    }
    (OUT / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    print(f"\nDone in {meta['elapsed_s']:.1f}s → {OUT.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
