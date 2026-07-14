#!/usr/bin/env python3
"""
Final open TODOs (Tier 1 scientific + figure regen).
Focusing OFF = paper main. Replicate seed policy matches Table 6.
"""
from __future__ import annotations

import csv
import math
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from microplastic_simulation import (
    AUTOCLAVE_WALL_WEAR,
    CONTAINER_PROFILES,
    DEFAULT_CLINICAL_REALISM,
    DEFAULT_MAX_TIME_S,
    EXPERIMENTAL_DIAMETERS_UM,
    ClinicalRealismConfig,
    ExperimentRunner,
    FluidProperties,
    ParticleState,
    TubeEnvironment,
    run_single_trial,
)

ROOT = Path(__file__).parent
OUT = ROOT / "Claude deliverables" / "final_open_todos"
FIGS = ROOT / "figs"
SEED = 42
N_REP = 20
REP_N = 25
W_PRISTINE = 1.0
W_AUTO = AUTOCLAVE_WALL_WEAR
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)


def cfg(**kw) -> ClinicalRealismConfig:
    base = dict(
        focus_enabled=False,
        use_brownian_kick=True,
        static_rate_bulk_per_s=0.0,
        size_snag_exponent=0.35,
        degraded_tubing_multiplier=1.0,
    )
    base.update(kw)
    return replace(DEFAULT_CLINICAL_REALISM, **base)


def _replicate_cell(args: tuple) -> dict:
    """Worker: one (alpha/state/d/band/flow) cell under Table-6 seeds."""
    (
        diameter_um,
        w,
        alpha,
        near_wall_band_m,
        u_avg,
        max_time,
        max_steps,
        state_label,
        extra,
    ) = args
    env = TubeEnvironment(
        length=0.100,
        width=0.004,
        fluid=FluidProperties(average_velocity=u_avg),
    )
    clinical = cfg(
        degraded_tubing_multiplier=w,
        size_snag_exponent=alpha,
        near_wall_band_m=near_wall_band_m,
    )
    fracs = []
    for rep in range(N_REP):
        seed = SEED + 10_000 * rep + int(diameter_um * 17)
        esc = 0
        for t in range(REP_N):
            rng = np.random.default_rng(seed + t)
            r = run_single_trial(
                env,
                diameter_um,
                t,
                rng,
                fast_mode=True,
                record_trajectory=False,
                max_time=max_time,
                max_physics_steps=max_steps,
                clinical_config=clinical,
                container_key="pristine",
            )
            if r.final_state == ParticleState.ESCAPED:
                esc += 1
        fracs.append(100.0 * esc / REP_N)
    mean = float(np.mean(fracs))
    sd = float(np.std(fracs, ddof=1))
    out = {
        "state": state_label,
        "d_um": float(diameter_um),
        "escape_pct_mean": round(mean, 1),
        "escape_pct_sd": round(sd, 1),
        "n_replicates": N_REP,
        "trials_per_rep": REP_N,
        "alpha": float(alpha),
        "near_wall_band_m": float(near_wall_band_m),
        "u_avg": float(u_avg),
    }
    out.update(extra)
    return out


def run_cells(jobs: List[tuple], label: str) -> List[dict]:
    print(f"\n=== {label}: {len(jobs)} cells × {N_REP}×{REP_N} on {N_WORKERS} workers ===", flush=True)
    rows: List[dict] = []
    t0 = time.perf_counter()
    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_replicate_cell, j): j for j in jobs}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            done += 1
            print(
                f"  [{done}/{len(jobs)}] {row.get('tag', '')} "
                f"{row['state']} {row['d_um']:.0f}µm → {row['escape_pct_mean']:.1f}±{row['escape_pct_sd']:.1f}% "
                f"({time.perf_counter()-t0:.0f}s)",
                flush=True,
            )
    return rows


def write_csv(path: Path, rows: List[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def load_table6() -> List[dict]:
    path = ROOT / "Claude deliverables" / "B_replicates.csv"
    rows = []
    with path.open() as f:
        for r in csv.DictReader(f):
            rows.append(
                {
                    "alpha": 0.35,
                    "state": r["state"],
                    "d_um": float(r["d_um"]),
                    "escape_pct_mean": float(r["escape_pct_mean"]),
                    "escape_pct_sd": float(r["escape_pct_sd"]),
                    "n_replicates": int(r["n_replicates"]),
                    "trials_per_rep": int(r["trials_per_rep"]),
                    "source": "B_replicates.csv (identical config)",
                }
            )
    return rows


def verify_table6_seed(env_u: float = 0.05) -> dict:
    """Recompute pristine 10 µm α=0.35 under Table-6 seeds; must be 87.0."""
    print("\n=== 1.1 VERIFY α=0.35 pristine 10 µm (must = 87.0) ===", flush=True)
    job = (
        10.0,
        W_PRISTINE,
        0.35,
        0.00030,
        env_u,
        DEFAULT_MAX_TIME_S,
        500_000,
        "pristine",
        {"tag": "VERIFY"},
    )
    row = _replicate_cell(job)
    ok = abs(row["escape_pct_mean"] - 87.0) < 0.05
    print(
        f"  VERIFY: {row['escape_pct_mean']}±{row['escape_pct_sd']}  "
        f"(Table 6 = 87.0±5.0) → {'PASS' if ok else 'FAIL'}",
        flush=True,
    )
    row["check_pass"] = ok
    return row


def tip_1_1_alpha_sweep() -> List[dict]:
    """α=0.35 from Table 6; rerun α∈{0,0.5,1.0} under same seeds."""
    verified = verify_table6_seed()
    if not verified.get("check_pass"):
        raise SystemExit("α=0.35 seed policy does not reproduce Table 6 — abort")

    base = load_table6()
    jobs = []
    for alpha in (0.0, 0.5, 1.0):
        for state, w in (("pristine", W_PRISTINE), ("autoclaved", W_AUTO)):
            for d in EXPERIMENTAL_DIAMETERS_UM:
                jobs.append(
                    (
                        d,
                        w,
                        alpha,
                        0.00030,
                        0.05,
                        DEFAULT_MAX_TIME_S,
                        500_000,
                        state,
                        {"tag": f"α={alpha}", "source": "rerun"},
                    )
                )
    new_rows = run_cells(jobs, "1.1 α-sweep (new α)")
    # Normalize columns
    out = []
    for r in base:
        out.append(
            {
                "alpha": 0.35,
                "state": r["state"],
                "d_um": r["d_um"],
                "escape_pct_mean": r["escape_pct_mean"],
                "escape_pct_sd": r["escape_pct_sd"],
                "n_replicates": r["n_replicates"],
                "trials_per_rep": r["trials_per_rep"],
                "source": r["source"],
            }
        )
    for r in new_rows:
        out.append(
            {
                "alpha": r["alpha"],
                "state": r["state"],
                "d_um": r["d_um"],
                "escape_pct_mean": r["escape_pct_mean"],
                "escape_pct_sd": r["escape_pct_sd"],
                "n_replicates": r["n_replicates"],
                "trials_per_rep": r["trials_per_rep"],
                "source": "rerun Table-6 seeds",
            }
        )
    out.sort(key=lambda x: (x["alpha"], x["state"], x["d_um"]))
    write_csv(OUT / "tables" / "E_alpha_sweep_table6_seeds.csv", out)

    # Also overwrite Claude deliverables E_sensitivity_alpha.csv with corrected
    write_csv(ROOT / "Claude deliverables" / "E_sensitivity_alpha.csv", out)

    lines = [
        "CORRECTED α-SWEEP (Table-6 replicate seed policy: 20×25)",
        "seed = 42 + 10000*rep + int(d*17); focusing OFF",
        f"VERIFY α=0.35 pristine 10 µm recompute = {verified['escape_pct_mean']}±{verified['escape_pct_sd']} (PASS)",
        "α=0.35 rows = B_replicates.csv (identical configuration)",
        "",
        "alpha   pristine_escape% (3/10/50)           autoclaved_escape% (3/10/50)",
    ]
    for alpha in (0.0, 0.35, 0.5, 1.0):
        p = {r["d_um"]: r for r in out if r["alpha"] == alpha and r["state"] == "pristine"}
        a = {r["d_um"]: r for r in out if r["alpha"] == alpha and r["state"] == "autoclaved"}
        lines.append(
            f"{alpha:<7.2f} {p[3.0]['escape_pct_mean']:.1f}±{p[3.0]['escape_pct_sd']:.1f} / "
            f"{p[10.0]['escape_pct_mean']:.1f}±{p[10.0]['escape_pct_sd']:.1f} / "
            f"{p[50.0]['escape_pct_mean']:.1f}±{p[50.0]['escape_pct_sd']:.1f}"
            f"     "
            f"{a[3.0]['escape_pct_mean']:.1f}±{a[3.0]['escape_pct_sd']:.1f} / "
            f"{a[10.0]['escape_pct_mean']:.1f}±{a[10.0]['escape_pct_sd']:.1f} / "
            f"{a[50.0]['escape_pct_mean']:.1f}±{a[50.0]['escape_pct_sd']:.1f}"
        )
    (OUT / "docs").mkdir(parents=True, exist_ok=True)
    (OUT / "docs" / "ALPHA_SWEEP_CORRECTED.txt").write_text("\n".join(lines) + "\n")
    (ROOT / "Claude deliverables" / "E_sensitivity.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    return out


def tip_1_2_band_sweep() -> List[dict]:
    print("\n=== 1.2 δ_w hazard-band sweep ===", flush=True)
    W = 0.004
    jobs = []
    for dw_um in (100, 200, 300, 500, 1000):
        dw = dw_um * 1e-6
        for d in EXPERIMENTAL_DIAMETERS_UM:
            # δ_w=300µm × α=0.35 × pristine already in Table 6 — skip recompute for speed? Still need all for consistency under same worker; for 300 use table6 when d known
            if dw_um == 300:
                continue
            jobs.append(
                (
                    d,
                    W_PRISTINE,
                    0.35,
                    dw,
                    0.05,
                    DEFAULT_MAX_TIME_S,
                    500_000,
                    "pristine",
                    {"tag": f"δw={dw_um}", "delta_w_um": dw_um},
                )
            )
    new_rows = run_cells(jobs, "1.2 δ_w sweep") if jobs else []

    # Inject Table-6 at δ_w=300
    t6 = load_table6()
    rows = []
    for dw_um in (100, 200, 300, 500, 1000):
        dw = dw_um * 1e-6
        geom = 2.0 * dw / W
        for d in EXPERIMENTAL_DIAMETERS_UM:
            if dw_um == 300:
                r = next(x for x in t6 if x["state"] == "pristine" and x["d_um"] == d)
                me, sd = r["escape_pct_mean"], r["escape_pct_sd"]
            else:
                r = next(
                    x
                    for x in new_rows
                    if x["d_um"] == d and int(x.get("delta_w_um", -1)) == dw_um
                )
                me, sd = r["escape_pct_mean"], r["escape_pct_sd"]
            ret = 100.0 - me
            rows.append(
                {
                    "delta_w_um": dw_um,
                    "geom_fraction_2dw_over_W": round(geom, 4),
                    "d_um": d,
                    "escape_pct_mean": round(me, 1),
                    "escape_pct_sd": round(sd, 1),
                    "retention_pct_mean": round(ret, 1),
                    "prediction_pct_if_Pcap1": round(100.0 * geom, 1),
                }
            )
            print(
                f"  δ_w={dw_um}µm d={d:.0f}: esc={me:.1f}% ret={ret:.1f}% pred={100*geom:.1f}%",
                flush=True,
            )

    write_csv(OUT / "tables" / "delta_w_band_sweep.csv", rows)

    import matplotlib.pyplot as plt

    by_dw = {r["delta_w_um"]: r for r in rows if r["d_um"] == 3.0}
    xs = sorted(by_dw.keys())
    measured = [by_dw[x]["retention_pct_mean"] for x in xs]
    predicted = [by_dw[x]["prediction_pct_if_Pcap1"] for x in xs]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(xs, predicted, "k--", lw=2, label=r"theory $2\delta_w/W$ ($\langle P_{\mathrm{cap}}\rangle=1$)")
    ax.plot(xs, measured, "o-", color="#c44e52", lw=2, ms=8, label="measured retention (3 µm, 20×25)")
    ax.set_xlabel(r"Hazard-band half-width $\delta_w$ (µm)")
    ax.set_ylabel("Retention (%)")
    ax.set_title("Retention tracks geometric band fraction")
    ax.legend()
    ax.grid(True, alpha=0.4)
    fig.tight_layout()
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "figures" / "fig_delta_w_retention.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGS / "fig_delta_w_retention.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return rows


def tip_1_3_slow_flow() -> List[dict]:
    print("\n=== 1.3 Routine infusion 125 mL/hr ===", flush=True)
    u = 125e-6 / 3600.0 / (math.pi * (0.002) ** 2)
    print(f"  ū = {u*1e3:.3f} mm/s", flush=True)
    jobs = []
    for state, w in (("pristine", W_PRISTINE), ("autoclaved", W_AUTO)):
        for d in EXPERIMENTAL_DIAMETERS_UM:
            jobs.append(
                (
                    d,
                    w,
                    0.35,
                    0.00030,
                    u,
                    120.0,
                    2_500_000,
                    state,
                    {"tag": "125mL/hr", "flow": "125_mL_hr"},
                )
            )
    slow = run_cells(jobs, "1.3 slow flow")
    rapid = [
        {
            "flow": "2260_mL_hr",
            "u_avg_mm_s": 50.0,
            "state": "pristine",
            "d_um": 3.0,
            "escape_pct_mean": 87.6,
            "escape_pct_sd": 7.4,
        },
        {
            "flow": "2260_mL_hr",
            "u_avg_mm_s": 50.0,
            "state": "pristine",
            "d_um": 10.0,
            "escape_pct_mean": 87.0,
            "escape_pct_sd": 5.0,
        },
        {
            "flow": "2260_mL_hr",
            "u_avg_mm_s": 50.0,
            "state": "pristine",
            "d_um": 50.0,
            "escape_pct_mean": 77.0,
            "escape_pct_sd": 9.9,
        },
        {
            "flow": "2260_mL_hr",
            "u_avg_mm_s": 50.0,
            "state": "autoclaved",
            "d_um": 3.0,
            "escape_pct_mean": 86.2,
            "escape_pct_sd": 8.5,
        },
        {
            "flow": "2260_mL_hr",
            "u_avg_mm_s": 50.0,
            "state": "autoclaved",
            "d_um": 10.0,
            "escape_pct_mean": 87.0,
            "escape_pct_sd": 5.0,
        },
        {
            "flow": "2260_mL_hr",
            "u_avg_mm_s": 50.0,
            "state": "autoclaved",
            "d_um": 50.0,
            "escape_pct_mean": 75.6,
            "escape_pct_sd": 10.4,
        },
    ]
    out = list(rapid)
    for r in slow:
        out.append(
            {
                "flow": "125_mL_hr",
                "u_avg_mm_s": round(u * 1e3, 3),
                "state": r["state"],
                "d_um": r["d_um"],
                "escape_pct_mean": r["escape_pct_mean"],
                "escape_pct_sd": r["escape_pct_sd"],
            }
        )
    write_csv(OUT / "tables" / "flow_rate_comparison.csv", out)
    return out


def tip_1_4_figures() -> None:
    print("\n=== 1.4 Regenerate figures (focusing OFF) ===", flush=True)
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)

    phi_path = ROOT / "Claude deliverables" / "C_phi_sweep.csv"
    grid: Dict[Tuple[float, float], float] = {}
    phi_degs: List[float] = []
    phi_walls: List[float] = []
    with phi_path.open() as f:
        for row in csv.DictReader(f):
            pd, pw = float(row["phi_deg"]), float(row["phi_wall"])
            grid[(pd, pw)] = float(row["ratio"])
            if pd not in phi_degs:
                phi_degs.append(pd)
            if pw not in phi_walls:
                phi_walls.append(pw)
    phi_degs = sorted(phi_degs)
    phi_walls = sorted(phi_walls)
    Z = np.array([[grid[(pd, pw)] for pd in phi_degs] for pw in phi_walls])
    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(
        Z,
        origin="lower",
        cmap="RdYlGn",
        aspect="auto",
        vmin=0.8,
        vmax=max(3.0, float(Z.max())),
    )
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
            ax.text(j, i, f"{Z[i, j]:.2f}", ha="center", va="center", fontsize=9)
    op_j, op_i = phi_degs.index(2.5), phi_walls.index(1.75)
    ax.scatter([op_j], [op_i], s=200, facecolors="none", edgecolors="k", lw=2)
    ax.contour(Z, levels=[1.0], colors="k", linewidths=1.2, origin="lower")
    fig.colorbar(im, ax=ax, label="escaped-count ratio A/P")
    fig.tight_layout()
    for dest in (FIGS / "fig_phi_sweep.png", OUT / "figures" / "fig_phi_sweep.png"):
        fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig_phi_sweep.png", flush=True)

    p_sub = 0.876 * 50 + 0.870 * 50
    a_sub = 0.862 * 196 + 0.870 * 154
    fig, ax = plt.subplots(figsize=(6, 5))
    bars = ax.bar(
        ["Pristine", "Autoclaved"],
        [p_sub, a_sub],
        color=["#4c72b0", "#c44e52"],
        edgecolor="k",
    )
    ax.bar_label(bars, fmt="%.1f", padding=3, fontsize=12)
    ax.set_ylabel("Absolute escaped count (sub-20 µm)")
    ax.set_title("Sub-20 µm delivered to patient\n(from replicate means — focusing OFF)")
    ax.text(0.5, max(p_sub, a_sub) * 0.88, f"A/P = {a_sub/p_sub:.2f}×", ha="center")
    ax.grid(axis="y", alpha=0.4)
    fig.tight_layout()
    for dest in (
        FIGS / "fig_sub20_escaped_counts.png",
        OUT / "figures" / "fig_sub20_escaped_counts.png",
    ):
        fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig_sub20_escaped_counts.png  ({p_sub:.1f} → {a_sub:.1f})", flush=True)

    env = TubeEnvironment(
        length=0.100,
        width=0.004,
        fluid=FluidProperties(average_velocity=0.05),
    )
    for key, w in (("pristine", W_PRISTINE), ("autoclaved", W_AUTO)):
        clinical = cfg(degraded_tubing_multiplier=w)
        runner = ExperimentRunner(
            environment=env,
            container_profile=CONTAINER_PROFILES[key],
            base_trials_per_size=50,
            seed=SEED,
            fast_mode=True,
            verbose=False,
            clinical_realism=True,
            clinical_config=clinical,
            trajectory_trials_per_group=12,
        )
        groups = runner.run_all()
        plot_traj_exaggerated(
            env,
            groups,
            save_path=FIGS / f"trajectory_exaggerated_{key}.png",
            title=f"Focusing OFF — {key}",
        )
        shutil.copy(
            FIGS / f"trajectory_exaggerated_{key}.png",
            OUT / "figures" / f"trajectory_exaggerated_{key}.png",
        )
        print(f"  wrote trajectory_exaggerated_{key}.png", flush=True)


def plot_traj_exaggerated(env, groups, save_path: Path, title: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 5))
    Lmm, Wmm = env.length * 1e3, env.width * 1e3
    ax.axhline(0, color="k", lw=1.5)
    ax.axhline(Wmm, color="k", lw=1.5)
    colors = {3.0: "#2ca02c", 10.0: "#ff7f0e", 50.0: "#d62728"}
    for gr in groups:
        c = colors.get(gr.diameter_um, "gray")
        for trial in gr.trials[:12]:
            if trial.trajectory.n_steps < 2:
                continue
            x, y = trial.trajectory.x_mm, trial.trajectory.y_mm
            if len(x) > 2500:
                idx = np.linspace(0, len(x) - 1, 2500, dtype=int)
                x, y = x[idx], y[idx]
            escaped = trial.final_state == ParticleState.ESCAPED
            ax.plot(
                x,
                y,
                color=c,
                lw=1.0,
                ls="-" if escaped else "--",
                alpha=0.55 if escaped else 0.9,
            )
            ax.scatter(
                x[-1],
                y[-1],
                marker="D" if escaped else "X",
                s=40,
                c=c,
                edgecolors="k",
                zorder=5,
            )
    ax.set_xlabel("Axial x (mm)")
    ax.set_ylabel("Lateral y (mm) [exaggerated]")
    ax.set_title(title + " (solid=escaped, dashed=stuck)")
    ax.set_xlim(-2, Lmm + 2)
    ax.set_ylim(-0.15, Wmm + 0.15)
    ax.set_aspect(8.0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def write_static_docs() -> None:
    docs = OUT / "docs"
    docs.mkdir(parents=True, exist_ok=True)
    (docs / "LATERAL_KICK.txt").write_text(
        "Lateral kick: REPLACED with Stokes–Einstein sqrt(2*D*dt), "
        "D = kB*T/(6*pi*mu*r), T=310 K, size-dependent "
        "(ClinicalRealismConfig.use_brownian_kick=True). "
        "The old size-independent lateral_kick_std_m=4e-7 was removed.\n"
    )
    (docs / "AI_DISCLOSURE.txt").write_text(
        "AI-use disclosure (paste near main.tex line 192):\n\n"
        "Portions of simulation code, figure-generation scripts, and manuscript "
        "editing were assisted by the Cursor AI coding assistant. All physical "
        "model choices, parameter values, experimental design, data, and scientific "
        "conclusions are the author's and were verified by the author.\n"
    )
    (docs / "REFERENCES_CHECKLIST.txt").write_text(
        "Reference verification (YOU must open each PDF):\n\n"
        "[5] USP ⟨788⟩ — cite the specific USP–NF edition year you opened\n"
        "[6] ISO 8536 — cite specific part number and year\n"
        "[7] Maxey & Riley 1983 — confirm year/volume/pages or cut\n"
        "[8] Crowe et al. 2011 — confirm or cut\n"
        "[9] Happel & Brenner 1983 — confirm or cut\n"
        "[10] Israelachvili — confirm edition/year or cut\n"
        "[11] Hawkins — replace with a specific recent paper on steam\n"
        "     sterilisation of polypropylene medical devices if one exists\n\n"
        "Dr. Barge city (line 87): ASK HIM — Deerfield, Illinois is flagged, not confirmed.\n"
    )
    (docs / "SOFTWARE_SECTION.txt").write_text(
        "§Software (replace red block ~line 968):\n\n"
        "All simulations use a single public repository containing\n"
        "microplastic_simulation.py, parameters.yaml, regenerate_paper_figures.py,\n"
        "requirements.txt, and README.md. Clone and run:\n\n"
        "  pip install -r requirements.txt\n"
        "  python regenerate_paper_figures.py\n\n"
        "URL: https://github.com/<YOUR_USERNAME>/iv-microplastic-transport\n"
        "(push this repo, then paste the live URL and commit SHA).\n\n"
        "Contradictory log C_run_log.txt (base_trials=20) was deleted; keep\n"
        "C_phi_sweep.txt / C_phi_sweep.csv (base_trials=30, N_esc_P=80).\n"
    )


def write_summary(alpha_rows, band_rows, flow_rows) -> None:
    lines = [
        "FINAL OPEN TODOS — OUTPUT SUMMARY",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "1.1 α-sweep: docs/ALPHA_SWEEP_CORRECTED.txt + tables/E_alpha_sweep_table6_seeds.csv",
        "1.2 δ_w sweep: tables/delta_w_band_sweep.csv + figs/fig_delta_w_retention.png",
        "1.3 Slow flow: tables/flow_rate_comparison.csv",
        "1.4 Figures in figs/: fig_phi_sweep.png, fig_sub20_escaped_counts.png,",
        "    trajectory_exaggerated_pristine.png, trajectory_exaggerated_autoclaved.png",
        "1.5 Lateral kick: REPLACED with Stokes–Einstein √(2DΔt)",
        "2.1 C_run_log.txt deleted; push GitHub yourself",
        "3.x Paste AI_DISCLOSURE.txt; ask Dr. Barge; verify REFERENCES_CHECKLIST.txt",
        "",
        "=== QUICK α TABLE ===",
    ]
    for alpha in (0.0, 0.35, 0.5, 1.0):
        p = {r["d_um"]: r for r in alpha_rows if r["alpha"] == alpha and r["state"] == "pristine"}
        if not p:
            continue
        lines.append(
            f"α={alpha}: pristine {p[3.0]['escape_pct_mean']}/{p[10.0]['escape_pct_mean']}/{p[50.0]['escape_pct_mean']}"
        )
    lines.append("")
    lines.append("=== FLOW RATE (pristine 10 µm) ===")
    for r in flow_rows:
        if r["state"] == "pristine" and r["d_um"] == 10.0:
            lines.append(
                f"{r['flow']}: escape {r['escape_pct_mean']}±{r['escape_pct_sd']}%"
            )
    (OUT / "SUMMARY.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


def clean_contradictory_logs() -> None:
    bad = ROOT / "Claude deliverables" / "C_run_log.txt"
    if bad.exists():
        bad.unlink()
        print(
            "Deleted contradictory C_run_log.txt (kept C_phi_sweep.* base_trials=30)",
            flush=True,
        )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    (OUT / "docs").mkdir(exist_ok=True)
    clean_contradictory_logs()
    write_static_docs()

    t0 = time.perf_counter()
    tip_1_4_figures()
    alpha_rows = tip_1_1_alpha_sweep()
    band_rows = tip_1_2_band_sweep()
    flow_rows = tip_1_3_slow_flow()
    write_summary(alpha_rows, band_rows, flow_rows)

    print(f"\nDone in {time.perf_counter()-t0:.1f}s → {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
