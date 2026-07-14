#!/usr/bin/env python3
"""
125 mL/hr slow-flow campaign at full Table-6 fidelity: 20×25.

Exact Stokes → Δt scaled ∝ 1/û in fast_mode (5e-5 → ~9e-4 s), so wall-clock
matches the rapid campaign. Seed: 42 + 10000*rep + int(d*17).
"""
from __future__ import annotations

import csv
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path

import numpy as np

from microplastic_simulation import (
    AUTOCLAVE_WALL_WEAR,
    DEFAULT_CLINICAL_REALISM,
    EXPERIMENTAL_DIAMETERS_UM,
    FluidProperties,
    ParticleState,
    TubeEnvironment,
    fast_mode_dt_cap,
    run_single_trial,
)

ROOT = Path(__file__).parent
OUT = ROOT / "Claude deliverables" / "final_open_todos"
ATTACH = ROOT / "PAPER_ATTACHMENTS"
SEED = 42
N_REP = 20
REP_N = 25
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)
U_SLOW = 125e-6 / 3600.0 / (math.pi * (0.002) ** 2)


def _cell(args: tuple) -> dict:
    d_um, w, state = args
    env = TubeEnvironment(
        length=0.100,
        width=0.004,
        fluid=FluidProperties(average_velocity=U_SLOW),
    )
    clinical = replace(
        DEFAULT_CLINICAL_REALISM,
        focus_enabled=False,
        use_brownian_kick=True,
        static_rate_bulk_per_s=0.0,
        size_snag_exponent=0.35,
        degraded_tubing_multiplier=w,
        near_wall_band_m=0.00030,
    )
    fracs = []
    for rep in range(N_REP):
        seed = SEED + 10_000 * rep + int(d_um * 17)
        esc = 0
        for t in range(REP_N):
            rng = np.random.default_rng(seed + t)
            r = run_single_trial(
                env,
                d_um,
                t,
                rng,
                fast_mode=True,
                record_trajectory=False,
                max_time=120.0,
                max_physics_steps=500_000,
                clinical_config=clinical,
                container_key="pristine",
            )
            if r.final_state == ParticleState.ESCAPED:
                esc += 1
        fracs.append(100.0 * esc / REP_N)
    return {
        "state": state,
        "d_um": float(d_um),
        "escape_pct_mean": round(float(np.mean(fracs)), 1),
        "escape_pct_sd": round(float(np.std(fracs, ddof=1)), 1),
        "n_replicates": N_REP,
        "trials_per_rep": REP_N,
    }


def smoke_dt() -> None:
    env = TubeEnvironment(
        length=0.100, width=0.004, fluid=FluidProperties(average_velocity=U_SLOW)
    )
    cap = fast_mode_dt_cap(env)
    print(f"ū={U_SLOW*1e3:.3f} mm/s  fast_mode Δt cap={cap:.3e} s  (want ~9e-4)", flush=True)
    clinical = replace(
        DEFAULT_CLINICAL_REALISM,
        focus_enabled=False,
        use_brownian_kick=True,
        static_rate_bulk_per_s=0.0,
        degraded_tubing_multiplier=1.0,
    )
    t0 = time.perf_counter()
    rng = np.random.default_rng(0)
    r = run_single_trial(
        env, 10.0, 0, rng, fast_mode=True, record_trajectory=False,
        max_time=120.0, clinical_config=clinical, container_key="pristine",
    )
    print(
        f"smoke 1 trial: {time.perf_counter()-t0:.2f}s  state={r.final_state.name}  "
        f"T={r.transit_time_s:.2f}s",
        flush=True,
    )


def write_outputs(slow_rows: list) -> None:
    rapid = [
        {"flow": "2260_mL_hr", "u_avg_mm_s": 50.0, "state": "pristine", "d_um": 3.0,
         "escape_pct_mean": 87.6, "escape_pct_sd": 7.4, "n_replicates": 20, "trials_per_rep": 25,
         "note": "Table 6"},
        {"flow": "2260_mL_hr", "u_avg_mm_s": 50.0, "state": "pristine", "d_um": 10.0,
         "escape_pct_mean": 87.0, "escape_pct_sd": 5.0, "n_replicates": 20, "trials_per_rep": 25,
         "note": "Table 6"},
        {"flow": "2260_mL_hr", "u_avg_mm_s": 50.0, "state": "pristine", "d_um": 50.0,
         "escape_pct_mean": 77.0, "escape_pct_sd": 9.9, "n_replicates": 20, "trials_per_rep": 25,
         "note": "Table 6"},
        {"flow": "2260_mL_hr", "u_avg_mm_s": 50.0, "state": "autoclaved", "d_um": 3.0,
         "escape_pct_mean": 86.2, "escape_pct_sd": 8.5, "n_replicates": 20, "trials_per_rep": 25,
         "note": "Table 6"},
        {"flow": "2260_mL_hr", "u_avg_mm_s": 50.0, "state": "autoclaved", "d_um": 10.0,
         "escape_pct_mean": 87.0, "escape_pct_sd": 5.0, "n_replicates": 20, "trials_per_rep": 25,
         "note": "Table 6"},
        {"flow": "2260_mL_hr", "u_avg_mm_s": 50.0, "state": "autoclaved", "d_um": 50.0,
         "escape_pct_mean": 75.6, "escape_pct_sd": 10.4, "n_replicates": 20, "trials_per_rep": 25,
         "note": "Table 6"},
    ]
    out_rows = list(rapid)
    for r in slow_rows:
        out_rows.append({
            "flow": "125_mL_hr",
            "u_avg_mm_s": round(U_SLOW * 1e3, 3),
            "state": r["state"],
            "d_um": r["d_um"],
            "escape_pct_mean": r["escape_pct_mean"],
            "escape_pct_sd": r["escape_pct_sd"],
            "n_replicates": r["n_replicates"],
            "trials_per_rep": r["trials_per_rep"],
            "note": "20×25 Table-6 seeds; fast_mode Δt∝1/û (~9e-4 s)",
        })

    for dest_dir in (OUT / "tables", ATTACH, ROOT / "Claude deliverables" / "final_open_todos" / "tables"):
        dest_dir.mkdir(parents=True, exist_ok=True)
        path = dest_dir / "flow_rate_comparison.csv"
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)

    lines = [
        "FLOW RATE COMPARISON — 20×25 AUTHORITATIVE SLOW RUN",
        f"ū_slow = {U_SLOW*1e3:.3f} mm/s (125 mL/hr); Δt cap scaled ∝ 1/û ≈ 9e-4 s",
        "seed = 42 + 10000*rep + int(d*17); focusing OFF",
        "",
        "flow          state       3 µm           10 µm          50 µm",
    ]
    for flow in ("2260_mL_hr", "125_mL_hr"):
        for state in ("pristine", "autoclaved"):
            cells = {
                r["d_um"]: r
                for r in out_rows
                if r["flow"] == flow and r["state"] == state
            }
            lines.append(
                f"{flow:<13} {state:<11} "
                f"{cells[3.0]['escape_pct_mean']:.1f}±{cells[3.0]['escape_pct_sd']:.1f}   "
                f"{cells[10.0]['escape_pct_mean']:.1f}±{cells[10.0]['escape_pct_sd']:.1f}   "
                f"{cells[50.0]['escape_pct_mean']:.1f}±{cells[50.0]['escape_pct_sd']:.1f}"
            )
    text = "\n".join(lines) + "\n"
    (OUT / "docs").mkdir(parents=True, exist_ok=True)
    (OUT / "docs" / "FLOW_RATE_COMPARISON.txt").write_text(text)
    (ATTACH / "FLOW_RATE_COMPARISON.txt").write_text(text)
    print(text, flush=True)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ATTACH.mkdir(parents=True, exist_ok=True)
    smoke_dt()

    jobs = []
    for state, w in (("pristine", 1.0), ("autoclaved", AUTOCLAVE_WALL_WEAR)):
        for d in EXPERIMENTAL_DIAMETERS_UM:
            jobs.append((d, w, state))

    print(f"\n=== 20×25 slow flow: {len(jobs)} cells on {N_WORKERS} workers ===", flush=True)
    t0 = time.perf_counter()
    rows = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_cell, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            done += 1
            print(
                f"  [{done}/{len(jobs)}] {row['state']} {row['d_um']:.0f}µm → "
                f"{row['escape_pct_mean']:.1f}±{row['escape_pct_sd']:.1f}% "
                f"({time.perf_counter()-t0:.0f}s)",
                flush=True,
            )
    rows.sort(key=lambda r: (r["state"], r["d_um"]))
    write_outputs(rows)
    print(f"Done in {time.perf_counter()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
