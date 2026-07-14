# IV Microplastic Transport Simulation

Single public codebase for the paper *Mitigating microplastic delivery from IV tubing*
(Aarush Goel).

## Reproduce every figure

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python regenerate_paper_figures.py
```

Outputs land in `paper_deliverables/` (tables, figures, docs).

Paper-facing plot filenames also live in `figs/`
(`fig_phi_sweep.png`, `fig_sub20_escaped_counts.png`,
`trajectory_exaggerated_pristine.png`, `trajectory_exaggerated_autoclaved.png`).
Authoritative focusing-OFF campaign tables: `Claude deliverables/`.

## Source of truth

| File | Role |
|------|------|
| `microplastic_simulation.py` | Physics + clinical hazard model |
| `parameters.yaml` | All constants (flow scenario, hazard rates) |
| `regenerate_paper_figures.py` | Regenerates every paper figure/table |
| `Claude deliverables/B_replicates.csv` | Table 6 (20×25, focusing OFF) |
| `archive/` | Deprecated Colab (do **not** use for reported numbers) |

## Flow scenario

Mean velocity ū = 50 mm/s through a 4 mm bore ≈ **2260 mL/hr**
(**rapid / gravity-bolus infusion**, not 100–125 mL/hr maintenance).
Routine-rate check: 125 mL/hr → ū ≈ 2.8 mm/s (see Tier-1 deliverables).

## Lateral Brownian kick

`ClinicalRealismConfig.use_brownian_kick=True` applies Stokes–Einstein
`√(2 D Δt)` with `D = k_B T / (6 π μ r)` at T = 310 K. The old fixed
`lateral_kick_std_m = 4e-7` is gone.

## License

Research / competition use. Cite the repository commit SHA with the paper.
