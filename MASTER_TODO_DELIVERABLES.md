# Master TODO — Deliverables Pack

Companion to the paper TODO. Everything below is what to hand to your co-author / put in `main.tex`.

**Regeneration command (after the background run finishes):**

```bash
python regenerate_paper_figures.py
```

All products: `paper_deliverables/`

---

## TIER 0 — DONE IN CODE

### 0.1 One codebase
- Canonical solver: `microplastic_simulation.py`
- Parameters: `parameters.yaml`
- One-command figures: `regenerate_paper_figures.py`
- Colab archived: `archive/Microplastic_Simulation_Colab.ipynb.DEPRECATED`
- Live Colab at repo root is **deprecated** (also in `.gitignore` for future commits)
- Git repo initialized locally (commit `d5ac53d`). `gh` is not installed on this machine — push yourself:

```bash
cd "/Users/aarushgoel/Microplastics Research"
# create empty public repo on GitHub, then:
git remote add origin https://github.com/<YOU>/iv-microplastic-transport.git
git push -u origin main
```

Paste the URL into §Software using `paper_deliverables/docs/SOFTWARE_SECTION.md`.

### 0.2 Defect capture gated on y
Snag only if `min_gap_to_wall ≤ near_wall_band_m` (300 µm). Centerline particles cannot hit wall defects.

### 0.3 Bulk static deleted
`static_rate_bulk_per_s = 0.0`. Remove “Static (Bulk Fluid)” column from Table 7.

### 0.4 Focusing ablation
Side-by-side table: `paper_deliverables/tables/table_focusing_ablation.csv`

**Smoke-test preview (pristine, n=60 / n=40):**

| Condition | 3 µm escape | 10 µm | 50 µm |
|-----------|-------------|-------|-------|
| Focusing ON (default) | 100% | 100% | 77% |
| Focusing OFF | 90% | — | 75% |

Matches the predicted pattern: with both Tier-0 fixes, small-particle size gradient collapses; tubing does not filter. **Say so in the paper.**

---

## TIER 1 — DONE IN CODE

| Item | Deliverable |
|------|-------------|
| 1.1 Replicates mean±s.d. | `tables/table6_replicates_mean_sd.csv` (+ Wilson CIs in `table6_outcomes.csv`) |
| 1.2 Paired seeding | Same defect layout for P/A; wear `w` only. See `docs/SEEDING.md` |
| 1.3 Exponent sensitivity | `tables/sensitivity_size_snag_exponent.csv` — α ∈ {0, 0.35, 0.5, 1} |
| 1.4 Brownian kick | Stokes–Einstein `√(2DΔt)`, T=310 K, size-dependent |
| 1.5 focus_efficiency | Documented: multiplies focusing blend α. See `docs/HAZARD_MODEL.md` |
| 1.6 Flow rate | Labeled **rapid / gravity-bolus ≈2260 mL/hr**. Copy `docs/FLOW_RATE.md` into Methods |

---

## TIER 2 — FIGURES

| Item | File |
|------|------|
| 2.1 Φ sweep contour | `figures/fig_phi_sweep_contour.png` + `tables/phi_sweep.csv` |
| 2.2 Sub-20 µm absolute escaped bar | `figures/fig_sub20_escaped_counts.png` |
| 2.3 Exaggerated trajectories | `figures/*/trajectory_exaggerated_*.png` |

---

## TIER 3 — PASTE-IN TEXTS (you must confirm)

| Item | File | Your action |
|------|------|-------------|
| AI disclosure | `docs/AI_DISCLOSURE.md` | Paste ~line 181; match competition rules |
| Dr. Barge city | — | **Ask him** (Deerfield IL flagged) |
| References [5][6][11] | `docs/REFERENCES_CHECKLIST.md` | Open each PDF; fix citations |

---

## TIER 4 — AFTER THE RUN FINISHES

1. Open `docs/NUMBERS_FOR_ABSTRACT.md`
2. Rewrite Abstract + Conclusion with new ≈ numbers (do not keep 94% / 2.26×)
3. Proofread, compile, ship

---

## Hazard equations for the empty § (reviewer TODO)

```
Near-wall gate:     near ⇔ min_gap ≤ δ_wall = 300 µm

Defect layout:      N_def ~ Poisson(ρ_def L)          # paired: wear does not change count
                    x_def ~ Uniform(0, L)

Defect snag:        if near ∧ |x−x_def|≤δ ∧ first visit:
                      P = p0 · w · (d/3)^α
                      (α=0.35 empirical; sensitivity shown)

Static:             if near:
                      λ = λ_near · w · (d/3)^α
                      · M_conn if in connector zone
                      P = 1 − exp(−λ Δt)
                    λ_bulk ≡ 0

Hydrodynamic:       stick if y≤r or y≥W−r

focus_efficiency:   α_eff = 0.88 · (1 − e^{−λ_f Δt})
Brownian:           D=kT/(6πμr),  Δy~N(0,√(2DΔt))
```

Full write-up: `docs/HAZARD_MODEL.md`

---

## Error bars (one line for the paper)

Single campaign: binomial Wilson 95% CI.  
n=25 ⇒ SE ≤ ±10 percentage points at p=0.5.  
Prefer `table6_replicates_mean_sd.csv` (mean ± s.d., 20 seeds).

---

## What you still do manually

1. `git push` to public GitHub; put URL in §Software  
2. Ask Dr. Barge his city  
3. Verify refs [5]–[11] yourself  
4. Rewrite Abstract after `NUMBERS_FOR_ABSTRACT.md` exists from the finished run  
5. Confirm competition AI disclosure wording  
