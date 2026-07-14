# FINAL OPEN TODOS — USER PACKAGE

Open this file first. Paths relative to workspace root.

---

## 1.1 Corrected α-sweep — DONE

VERIFY: α = 0.35 pristine 10 µm = **87.0 ± 5.0** ≡ Table 6 — **PASS**.

| α | pristine escape % (3 / 10 / 50 µm) | autoclaved escape % (3 / 10 / 50 µm) |
|---|-----------------------------------|-------------------------------------|
| 0.00 | 87.6±7.4 / 87.6±5.0 / 79.4±9.8 | 86.2±8.5 / 87.0±5.0 / 78.2±9.7 |
| **0.35** | **87.6±7.4 / 87.0±5.0 / 77.0±9.9** | **86.2±8.5 / 87.0±5.0 / 75.6±10.4** |
| 0.50 | 87.6±7.4 / 87.0±5.0 / 76.2±10.7 | 86.2±8.5 / 86.6±5.4 / 74.2±10.4 |
| 1.00 | 87.6±7.4 / 86.6±5.4 / 73.4±11.0 | 86.2±8.5 / 86.6±5.4 / 72.6±10.9 |

Files: `Claude deliverables/final_open_todos/docs/ALPHA_SWEEP_CORRECTED.txt`,
`Claude deliverables/E_sensitivity_alpha.csv`

---

## 1.2 δ_w sweep — DONE

Retention vs `2δ_w/W` (⟨P_cap⟩=1). Pristine, focusing OFF, 20×25.

| δ_w (µm) | pred % | ret 3 µm | ret 10 µm | ret 50 µm |
|----------|--------|----------|-----------|-----------|
| 100 | 5.0 | 3.2 | 4.6 | 18.8 |
| 200 | 10.0 | 7.6 | 9.6 | 20.8 |
| 300 | 15.0 | 12.4 | 13.0 | 23.0 |
| 500 | 25.0 | 17.8 | 21.6 | 30.6 |
| 1000 | 50.0 | 30.4 | 34.8 | 47.6 |

**Verdict:** For 3–10 µm, measured retention **tracks** `2δ_w/W` (monotone / near-linear), at ~0.6–0.85× perfect-capture — geometry sets the scale; ⟨P_cap⟩ high but not exactly 1. 50 µm sits above the small-size line (buoyancy into the top wall). Plot: `figs/fig_delta_w_retention.png`.

---

## 1.3 Rapid vs 125 mL/hr — DONE (authoritative 20×25)

Δt scaled ∝ 1/û in `fast_mode` (5×10⁻⁵ → ≈9×10⁻⁴ s). Wall-clock ≈ rapid campaign (~11 min).

| flow | state | 3 µm | 10 µm | 50 µm | n |
|------|-------|------|-------|-------|---|
| 2260 mL/hr (ū=50) | pristine | 87.6±7.4 | 87.0±5.0 | 77.0±9.9 | **20×25** Table 6 |
| 2260 mL/hr | autoclaved | 86.2±8.5 | 87.0±5.0 | 75.6±10.4 | **20×25** |
| **125 mL/hr** (ū=2.763) | pristine | **84.4±8.5** | **80.6±7.0** | **0.0±0.0** | **20×25** |
| 125 mL/hr | autoclaved | 84.4±8.5 | 80.6±7.0 | 0.0±0.0 | **20×25** |

- Seed: `42+10000*rep+int(d*17)`. Supersedes the emergency 1×25 (76/68).
- pristine 10 µm: **87.0 → 80.6 (−6.4 pp)**; 3 µm: **87.6 → 84.4 (−3.2 pp)**. Longer residence → higher retention; conclusion survives.
- 50 µm: **0% escape** at 125 mL/hr (Λ-driven top-wall sweep).
- P and A escape % match at slow flow (paired RNG; saturation).

Files: `Claude deliverables/final_open_todos/tables/flow_rate_comparison.csv`,
`PAPER_ATTACHMENTS/flow_rate_comparison.csv`

---

## 1.4 Figures — DONE

All under `figs/`:

- `figs/fig_phi_sweep.png` — vertical contours (Φ_wall inert)
- `figs/fig_sub20_escaped_counts.png` — **87.3 → 302.9**
- `figs/trajectory_exaggerated_pristine.png`
- `figs/trajectory_exaggerated_autoclaved.png` (distinct; no centreline snap)
- `figs/fig_delta_w_retention.png`

---

## 1.5 Lateral kick — DONE

**Replaced** with Stokes–Einstein `√(2 D Δt)`, `D = k_B T / (6 π μ r)` (T=310 K); old fixed `lateral_kick_std_m=4e-7` removed.

---

## 2.1 Reproducibility — LOCAL DONE; you push

- `C_run_log.txt` **gone**. Keep `C_phi_sweep.txt` / `.csv` (base_trials=30, N_esc_P=80).
- README OK. `gh` not installed; no git remote.
- Push template + placeholder URL: `Claude deliverables/final_open_todos/docs/SOFTWARE_SECTION.txt`

---

## 3 Admin drafts — READY; human actions remain

| Item | File / action |
|------|----------------|
| AI disclosure | Paste `docs/AI_DISCLOSURE.txt` ~line 192 |
| Dr. Barge city | **Ask him** — Deerfield IL flagged only |
| Refs [5]–[11] | Open PDFs; `docs/REFERENCES_CHECKLIST.txt` |
| Red TODOs | Delete every `\textcolor{red}{...}` |
| Overleaf | Compile clean and ship |

---

## Still human-only

1. Push public GitHub → paste URL + SHA into §Software  
2. Ask Dr. Barge his city  
3. Verify refs [5]–[11] yourself  
4. Clear red blocks + Overleaf compile  

Do **not** rewrite the abstract for these checks.
