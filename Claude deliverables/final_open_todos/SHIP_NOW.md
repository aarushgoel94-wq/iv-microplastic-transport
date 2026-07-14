# SHIP NOW — paste pack for `main.tex`

Open this file. Paste blocks in order. Do **not** invent USP edition numbers. Do **not** guess Dr Barge’s city.

Live public repo confirmed HTTP 200:
`https://github.com/aarushgoel94-wq/iv-microplastic-transport`
Tip SHA (main): `2e44a7ec5afcd79f8a0fb24922fab29f619f617a`
(If `git rev-parse origin/main` later differs by a docs-only tip-SHA refresh, either SHA is fine — code + figs unchanged.)

Five Overleaf PNGs confirmed bit-identical in:
- `~/Downloads/overleaf_project/figs/`
- `Claude deliverables/final_open_todos/figures/`
- repo `figs/`
- Desktop zip `~/Desktop/OVERLEAF_UPLOAD_figs.zip` (Finder revealed)

**Overleaf API upload: NOT possible** — no Overleaf credentials / env tokens on this machine. Drag the Desktop zip into Overleaf `figs/` yourself.

---

## Checklist

### Agent-done

- [x] Public GitHub URL live + tip SHA recorded
- [x] §Software paste block (URL + SHA)
- [x] AI disclosure paste block
- [x] Lateral-kick one-liner paste block
- [x] δ_w = 300 µm origin one-sentence (modelling guess, honest)
- [x] Refs [5][6][11] candidates written (verify-by-opening flagged)
- [x] Dr Barge: ASK HIM — Deerfield IL flagged only (not pasted as fact)
- [x] Flow 20×25 table numbers ready
- [x] Five PNGs confirmed in Downloads/overleaf_project/figs/ + Desktop zip
- [x] `OVERLEAF_UPLOAD_figs.zip` on Desktop refreshed; Finder opened on zip
- [x] Overleaf API checked — no credentials → user must drag-upload
- [x] This pack (`SHIP_NOW.md`) written

### Human-only remaining

- [ ] Paste §Software / AI / lateral kick / δ_w into `main.tex`; delete every `\textcolor{red}{...}`
- [ ] Ask Dr Barge his city (Deerfield, Illinois = Baxter HQ flag only)
- [ ] Open USP ⟨788⟩ / ISO 8536 / [11] candidate PDFs; paste only after you verify
- [ ] Drag `~/Desktop/OVERLEAF_UPLOAD_figs.zip` into Overleaf project `figs/` (or unzip + upload the 5 PNGs)
- [ ] Overleaf compile clean and ship

---

## 1. §Software — URL + SHA

Paste over the red §Software block (~line 1019 / ~968 depending on draft):

```
All simulations use a single public repository:

  https://github.com/aarushgoel94-wq/iv-microplastic-transport

Tip commit (main): 2e44a7ec5afcd79f8a0fb24922fab29f619f617a

Contents: microplastic_simulation.py, parameters.yaml,
regenerate_paper_figures.py, requirements.txt, README.md, and figs/
(fig_phi_sweep.png, fig_sub20_escaped_counts.png, fig_delta_w_retention.png,
trajectory_exaggerated_pristine.png, trajectory_exaggerated_autoclaved.png).

Clone and regenerate figures:

  pip install -r requirements.txt
  python regenerate_paper_figures.py

Contradictory log C_run_log.txt was deleted; keep C_phi_sweep.txt / .csv
(base_trials=30, N_esc_P=80).
```

Source: `docs/SOFTWARE_SECTION.txt` (SHA refreshed to live tip).

---

## 2. AI disclosure

Paste near `main.tex` ~line 192:

```
Portions of simulation code, figure-generation scripts, and manuscript editing were assisted by the Cursor AI coding assistant. All physical model choices, parameter values, experimental design, data, and scientific conclusions are the author's and were verified by the author.
```

---

## 3. Lateral kick one-liner

Paste ~line 757:

```
Lateral Brownian kicks use the Stokes–Einstein step √(2 D Δt) with D = k_B T / (6 π μ r) at T = 310 K; the earlier size-independent lateral_kick_std_m = 4×10⁻⁷ m was removed.
```

---

## 4. δ_w = 300 µm origin (honest: modelling guess)

Paste ~line 1310 / ~1550 (near-wall hazard text):

```
The near-wall hazard half-width δ_w = 300 µm is a modelling choice (≈ 7.5% of the 4 mm bore per wall), not a measured adhesion-layer thickness; the δ_w sweep shows retention tracks 2δ_w/W for sub-20 µm fragments, so the operating-point retention is set by that geometric band rather than by an independent stickiness parameter.
```

---

## 5. Refs [5][6][11] — candidates only (VERIFY BY OPENING)

**Do not paste until you have opened the source.** No invented USP edition numbers.

### [5] USP ⟨788⟩ Particulate Matter in Injections

Official text lives in USP–NF Online (paywalled). Cite the edition / revision date **printed on the page you opened**. Defensible skeleton if you used the current online official text:

```
United States Pharmacopeia. ⟨788⟩ Particulate Matter in Injections.
In: USP–NF [online]. Rockville, MD: United States Pharmacopeial Convention.
```

Then append the edition/revision date from your screen. There is also a public revision PDF (`revisionGeneralChapter788.pdf`) — if that is what you used, cite that revision’s printed year.

### [6] ISO 8536 — Infusion equipment for medical use

If you meant gravity-feed single-use infusion sets, candidate (confirm on iso.org/standard/70730.html against the copy you opened):

```
ISO 8536-4:2019, Infusion equipment for medical use — Part 4:
Infusion sets for single use, gravity feed. Sixth edition, 2019-09.
```

If you meant filters / another part, cite that part instead.

### [11] Replace Hawkins (or CUT)

Strongest match for steam disinfection of disposable medical devices → MPs (incl. PP):

```
Wang, … et al. (2023). Microplastics released from disposable medical devices
and their toxic responses in Caenorhabditis elegans. Environmental Research
238, 117345. https://doi.org/10.1016/j.envres.2023.117345
```

Alternative if the claim is only “PP survives autoclave,” not shedding:

```
Fischer, K.M. & Howell, A.P. (2021). Reusability of autoclaved 3D printed
polypropylene… 3D Printing in Medicine 7, 20.
https://doi.org/10.1186/s41205-021-00111-x
```

Open the PDF, check authors/volume/pages, then replace or delete Hawkins.

Full notes: `docs/REFERENCES_FILL_IN.txt`, `docs/REFERENCES_CHECKLIST.txt`.

---

## 6. Dr Barge city — ASK HIM (do not invent)

`main.tex` ~line 87:

**Ask Dr Barge.** Deerfield, Illinois (Baxter HQ) is flagged only — **not confirmed**. Leave the placeholder until he answers.

---

## 7. Flow 20×25 table numbers (authoritative)

ū_slow = 2.763 mm/s (125 mL/hr); Δt cap scaled ∝ 1/û ≈ 9×10⁻⁴ s; seed = `42 + 10000*rep + int(d*17)`; focusing OFF; n = 20×25.

| flow | state | 3 µm | 10 µm | 50 µm |
|------|-------|------|-------|-------|
| 2260 mL/hr | pristine | 87.6±7.4 | 87.0±5.0 | 77.0±9.9 |
| 2260 mL/hr | autoclaved | 86.2±8.5 | 87.0±5.0 | 75.6±10.4 |
| 125 mL/hr | pristine | 84.4±8.5 | 80.6±7.0 | 0.0±0.0 |
| 125 mL/hr | autoclaved | 84.4±8.5 | 80.6±7.0 | 0.0±0.0 |

Plain-text paste:

```
flow          state       3 µm           10 µm          50 µm
2260_mL_hr    pristine    87.6±7.4   87.0±5.0   77.0±9.9
2260_mL_hr    autoclaved  86.2±8.5   87.0±5.0   75.6±10.4
125_mL_hr     pristine    84.4±8.5   80.6±7.0   0.0±0.0
125_mL_hr     autoclaved  84.4±8.5   80.6±7.0   0.0±0.0
```

Files: `docs/FLOW_RATE_COMPARISON.txt`, `tables/flow_rate_comparison.csv`.

---

## 8. Overleaf figures — exact 5 filenames + paths

Exact filenames (only these five for the new figs drop-in):

1. `fig_phi_sweep.png`
2. `fig_sub20_escaped_counts.png`
3. `fig_delta_w_retention.png`
4. `trajectory_exaggerated_pristine.png`
5. `trajectory_exaggerated_autoclaved.png`

| Where | Path |
|-------|------|
| Desktop zip (upload this) | `/Users/aarushgoel/Desktop/OVERLEAF_UPLOAD_figs.zip` |
| Local Overleaf project figs | `/Users/aarushgoel/Downloads/overleaf_project/figs/` |
| Workspace copy | `Claude deliverables/final_open_todos/figures/` |
| Public repo | `figs/` on GitHub tip |

**Action:** In Overleaf → project → `figs/` folder → Upload / drag `OVERLEAF_UPLOAD_figs.zip` (or the five PNGs). Agent cannot API-upload without Overleaf login credentials.

---

## End

Science Tier 1 is done. Everything above that an agent can finish is in this pack. Remaining work is human paste + verify + Overleaf drag-upload + compile.
