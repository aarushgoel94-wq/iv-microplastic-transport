# LAST PASTES — refs, AI tool name, paired note, Overleaf

Science is closed. Use this file for the four red blocks + AI name.

---

## Five figures (in this chat AND on disk)

Valid PNGs — confirmed rendered above in Cursor chat:

1. `fig_phi_sweep.png` — vertical contours, operating point 2.44 at (2.5, 1.75)
2. `fig_sub20_escaped_counts.png` — 87.3 → 302.9, A/P = 3.47×
3. `fig_delta_w_retention.png` — measured < 2δ_w/W but tracks it
4. `trajectory_exaggerated_pristine.png` — Focusing OFF; stuck vs escaped markers
5. `trajectory_exaggerated_autoclaved.png` — Focusing OFF; distinct from pristine

**Drag into Overleaf `figs/`:**
`/Users/aarushgoel/Desktop/OVERLEAF_UPLOAD_figs.zip`

Also mirrored at:
`~/Downloads/overleaf_project/figs/`
`Microplastics Research/figs/`
https://github.com/aarushgoel94-wq/iv-microplastic-transport/tree/main/figs

I cannot log into your Overleaf account from this machine.

---

## AI disclosure (use BOTH tools honestly)

Paste near the AUTHOR comment:

```
Manuscript drafting and editing were assisted by Claude (Anthropic); simulation
code, figure generation, and manuscript engineering were assisted by the Cursor
AI coding assistant. All physical model choices, parameter values, experimental
design, data, and scientific conclusions are the author's and were verified by
the author.
```

---

## Paired-design aside (~line 1046) — decision: KEEP as plain text

Drop the `\textcolor{red}{...}` wrapper; keep the scientific content as a normal
sentence. The paired inlet/defect layout is already implemented and is a strength,
not a “future improvement” to delete.

Suggested plain wording:

```
Defect sites and inlet seeds are paired across pristine and autoclaved runs
(shared spatial layout RNG, independent of wall wear), so differences in escape
counts are attributable to injection multiplicity and wear rather than resampling
noise in the defect field.
```

---

## [5] USP ⟨788⟩ — WHAT I OPENED

I opened USP’s public revision PDF of ⟨788⟩ Particulate Matter in Injections:
https://www.uspnf.com/sites/default/files/usp_pdf/EN/USPNF/revisionGeneralChapter788.pdf

That PDF is a chapter text **without a printed USP–NF edition number on the title
block** (no “USP 43–NF 38” etc.). Do **not** invent an edition.

If this public revision is what you cite, a honest form is:

```
United States Pharmacopeia. ⟨788⟩ Particulate Matter in Injections.
USP–NF revision text (public PDF). Rockville, MD: United States Pharmacopeial
Convention. Available at:
https://www.uspnf.com/sites/default/files/usp_pdf/EN/USPNF/revisionGeneralChapter788.pdf
(accessed 14 July 2026).
```

If you have USP–NF Online access, prefer the edition banner on *your* page and
replace the URL form with that edition string.

---

## [6] ISO 8536 — VERIFIED on iso.org

Opened https://www.iso.org/standard/70730.html :

```
ISO 8536-4:2019. Infusion equipment for medical use — Part 4: Infusion sets for
single use, gravity feed. Sixth edition, 2019-09. Geneva: International
Organization for Standardization.
```

(Confirmed current; last reviewed 2025; edition 6.)

---

## [11] Replace Hawkins — VERIFIED via Crossref (DOI)

Crossref for `10.1016/j.envres.2023.117345` returns **volume 239**, article
**117345** (not 238 — use Crossref):

```
Zhou, T., Wu, J., Hu, X., Cao, Z., Yang, B., Li, Y., Zhao, Y., Ding, Y.,
Liu, Y. & Xu, A. (2023). Microplastics released from disposable medical devices
and their toxic responses in Caenorhabditis elegans. Environmental Research,
239, 117345. https://doi.org/10.1016/j.envres.2023.117345
```

---

## Repo consistency (Φ-sweep N_esc_P=80)

- `C_run_log.txt` (N_esc_P=51) is **deleted**.
- Authoritative: `Claude deliverables/C_phi_sweep.txt` — base_trials=30, **N_esc_P=80**,
  ratio/Φ_deg = **0.9721 ± 0.0175**.
- `regenerate_paper_figures.py` was **fixed**: Φ-sweep now loads that CSV (refuses
  if N_esc_P ≠ 80) instead of silently re-running at base=20. Paper figs are
  shipped to `figs/` under the five Overleaf names.

Verify locally:

```bash
git clone https://github.com/aarushgoel94-wq/iv-microplastic-transport /tmp/iv-check
cd /tmp/iv-check
pip install -r requirements.txt
python -c "from pathlib import Path; import csv
rows=list(csv.DictReader(open('Claude deliverables/C_phi_sweep.csv')))
assert all(int(float(r['N_esc_P']))==80 for r in rows)
print('N_esc_P=80 OK; n_rows', len(rows))"
python regenerate_paper_figures.py   # full regen is long; phi step uses CSV
```
