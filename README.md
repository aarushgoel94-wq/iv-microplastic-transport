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

## Source of truth

| File | Role |
|------|------|
| `microplastic_simulation.py` | Physics + clinical hazard model |
| `parameters.yaml` | All constants (flow scenario, hazard rates) |
| `regenerate_paper_figures.py` | Regenerates every paper figure/table |
| `archive/` | Deprecated Colab (do **not** use for reported numbers) |

## Flow scenario

Mean velocity ū = 50 mm/s through a 4 mm bore ≈ **2260 mL/hr**
(**rapid / gravity-bolus infusion**, not 100–125 mL/hr maintenance).

## License

Research / competition use. Cite the repository commit SHA with the paper.
