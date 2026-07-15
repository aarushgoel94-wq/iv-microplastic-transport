# Verification vs validation

These are different claims. Mixing them is how numerical papers get rejected.

## Verification — *does the code solve the equations it claims to solve?*

**Yes, for the Stokes–drag integrator.**

Simplified scenario with a closed-form solution (still saline, no walls, no
clinical adhesion, no focusing, no Brownian motion):

$$
\dot{v}_y = -\beta\,v_y + a_{\mathrm{body}},\qquad
\beta = \frac{9\mu}{2\rho_p r^2},\qquad
a_{\mathrm{body}} = \frac{\rho_f-\rho_p}{\rho_p}\,g.
$$

Exact solution (zero initial velocity):

$$
v_\infty = \frac{a_{\mathrm{body}}}{\beta}
  = \frac{2 r^2 (\rho_f-\rho_p) g}{9\mu},\qquad
v(t)=v_\infty\bigl(1-e^{-\beta t}\bigr),\qquad
y(t)=v_\infty\,t - \frac{v_\infty}{\beta}\bigl(1-e^{-\beta t}\bigr).
$$

The code uses the **exact** one-step map of this ODE (exponential Stokes update),
not forward Euler. Against the closed form, the maximum relative trajectory error
across Δt ∈ [10⁻⁵, 2×10⁻²] s is **1.012e-15**
(PASS — numerical trajectory matches the closed-form Stokes solution to machine precision for every tested Δt (exact update, not Euler).).

Figures: `fig_verification_stokes.png`.
Table: `tables/stokes_verification_errors.csv`.

### Time-step independence (same 25 trials)

Holding inlet seeds fixed (seed = 42 + ⌊d·17⌋ + trial_id) and varying only Δt:

| Δt (s) | escape % (10 µm) | mean transit (s) |
|--------|--------------------------|------------------|
| 1.0e-05 | 92.0 | 2.783 |
| 2.5e-05 | 92.0 | 2.783 |
| 5.0e-05 | 92.0 | 2.783 |
| 1.0e-04 | 92.0 | 2.783 |
| 2.5e-04 | 92.0 | 2.783 |
| 5.0e-04 | 92.0 | 2.783 |
| 1.0e-03 | 92.0 | 2.783 |

Paper fast_mode cap Δt = 5×10⁻⁵ s → **92.0%** escape.
Across Δt ≤ 5×10⁻⁴ s the escape fraction spans only **0.0 pp**.
Outcomes are therefore timestep-independent in the geometric-accuracy regime
used for the paper; coarser Δt eventually drifts when wall-clearance sampling
degrades (expected — geometric accuracy, not stability, sets the Δt ceiling for
an exact Stokes integrator).

Figure: `fig_timestep_independence.png`.

## Validation — *does the model represent clinical IV transport?*

**Not demonstrated by this study, and we do not claim it.**

Validation would require measured fragment escape / retention in real tubing
(counts vs size, pristine vs autoclaved) against which to compare Table 6.
We do not have those data here.

What we *do* have is **structural** evidence inside the model:

- near-wall hazard saturation ⇒ Φ_wall and α are weak levers;
- retention tracks 2δ_w/W for sub-20 µm fragments;
- paired seeding isolates injection multiplicity from wall wear.

Those are consistency checks and sensitivity tests — useful, but they sit on the
**verification / sensitivity** side of the boundary. Calling them “validation”
would overclaim.

## One sentence for the paper

> The Stokes integrator was verified against the closed-form buoyant-sphere
> solution (machine-precision agreement for all tested Δt); clinical escape
> fractions were shown to be independent of Δt across the geometric-accuracy
> window used in the campaign; comparison to measured patient-line particle
> counts (external validation) is left for future work.

## Pastes for figure captions

**fig_verification_stokes.png** — Verification against the exact Stokes solution
for a 50 µm polypropylene sphere in still saline
(no walls, no clinical adhesion). Left: exact $y(t)$ vs the production
integrator. Right: global trajectory error vs Δt (floor at rounding precision).

**fig_timestep_independence.png** — Time-step independence for the *full* clinical
tube model (focusing OFF). Same 25 inlet seeds; only Δt varies. Escape fraction
and mean transit time are flat through the paper’s Δt = 5×10⁻⁵ s operating point.
