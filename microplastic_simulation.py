"""
2D microplastic particle simulation in medical tubing.

Models laminar (Poiseuille) flow, Stokes drag, gravity/buoyancy, and wall collisions.
Particle states: 'flowing', 'escaped', 'stuck'.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field, replace
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

EXPERIMENTAL_DIAMETERS_UM: Tuple[float, ...] = (3.0, 10.0, 50.0)
TRIALS_PER_GROUP: int = 100
SUB_DANGEROUS_SIZE_THRESHOLD_UM: float = 20.0


# ---------------------------------------------------------------------------
# Physical constants (SI units throughout)
# ---------------------------------------------------------------------------

GRAVITY = 9.81  # m/s²

# Polypropylene (PP) — typical injection-moulded grade
POLYPROPYLENE_DENSITY = 905.0  # kg/m³

# Saline / blood-like carrier fluid (defaults toward physiological saline)
SALINE_DENSITY = 1000.0  # kg/m³
SALINE_DYNAMIC_VISCOSITY = 1.0e-3  # Pa·s (~1 mPa·s at ~37 °C)

BLOOD_DENSITY = 1060.0  # kg/m³
BLOOD_DYNAMIC_VISCOSITY = 3.5e-3  # Pa·s (shear-thinning; representative bulk value)

# Euler integration limits (dt in seconds — NOT milliseconds)
EULER_RELAXATION_STEPS: int = 50
EULER_CFL_FRACTION: float = 0.1
EULER_WALL_SAFETY: float = 0.25
DT_MIN: float = 1e-9
DT_MAX: float = 1e-5
# fast_mode geometric dt cap at the paper's rapid û = 50 mm/s. Exact Stokes
# integration is unconditionally stable, so at slower û the same wall-clearance
# accuracy allows dt ∝ 1/û (e.g. 125 mL/hr ≈ 18× slower → cap ≈ 9×10⁻⁴ s).
FAST_MODE_DT_CAP_AT_U_REF: float = 5.0e-5
FAST_MODE_U_REF_M_S: float = 0.05
FAST_MODE_DT_CAP_MAX: float = 2.0e-3
MAX_TRAJECTORY_POINTS: int = 10_000  # stored points per trial (plots need far fewer)
MAX_PHYSICS_STEPS: int = 500_000    # hard cap per trial — prevents Colab timeouts
DEFAULT_MAX_TIME_S: float = 20.0    # most particles finish within ~5 s simulated time
COLAB_TRIALS_PER_GROUP: int = 25    # use 100 only on a machine with hours to spare

# Drag-dominated trajectory coupling (Stokes + Poiseuille)
DRAG_LOCK_BETA_DT: float = 5.0
"""When β·dt exceeds this, axial speed snaps to the local parabolic stream u(y)."""

TUBULE_LATERAL_RELAX_RATE: float = 100.0
"""1/s — exponential pull toward centerline, scaled by log(β/β_ref) for small particles."""

REFERENCE_FOCUS_DIAMETER_UM: float = 50.0
"""Particles at or above this size receive no centerline focusing — buoyancy drives wall drift."""


# ---------------------------------------------------------------------------
# Clinical / hospital-line realism (beyond ideal laminar physics)
# ---------------------------------------------------------------------------

BOLTZMANN_K = 1.380649e-23  # J/K
BODY_TEMPERATURE_K = 310.0  # ~37 °C saline in IV line


@dataclass(frozen=True)
class ClinicalRealismConfig:
    """
    Stochastic retention mechanisms observed in real IV lines:

    - Electrostatic adhesion on polymer tubing (near-wall charge patches only)
    - Micro-asperities / jagged edges (near-wall gated; cannot snag on centerline)
    - Elevated capture at connectors, clamps, and injection ports
    - Optional tubule-focusing heuristic (focus_enabled) scaled by focus_efficiency
    - Thermal Brownian lateral kicks (size-dependent Stokes–Einstein)
    """

    enabled: bool = True
    static_rate_near_wall_per_s: float = 0.18
    """Hazard rate when particle is in the slow-flow band near a wall."""
    static_rate_bulk_per_s: float = 0.0
    """MUST stay 0 — no adhesive surface exists in the bulk lumen."""
    near_wall_band_m: float = 0.00030
    """Within 300 µm of wall surface → static risk AND defect snag eligibility."""
    defect_density_per_m: float = 28.0
    """Mean count of roughness / flash defect sites per meter (layout; wear scales rate)."""
    defect_capture_prob: float = 0.028
    """One-shot snagging probability when a near-wall particle first enters a defect zone."""
    defect_capture_width_m: float = 0.0006
    """Axial window (0.6 mm) around each defect where snagging is evaluated once."""
    connector_zone_fractions: Tuple[Tuple[float, float], ...] = (
        (0.05, 0.10),   # proximal injection port / spike
        (0.30, 0.37),   # mid-line connector or clamp
        (0.60, 0.68),   # distal roller / access port
    )
    connector_rate_multiplier: float = 3.0
    focus_enabled: bool = False
    """Paper main = OFF. Set True only for Appendix A focusing ablation."""
    focus_efficiency: float = 0.88
    """
    Multiplier on the focusing blend α each step: α_eff = focus_efficiency × (1 − e^{−λdt}).
    Models imperfect wall-repulsion; 1.0 = full heuristic, 0.0 = no focusing.
    Inert when focus_enabled is False.
    """
    use_brownian_kick: bool = True
    """If True, lateral kick uses Stokes–Einstein √(2 D Δt); if False, no lateral kick."""
    temperature_k: float = BODY_TEMPERATURE_K
    size_snag_exponent: float = 0.35
    """Empirical size weight (d/3)^α. Sensitivity: vary α ∈ [0, 1]."""
    degraded_tubing_multiplier: float = 1.0
    """
    Wall / adhesion wear factor w. Applied directly to hazard rates (paired layout).
    Pristine: w=1.0. Autoclaved operating point: w=1.75. Set explicitly in config —
    do NOT gate on container_key (that bug made Φ_wall sweeps inert).
    """
    paired_defect_layout: bool = True
    """
    If True, defect x-positions are drawn with density ρ_def·L (shared layout across P/A);
    only the adhesion multiplier w differs.
    """


AUTOCLAVE_WALL_WEAR = 1.75  # operating-point Φ_wall for autoclaved containers


DEFAULT_CLINICAL_REALISM = ClinicalRealismConfig()


@dataclass
class ClinicalRealismContext:
    """Per-trial hospital-line hazard context (defect sites drawn once per trial)."""

    config: ClinicalRealismConfig
    environment: TubeEnvironment
    defect_x_m: np.ndarray
    tubing_wear_multiplier: float
    rng: np.random.Generator
    _triggered_defects: set = field(default_factory=set, repr=False)

    @classmethod
    def for_trial(
        cls,
        environment: TubeEnvironment,
        rng: np.random.Generator,
        *,
        config: ClinicalRealismConfig = DEFAULT_CLINICAL_REALISM,
        container_key: str = "pristine",
        shared_defect_x_m: Optional[np.ndarray] = None,
    ) -> ClinicalRealismContext:
        # Tier-0 bugfix: w MUST come from the config field itself, never from container_key alone.
        wear = float(config.degraded_tubing_multiplier)
        if shared_defect_x_m is not None:
            defect_x = np.asarray(shared_defect_x_m, dtype=float)
        elif config.paired_defect_layout:
            # Paired design: same spatial defect field for pristine and autoclaved.
            expected_defects = config.defect_density_per_m * environment.length  # wear=1 for layout
            n_defects = int(rng.poisson(expected_defects))
            defect_x = (
                np.sort(rng.uniform(0.0, environment.length, size=n_defects))
                if n_defects > 0
                else np.array([], dtype=float)
            )
        else:
            expected_defects = config.defect_density_per_m * environment.length * wear
            n_defects = int(rng.poisson(expected_defects))
            defect_x = (
                np.sort(rng.uniform(0.0, environment.length, size=n_defects))
                if n_defects > 0
                else np.array([], dtype=float)
            )
        return cls(
            config=config,
            environment=environment,
            defect_x_m=defect_x,
            tubing_wear_multiplier=wear,
            rng=rng,
        )

    def _in_connector_zone(self, x_m: float) -> bool:
        length = self.environment.length
        for start_frac, end_frac in self.config.connector_zone_fractions:
            if start_frac * length <= x_m <= end_frac * length:
                return True
        return False

    def adjust_focus_blend(self, focus_blend: float) -> float:
        if not self.config.focus_enabled:
            return 0.0
        return focus_blend * self.config.focus_efficiency

    def brownian_kick_std(self, radius_m: float, dynamic_viscosity: float, dt: float) -> float:
        """Stokes–Einstein lateral step std: √(2 D Δt), D = kT / (6 π μ r)."""
        if not self.config.use_brownian_kick or dt <= 0.0 or radius_m <= 0.0:
            return 0.0
        D = BOLTZMANN_K * self.config.temperature_k / (6.0 * np.pi * dynamic_viscosity * radius_m)
        return float(np.sqrt(2.0 * D * dt))

    def lateral_imperfection_kick(self, radius_m: float = 1.5e-6, dynamic_viscosity: float = 1e-3, dt: float = 5e-5) -> float:
        std = self.brownian_kick_std(radius_m, dynamic_viscosity, dt)
        if std <= 0.0:
            return 0.0
        return float(self.rng.normal(0.0, std))

    def size_snag_factor(self, diameter_um: float) -> float:
        return (float(diameter_um) / 3.0) ** self.config.size_snag_exponent

    def try_clinical_retention(self, particle: "Particle", dt: float) -> Optional[str]:
        """
        Stochastic retention check after hydrodynamic step.

        Returns retention mechanism label if particle sticks, else None.
        Defect snags and electrostatic adhesion require near-wall proximity.
        """
        if not self.config.enabled or particle.state != ParticleState.FLOWING:
            return None

        cfg = self.config
        x_m = float(particle.position[0])
        clearance_bottom, clearance_top = self.environment.wall_clearances(
            particle.position, particle.radius
        )
        min_gap = min(clearance_bottom, clearance_top)
        near_wall = min_gap <= cfg.near_wall_band_m
        size_snag = self.size_snag_factor(particle.diameter_um)

        # Jagged-edge defect — near-wall gate + one evaluation per defect per trial
        if near_wall:
            for defect_idx, defect_x in enumerate(self.defect_x_m):
                if defect_idx in self._triggered_defects:
                    continue
                if abs(x_m - defect_x) <= cfg.defect_capture_width_m:
                    self._triggered_defects.add(defect_idx)
                    p_snag = cfg.defect_capture_prob * self.tubing_wear_multiplier * size_snag
                    if self.rng.random() < p_snag:
                        return "jagged_edge_defect"

        # Electrostatic adhesion — near-wall only (bulk rate is identically zero)
        if not near_wall:
            return None

        rate = cfg.static_rate_near_wall_per_s
        mechanism = "electrostatic_near_wall"

        if self._in_connector_zone(x_m):
            rate *= cfg.connector_rate_multiplier
            mechanism = "connector_static_clamp"

        rate *= self.tubing_wear_multiplier
        rate *= size_snag

        if rate <= 0.0:
            return None

        hazard = 1.0 - np.exp(-rate * dt)
        if self.rng.random() < hazard:
            return mechanism
        return None


class ParticleState(str, Enum):
    FLOWING = "flowing"
    ESCAPED = "escaped"
    STUCK = "stuck"


def diameter_um_to_radius_m(diameter_um: float) -> float:
    """Convert particle diameter (µm) to radius (m)."""
    return (float(diameter_um) * 1e-6) / 2.0


def normalize_diameter_um(diameter_um: float) -> float:
    """Canonical float key for diameter dict lookups across groups."""
    return float(diameter_um)


# ---------------------------------------------------------------------------
# IV / medical container profiles (sterilization history)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ContainerProfile:
    """
    Models how container sterilization history affects microplastic release.

    Attributes:
        degradation_factor: Multiplier on total particle injection vs. baseline.
        size_weights: Relative release probability per diameter (µm). Normalized
            internally to allocate trials across size groups.
    """

    key: str
    display_name: str
    comparison_label: str
    degradation_factor: float
    size_weights: Dict[float, float]

    def normalized_weights(self, diameters_um: Tuple[float, ...]) -> Dict[float, float]:
        weights = {d: self.size_weights.get(d, 0.0) for d in diameters_um}
        total = sum(weights.values())
        if total <= 0.0:
            raise ValueError(f"Container profile '{self.key}' has no positive size weights.")
        return {d: w / total for d, w in weights.items()}

    def trial_counts(
        self,
        diameters_um: Tuple[float, ...] = EXPERIMENTAL_DIAMETERS_UM,
        base_trials_per_size: int = TRIALS_PER_GROUP,
    ) -> Dict[float, int]:
        """
        Allocate simulation trials per size group.

        Baseline budget = base_trials_per_size × number of size groups.
        Degraded containers scale total injections and skew toward smaller fragments.
        """
        normalized = self.normalized_weights(diameters_um)
        total_budget = int(
            round(base_trials_per_size * len(diameters_um) * self.degradation_factor)
        )
        raw = {d: total_budget * normalized[d] for d in diameters_um}
        counts = {d: max(1, int(round(n))) for d, n in raw.items()}

        # Preserve total budget after integer rounding
        diff = total_budget - sum(counts.values())
        if diff != 0:
            adjust_d = max(diameters_um, key=lambda d: raw[d])
            counts[adjust_d] = max(1, counts[adjust_d] + diff)

        return counts

    @property
    def description(self) -> str:
        return (
            f"{self.display_name} — degradation ×{self.degradation_factor:.1f}, "
            f"label: {self.comparison_label}"
        )

    def injection_summary(
        self,
        diameters_um: Tuple[float, ...] = EXPERIMENTAL_DIAMETERS_UM,
        base_trials_per_size: int = TRIALS_PER_GROUP,
    ) -> dict:
        counts = self.trial_counts(diameters_um, base_trials_per_size)
        total = sum(counts.values())
        sub_threshold = sum(
            n for d, n in counts.items() if d < SUB_DANGEROUS_SIZE_THRESHOLD_UM
        )
        return {
            "container_key": self.key,
            "container_name": self.display_name,
            "comparison_label": self.comparison_label,
            "degradation_factor": self.degradation_factor,
            "trial_counts": counts,
            "total_particles": total,
            "sub_20um_particles": sub_threshold,
            "sub_20um_fraction": sub_threshold / total if total else 0.0,
            "baseline_total": base_trials_per_size * len(diameters_um),
        }


CONTAINER_PROFILES: Dict[str, ContainerProfile] = {
    "pristine": ContainerProfile(
        key="pristine",
        display_name="Pristine / Standard Container",
        comparison_label="Safe",
        degradation_factor=1.0,
        size_weights={3.0: 1.0, 10.0: 1.0, 50.0: 1.0},
    ),
    "autoclaved": ContainerProfile(
        key="autoclaved",
        display_name="Autoclaved / Harsh Sterilization Container",
        comparison_label="Sterilization-Degraded",
        degradation_factor=2.5,
        size_weights={3.0: 2.8, 10.0: 2.2, 50.0: 0.35},
    ),
}


def list_container_profiles() -> List[ContainerProfile]:
    return list(CONTAINER_PROFILES.values())


def get_container_profile(key: str) -> ContainerProfile:
    normalized = key.strip().lower().replace(" ", "_")
    aliases = {
        "standard": "pristine",
        "safe": "pristine",
        "pristine_standard": "pristine",
        "harsh": "autoclaved",
        "sterilization_degraded": "autoclaved",
        "gamma": "autoclaved",
    }
    resolved = aliases.get(normalized, normalized)
    if resolved not in CONTAINER_PROFILES:
        valid = ", ".join(sorted(CONTAINER_PROFILES))
        raise ValueError(f"Unknown container type '{key}'. Valid options: {valid}")
    return CONTAINER_PROFILES[resolved]


def prompt_container_selection() -> ContainerProfile:
    """Interactive menu for selecting the IV container sterilization profile."""
    profiles = list_container_profiles()
    print()
    print("=" * 62)
    print("  SELECT IV / MEDICAL CONTAINER TYPE")
    print("=" * 62)
    for index, profile in enumerate(profiles, start=1):
        summary = profile.injection_summary()
        counts = summary["trial_counts"]
        count_str = ", ".join(f"{d:.0f} µm: {counts[d]}" for d in sorted(counts))
        print(f"\n  [{index}] {profile.display_name}")
        print(f"      Comparison label: {profile.comparison_label}")
        print(f"      Degradation factor: ×{profile.degradation_factor:.1f}")
        print(f"      Planned injections ({summary['total_particles']} particles): {count_str}")
        print(
            f"      Sub-{SUB_DANGEROUS_SIZE_THRESHOLD_UM:.0f} µm share: "
            f"{summary['sub_20um_fraction']:.1%}"
        )

    while True:
        try:
            choice = input("\nEnter container number (1–2) or key [pristine]: ").strip()
        except EOFError:
            print("\nNo input detected; defaulting to Pristine / Standard Container.")
            return CONTAINER_PROFILES["pristine"]

        if not choice:
            return CONTAINER_PROFILES["pristine"]

        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(profiles):
                return profiles[idx - 1]
            print(f"  Please enter a number between 1 and {len(profiles)}.")
            continue

        try:
            return get_container_profile(choice)
        except ValueError as exc:
            print(f"  {exc}")


def resolve_container_profile(argv: Optional[List[str]] = None) -> ContainerProfile:
    """Resolve container from CLI flag or interactive prompt."""
    parser = argparse.ArgumentParser(description="Microplastic medical tubing simulation")
    parser.add_argument(
        "--container",
        "-c",
        choices=list(CONTAINER_PROFILES.keys()),
        help="Container profile: pristine (safe baseline) or autoclaved (degraded)",
    )
    args, _ = parser.parse_known_args(argv)

    if args.container:
        profile = get_container_profile(args.container)
        print(f"\nContainer selected via CLI: {profile.display_name} [{profile.comparison_label}]")
        return profile

    if sys.stdin.isatty():
        return prompt_container_selection()

    print("Non-interactive mode: defaulting to Pristine / Standard Container.")
    return CONTAINER_PROFILES["pristine"]


@dataclass
class FluidProperties:
    """
    Carrier fluid parameters.

    Default average_velocity = 0.05 m/s corresponds to ~2260 mL/hr through a 4 mm
    bore — a *rapid / gravity-bolus infusion* scenario (not 100–125 mL/hr maintenance).
    See parameters.yaml and §Methods for the explicit scenario label (Tier 1.6).
    """

    density: float = SALINE_DENSITY
    dynamic_viscosity: float = SALINE_DYNAMIC_VISCOSITY
    average_velocity: float = 0.05  # m/s mean axial speed (rapid infusion)


@dataclass
class TubeEnvironment:
    """
    2D cross-section of medical tubing: axial x ∈ [0, length], lateral y ∈ [0, width].

    Flow is laminar Poiseuille (parabolic) between parallel walls at y = 0 and y = width.
    """

    length: float = 0.100  # m (100 mm)
    width: float = 0.004  # m (4 mm diameter)
    fluid: FluidProperties = field(default_factory=FluidProperties)

    @property
    def centerline_y(self) -> float:
        return self.width / 2.0

    @property
    def max_centerline_velocity(self) -> float:
        """Peak axial speed at the channel center (u_max = 1.5 × u_avg for parabolic profile)."""
        return 1.5 * self.fluid.average_velocity

    def axial_velocity(self, y: float) -> float:
        """
        Parabolic Poiseuille profile: u(y) = u_max [1 − ((y − y_c) / (W/2))²].

        Zero at both walls, maximum at the centerline.
        """
        half_width = self.width / 2.0
        normalized = (y - self.centerline_y) / half_width
        return self.max_centerline_velocity * (1.0 - normalized**2)

    def fluid_velocity_at(self, position: np.ndarray) -> np.ndarray:
        """Return fluid velocity vector [u_x, u_y] at particle position."""
        return np.array([self.axial_velocity(float(position[1])), 0.0])

    def verify_parabolic_profile(self) -> None:
        """Debug: confirm u=0 at walls and u=u_max at centerline."""
        u_wall = self.axial_velocity(0.0)
        u_top = self.axial_velocity(self.width)
        u_center = self.axial_velocity(self.centerline_y)
        assert abs(u_wall) < 1e-15 and abs(u_top) < 1e-15, "Wall velocity must be zero"
        assert u_center == self.max_centerline_velocity, "Centerline must equal u_max"
        half = self.width / 2.0
        u_quarter = self.axial_velocity(self.centerline_y + 0.5 * half)
        expected = self.max_centerline_velocity * 0.75
        assert abs(u_quarter - expected) < 1e-12, "Parabolic profile shape incorrect"

    def contains(self, position: np.ndarray, radius: float = 0.0) -> bool:
        """True if the particle (including radius) lies inside the tube interior."""
        x, y = position
        return (
            radius <= x <= self.length - radius
            and radius <= y <= self.width - radius
        )

    def sphere_bottom_y(self, center_y: float, radius: float) -> float:
        """Axial cross-section: lowest point of the sphere (y − r)."""
        return center_y - radius

    def sphere_top_y(self, center_y: float, radius: float) -> float:
        """Axial cross-section: highest point of the sphere (y + r)."""
        return center_y + radius

    def wall_clearances(self, position: np.ndarray, radius: float) -> Tuple[float, float]:
        """
        Gap from sphere surface to bottom wall (y=0) and top wall (y=width).

        Equivalent to:
          bottom contact when center_y <= radius   (surface_y_bottom <= 0)
          top contact    when center_y >= width − radius (surface_y_top >= width)

        Larger radius shrinks both clearances at the same center position.
        Returns (clearance_bottom, clearance_top) in meters.
        """
        _, y = position
        clearance_bottom = y - radius
        clearance_top = self.width - (y + radius)
        return clearance_bottom, clearance_top

    def wall_contact_thresholds(self, radius: float) -> dict:
        """
        Center-y positions where a sphere of this radius contacts each wall.

        Larger particles contact the top/bottom walls at lower center heights
        (e.g. 50 µm contacts top at y = width − 25 µm; 3 µm at y = width − 1.5 µm).
        """
        return {
            "bottom_contact_y_m": radius,
            "top_contact_y_m": self.width - radius,
            "interior_band_m": self.width - 2.0 * radius,
            "radius_m": radius,
        }

    def check_wall_collision(self, position: np.ndarray, radius: float) -> bool:
        """
        True when the sphere surface touches or crosses the top or bottom wall.

        Uses physical radius: bottom if y <= radius, top if y >= width − radius.
        """
        _, y = position
        return self.sphere_bottom_y(y, radius) <= 0.0 or self.sphere_top_y(y, radius) >= self.width

    def snap_to_wall_contact(self, position: np.ndarray, radius: float) -> np.ndarray:
        """Clamp center position so the sphere rests against the contacted wall."""
        pos = position.copy()
        _, y = pos
        penetration_bottom = radius - y
        penetration_top = (y + radius) - self.width
        if penetration_bottom > 0.0 and penetration_top > 0.0:
            if penetration_top >= penetration_bottom:
                pos[1] = self.width - radius
            else:
                pos[1] = radius
        elif penetration_bottom > 0.0:
            pos[1] = radius
        elif penetration_top > 0.0:
            pos[1] = self.width - radius
        return pos


@dataclass
class Particle:
    """
    Spherical microplastic particle with polypropylene density.

    Size-dependent physics (all derived from radius r):
      volume  V = (4/3) π r³
      mass    m = ρ_p V
      drag    F_d = 6 π μ r (u_fluid − u_particle)   [Stokes, per axis]
      weight  F_g = −ρ_p V g  (downward)
      buoyancy F_b = +ρ_f V g  (upward)
    """

    position: np.ndarray
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(2))
    diameter_um: float = 0.0
    density: float = POLYPROPYLENE_DENSITY
    state: ParticleState = ParticleState.FLOWING
    id: int = 0
    retention_mechanism: Optional[str] = None

    def __post_init__(self) -> None:
        self.position = np.asarray(self.position, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)
        if self.diameter_um <= 0.0:
            raise ValueError("Particle requires a positive diameter_um (micrometers).")
        self._sync_size_properties()

    def _sync_size_properties(self) -> None:
        """Cache volume and mass from the current diameter."""
        r = self.radius
        self._volume_m3 = (4.0 / 3.0) * np.pi * r**3
        self._mass_kg = self.density * self._volume_m3

    @property
    def radius(self) -> float:
        """Radius in meters, always computed from diameter_um."""
        return diameter_um_to_radius_m(self.diameter_um)

    @property
    def mass(self) -> float:
        return self._mass_kg

    @property
    def volume(self) -> float:
        return self._volume_m3

    @property
    def projected_area(self) -> float:
        return np.pi * self.radius**2

    @classmethod
    def from_diameter_um(
        cls,
        diameter_um: float,
        position: np.ndarray,
        *,
        trial_id: int = 0,
        density: float = POLYPROPYLENE_DENSITY,
    ) -> Particle:
        """Create a particle from experimental diameter (µm) and inlet position."""
        return cls(
            position=position,
            diameter_um=float(diameter_um),
            density=density,
            id=trial_id,
        )

    def stokes_drag_coefficient(self, dynamic_viscosity: float) -> float:
        """6 π μ r — multiplies relative velocity to yield drag force (N)."""
        return 6.0 * np.pi * dynamic_viscosity * self.radius

    def reynolds_number(self, relative_speed: float, dynamic_viscosity: float, fluid_density: float) -> float:
        if relative_speed < 1e-20:
            return 0.0
        return 2.0 * self.radius * fluid_density * relative_speed / dynamic_viscosity

    def compute_forces(self, environment: TubeEnvironment) -> Tuple[np.ndarray, dict]:
        """
        Size-dependent force summation for diagnostics and legacy Euler paths.

        Returns total force vector and a breakdown dict for diagnostics.
        """
        mu = environment.fluid.dynamic_viscosity
        rho_f = environment.fluid.density
        r = self.radius
        volume = self.volume
        mass = self.mass

        fluid_velocity = environment.fluid_velocity_at(self.position)
        relative_velocity = fluid_velocity - self.velocity

        # Stokes drag: F_drag = 6 π μ r (u_fluid − u_particle)
        drag_coefficient = self.stokes_drag_coefficient(mu)
        drag_force = drag_coefficient * relative_velocity

        # Gravity (down) and buoyancy (up) both scale with volume ∝ r³
        gravity_force = np.array([0.0, -mass * GRAVITY])
        buoyancy_force = np.array([0.0, rho_f * volume * GRAVITY])
        body_force = gravity_force + buoyancy_force

        total_force = drag_force + body_force

        breakdown = {
            "radius_m": r,
            "diameter_um": self.diameter_um,
            "mass_kg": mass,
            "volume_m3": volume,
            "drag_force": drag_force.copy(),
            "gravity_force": gravity_force.copy(),
            "buoyancy_force": buoyancy_force.copy(),
            "body_force": body_force.copy(),
            "total_force": total_force.copy(),
        }
        return total_force, breakdown

    def drag_force(
        self,
        fluid_velocity: np.ndarray,
        fluid_density: float,
        dynamic_viscosity: float,
    ) -> np.ndarray:
        relative_velocity = fluid_velocity - self.velocity
        return self.stokes_drag_coefficient(dynamic_viscosity) * relative_velocity

    def buoyancy_gravity_force(self, fluid_density: float) -> np.ndarray:
        net_vertical = (fluid_density - self.density) * self.volume * GRAVITY
        return np.array([0.0, net_vertical])

    def total_force(self, environment: TubeEnvironment) -> np.ndarray:
        return self.compute_forces(environment)[0]

    def acceleration(self, environment: TubeEnvironment) -> np.ndarray:
        """a = F_total / m — drag term scales ∝ 1/r², buoyancy/gravity ∝ 1 (vertical)."""
        force, _ = self.compute_forces(environment)
        return force / self.mass

    def terminal_buoyancy_rise_velocity(self, fluid_density: float, dynamic_viscosity: float) -> float:
        """Stokes terminal rise speed for buoyant PP sphere: v_t = 2 r² Δρ g / (9 μ)."""
        delta_rho = fluid_density - self.density
        return 2.0 * self.radius**2 * delta_rho * GRAVITY / (9.0 * dynamic_viscosity)

    def stokes_drag_rate(self, dynamic_viscosity: float) -> float:
        """
        Drag coupling rate β = (6 π μ r) / m = 9 μ / (2 ρ_p r²).

        Scales as 1/r²: small particles lock onto the local fluid velocity;
        large particles lag the stream and respond more to body forces.
        """
        return self.stokes_drag_coefficient(dynamic_viscosity) / self.mass

    def net_body_acceleration_y(self, fluid_density: float) -> float:
        """Vertical acceleration from weight + buoyancy: a_y = (ρ_f − ρ_p)/ρ_p × g."""
        return (fluid_density - self.density) / self.density * GRAVITY

    def apply_tubule_centerline_focus(
        self,
        environment: TubeEnvironment,
        y_m: float,
        dynamic_viscosity: float,
        dt: float,
        clinical: Optional[ClinicalRealismContext] = None,
    ) -> float:
        """
        Drag-dominated tubule focusing: small particles relax toward y = centerline
        where parabolic flow is fastest. No effect at or below the 50 µm reference size.
        focus_efficiency scales the blend; focus_enabled=False disables entirely.
        """
        y_out = y_m
        beta = self.stokes_drag_rate(dynamic_viscosity)
        beta_ref = stokes_drag_rate_from_radius(
            diameter_um_to_radius_m(REFERENCE_FOCUS_DIAMETER_UM),
            dynamic_viscosity,
            self.density,
        )
        focusing_on = True
        if clinical is not None and not clinical.config.focus_enabled:
            focusing_on = False

        if focusing_on and beta > beta_ref:
            y_c = environment.centerline_y
            lateral_rate = TUBULE_LATERAL_RELAX_RATE * np.log1p(beta / beta_ref)
            focus_blend = 1.0 - np.exp(-lateral_rate * dt)
            if clinical is not None:
                focus_blend = clinical.adjust_focus_blend(focus_blend)
            y_out = y_m + focus_blend * (y_c - y_m)

        if clinical is not None:
            y_out += clinical.lateral_imperfection_kick(
                radius_m=self.radius,
                dynamic_viscosity=dynamic_viscosity,
                dt=dt,
            )
        # Keep center inside tube interior for this radius
        y_out = float(np.clip(y_out, self.radius, environment.width - self.radius))
        return y_out

    def _stick_to_wall(
        self,
        environment: TubeEnvironment,
        mechanism: str = "hydrodynamic_wall",
    ) -> None:
        """Mark stuck and snap using this particle's physical radius (always from diameter_um)."""
        r = self.radius
        self.state = ParticleState.STUCK
        self.retention_mechanism = mechanism
        self.position[:] = environment.snap_to_wall_contact(self.position, r)
        self.velocity[:] = 0.0

    def _stick_clinical(
        self,
        environment: TubeEnvironment,
        mechanism: str,
    ) -> None:
        """Clinical retention — snap to nearest wall (static / defect capture on tubing surface)."""
        r = self.radius
        self.state = ParticleState.STUCK
        self.retention_mechanism = mechanism
        if environment.check_wall_collision(self.position, r):
            self.position[:] = environment.snap_to_wall_contact(self.position, r)
        else:
            clearance_bottom, clearance_top = environment.wall_clearances(self.position, r)
            self.position[1] = r if clearance_bottom <= clearance_top else environment.width - r
        self.velocity[:] = 0.0

    def update(
        self,
        environment: TubeEnvironment,
        dt: float,
        clinical: Optional[ClinicalRealismContext] = None,
    ) -> None:
        """
        Advance one timestep with exact Stokes-drag integration.

        Vertical motion: buoyancy (all sizes) + tubule focusing (small particles only).
        Axial motion: coupled to local parabolic stream u(y) at the updated height;
        drag-dominated particles (β·dt large) ride the center fast lane.
        """
        if self.state != ParticleState.FLOWING:
            return

        if environment.check_wall_collision(self.position, self.radius):
            self._stick_to_wall(environment)
            return

        mu = environment.fluid.dynamic_viscosity
        rho_f = environment.fluid.density
        beta = self.stokes_drag_rate(mu)
        a_body_y = self.net_body_acceleration_y(rho_f)
        v_term_y = a_body_y / beta

        v0 = self.velocity.copy()
        decay = np.exp(-beta * dt)
        if beta > 1.0e-30:
            drift_factor = (1.0 - decay) / beta
        else:
            drift_factor = dt

        y0 = float(self.position[1])
        self.velocity[1] = v_term_y + (v0[1] - v_term_y) * decay
        dy = v_term_y * dt + (v0[1] - v_term_y) * drift_factor
        y1 = self.apply_tubule_centerline_focus(environment, y0 + dy, mu, dt, clinical=clinical)
        self.position[1] = y1

        u_stream = environment.axial_velocity(0.5 * (y0 + y1))
        if beta * dt >= DRAG_LOCK_BETA_DT:
            self.velocity[0] = u_stream
        else:
            self.velocity[0] = u_stream + (v0[0] - u_stream) * decay
        self.position[0] += u_stream * dt + (v0[0] - u_stream) * drift_factor

        if environment.check_wall_collision(self.position, self.radius):
            self._stick_to_wall(environment)
            return

        if self.position[0] - self.radius >= environment.length:
            self.state = ParticleState.ESCAPED
            return

        if self.position[0] < self.radius:
            self.position[0] = self.radius
            self.velocity[0] = max(self.velocity[0], 0.0)


def debug_print_particle_physics(
    environment: TubeEnvironment,
    diameter_um: float,
    *,
    label: str = "",
) -> Particle:
    """
    Debug helper: print mass and drag for a probe particle at this diameter.
    Called before each group's trials to verify size-dependent physics.
    """
    diameter_um = normalize_diameter_um(diameter_um)
    radius = diameter_um_to_radius_m(diameter_um)
    probe = Particle.from_diameter_um(
        diameter_um,
        np.array([radius, environment.centerline_y]),
        trial_id=-1,
    )
    mu = environment.fluid.dynamic_viscosity
    test_fluid_velocity = np.array([0.05, 0.0])
    drag_force = probe.stokes_drag_coefficient(mu) * (test_fluid_velocity - probe.velocity)
    _, breakdown = probe.compute_forces(environment)

    tag = f"{label} " if label else ""
    print(
        f"  [DEBUG {tag}d={diameter_um:.0f} µm] "
        f"particle.diameter_um={probe.diameter_um:.1f}  "
        f"r={probe.radius * 1e6:.3f} µm  "
        f"mass={probe.mass:.4e} kg  "
        f"|F_drag|@0.05 m/s={np.linalg.norm(drag_force):.4e} N  "
        f"F_body_y={breakdown['body_force'][1]:.4e} N"
    )
    return probe


def debug_compare_reference_particles(environment: TubeEnvironment) -> None:
    """Print side-by-side physics for 3 µm vs 50 µm before trials start."""
    print("\n  [DEBUG] Physics check — 3 µm vs 50 µm (must differ before trials run):")
    p_small = debug_print_particle_physics(environment, 3.0, label="small")
    p_large = debug_print_particle_physics(environment, 50.0, label="large")
    mu = environment.fluid.dynamic_viscosity
    rho_f = environment.fluid.density
    mass_ratio = p_large.mass / p_small.mass
    drag_ratio = p_large.stokes_drag_coefficient(mu) / p_small.stokes_drag_coefficient(mu)
    beta_ratio = p_small.stokes_drag_rate(mu) / p_large.stokes_drag_rate(mu)
    v_term_small = p_small.terminal_buoyancy_rise_velocity(rho_f, mu)
    v_term_large = p_large.terminal_buoyancy_rise_velocity(rho_f, mu)
    print(
        f"  [DEBUG] Ratios (50 µm / 3 µm): mass={mass_ratio:.1f}×  "
        f"drag_coeff={drag_ratio:.1f}×  (mass should ≈ (50/3)³ ≈ 4630×)"
    )
    print(
        f"  [DEBUG] Drag coupling β (3/50 µm)={beta_ratio:.0f}× — small particles track fluid faster"
    )
    print(
        f"  [DEBUG] Terminal wall drift v_y: 3 µm={v_term_small*1e6:.2f} µm/s  "
        f"50 µm={v_term_large*1e6:.1f} µm/s  (large particles drift to walls faster)"
    )


def verify_particle_size_scaling(diameters_um: Tuple[float, ...] = EXPERIMENTAL_DIAMETERS_UM) -> None:
    """Print mass/drag/terminal-velocity ratios to confirm size-dependent physics."""
    env = TubeEnvironment()
    ref = diameters_um[0]
    p_ref = Particle.from_diameter_um(ref, np.array([0.0, env.centerline_y]))
    print("Particle size scaling check (polypropylene in saline):")
    print(f"  {'d (µm)':>8} {'mass (kg)':>14} {'drag coeff':>14} {'v_rise (µm/s)':>14}")
    for d in diameters_um:
        p = Particle.from_diameter_um(d, np.array([0.0, env.centerline_y]))
        drag_c = p.stokes_drag_coefficient(env.fluid.dynamic_viscosity)
        v_rise = p.terminal_buoyancy_rise_velocity(env.fluid.density, env.fluid.dynamic_viscosity)
        print(f"  {d:8.1f} {p.mass:14.3e} {drag_c:14.3e} {v_rise * 1e6:14.1f}")
    p_large = Particle.from_diameter_um(diameters_um[-1], np.array([0.0, env.centerline_y]))
    mass_ratio = p_large.mass / p_ref.mass
    radius_ratio = p_large.radius / p_ref.radius
    print(
        f"  Mass ratio ({diameters_um[-1]:.0f}/{ref:.0f} µm): {mass_ratio:.1f}×  "
        f"(expected (r_large/r_small)³ = {radius_ratio**3:.1f}×)"
    )


class Simulation:
    """Orchestrates the time-stepping loop over many particles."""

    def __init__(
        self,
        environment: TubeEnvironment,
        particles: List[Particle],
        dt: float = 1.0e-5,
        max_time: float = 10.0,
    ) -> None:
        self.environment = environment
        self.particles = particles
        self.dt = dt
        self.max_time = max_time
        self.time = 0.0
        self.history: List[dict] = []

    def step(self) -> None:
        for particle in self.particles:
            particle.update(self.environment, self.dt)
        self.time += self.dt

    def run(self, record_interval: Optional[int] = None) -> None:
        """
        Run until all particles are terminal or max_time is reached.

        If record_interval is set (in steps), snapshot state every N steps.
        """
        step_count = 0
        while self.time < self.max_time:
            if all(p.state != ParticleState.FLOWING for p in self.particles):
                break

            self.step()
            step_count += 1

            if record_interval and step_count % record_interval == 0:
                self._record_snapshot()

        self._record_snapshot()

    def _record_snapshot(self) -> None:
        self.history.append(
            {
                "time": self.time,
                "particles": [
                    {
                        "id": p.id,
                        "position": p.position.copy(),
                        "velocity": p.velocity.copy(),
                        "state": p.state.value,
                    }
                    for p in self.particles
                ],
            }
        )

    def summary(self) -> dict:
        counts = {state.value: 0 for state in ParticleState}
        for p in self.particles:
            counts[p.state.value] += 1
        return {
            "simulation_time_s": self.time,
            "tube_length_mm": self.environment.length * 1e3,
            "tube_width_mm": self.environment.width * 1e3,
            "particle_counts": counts,
        }


def random_inlet_y(environment: TubeEnvironment, radius: float, rng: np.random.Generator) -> float:
    """
    Randomized vertical entry position at the tube inlet.

    Uniform over the interior channel so the full sphere fits between walls.
    """
    margin = radius
    return float(rng.uniform(margin, environment.width - margin))


def create_inlet_particle(
    environment: TubeEnvironment,
    diameter_um: float,
    rng: np.random.Generator,
    trial_id: int = 0,
) -> Particle:
    """
    Instantiate one particle at the tube inlet (x = 0 plane).

    The particle center is placed at x = radius so the leading hemisphere lies
    at the inlet plane while the sphere remains inside the tube.
    """
    diameter_um = normalize_diameter_um(diameter_um)
    radius = diameter_um_to_radius_m(diameter_um)
    y = random_inlet_y(environment, radius, rng)
    particle = Particle.from_diameter_um(
        diameter_um,
        position=np.array([radius, y]),
        trial_id=trial_id,
    )
    if normalize_diameter_um(particle.diameter_um) != diameter_um:
        raise RuntimeError(
            f"Particle diameter mismatch: requested {diameter_um} µm, "
            f"got {particle.diameter_um} µm"
        )
    return particle


@dataclass
class TrajectoryLog:
    """Continuous dependent variable: particle (x, y) at every Euler timestep."""

    time_s: np.ndarray
    x_m: np.ndarray
    y_m: np.ndarray

    def __post_init__(self) -> None:
        self.time_s = np.asarray(self.time_s, dtype=float)
        self.x_m = np.asarray(self.x_m, dtype=float)
        self.y_m = np.asarray(self.y_m, dtype=float)

    @property
    def n_steps(self) -> int:
        return len(self.time_s)

    @property
    def x_mm(self) -> np.ndarray:
        return self.x_m * 1e3

    @property
    def y_mm(self) -> np.ndarray:
        return self.y_m * 1e3

    @classmethod
    def from_samples(
        cls,
        times: List[float],
        x_coords: List[float],
        y_coords: List[float],
        max_points: int = MAX_TRAJECTORY_POINTS,
    ) -> TrajectoryLog:
        n = len(times)
        if n > max_points:
            indices = np.linspace(0, n - 1, max_points, dtype=int)
            indices[-1] = n - 1
            times = [times[i] for i in indices]
            x_coords = [x_coords[i] for i in indices]
            y_coords = [y_coords[i] for i in indices]
        return cls(time_s=np.array(times), x_m=np.array(x_coords), y_m=np.array(y_coords))


@dataclass
class DiscreteOutcomeTally:
    """Running discrete counts of terminal outcomes for one particle-size group."""

    diameter_um: float
    escaped: int = 0
    stuck: int = 0
    incomplete: int = 0

    @property
    def total(self) -> int:
        return self.escaped + self.stuck + self.incomplete

    def record(self, state: ParticleState) -> None:
        if state == ParticleState.ESCAPED:
            self.escaped += 1
        elif state == ParticleState.STUCK:
            self.stuck += 1
        else:
            self.incomplete += 1

    @property
    def escaped_percent(self) -> float:
        return 100.0 * self.escaped / self.total if self.total else 0.0

    @property
    def stuck_percent(self) -> float:
        return 100.0 * self.stuck / self.total if self.total else 0.0

    @property
    def incomplete_percent(self) -> float:
        return 100.0 * self.incomplete / self.total if self.total else 0.0

    def to_dict(self) -> dict:
        return {
            "diameter_um": self.diameter_um,
            "total_trials": self.total,
            "escaped_count": self.escaped,
            "stuck_count": self.stuck,
            "incomplete_count": self.incomplete,
            "escaped_percent": self.escaped_percent,
            "stuck_percent": self.stuck_percent,
            "incomplete_percent": self.incomplete_percent,
        }


@dataclass
class ExperimentDataLog:
    """
    Central log for discrete outcomes and continuous trajectory data.

    Discrete: running escape/stuck tallies per size group.
    Continuous: full (t, x, y) path for every trial.
    """

    tallies: Dict[float, DiscreteOutcomeTally] = field(default_factory=dict)
    trajectories: Dict[float, List[TrajectoryLog]] = field(default_factory=dict)

    def register_group(self, diameter_um: float) -> None:
        diameter_um = normalize_diameter_um(diameter_um)
        if diameter_um not in self.tallies:
            self.tallies[diameter_um] = DiscreteOutcomeTally(diameter_um=diameter_um)
            self.trajectories[diameter_um] = []

    def reset_group(self, diameter_um: float) -> None:
        """Clear prior trial data for this size group before a fresh run."""
        diameter_um = normalize_diameter_um(diameter_um)
        self.tallies[diameter_um] = DiscreteOutcomeTally(diameter_um=diameter_um)
        self.trajectories[diameter_um] = []

    def reset_all(self) -> None:
        """Clear all logged trial data (call before run_all)."""
        self.tallies.clear()
        self.trajectories.clear()

    def record_trial(
        self,
        diameter_um: float,
        final_state: ParticleState,
        trajectory: TrajectoryLog,
    ) -> None:
        diameter_um = normalize_diameter_um(diameter_um)
        self.register_group(diameter_um)
        self.tallies[diameter_um].record(final_state)
        self.trajectories[diameter_um].append(trajectory)

    def get_tally(self, diameter_um: float) -> DiscreteOutcomeTally:
        diameter_um = normalize_diameter_um(diameter_um)
        self.register_group(diameter_um)
        return self.tallies[diameter_um]

    def total_trajectory_points(self) -> int:
        return sum(traj.n_steps for logs in self.trajectories.values() for traj in logs)

    def save_trajectories_npz(self, path: Path) -> None:
        """Persist trajectory arrays for offline plotting."""
        payload: dict = {}
        for diameter_um, logs in self.trajectories.items():
            key = f"{diameter_um:.0f}um"
            for trial_idx, traj in enumerate(logs):
                prefix = f"{key}_trial{trial_idx:03d}"
                payload[f"{prefix}_time_s"] = traj.time_s
                payload[f"{prefix}_x_m"] = traj.x_m
                payload[f"{prefix}_y_m"] = traj.y_m
        np.savez_compressed(path, **payload)

@dataclass
class TrialResult:
    """Outcome of a single independent trial."""

    trial_id: int
    diameter_um: float
    initial_y_m: float
    final_state: ParticleState
    transit_time_s: float
    final_position: np.ndarray
    final_velocity: np.ndarray
    trajectory: TrajectoryLog
    retention_mechanism: Optional[str] = None


@dataclass
class GroupResult:
    """Aggregated outcomes for one particle-size experimental group."""

    diameter_um: float
    trials: List[TrialResult]
    planned_trials: Optional[int] = None

    @property
    def n_trials(self) -> int:
        return len(self.trials)

    @property
    def escape_count(self) -> int:
        return sum(1 for t in self.trials if t.final_state == ParticleState.ESCAPED)

    @property
    def stuck_count(self) -> int:
        return sum(1 for t in self.trials if t.final_state == ParticleState.STUCK)

    @property
    def incomplete_count(self) -> int:
        """Trials still flowing when the time limit was reached."""
        return sum(1 for t in self.trials if t.final_state == ParticleState.FLOWING)

    @property
    def escape_fraction(self) -> float:
        return self.escape_count / self.n_trials if self.n_trials else 0.0

    @property
    def stuck_fraction(self) -> float:
        return self.stuck_count / self.n_trials if self.n_trials else 0.0

    def mean_transit_time_escaped_s(self) -> Optional[float]:
        escaped = [t.transit_time_s for t in self.trials if t.final_state == ParticleState.ESCAPED]
        return float(np.mean(escaped)) if escaped else None

    def summary(self) -> dict:
        return {
            "diameter_um": self.diameter_um,
            "n_trials": self.n_trials,
            "escaped": self.escape_count,
            "stuck": self.stuck_count,
            "incomplete": self.incomplete_count,
            "escape_fraction": self.escape_fraction,
            "stuck_fraction": self.stuck_fraction,
            "mean_transit_time_escaped_s": self.mean_transit_time_escaped_s(),
        }


def stokes_relaxation_time(
    radius: float,
    dynamic_viscosity: float,
    density: float = POLYPROPYLENE_DENSITY,
) -> float:
    """Momentum relaxation time τ = 2 ρ r² / (9 μ) for Stokes drag on a sphere."""
    return 2.0 * density * radius**2 / (9.0 * dynamic_viscosity)


def stokes_drag_rate_from_radius(
    radius: float,
    dynamic_viscosity: float,
    density: float = POLYPROPYLENE_DENSITY,
) -> float:
    """Drag coupling rate β = 9 μ / (2 ρ r²) from radius alone."""
    return 9.0 * dynamic_viscosity / (2.0 * density * radius**2)


def verify_parabolic_flow_profile(environment: TubeEnvironment) -> None:
    """Confirm Poiseuille profile: zero at walls, maximum at centerline."""
    environment.verify_parabolic_profile()
    print("\n  [DEBUG] Parabolic flow profile u(y):")
    print(f"    u_wall={environment.axial_velocity(0.0)*1e3:.2f} mm/s  "
          f"u_center={environment.axial_velocity(environment.centerline_y)*1e3:.2f} mm/s  "
          f"u_avg={environment.fluid.average_velocity*1e3:.2f} mm/s")
    for label, y in [("center", environment.centerline_y), ("mid-low", environment.width * 0.25),
                     ("mid-high", environment.width * 0.75)]:
        print(f"    y={label}: u={environment.axial_velocity(y)*1e3:.2f} mm/s")


def fast_mode_dt_cap(environment: TubeEnvironment, dt_max: float = DT_MAX) -> float:
    """
    Geometric fast_mode Δt ceiling, scaled so wall-clock ∝ campaign length not û.

    At û_ref = 50 mm/s the cap is FAST_MODE_DT_CAP_AT_U_REF (5×10⁻⁵ s). At slower
    mean speed the exact Stokes integrator stays stable, so the same clearance
    accuracy permits Δt × (û_ref / û).
    """
    base = min(dt_max * 5.0, FAST_MODE_DT_CAP_AT_U_REF)
    u = max(float(environment.fluid.average_velocity), 1e-12)
    return float(min(FAST_MODE_DT_CAP_MAX, base * (FAST_MODE_U_REF_M_S / u)))


def compute_euler_timestep(
    particle: Particle,
    environment: TubeEnvironment,
    *,
    fixed_dt: Optional[float] = None,
    relaxation_steps: int = EULER_RELAXATION_STEPS,
    cfl_fraction: float = EULER_CFL_FRACTION,
    wall_safety: float = EULER_WALL_SAFETY,
    dt_min: float = DT_MIN,
    dt_max: float = DT_MAX,
    fast_mode: bool = False,
    analytic_stokes: bool = True,
) -> float:
    """
    Adaptive timestep for the particle integrator.

    With analytic Stokes drag (default), stability no longer requires dt ≪ τ.
    Criteria (minimum):
      1. CFL displacement:   dt ≤ cfl_fraction × r / |v|
      2. Wall approach:      dt ≤ wall_safety × gap_to_wall / |v_y|
      3. Legacy Euler only:    dt ≤ τ / relaxation_steps

    fast_mode uses a uniform geometric cap scaled ∝ 1/û (exact Stokes).
    """
    cap = fast_mode_dt_cap(environment, dt_max) if fast_mode else dt_max

    if fixed_dt is not None:
        return float(np.clip(fixed_dt, dt_min, cap))

    mu = environment.fluid.dynamic_viscosity
    radius = particle.radius
    tau = stokes_relaxation_time(radius, mu, particle.density)

    particle_speed = float(np.linalg.norm(particle.velocity))
    fluid_speed = float(np.linalg.norm(environment.fluid_velocity_at(particle.position)))
    ref_speed = max(particle_speed, fluid_speed, 1e-12)
    dt_cfl = cfl_fraction * radius / ref_speed

    clearance_bottom, clearance_top = environment.wall_clearances(particle.position, radius)
    nearest_wall_gap = max(min(clearance_bottom, clearance_top), 1e-15)
    vertical_speed = max(abs(float(particle.velocity[1])), 1e-12)
    dt_wall = wall_safety * nearest_wall_gap / vertical_speed

    if analytic_stokes:
        if fast_mode:
            # Drag is integrated exactly — use a uniform step size for all diameters.
            # (Particle-radius CFL would shrink dt for 3 µm and starve them of steps.)
            return float(np.clip(min(dt_wall, cap), dt_min, cap))
        u_ref = max(particle_speed, fluid_speed, 1e-12)
        dt_flow = cfl_fraction * environment.length / (100.0 * u_ref)
        return float(np.clip(min(dt_flow, dt_wall), dt_min, cap))

    candidates = [dt_cfl, dt_wall, tau / relaxation_steps]
    return float(np.clip(min(candidates), dt_min, cap))


def stable_timestep(
    radius: float,
    dynamic_viscosity: float,
    density: float = POLYPROPYLENE_DENSITY,
) -> float:
    """Legacy helper — prefer compute_euler_timestep for full adaptive logic."""
    tau = stokes_relaxation_time(radius, dynamic_viscosity, density)
    return float(np.clip(tau / EULER_RELAXATION_STEPS, DT_MIN, DT_MAX))


def verify_wall_collision_scaling(environment: TubeEnvironment) -> None:
    """Confirm larger particles contact walls at lower center-y than smaller ones."""
    print("\n  [DEBUG] Wall contact thresholds (center-y at contact):")
    for d in EXPERIMENTAL_DIAMETERS_UM:
        r = diameter_um_to_radius_m(d)
        t = environment.wall_contact_thresholds(r)
        print(
            f"    {d:4.0f} µm (r={r*1e6:.1f} µm): bottom y≥{t['bottom_contact_y_m']*1e3:.4f} mm  "
            f"top y≤{t['top_contact_y_m']*1e3:.4f} mm  "
            f"interior band={t['interior_band_m']*1e3:.4f} mm"
        )
    r_small = diameter_um_to_radius_m(3.0)
    r_large = diameter_um_to_radius_m(50.0)
    t_small = environment.wall_contact_thresholds(r_small)
    t_large = environment.wall_contact_thresholds(r_large)
    assert t_large["top_contact_y_m"] < t_small["top_contact_y_m"], (
        "Wall scaling error: larger radius must contact top wall at lower center-y"
    )
    assert t_large["bottom_contact_y_m"] > t_small["bottom_contact_y_m"], (
        "Wall scaling error: larger radius must contact bottom wall at higher center-y"
    )
    print(
        f"  [DEBUG] 50 µm contacts top wall {((t_small['top_contact_y_m'] - t_large['top_contact_y_m'])*1e6):.1f} µm "
        f"sooner (lower y) than 3 µm"
    )

    # Same center-y trajectory: only the larger sphere should intersect first.
    shared_y = environment.width - 3.0e-6
    pos = np.array([0.05, shared_y])
    hit_small = environment.check_wall_collision(pos, r_small)
    hit_large = environment.check_wall_collision(pos, r_large)
    print(
        f"  [DEBUG] Shared center-y={shared_y*1e3:.4f} mm: "
        f"3 µm collision={hit_small}  50 µm collision={hit_large}"
    )
    if hit_large and not hit_small:
        print("  [DEBUG] ✓ Radius expands collision envelope — large particle clips wall first")
    assert hit_large or not hit_small, (
        "Wall scaling error: 3 µm must not collide before 50 µm on the same center-y path"
    )


def debug_print_timestep_scaling(environment: TubeEnvironment) -> None:
    """Print adaptive dt for reference particle sizes at inlet centerline."""
    print("\n  [DEBUG] Adaptive Euler dt at inlet centerline (v=0):")
    for d in EXPERIMENTAL_DIAMETERS_UM:
        p = Particle.from_diameter_um(d, np.array([diameter_um_to_radius_m(d), environment.centerline_y]))
        dt = compute_euler_timestep(p, environment)
        tau = stokes_relaxation_time(p.radius, environment.fluid.dynamic_viscosity, p.density)
        print(f"    {d:4.0f} µm: dt={dt:.3e} s  (τ={tau:.3e} s, dt/τ={dt/tau:.4f})")


def run_single_trial(
    environment: TubeEnvironment,
    diameter_um: float,
    trial_id: int,
    rng: np.random.Generator,
    dt: Optional[float] = None,
    max_time: float = DEFAULT_MAX_TIME_S,
    record_trajectory: bool = True,
    debug_first_trial: bool = False,
    fast_mode: bool = False,
    max_physics_steps: int = MAX_PHYSICS_STEPS,
    clinical: Optional[ClinicalRealismContext] = None,
    clinical_config: Optional[ClinicalRealismConfig] = None,
    container_key: str = "pristine",
) -> TrialResult:
    """
    Run one trial: single inlet particle until it escapes or sticks.

    Records (x, y) at the initial state and after every Euler timestep.
    Returns immediately when the particle reaches a terminal state.
    """
    diameter_um = normalize_diameter_um(diameter_um)
    particle = create_inlet_particle(environment, diameter_um, rng, trial_id=trial_id)

    cfg = clinical_config if clinical_config is not None else DEFAULT_CLINICAL_REALISM
    if clinical is None and cfg.enabled:
        # Paired layout RNG: independent of container wear (same defects for P and A)
        layout_seed = (
            int(diameter_um * 1000)
            + int(trial_id)
            + 777_001
        )
        layout_rng = np.random.default_rng(layout_seed)
        expected = cfg.defect_density_per_m * environment.length
        n_def = int(layout_rng.poisson(expected))
        shared_defects = (
            np.sort(layout_rng.uniform(0.0, environment.length, size=n_def))
            if n_def > 0
            else np.array([], dtype=float)
        )
        clinical = ClinicalRealismContext.for_trial(
            environment,
            rng,
            config=cfg,
            container_key=container_key,
            shared_defect_x_m=shared_defects if cfg.paired_defect_layout else None,
        )
    elif clinical is not None and not cfg.enabled:
        clinical = None

    if debug_first_trial:
        mu = environment.fluid.dynamic_viscosity
        test_velocity = np.array([0.05, 0.0])
        drag = particle.stokes_drag_coefficient(mu) * test_velocity
        print(
            f"    [DEBUG trial 0] group={diameter_um:.0f} µm  "
            f"particle.diameter_um={particle.diameter_um:.1f}  "
            f"mass={particle.mass:.4e} kg  "
            f"|F_drag|={np.linalg.norm(drag):.4e} N  "
            f"rng_sample_y={particle.position[1] * 1e3:.4f} mm"
        )

    initial_y = float(particle.position[1])
    fixed_dt = dt  # None → fully adaptive per step
    step_dt = compute_euler_timestep(particle, environment, fixed_dt=fixed_dt, fast_mode=fast_mode)

    times: List[float] = [0.0]
    x_coords: List[float] = [float(particle.position[0])]
    y_coords: List[float] = [float(particle.position[1])]
    time = 0.0

    if debug_first_trial:
        print(f"    [DEBUG trial 0] adaptive dt={step_dt:.3e} s  fast_mode={fast_mode}")

    est_steps = min(int(max_time / max(step_dt, DT_MIN)) + 2, max_physics_steps)
    record_stride = max(1, est_steps // MAX_TRAJECTORY_POINTS) if record_trajectory else 1
    if debug_first_trial and record_stride > 1:
        print(f"    [DEBUG trial 0] trajectory record_stride={record_stride}")

    step_index = 0
    while (
        time < max_time
        and particle.state == ParticleState.FLOWING
        and step_index < max_physics_steps
    ):
        step_dt = compute_euler_timestep(particle, environment, fixed_dt=fixed_dt, fast_mode=fast_mode)
        particle.update(environment, step_dt, clinical=clinical)
        if clinical is not None and particle.state == ParticleState.FLOWING:
            mechanism = clinical.try_clinical_retention(particle, step_dt)
            if mechanism is not None:
                particle._stick_clinical(environment, mechanism)
        time += step_dt
        step_index += 1
        if record_trajectory and (step_index % record_stride == 0):
            times.append(time)
            x_coords.append(float(particle.position[0]))
            y_coords.append(float(particle.position[1]))

    # Always retain the terminal point for outcome markers
    if record_trajectory and (not times or times[-1] != time):
        times.append(time)
        x_coords.append(float(particle.position[0]))
        y_coords.append(float(particle.position[1]))

    trajectory = TrajectoryLog.from_samples(times, x_coords, y_coords)

    return TrialResult(
        trial_id=trial_id,
        diameter_um=diameter_um,
        initial_y_m=initial_y,
        final_state=particle.state,
        transit_time_s=time,
        final_position=particle.position.copy(),
        final_velocity=particle.velocity.copy(),
        trajectory=trajectory,
        retention_mechanism=particle.retention_mechanism,
    )


class ExperimentRunner:
    """
    Execute repeated independent trials across experimental particle-size groups.

    Trial counts per diameter are set by the container profile (sterilization history).
    """

    def __init__(
        self,
        environment: TubeEnvironment,
        container_profile: ContainerProfile,
        diameters_um: Tuple[float, ...] = EXPERIMENTAL_DIAMETERS_UM,
        base_trials_per_size: int = TRIALS_PER_GROUP,
        dt: Optional[float] = None,
        max_time: float = DEFAULT_MAX_TIME_S,
        seed: int = 42,
        verbose: bool = True,
        debug: bool = False,
        fast_mode: bool = True,
        trajectory_trials_per_group: int = 2,
        max_physics_steps: int = MAX_PHYSICS_STEPS,
        data_log: Optional[ExperimentDataLog] = None,
        clinical_realism: bool = True,
        clinical_config: Optional[ClinicalRealismConfig] = None,
    ) -> None:
        self.environment = environment
        self.container_profile = container_profile
        self.diameters_um = diameters_um
        self.base_trials_per_size = base_trials_per_size
        self.trial_counts = container_profile.trial_counts(diameters_um, base_trials_per_size)
        self.injection_plan = container_profile.injection_summary(diameters_um, base_trials_per_size)
        self.dt = dt
        self.max_time = max_time
        self.seed = seed
        self.verbose = verbose
        self.debug = debug
        self.fast_mode = fast_mode
        self.trajectory_trials_per_group = trajectory_trials_per_group
        self.max_physics_steps = max_physics_steps
        self.data_log = data_log if data_log is not None else ExperimentDataLog()
        self.group_results: List[GroupResult] = []
        self.clinical_realism = clinical_realism
        if not clinical_realism:
            self.clinical_config = ClinicalRealismConfig(enabled=False)
        elif clinical_config is not None:
            self.clinical_config = clinical_config
        else:
            # Default: pristine w=1.0; autoclaved operating point w=1.75
            w = AUTOCLAVE_WALL_WEAR if container_profile.key == "autoclaved" else 1.0
            self.clinical_config = replace(
                DEFAULT_CLINICAL_REALISM,
                degraded_tubing_multiplier=w,
                enabled=True,
            )

    def run_group(
        self,
        diameter_um: float,
        n_trials: int,
        group_seed: int,
    ) -> GroupResult:
        diameter_um = normalize_diameter_um(diameter_um)
        self.data_log.reset_group(diameter_um)
        tally = self.data_log.get_tally(diameter_um)
        trials: List[TrialResult] = []

        if self.debug:
            debug_print_particle_physics(
                self.environment,
                diameter_um,
                label=f"group {diameter_um:.0f}µm pre-loop",
            )

        for trial_id in range(n_trials):
            trial_rng = np.random.default_rng(group_seed + trial_id)
            record_traj = trial_id < self.trajectory_trials_per_group
            result = run_single_trial(
                self.environment,
                diameter_um,
                trial_id,
                trial_rng,
                dt=self.dt,
                max_time=self.max_time,
                record_trajectory=record_traj,
                debug_first_trial=self.debug and trial_id == 0,
                fast_mode=self.fast_mode,
                max_physics_steps=self.max_physics_steps,
                clinical_config=self.clinical_config,
                container_key=self.container_profile.key,
            )
            if result.diameter_um != diameter_um:
                raise RuntimeError(
                    f"Trial result diameter mismatch: group={diameter_um}, result={result.diameter_um}"
                )
            trials.append(result)
            self.data_log.record_trial(diameter_um, result.final_state, result.trajectory)

            if self.verbose and n_trials >= 25 and (trial_id + 1) % max(1, n_trials // 4) == 0:
                print(
                    f"  {diameter_um:.0f} µm: {trial_id + 1}/{n_trials} trials  "
                    f"[escaped {tally.escaped_percent:.1f}%, stuck {tally.stuck_percent:.1f}%]"
                )

        return GroupResult(diameter_um=diameter_um, trials=trials, planned_trials=n_trials)

    def run_all(self) -> List[GroupResult]:
        self.group_results = []
        self.data_log.reset_all()

        if self.debug:
            debug_compare_reference_particles(self.environment)
            verify_parabolic_flow_profile(self.environment)
            verify_wall_collision_scaling(self.environment)
            debug_print_timestep_scaling(self.environment)

        if self.verbose:
            print()
            print("=" * 62)
            print(f"  CONTAINER: {self.container_profile.display_name}")
            print(f"  Profile: {self.container_profile.comparison_label}  "
                  f"(degradation ×{self.container_profile.degradation_factor:.1f})")
            print(f"  Total particles to inject: {self.injection_plan['total_particles']}")
            mode = "FAST (Colab-safe)" if self.fast_mode else "FULL precision"
            print(f"  Mode: {mode}  |  max_time={self.max_time}s  |  max_steps={self.max_physics_steps:,}")
            print("=" * 62)

        for group_index, diameter_um in enumerate(self.diameters_um):
            diameter_um = normalize_diameter_um(diameter_um)
            n_trials = self.trial_counts[diameter_um]
            if self.verbose:
                print(
                    f"\nRunning group: {diameter_um:.0f} µm diameter "
                    f"({n_trials} trials, injection-weighted)"
                )
            group_seed = (
                self.seed
                + int(diameter_um * 1000)
                + group_index * 10_000
            )
            # Note: degradation_factor is intentionally excluded so pristine and
            # autoclaved share the same trial RNG streams (paired design, Tier 1.2).
            if self.debug:
                print(f"  [DEBUG] group_seed={group_seed} (unique per size group; paired across containers)")
            self.group_results.append(self.run_group(diameter_um, n_trials, group_seed))
        return self.group_results

    def summary_by_group(self) -> List[dict]:
        return [group.summary() for group in self.group_results]


def print_experiment_report(
    group_results: List[GroupResult],
    data_log: Optional[ExperimentDataLog] = None,
    container_profile: Optional[ContainerProfile] = None,
    injection_plan: Optional[dict] = None,
) -> None:
    """Print final discrete-outcome summary for hypothesis evaluation."""
    total_trials = sum(g.n_trials for g in group_results)

    print()
    print("=" * 62)
    print("  EXPERIMENTAL RESULTS — DISCRETE OUTCOMES")
    print("=" * 62)
    print()

    if container_profile is not None:
        print(f"  Container type: {container_profile.display_name}")
        print(f"  Safety comparison: {container_profile.comparison_label}")
        print(f"  Degradation factor: ×{container_profile.degradation_factor:.1f}")
        if injection_plan:
            print(
                f"  Particles injected: {injection_plan['total_particles']} "
                f"(baseline would be {injection_plan['baseline_total']})"
            )
            print(
                f"  Sub-{SUB_DANGEROUS_SIZE_THRESHOLD_UM:.0f} µm injection share: "
                f"{injection_plan['sub_20um_fraction']:.1%} "
                f"({injection_plan['sub_20um_particles']} particles)"
            )
        print()

    print("  Dependent variables: escaped (reached tube exit) vs. stuck (wall)")
    print("  (Note: container injection share ≠ physics escape rate)")
    print(f"  Total trials completed: {total_trials}")
    print()
    print(f"  {'Size (µm)':<12} {'Trials':>8} {'Escaped':>10} {'Stuck':>10} {'Escaped %':>12} {'Stuck %':>10}")
    print("  " + "-" * 64)

    for group in group_results:
        tally = data_log.get_tally(group.diameter_um) if data_log else None
        if tally:
            escaped = tally.escaped
            stuck = tally.stuck
            incomplete = tally.incomplete
            escaped_pct = tally.escaped_percent
            stuck_pct = tally.stuck_percent
        else:
            escaped = group.escape_count
            stuck = group.stuck_count
            incomplete = group.incomplete_count
            escaped_pct = 100.0 * group.escape_fraction
            stuck_pct = 100.0 * group.stuck_fraction

        escaped_str = f"{escaped} ({escaped_pct:.1f}%)"
        stuck_str = f"{stuck} ({stuck_pct:.1f}%)"
        print(
            f"  {group.diameter_um:<12.0f} {group.n_trials:>8d} "
            f"{escaped_str:>10} {stuck_str:>10} {escaped_pct:>11.1f}% {stuck_pct:>9.1f}%"
        )
        if incomplete:
            print(f"             ({incomplete} trials still flowing at time limit)")

    print("  " + "-" * 64)
    print()

    if data_log is not None:
        print("  CONTINUOUS DATA LOGGED")
        print(f"    Trajectory records: {total_trials} trials")
        print(f"    Total Euler steps stored: {data_log.total_trajectory_points():,}")
        for diameter_um in sorted(data_log.trajectories):
            logs = data_log.trajectories[diameter_um]
            steps = sum(t.n_steps for t in logs)
            print(f"    {diameter_um:.0f} µm group: {len(logs)} paths, {steps:,} (t, x, y) samples")
        print()
        print("  Access trajectories via: runner.data_log.trajectories[diameter_um][trial_id]")
        print("  Coordinates in meters (x_m, y_m); use .x_mm / .y_mm for millimeters.")
        print()

    print("=" * 62)


# ---------------------------------------------------------------------------
# Visualization module
# ---------------------------------------------------------------------------

SIZE_COLORS: Dict[float, str] = {
    3.0: "#2ca02c",   # green
    10.0: "#ff7f0e",  # orange
    50.0: "#d62728",  # red
}

ESCAPE_MARKER = "D"  # diamond
STUCK_MARKER = "X"
FLOWING_MARKER = "o"


def _apply_plot_style() -> None:
    import matplotlib.pyplot as plt

    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except OSError:
        plt.style.use("ggplot")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
        }
    )


def _subsample_path(x: np.ndarray, y: np.ndarray, max_points: int = 3000) -> Tuple[np.ndarray, np.ndarray]:
    if len(x) <= max_points:
        return x, y
    indices = np.linspace(0, len(x) - 1, max_points, dtype=int)
    return x[indices], y[indices]


def _select_sample_trials(group: GroupResult, n_samples: int = 3) -> List[TrialResult]:
    """Pick a representative mix of escaped and stuck trials for plotting."""
    escaped = [t for t in group.trials if t.final_state == ParticleState.ESCAPED]
    stuck = [t for t in group.trials if t.final_state == ParticleState.STUCK]
    samples: List[TrialResult] = []

    if escaped:
        samples.append(escaped[0])
    if stuck:
        samples.append(stuck[0])

    for trial in group.trials:
        if trial not in samples:
            samples.append(trial)
        if len(samples) >= n_samples:
            break

    return samples[:n_samples]


def _draw_tube_outline(ax, length_mm: float, width_mm: float) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    tube = Rectangle(
        (0.0, 0.0),
        length_mm,
        width_mm,
        fill=False,
        edgecolor="black",
        linewidth=2.0,
        linestyle="-",
        zorder=1,
    )
    ax.add_patch(tube)
    ax.axhline(0.0, color="black", linewidth=1.5, zorder=1)
    ax.axhline(width_mm, color="black", linewidth=1.5, zorder=1)
    ax.axvline(0.0, color="black", linewidth=1.5, linestyle="--", alpha=0.5, zorder=1)
    ax.axvline(length_mm, color="black", linewidth=1.5, linestyle="--", alpha=0.5, zorder=1)


def plot_trajectory_figure(
    environment: TubeEnvironment,
    group_results: List[GroupResult],
    *,
    container_profile: Optional[ContainerProfile] = None,
    samples_per_group: int = 3,
    save_path: Optional[Path] = None,
    show: bool = False,
):
    """
    Figure 1: 2D tube cross-section with sample particle trajectories per size group.

    Colors: green (3 µm), orange (10 µm), red (50 µm).
    Terminal outcomes marked with diamonds (escaped) or X (stuck).
    """
    import matplotlib.pyplot as plt

    _apply_plot_style()
    length_mm = environment.length * 1e3
    width_mm = environment.width * 1e3

    fig, ax = plt.subplots(figsize=(12, 4))
    _draw_tube_outline(ax, length_mm, width_mm)

    for group in group_results:
        color = SIZE_COLORS.get(group.diameter_um, "gray")
        for trial in _select_sample_trials(group, n_samples=samples_per_group):
            x_mm, y_mm = _subsample_path(trial.trajectory.x_mm, trial.trajectory.y_mm)
            label = f"{group.diameter_um:.0f} µm (trial {trial.trial_id})"
            ax.plot(x_mm, y_mm, color=color, linewidth=1.4, alpha=0.85, label=label, zorder=2)

            end_x = trial.trajectory.x_mm[-1]
            end_y = trial.trajectory.y_mm[-1]
            if trial.final_state == ParticleState.ESCAPED:
                ax.scatter(
                    end_x, end_y,
                    marker=ESCAPE_MARKER,
                    s=70,
                    color=color,
                    edgecolors="black",
                    linewidths=0.6,
                    zorder=4,
                )
            elif trial.final_state == ParticleState.STUCK:
                ax.scatter(
                    end_x, end_y,
                    marker=STUCK_MARKER,
                    s=80,
                    color=color,
                    edgecolors="black",
                    linewidths=0.6,
                    zorder=4,
                )
            else:
                ax.scatter(
                    end_x, end_y,
                    marker=FLOWING_MARKER,
                    s=50,
                    color=color,
                    edgecolors="black",
                    linewidths=0.6,
                    zorder=4,
                )

    ax.scatter([], [], marker=ESCAPE_MARKER, s=70, color="gray", edgecolors="black", label="Escaped")
    ax.scatter([], [], marker=STUCK_MARKER, s=80, color="gray", edgecolors="black", label="Stuck on wall")

    ax.set_xlim(-2.0, length_mm + 5.0)
    ax.set_ylim(-0.3, width_mm + 0.3)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Axial position (mm)")
    ax.set_ylabel("Lateral position (mm)")
    title = "Microplastic Particle Trajectories in Medical Tubing (2D Cross-Section)"
    if container_profile is not None:
        title += f"\nContainer: {container_profile.display_name} [{container_profile.comparison_label}]"
    ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0, framealpha=0.95)
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def plot_outcome_bar_chart(
    group_results: List[GroupResult],
    data_log: Optional[ExperimentDataLog] = None,
    *,
    container_profile: Optional[ContainerProfile] = None,
    save_path: Optional[Path] = None,
    show: bool = False,
):
    """
    Figure 2: Grouped bar chart of escape rate vs. wall retention (stuck) by particle size.
    """
    import matplotlib.pyplot as plt

    _apply_plot_style()

    diameters = [g.diameter_um for g in group_results]
    x_labels = [f"{d:.0f} µm" for d in diameters]
    x_pos = np.arange(len(diameters))
    bar_width = 0.35

    escape_pcts: List[float] = []
    stuck_pcts: List[float] = []

    for group in group_results:
        if data_log is not None:
            tally = data_log.get_tally(group.diameter_um)
            escape_pcts.append(tally.escaped_percent)
            stuck_pcts.append(tally.stuck_percent)
        else:
            escape_pcts.append(100.0 * group.escape_fraction)
            stuck_pcts.append(100.0 * group.stuck_fraction)

    fig, ax = plt.subplots(figsize=(8, 5))
    bars_escape = ax.bar(
        x_pos - bar_width / 2,
        escape_pcts,
        bar_width,
        label="Escape rate",
        color="#4c72b0",
        edgecolor="black",
        linewidth=0.6,
    )
    bars_stuck = ax.bar(
        x_pos + bar_width / 2,
        stuck_pcts,
        bar_width,
        label="Wall retention (stuck)",
        color="#c44e52",
        edgecolor="black",
        linewidth=0.6,
    )

    ax.bar_label(bars_escape, fmt="%.1f%%", padding=3, fontsize=9)
    ax.bar_label(bars_stuck, fmt="%.1f%%", padding=3, fontsize=9)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel("Particle diameter (µm)")
    ax.set_ylabel("Percentage of trials (%)")
    title = "Escape Rate vs. Wall Retention by Particle Size"
    if container_profile is not None:
        title += f"\nContainer: {container_profile.comparison_label} — {container_profile.display_name}"
    ax.set_title(title)
    ax.set_ylim(0, 105)
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(axis="y", linestyle="--", alpha=0.6)
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def plot_injection_profile_chart(
    container_profile: ContainerProfile,
    injection_plan: dict,
    *,
    save_path: Optional[Path] = None,
    show: bool = False,
):
    """Bar chart of planned particle injections by size for the selected container."""
    import matplotlib.pyplot as plt

    _apply_plot_style()

    diameters = sorted(injection_plan["trial_counts"].keys())
    counts = [injection_plan["trial_counts"][d] for d in diameters]
    colors = [SIZE_COLORS.get(d, "gray") for d in diameters]
    labels = [f"{d:.0f} µm" for d in diameters]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, counts, color=colors, edgecolor="black", linewidth=0.6)
    ax.bar_label(bars, fmt="%d", padding=3, fontsize=9)

    sub_frac = injection_plan["sub_20um_fraction"]
    ax.set_xlabel("Particle diameter (µm)")
    ax.set_ylabel("Injected particles (count)")
    ax.set_title(
        "Container Release Profile by Particle Size\n"
        f"{container_profile.display_name} [{container_profile.comparison_label}]  "
        f"(degradation ×{container_profile.degradation_factor:.1f}, "
        f"sub-{SUB_DANGEROUS_SIZE_THRESHOLD_UM:.0f} µm: {sub_frac:.1%})"
    )
    ax.grid(axis="y", linestyle="--", alpha=0.6)
    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def generate_experiment_figures(
    environment: TubeEnvironment,
    group_results: List[GroupResult],
    data_log: Optional[ExperimentDataLog] = None,
    container_profile: Optional[ContainerProfile] = None,
    injection_plan: Optional[dict] = None,
    output_dir: Optional[Path] = None,
    *,
    samples_per_group: int = 3,
    show: bool = False,
) -> Tuple[Path, Path, Optional[Path]]:
    """Generate and save experiment figures. Returns paths to saved files."""
    if output_dir is None:
        output_dir = Path(__file__).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"_{container_profile.key}" if container_profile else ""
    trajectory_path = output_dir / f"trajectory_plot{suffix}.png"
    bar_chart_path = output_dir / f"outcome_bar_chart{suffix}.png"
    injection_path: Optional[Path] = None

    plot_trajectory_figure(
        environment,
        group_results,
        container_profile=container_profile,
        samples_per_group=samples_per_group,
        save_path=trajectory_path,
        show=show,
    )
    plot_outcome_bar_chart(
        group_results,
        data_log=data_log,
        container_profile=container_profile,
        save_path=bar_chart_path,
        show=show,
    )

    if container_profile is not None and injection_plan is not None:
        injection_path = output_dir / f"injection_profile{suffix}.png"
        plot_injection_profile_chart(
            container_profile,
            injection_plan,
            save_path=injection_path,
            show=show,
        )

    print()
    print("  FIGURES GENERATED")
    print(f"    Trajectory plot:   {trajectory_path}")
    print(f"    Outcome bar chart: {bar_chart_path}")
    if injection_path is not None:
        print(f"    Injection profile: {injection_path}")

    return trajectory_path, bar_chart_path, injection_path


def main() -> None:
    verify_particle_size_scaling()

    container = resolve_container_profile()

    env = TubeEnvironment(
        length=0.100,
        width=0.004,
        fluid=FluidProperties(
            density=SALINE_DENSITY,
            dynamic_viscosity=SALINE_DYNAMIC_VISCOSITY,
            average_velocity=0.05,
        ),
    )

    runner = ExperimentRunner(
        environment=env,
        container_profile=container,
        diameters_um=EXPERIMENTAL_DIAMETERS_UM,
        base_trials_per_size=COLAB_TRIALS_PER_GROUP,
        max_time=DEFAULT_MAX_TIME_S,
        fast_mode=True,
        debug=False,
        seed=42,
    )
    group_results = runner.run_all()
    print_experiment_report(
        group_results,
        data_log=runner.data_log,
        container_profile=container,
        injection_plan=runner.injection_plan,
    )

    output_path = Path(__file__).parent / f"trajectory_data_{container.key}.npz"
    runner.data_log.save_trajectories_npz(output_path)
    print(f"  Trajectory arrays saved to: {output_path}")

    generate_experiment_figures(
        env,
        group_results,
        data_log=runner.data_log,
        container_profile=container,
        injection_plan=runner.injection_plan,
        output_dir=Path(__file__).parent,
    )
    print("=" * 62)


if __name__ == "__main__":
    main()
