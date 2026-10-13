"""Validação do propagador numérico contra a solução analítica de Kepler.

Para o modelo de dois corpos existe solução exata (:func:`aurora.physics.kepler.propagate_kepler`).
Comparando cada amostra da trajetória integrada com essa solução, medimos o erro **global**
do integrador; a conservação da energia específica e do momento angular é uma verificação
adicional, independente da solução analítica.

Os cenários padrão cobrem órbitas circulares baixa e geoestacionária, elipses muito
excêntricas, uma órbita lunar e uma trajetória hiperbólica. Os parâmetros são ilustrativos
(não descrevem satélites reais específicos).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..core.validation import as_positive
from ..data.bodies import EARTH, MOON, CelestialBody
from ..errors import InvalidParameterError
from ..physics.gravitation import TwoBodyGravity
from ..physics.kepler import propagate_kepler, state_from_elements
from ..physics.orbits import circular_orbital_speed, escape_speed, orbital_period
from .propagator import IntegratorSettings, Trajectory, propagate


@dataclass(frozen=True)
class KeplerComparison:
    """Diferença entre a trajetória integrada e a solução analítica (em SI)."""

    max_position_error: float
    max_velocity_error: float
    final_position_error: float
    final_velocity_error: float


@dataclass(frozen=True)
class ConservationErrors:
    """Maiores variações relativas das constantes do movimento ao longo da trajetória."""

    energy: float
    """max |ε(t) − ε(0)| / |ε(0)| (ou absoluta, em J/kg, se ε(0) = 0)."""
    angular_momentum: float
    """max |h(t) − h(0)| / |h(0)|  (ou absoluta, em m²/s, se h(0) = 0)."""


def kepler_error_series(trajectory: Trajectory, mu: float) -> tuple[np.ndarray, np.ndarray]:
    """Erros |Δr| (m) e |Δv| (m/s) de cada amostra em relação à solução de Kepler.

    A referência parte do estado inicial da trajetória (o mesmo modelo de dois corpos), de
    modo que a diferença mede apenas o erro numérico da integração.
    """
    mu = as_positive("mu", mu)
    r0, v0 = trajectory.positions[0], trajectory.velocities[0]
    t0 = trajectory.times[0]
    dr = np.empty(trajectory.times.size)
    dv = np.empty(trajectory.times.size)
    for k, t in enumerate(trajectory.times):
        r_ref, v_ref = propagate_kepler(mu, r0, v0, float(t - t0))
        dr[k] = math.hypot(*(trajectory.positions[k] - r_ref))
        dv[k] = math.hypot(*(trajectory.velocities[k] - v_ref))
    return dr, dv


def compare_with_kepler(trajectory: Trajectory, mu: float) -> KeplerComparison:
    """Compara cada amostra com a solução de Kepler a partir do estado inicial da trajetória."""
    dr, dv = kepler_error_series(trajectory, mu)
    return KeplerComparison(
        max_position_error=float(dr.max()),
        max_velocity_error=float(dv.max()),
        final_position_error=float(dr[-1]),
        final_velocity_error=float(dv[-1]),
    )


def specific_energy_series(trajectory: Trajectory, mu: float) -> np.ndarray:
    """Energia específica ε = v²/2 − μ/r em cada amostra (J/kg)."""
    r = np.linalg.norm(trajectory.positions, axis=1)
    v2 = np.einsum("ij,ij->i", trajectory.velocities, trajectory.velocities)
    return 0.5 * v2 - mu / r


def angular_momentum_series(trajectory: Trajectory) -> np.ndarray:
    """Módulo do momento angular específico |r⃗ × v⃗| em cada amostra (m²/s)."""
    return np.linalg.norm(np.cross(trajectory.positions, trajectory.velocities), axis=1)


def _relative_drift(series: np.ndarray) -> float:
    reference = abs(float(series[0]))
    drift = float(np.max(np.abs(series - series[0])))
    return drift / reference if reference > 0.0 else drift


def conservation_errors(trajectory: Trajectory, mu: float) -> ConservationErrors:
    """Variação relativa máxima da energia e do momento angular (constantes no modelo de dois corpos)."""
    mu = as_positive("mu", mu)
    return ConservationErrors(
        energy=_relative_drift(specific_energy_series(trajectory, mu)),
        angular_momentum=_relative_drift(angular_momentum_series(trajectory)),
    )


@dataclass(frozen=True)
class ValidationCase:
    """Cenário de validação: estado inicial em torno de um corpo e duração (SI)."""

    key: str
    description: str
    body: CelestialBody
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    duration: float


@dataclass(frozen=True)
class ValidationResult:
    case: ValidationCase
    settings: IntegratorSettings
    comparison: KeplerComparison
    conservation: ConservationErrors
    evaluations: int
    samples: int


def _case_circular(key: str, description: str, body: CelestialBody, altitude: float, orbits: float) -> ValidationCase:
    r = body.radius + altitude
    v = circular_orbital_speed(body.mu, r)
    return ValidationCase(key, description, body, (r, 0.0, 0.0), (0.0, v, 0.0), orbits * orbital_period(body.mu, r))


def _case_ellipse(
    key: str, description: str, body: CelestialBody, h_peri: float, h_apo: float, inc_deg: float, orbits: float
) -> ValidationCase:
    rp, ra = body.radius + h_peri, body.radius + h_apo
    a, e = 0.5 * (rp + ra), (ra - rp) / (ra + rp)
    r, v = state_from_elements(body.mu, a, e, math.radians(inc_deg), 0.0, math.radians(270.0), 0.0)
    return ValidationCase(key, description, body, tuple(r), tuple(v), orbits * orbital_period(body.mu, a))


def standard_cases() -> tuple[ValidationCase, ...]:
    """Cenários padrão de validação (parâmetros ilustrativos)."""
    r_esc = EARTH.radius + 400e3
    v_hyp = 1.2 * escape_speed(EARTH.mu, r_esc)
    return (
        _case_circular("leo", "Terra: órbita circular de 400 km, 10 órbitas", EARTH, 400e3, 10),
        _case_circular("geo", "Terra: órbita geoestacionária (35 786 km), 2 órbitas", EARTH, 35_786e3, 2),
        _case_ellipse("gto", "Terra: elipse 200 × 35 786 km (e ≈ 0,73), 5 órbitas", EARTH, 200e3, 35_786e3, 27.0, 5),
        _case_ellipse(
            "molniya",
            "Terra: elipse tipo Molniya 500 × 39 800 km, i = 63,4°, 3 órbitas",
            EARTH,
            500e3,
            39_800e3,
            63.4,
            3,
        ),
        _case_circular("lua", "Lua: órbita circular de 100 km, 5 órbitas", MOON, 100e3, 5),
        ValidationCase(
            "hiperbolica",
            "Terra: hiperbólica, v = 1,2·v_escape a 400 km, 2 dias",
            EARTH,
            (r_esc, 0.0, 0.0),
            (0.0, v_hyp, 0.0),
            2 * 86_400.0,
        ),
    )


def run_case(case: ValidationCase, settings: IntegratorSettings | None = None, n_points: int = 400) -> ValidationResult:
    """Integra um cenário e mede o erro contra Kepler e a conservação das constantes."""
    settings = settings if settings is not None else IntegratorSettings()
    trajectory = propagate(
        TwoBodyGravity(case.body.mu), case.position, case.velocity, case.duration, settings=settings, n_points=n_points
    )
    return ValidationResult(
        case=case,
        settings=settings,
        comparison=compare_with_kepler(trajectory, case.body.mu),
        conservation=conservation_errors(trajectory, case.body.mu),
        evaluations=trajectory.evaluations,
        samples=trajectory.times.size,
    )


def run_validation(
    settings: IntegratorSettings | None = None, keys: tuple[str, ...] | None = None, n_points: int = 400
) -> list[ValidationResult]:
    """Executa os cenários padrão (ou só os de ``keys``)."""
    cases = standard_cases()
    if keys is not None:
        known = {c.key for c in cases}
        unknown = [k for k in keys if k not in known]
        if unknown:
            raise InvalidParameterError(
                f"Cenário(s) desconhecido(s): {', '.join(unknown)}. Disponíveis: {', '.join(sorted(known))}."
            )
        cases = tuple(c for c in cases if c.key in keys)
    return [run_case(case, settings, n_points) for case in cases]
