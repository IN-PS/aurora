"""Mecânica orbital analítica do problema de dois corpos (órbitas circulares e cônicas).

Fórmulas (SI; μ = G·M do corpo central)::

    v_c   = √(μ/r)                      velocidade orbital circular
    T     = 2π·√(a³/μ)                  período (3ª lei de Kepler; a = r na órbita circular)
    r     = ∛(μ·T²/(4π²))               raio da órbita circular de período T
    ε     = v²/2 − μ/r                  energia orbital específica
    v     = √(μ·(2/r − 1/a))            equação vis-viva
    v_esc = √(2μ/r)                     velocidade de escape
    a     = −μ/(2ε)                     semieixo maior a partir da energia
    r_p   = a·(1 − e),  r_a = a·(1 + e) raios de periapse e apoapse

Convenções: ``a > 0`` para órbitas elípticas e ``a < 0`` para hiperbólicas. Grandezas
ilimitadas de trajetórias abertas (apoapse, período) são representadas por ``math.inf``.

Estes resultados são analíticos e idealizados; não constituem planos de voo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

import numpy as np

from ..core.formula import FormulaResult, get_formula, register_formula
from ..core.validation import as_non_negative, as_positive, as_real, as_vector
from ..data.bodies import CelestialBody
from ..errors import DomainError, InvalidParameterError
from . import hypotheses as H

_TWO_PI = 2.0 * math.pi
_ORBIT_HYPOTHESES = (H.TWO_BODY, H.NO_PERTURBATIONS, H.CLASSICAL, H.MU_GIVEN, H.RADIUS_FROM_CENTER)
_CIRCULAR = "Órbita perfeitamente circular (excentricidade nula)."
_MU = ("Parâmetro gravitacional do corpo central (G*M)", "m^3/s^2")
_R = ("Distância ao centro do corpo central", "m")


@register_formula(
    key="velocidade_orbital_circular",
    title="Velocidade orbital circular",
    expression="v = sqrt(mu/r)",
    unit="m/s",
    parameters={"mu": _MU, "r": ("Raio da órbita circular", "m")},
    hypotheses=(_CIRCULAR, *_ORBIT_HYPOTHESES),
)
def circular_orbital_speed(mu: float, r: float) -> float:
    """Velocidade de uma órbita circular de raio ``r``, em m/s."""
    mu, r = as_positive("mu", mu), as_positive("r", r)
    return math.sqrt(mu / r)


@register_formula(
    key="periodo_orbital_eliptico",
    title="Período orbital (3ª lei de Kepler)",
    expression="T = 2*pi*sqrt(a^3/mu)",
    unit="s",
    parameters={"mu": _MU, "a": ("Semieixo maior da elipse (positivo)", "m")},
    hypotheses=("Órbita fechada (elíptica ou circular).", *_ORBIT_HYPOTHESES),
)
def orbital_period(mu: float, a: float) -> float:
    """Período de uma órbita elíptica de semieixo maior ``a``, em s."""
    mu, a = as_positive("mu", mu), as_positive("a", a)
    return _TWO_PI * a * math.sqrt(a / mu)  # = 2π√(a³/μ), sem calcular a³ (evita overflow)


@register_formula(
    key="periodo_orbital_circular",
    title="Período orbital circular",
    expression="T = 2*pi*sqrt(r^3/mu)",
    unit="s",
    parameters={"mu": _MU, "r": ("Raio da órbita circular", "m")},
    hypotheses=(_CIRCULAR, *_ORBIT_HYPOTHESES),
)
def circular_orbital_period(mu: float, r: float) -> float:
    """Período de uma órbita circular de raio ``r``, em s."""
    mu, r = as_positive("mu", mu), as_positive("r", r)
    return _TWO_PI * r * math.sqrt(r / mu)


@register_formula(
    key="raio_orbita_circular",
    title="Raio da órbita circular com período dado",
    expression="r = (mu*T^2/(4*pi^2))^(1/3)",
    unit="m",
    parameters={"mu": _MU, "T": ("Período orbital desejado", "s")},
    hypotheses=(_CIRCULAR, *_ORBIT_HYPOTHESES),
)
def circular_orbit_radius(mu: float, T: float) -> float:  # noqa: N803 (símbolo da fórmula)
    """Raio (m) da órbita circular cujo período é ``T``."""
    mu, T = as_positive("mu", mu), as_positive("T", T)
    return mu ** (1.0 / 3.0) * (T / _TWO_PI) ** (2.0 / 3.0)


@register_formula(
    key="energia_especifica",
    title="Energia orbital específica (dois corpos)",
    expression="epsilon = v^2/2 - mu/r",
    unit="J/kg",
    parameters={"mu": _MU, "r": _R, "v": ("Módulo da velocidade relativa ao corpo central", "m/s")},
    hypotheses=(H.TWO_BODY, H.POTENTIAL_ZERO_AT_INFINITY, H.CLASSICAL, H.MU_GIVEN, H.RADIUS_FROM_CENTER),
)
def specific_orbital_energy(mu: float, r: float, v: float) -> float:
    """Energia mecânica por unidade de massa do veículo, em J/kg (= m²/s²).

    ε < 0: órbita fechada (elíptica); ε = 0: parabólica; ε > 0: hiperbólica.
    """
    mu, r, v = as_positive("mu", mu), as_positive("r", r), as_non_negative("v", v)
    return 0.5 * v * v - mu / r


@register_formula(
    key="velocidade_vis_viva",
    title="Velocidade pela equação vis-viva",
    expression="v = sqrt(mu*(2/r - 1/a))",
    unit="m/s",
    parameters={
        "mu": _MU,
        "r": _R,
        "a": ("Semieixo maior (positivo: elipse; negativo: hipérbole)", "m"),
    },
    hypotheses=_ORBIT_HYPOTHESES,
)
def vis_viva_speed(mu: float, r: float, a: float) -> float:
    """Velocidade (m/s) na distância ``r`` de uma órbita de semieixo maior ``a``."""
    mu, r = as_positive("mu", mu), as_positive("r", r)
    a = as_real("a", a)
    if a == 0.0:
        raise InvalidParameterError("O semieixo maior 'a' não pode ser zero.")
    term = 2.0 / r - 1.0 / a
    if term < 0.0:
        raise DomainError(
            "r > 2a: essa distância não é alcançada por uma órbita elíptica com esse semieixo maior."
        )
    return math.sqrt(mu * term)


@register_formula(
    key="velocidade_escape",
    title="Velocidade de escape",
    expression="v_esc = sqrt(2*mu/r)",
    unit="m/s",
    parameters={"mu": _MU, "r": _R},
    hypotheses=(*_ORBIT_HYPOTHESES, H.POTENTIAL_ZERO_AT_INFINITY),
)
def escape_speed(mu: float, r: float) -> float:
    """Menor velocidade (m/s) que leva a uma trajetória aberta (ε = 0) a partir da distância ``r``."""
    mu, r = as_positive("mu", mu), as_positive("r", r)
    return math.sqrt(2.0 * mu / r)


@register_formula(
    key="semieixo_maior_por_energia",
    title="Semieixo maior a partir da energia específica",
    expression="a = -mu/(2*epsilon)",
    unit="m",
    parameters={"mu": _MU, "epsilon": ("Energia orbital específica (não nula)", "J/kg")},
    hypotheses=(*_ORBIT_HYPOTHESES, "Resultado negativo indica trajetória hiperbólica (ε > 0)."),
)
def semi_major_axis_from_energy(mu: float, epsilon: float) -> float:
    """Semieixo maior (m); positivo para ε < 0 (elipse) e negativo para ε > 0 (hipérbole)."""
    mu, epsilon = as_positive("mu", mu), as_real("epsilon", epsilon)
    if epsilon == 0.0:
        raise DomainError("ε = 0 corresponde a uma trajetória parabólica: semieixo maior infinito.")
    return -mu / (2.0 * epsilon)


def _ellipse_params(a: float, e: float) -> tuple[float, float]:
    a = as_positive("a", a)
    e = as_non_negative("e", e)
    if e >= 1.0:
        raise DomainError("Para órbitas elípticas a excentricidade deve satisfazer 0 ≤ e < 1.")
    return a, e


@register_formula(
    key="raio_periapse",
    title="Raio de periapse (menor distância ao centro)",
    expression="r_p = a*(1 - e)",
    unit="m",
    parameters={"a": ("Semieixo maior", "m"), "e": ("Excentricidade (0 ≤ e < 1)", "1")},
    hypotheses=("Órbita elíptica ou circular.", *_ORBIT_HYPOTHESES),
)
def periapsis_radius(a: float, e: float) -> float:
    """Raio de periapse (m) de uma órbita elíptica."""
    a, e = _ellipse_params(a, e)
    return a * (1.0 - e)


@register_formula(
    key="raio_apoapse",
    title="Raio de apoapse (maior distância ao centro)",
    expression="r_a = a*(1 + e)",
    unit="m",
    parameters={"a": ("Semieixo maior", "m"), "e": ("Excentricidade (0 ≤ e < 1)", "1")},
    hypotheses=("Órbita elíptica ou circular.", *_ORBIT_HYPOTHESES),
)
def apoapsis_radius(a: float, e: float) -> float:
    """Raio de apoapse (m) de uma órbita elíptica."""
    a, e = _ellipse_params(a, e)
    return a * (1.0 + e)


# --------------------------------------------------------------------------------------
# Elementos orbitais a partir do vetor de estado
# --------------------------------------------------------------------------------------


class OrbitType(Enum):
    """Classificação da trajetória no problema de dois corpos."""

    CIRCULAR = "circular"
    ELLIPTIC = "elíptica"
    PARABOLIC = "parabólica"
    HYPERBOLIC = "hiperbólica"
    RADIAL = "radial (degenerada: momento angular nulo)"


@dataclass(frozen=True)
class TwoBodyOrbit:
    """Grandezas da cônica definida por um vetor de estado (todas em SI)."""

    orbit_type: OrbitType
    specific_energy: float
    """ε em J/kg."""
    angular_momentum: np.ndarray
    """h⃗ = r⃗ × v⃗ em m²/s."""
    eccentricity_vector: np.ndarray
    eccentricity: float
    semi_major_axis: float
    """m; ``inf`` para a parábola; negativo para a hipérbole."""
    semi_latus_rectum: float
    """p = h²/μ em m."""
    periapsis_radius: float
    apoapsis_radius: float
    """m; ``inf`` para trajetórias abertas."""
    period: float
    """s; ``inf`` para trajetórias abertas."""
    inclination: float
    """rad em relação ao plano xy do referencial; ``nan`` se h⃗ = 0."""


def orbit_from_state(mu: float, position: object, velocity: object, tol: float = 1e-10) -> TwoBodyOrbit:
    """Calcula a cônica de dois corpos a partir de posição (m) e velocidade (m/s) relativas.

    A classificação usa a excentricidade com tolerância ``tol``: ``e < tol`` é circular e
    ``|e − 1| ≤ tol`` é parabólica. Trajetórias com momento angular desprezível
    (``|h| ≤ tol·|r|·|v|``) são classificadas como radiais.
    """
    mu = as_positive("mu", mu)
    tol = as_positive("tol", tol)
    r_vec = as_vector("position", position)
    v_vec = as_vector("velocity", velocity)
    r = float(np.linalg.norm(r_vec))
    if r == 0.0:
        raise DomainError("A posição não pode coincidir com o centro do corpo central.")
    v2 = float(np.dot(v_vec, v_vec))
    v = math.sqrt(v2)

    energy = 0.5 * v2 - mu / r
    h_vec = np.cross(r_vec, v_vec)
    h = float(np.linalg.norm(h_vec))
    e_vec = ((v2 - mu / r) * r_vec - float(np.dot(r_vec, v_vec)) * v_vec) / mu
    e = float(np.linalg.norm(e_vec))
    p = h * h / mu

    if h <= tol * r * v:
        kind = OrbitType.RADIAL
    elif e < tol:
        kind = OrbitType.CIRCULAR
    elif abs(e - 1.0) <= tol:
        kind = OrbitType.PARABOLIC
    elif e < 1.0:
        kind = OrbitType.ELLIPTIC
    else:
        kind = OrbitType.HYPERBOLIC

    if kind is OrbitType.PARABOLIC or energy == 0.0:
        a = math.inf
    else:
        a = -mu / (2.0 * energy)

    if kind is OrbitType.RADIAL:
        r_p = 0.0
        r_a = 2.0 * a if energy < 0.0 else math.inf
    else:
        r_p = p / (1.0 + e)
        r_a = p / (1.0 - e) if kind in (OrbitType.CIRCULAR, OrbitType.ELLIPTIC) else math.inf

    period = _TWO_PI * a * math.sqrt(a / mu) if energy < 0.0 and math.isfinite(a) else math.inf
    inclination = math.acos(max(-1.0, min(1.0, float(h_vec[2]) / h))) if h > 0.0 else math.nan

    return TwoBodyOrbit(
        orbit_type=kind,
        specific_energy=energy,
        angular_momentum=h_vec,
        eccentricity_vector=e_vec,
        eccentricity=e,
        semi_major_axis=a,
        semi_latus_rectum=p,
        periapsis_radius=r_p,
        apoapsis_radius=r_a,
        period=period,
        inclination=inclination,
    )


# --------------------------------------------------------------------------------------
# Resumo de órbita circular em torno de um corpo de referência
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class CircularOrbitSummary:
    """Grandezas de uma órbita circular em torno de um corpo celeste de referência."""

    body: CelestialBody
    altitude: float
    """m acima do raio de referência do corpo."""
    radius: float
    """m a partir do centro do corpo."""
    results: tuple[FormulaResult, ...]


def radius_from_altitude(body: CelestialBody, altitude: float) -> float:
    """Raio orbital ``r = R + h`` (m), com ``h`` medida a partir do raio de referência do corpo."""
    return body.radius + as_real("altitude", altitude)


def circular_orbit_summary(body: CelestialBody, radius: float) -> CircularOrbitSummary:
    """Velocidade, período, energia, velocidade de escape e gravidade local de uma órbita circular.

    ``radius`` (m) é medido a partir do centro do corpo, modelado como esfera de raio
    ``body.radius``. Raios menores que o do corpo são rejeitados, pois a órbita o
    intersectaria.
    """
    r = as_positive("radius", radius)
    if r < body.radius:
        raise DomainError(
            f"Raio orbital ({r:.10g} m) menor que o {body.radius_description} do corpo {body.name} "
            f"({body.radius:.10g} m): a órbita intersectaria o corpo (altitude negativa)."
        )
    altitude = r - body.radius
    speed = get_formula("velocidade_orbital_circular").evaluate(mu=body.mu, r=r)
    results = (
        speed,
        get_formula("periodo_orbital_circular").evaluate(mu=body.mu, r=r),
        get_formula("energia_especifica").evaluate(mu=body.mu, r=r, v=speed.value),
        get_formula("velocidade_escape").evaluate(mu=body.mu, r=r),
        get_formula("aceleracao_gravitacional").evaluate(mu=body.mu, r=r),
    )
    return CircularOrbitSummary(body=body, altitude=altitude, radius=r, results=results)
