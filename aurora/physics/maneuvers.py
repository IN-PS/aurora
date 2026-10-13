"""Manobras impulsivas: equação de Tsiolkovsky e transferência de Hohmann.

Fórmulas (SI; g0 = 9,80665 m/s²)::

    Δv  = Isp·g0·ln(m0/mf)                          equação do foguete (Tsiolkovsky)
    m_p = m0·(1 − exp(−Δv/(Isp·g0)))                 propelente para um dado Δv
    a_t = (r1 + r2)/2                                semieixo maior da elipse de transferência
    Δv1 = √(μ(2/r1 − 1/a_t)) − √(μ/r1)               1ª queima (no raio r1)
    Δv2 = √(μ/r2) − √(μ(2/r2 − 1/a_t))               2ª queima (no raio r2)
    t_H = π·√(a_t³/μ)                                tempo de transferência (meia elipse)

A transferência de Hohmann é a transferência impulsiva de duas queimas entre órbitas
circulares e coplanares. Δv1 e Δv2 negativos indicam queimas retrógradas (redução de órbita).

Hipóteses: queimas impulsivas (instantâneas, sem perdas gravitacionais), Isp constante,
problema de dois corpos. São resultados idealizados, não planos de voo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..core.formula import register_formula
from ..core.validation import as_non_negative, as_positive
from ..data.constants import G0
from ..errors import DomainError
from . import hypotheses as H
from .orbits import _REFERENCES, circular_orbital_speed, vis_viva_speed

_IMPULSIVE = "Queima impulsiva: variação instantânea de velocidade, sem perdas gravitacionais nem de arrasto."
_CONSTANT_ISP = "Impulso específico constante durante a queima; g0 = 9,80665 m/s² (valor convencional)."
_HOHMANN = "Órbitas inicial e final circulares e coplanares; duas queimas tangenciais."
_ROCKET = "Sem forças externas durante a queima; a massa varia apenas pela ejeção de propelente."


@register_formula(
    key="delta_v_tsiolkovsky",
    title="Variação de velocidade (equação de Tsiolkovsky)",
    expression="dv = Isp*g0*ln(m0/mf)",
    unit="m/s",
    parameters={
        "isp": ("Impulso específico", "s"),
        "m0": ("Massa inicial (veículo + propelente)", "kg"),
        "mf": ("Massa final (após consumir o propelente)", "kg"),
    },
    hypotheses=(_IMPULSIVE, _CONSTANT_ISP, _ROCKET, H.CLASSICAL),
)
def rocket_delta_v(isp: float, m0: float, mf: float) -> float:
    """Δv (m/s) obtido ao queimar propelente de ``m0`` até ``mf`` com impulso específico ``isp`` (s)."""
    isp, m0, mf = as_positive("isp", isp), as_positive("m0", m0), as_positive("mf", mf)
    if mf > m0:
        raise DomainError("A massa final não pode ser maior que a massa inicial.")
    return isp * G0 * math.log(m0 / mf)


@register_formula(
    key="massa_propelente",
    title="Massa de propelente para um dado Δv (Tsiolkovsky)",
    expression="m_p = m0*(1 - exp(-dv/(Isp*g0)))",
    unit="kg",
    parameters={
        "isp": ("Impulso específico", "s"),
        "m0": ("Massa inicial (veículo + propelente)", "kg"),
        "dv": ("Variação de velocidade desejada (módulo)", "m/s"),
    },
    hypotheses=(_IMPULSIVE, _CONSTANT_ISP, _ROCKET, H.CLASSICAL),
)
def propellant_for_delta_v(isp: float, m0: float, dv: float) -> float:
    """Massa de propelente (kg) consumida para produzir ``dv`` (m/s) partindo da massa ``m0``."""
    isp, m0, dv = as_positive("isp", isp), as_positive("m0", m0), as_non_negative("dv", dv)
    return m0 * -math.expm1(-dv / (isp * G0))  # expm1: preciso também para Δv pequeno


@dataclass(frozen=True)
class HohmannTransfer:
    """Transferência de Hohmann entre órbitas circulares de raios ``r1`` e ``r2`` (SI)."""

    mu: float
    r1: float
    r2: float
    transfer_semi_major_axis: float
    speed_initial_circular: float
    speed_transfer_departure: float
    delta_v1: float
    speed_transfer_arrival: float
    speed_final_circular: float
    delta_v2: float
    delta_v_total: float
    transfer_time: float


def hohmann_transfer(mu: float, r1: float, r2: float) -> HohmannTransfer:
    """Calcula a transferência de Hohmann de ``r1`` para ``r2`` (raios a partir do centro, m)."""
    mu, r1, r2 = as_positive("mu", mu), as_positive("r1", r1), as_positive("r2", r2)
    a_t = 0.5 * (r1 + r2)
    v1 = circular_orbital_speed(mu, r1)
    v2 = circular_orbital_speed(mu, r2)
    vt1 = vis_viva_speed(mu, r1, a_t)
    vt2 = vis_viva_speed(mu, r2, a_t)
    dv1 = vt1 - v1
    dv2 = v2 - vt2
    return HohmannTransfer(
        mu=mu,
        r1=r1,
        r2=r2,
        transfer_semi_major_axis=a_t,
        speed_initial_circular=v1,
        speed_transfer_departure=vt1,
        delta_v1=dv1,
        speed_transfer_arrival=vt2,
        speed_final_circular=v2,
        delta_v2=dv2,
        delta_v_total=abs(dv1) + abs(dv2),
        transfer_time=math.pi * a_t * math.sqrt(a_t / mu),
    )


_MU = ("Parâmetro gravitacional do corpo central (G*M)", "m^3/s^2")
_R1 = ("Raio da órbita circular inicial (a partir do centro)", "m")
_R2 = ("Raio da órbita circular final (a partir do centro)", "m")
_HOHMANN_HYPOTHESES = (_HOHMANN, _IMPULSIVE, H.TWO_BODY, H.NO_PERTURBATIONS, H.CLASSICAL, H.MU_GIVEN)


@register_formula(
    key="delta_v_hohmann",
    title="Δv total da transferência de Hohmann",
    expression="dv = |dv1| + |dv2|",
    unit="m/s",
    parameters={"mu": _MU, "r1": _R1, "r2": _R2},
    hypotheses=_HOHMANN_HYPOTHESES,
    references=_REFERENCES,
)
def hohmann_delta_v(mu: float, r1: float, r2: float) -> float:
    """Δv total (m/s) da transferência de Hohmann entre órbitas circulares de raios r1 e r2."""
    return hohmann_transfer(mu, r1, r2).delta_v_total


@register_formula(
    key="tempo_hohmann",
    title="Tempo de transferência de Hohmann",
    expression="t = pi*sqrt(((r1 + r2)/2)^3/mu)",
    unit="s",
    parameters={"mu": _MU, "r1": _R1, "r2": _R2},
    hypotheses=_HOHMANN_HYPOTHESES,
    references=_REFERENCES,
)
def hohmann_time(mu: float, r1: float, r2: float) -> float:
    """Duração (s) da transferência: metade do período da elipse de transferência."""
    return hohmann_transfer(mu, r1, r2).transfer_time
