"""Gravitação universal de Newton.

Fórmulas (SI)::

    F  = G·m1·m2 / r²          força entre duas massas
    μ  = G·M                   parâmetro gravitacional
    g  = μ / r²                módulo da aceleração gravitacional
    a⃗  = −μ·r⃗ / |r⃗|³           aceleração vetorial devida a uma massa pontual

Ponto de extensão para modelos futuros
--------------------------------------
Qualquer objeto com o método ``acceleration(position, velocity, time)`` satisfaz o protocolo
:class:`GravityModel`. :class:`NewtonianGravity` soma a atração de várias massas pontuais
(planetas, luas). Modelos relativísticos (por exemplo, a correção pós-newtoniana de primeira
ordem, que depende da velocidade) ou fontes móveis lidas de efemérides poderão ser
implementados com a mesma assinatura e usados pelo integrador orbital da etapa 2 sem
alterá-lo — por isso ``velocity`` e ``time`` fazem parte da interface, embora o modelo
newtoniano estático não os utilize.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from ..core.formula import register_formula
from ..core.validation import as_non_negative, as_positive, as_real, as_vector
from ..data.bodies import CelestialBody
from ..data.constants import G
from ..errors import DomainError, InvalidParameterError
from . import hypotheses as H

_UNDERFLOW_MESSAGE = (
    "O resultado é menor, em módulo, que o menor valor representável em ponto flutuante de 64 bits."
)


@register_formula(
    key="forca_gravitacional",
    title="Força gravitacional entre dois corpos (lei de Newton)",
    expression="F = G*m1*m2/r^2",
    unit="N",
    parameters={
        "m1": ("Massa do corpo 1", "kg"),
        "m2": ("Massa do corpo 2", "kg"),
        "r": ("Distância entre os centros de massa", "m"),
    },
    hypotheses=(H.POINT_MASSES, H.CLASSICAL, H.G_VALUE),
)
def gravitational_force(m1: float, m2: float, r: float) -> float:
    """Módulo da força de atração gravitacional, em N."""
    m1, m2, r = as_non_negative("m1", m1), as_non_negative("m2", m2), as_positive("r", r)
    force = G * (m1 / r) * (m2 / r)  # intercalado: nem r² nem G·m1·m2 sofrem under/overflow à toa
    if force == 0.0 and m1 > 0.0 and m2 > 0.0:
        raise DomainError(_UNDERFLOW_MESSAGE)
    return force


@register_formula(
    key="parametro_gravitacional",
    title="Parâmetro gravitacional",
    expression="mu = G*M",
    unit="m^3/s^2",
    parameters={"M": ("Massa do corpo central", "kg")},
    hypotheses=(H.G_VALUE,),
)
def gravitational_parameter(M: float) -> float:  # noqa: N803 (símbolo da fórmula)
    """``μ = G·M`` em m³/s²."""
    return G * as_non_negative("M", M)


@register_formula(
    key="aceleracao_gravitacional",
    title="Aceleração gravitacional a uma distância r",
    expression="g = mu/r^2",
    unit="m/s^2",
    parameters={
        "mu": ("Parâmetro gravitacional do corpo (G*M)", "m^3/s^2"),
        "r": ("Distância ao centro do corpo", "m"),
    },
    hypotheses=(H.POINT_MASSES, H.CLASSICAL, H.MU_GIVEN, "Não inclui o efeito centrífugo da rotação do corpo."),
)
def gravitational_acceleration(mu: float, r: float) -> float:
    """Módulo da aceleração gravitacional ``μ/r²`` em m/s²."""
    mu, r = as_non_negative("mu", mu), as_positive("r", r)
    g = mu / r / r
    if g == 0.0 and mu > 0.0:
        raise DomainError(_UNDERFLOW_MESSAGE)
    return g


def point_mass_acceleration(mu: float, position: object, source_position: object = (0.0, 0.0, 0.0)) -> np.ndarray:
    """Aceleração vetorial (m/s²) em ``position`` devida a uma massa pontual em ``source_position``."""
    mu = as_non_negative("mu", mu)
    d = as_vector("position", position) - as_vector("source_position", source_position)
    r = math.hypot(*d)  # norma com reescala interna: sem overflow/underflow de |d|²
    if r == 0.0:
        raise DomainError("A posição coincide com a da massa pontual: aceleração indefinida (singularidade).")
    with np.errstate(over="ignore", invalid="ignore"):
        acc = -(mu / r / r) * (d / r)
    if not np.all(np.isfinite(acc)):
        raise DomainError("A aceleração calculada não é representável em ponto flutuante de 64 bits.")
    return acc


@runtime_checkable
class GravityModel(Protocol):
    """Interface de um modelo de campo gravitacional (ver docstring do módulo)."""

    def acceleration(
        self, position: np.ndarray, velocity: np.ndarray | None = None, time: float = 0.0
    ) -> np.ndarray:
        """Aceleração (m/s²) sobre uma partícula de teste no estado dado."""
        ...


@dataclass(frozen=True, eq=False)  # igualdade por identidade: ndarray não define == escalar
class PointMass:
    """Massa pontual com parâmetro gravitacional ``mu`` (m³/s²) fixa em ``position`` (m)."""

    name: str
    mu: float
    position: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "mu", as_non_negative("mu", self.mu))
        pos = as_vector("position", self.position)
        pos.flags.writeable = False
        object.__setattr__(self, "position", pos)

    @classmethod
    def from_body(cls, body: CelestialBody, position: object = (0.0, 0.0, 0.0)) -> PointMass:
        """Cria a massa pontual equivalente a um corpo celeste de referência."""
        return cls(name=body.name, mu=body.mu, position=np.asarray(position))


class NewtonianGravity:
    """Campo newtoniano produzido por massas pontuais estáticas (superposição linear).

    Hipóteses: fontes fixas no referencial inercial adotado; sem correções relativísticas;
    sem efeitos de maré ou de forma (cada fonte é esfericamente simétrica).
    """

    def __init__(self, sources: Iterable[PointMass]) -> None:
        self.sources: tuple[PointMass, ...] = tuple(sources)
        if not self.sources:
            raise InvalidParameterError("É necessário ao menos uma fonte gravitacional.")
        if not all(isinstance(source, PointMass) for source in self.sources):
            raise InvalidParameterError("Todas as fontes gravitacionais devem ser objetos PointMass.")

    def acceleration(
        self, position: np.ndarray, velocity: np.ndarray | None = None, time: float = 0.0
    ) -> np.ndarray:
        as_real("time", time)
        if velocity is not None:
            as_vector("velocity", velocity)  # não usada pelo modelo newtoniano estático, mas validada
        total = np.zeros(3)
        for source in self.sources:
            total += point_mass_acceleration(source.mu, position, source.position)
        return total
