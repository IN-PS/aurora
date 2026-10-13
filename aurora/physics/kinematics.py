"""Cinemática do movimento com aceleração constante (movimento uniformemente variado — MUV).

Fórmulas (SI)::

    x(t)  = x0 + v0·t + a·t²/2         posição
    v(t)  = v0 + a·t                   velocidade
    Δx    = v0·t + a·t²/2              deslocamento
    a_m   = (v − v0)/Δt                aceleração média
    |v|   = √(v0² + 2·a·Δx)            equação de Torricelli

As fórmulas escalares (uma dimensão) estão registradas para uso pela interface. A função
vetorial :func:`constant_acceleration_state` aplica as mesmas relações em três dimensões.
"""

from __future__ import annotations

import math

import numpy as np

from ..core.formula import register_formula
from ..core.validation import as_nonzero, as_real, as_vector
from ..errors import DomainError
from . import hypotheses as H
from .quadratic import solve_quadratic

_MUV_HYPOTHESES = (H.CONSTANT_ACCELERATION, H.ONE_DIMENSIONAL, H.INERTIAL_FRAME, H.CLASSICAL)


@register_formula(
    key="posicao_muv",
    title="Posição no movimento uniformemente variado",
    expression="x = x0 + v0*t + a*t^2/2",
    unit="m",
    parameters={
        "x0": ("Posição inicial", "m"),
        "v0": ("Velocidade inicial", "m/s"),
        "a": ("Aceleração (constante)", "m/s^2"),
        "t": ("Tempo decorrido desde o instante inicial", "s"),
    },
    hypotheses=_MUV_HYPOTHESES,
)
def position_constant_acceleration(x0: float, v0: float, a: float, t: float) -> float:
    """Posição ``x(t)`` em m."""
    x0, v0, a, t = as_real("x0", x0), as_real("v0", v0), as_real("a", a), as_real("t", t)
    return x0 + (v0 * t + 0.5 * a * t * t)


@register_formula(
    key="velocidade_muv",
    title="Velocidade no movimento uniformemente variado",
    expression="v = v0 + a*t",
    unit="m/s",
    parameters={
        "v0": ("Velocidade inicial", "m/s"),
        "a": ("Aceleração (constante)", "m/s^2"),
        "t": ("Tempo decorrido desde o instante inicial", "s"),
    },
    hypotheses=_MUV_HYPOTHESES,
)
def velocity_constant_acceleration(v0: float, a: float, t: float) -> float:
    """Velocidade ``v(t)`` em m/s."""
    v0, a, t = as_real("v0", v0), as_real("a", a), as_real("t", t)
    return v0 + a * t


@register_formula(
    key="deslocamento_muv",
    title="Deslocamento no movimento uniformemente variado",
    expression="dx = v0*t + a*t^2/2",
    unit="m",
    parameters={
        "v0": ("Velocidade inicial", "m/s"),
        "a": ("Aceleração (constante)", "m/s^2"),
        "t": ("Tempo decorrido desde o instante inicial", "s"),
    },
    hypotheses=_MUV_HYPOTHESES,
)
def displacement_constant_acceleration(v0: float, a: float, t: float) -> float:
    """Deslocamento ``Δx = x(t) − x0`` em m."""
    v0, a, t = as_real("v0", v0), as_real("a", a), as_real("t", t)
    return v0 * t + 0.5 * a * t * t


@register_formula(
    key="aceleracao_media",
    title="Aceleração média",
    expression="a_m = (v - v0)/dt",
    unit="m/s^2",
    parameters={
        "v0": ("Velocidade no início do intervalo", "m/s"),
        "v": ("Velocidade no fim do intervalo", "m/s"),
        "dt": ("Duração do intervalo (não nula)", "s"),
    },
    hypotheses=(H.ONE_DIMENSIONAL, H.INERTIAL_FRAME),
)
def average_acceleration(v0: float, v: float, dt: float) -> float:
    """Aceleração média ``(v − v0)/Δt`` em m/s²."""
    v0, v, dt = as_real("v0", v0), as_real("v", v), as_nonzero("dt", dt)
    return (v - v0) / dt


@register_formula(
    key="velocidade_torricelli",
    title="Módulo da velocidade pela equação de Torricelli",
    expression="|v| = sqrt(v0^2 + 2*a*dx)",
    unit="m/s",
    parameters={
        "v0": ("Velocidade inicial", "m/s"),
        "a": ("Aceleração (constante)", "m/s^2"),
        "dx": ("Deslocamento", "m"),
    },
    hypotheses=_MUV_HYPOTHESES,
)
def torricelli_speed(v0: float, a: float, dx: float) -> float:
    """Módulo da velocidade após o deslocamento ``dx``, em m/s.

    Se ``v0² + 2·a·Δx < 0`` o deslocamento nunca é atingido (o móvel inverte o sentido
    antes) e é levantado :class:`~aurora.errors.DomainError`.
    """
    v0, a, dx = as_real("v0", v0), as_real("a", a), as_real("dx", dx)
    radicand = v0 * v0 + 2.0 * a * dx
    if radicand < 0.0:
        raise DomainError(
            "v0² + 2·a·Δx < 0: com essa velocidade inicial e aceleração o deslocamento informado nunca é atingido."
        )
    return math.sqrt(radicand)


def constant_acceleration_state(
    r0: object, v0: object, a: object, t: float
) -> tuple[np.ndarray, np.ndarray]:
    """Posição e velocidade em 3D após o tempo ``t`` com aceleração constante.

    ``r0`` (m), ``v0`` (m/s) e ``a`` (m/s²) são vetores de 3 componentes; ``t`` em s.
    Retorna ``(r, v)`` como vetores ``numpy`` de ``float64``.
    """
    r0v, v0v, av = as_vector("r0", r0), as_vector("v0", v0), as_vector("a", a)
    t = as_real("t", t)
    r = r0v + (v0v * t + 0.5 * av * t * t)
    v = v0v + av * t
    return r, v


def time_to_reach_position(x0: float, v0: float, a: float, x: float) -> tuple[float, ...]:
    """Instantes ``t ≥ 0`` em que o móvel em MUV passa pela posição ``x``.

    Resolve ``a/2·t² + v0·t + (x0 − x) = 0`` com a calculadora de segundo grau. Retorna
    uma tupla ordenada com zero, um ou dois instantes (s). Se ``a = v0 = 0`` e ``x = x0`` o
    móvel permanece sempre na posição e há infinitas soluções: levanta ``DomainError``.
    """
    x0, v0, a, x = as_real("x0", x0), as_real("v0", v0), as_real("a", a), as_real("x", x)
    gap = x0 - x
    if a == 0.0:
        if v0 == 0.0:
            if gap == 0.0:
                raise DomainError("Móvel em repouso na própria posição alvo: infinitas soluções.")
            return ()
        t = -gap / v0 + 0.0
        return (t,) if t >= 0.0 else ()
    solution = solve_quadratic(0.5 * a, v0, gap)
    return tuple(t for t in solution.real_roots if t >= 0.0)
