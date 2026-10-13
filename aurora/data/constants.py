"""Constantes físicas e astronômicas usadas pelos modelos, em unidades do SI.

Cada constante registra a fonte do valor e se ele é exato por definição. Valores com
incerteza experimental (como ``G``) são usados pelo seu valor central; a incerteza não é
propagada aos resultados nesta versão.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PhysicalConstant:
    """Constante física com unidade e procedência."""

    symbol: str
    name: str
    value: float
    unit: str
    source: str
    exact: bool = False
    relative_uncertainty: float | None = None


GRAVITATIONAL_CONSTANT = PhysicalConstant(
    symbol="G",
    name="Constante da gravitação universal",
    value=6.67430e-11,
    unit="m^3 kg^-1 s^-2",
    source="CODATA 2018 (valor recomendado 6,67430(15) × 10⁻¹¹)",
    relative_uncertainty=2.2e-5,
)

SPEED_OF_LIGHT = PhysicalConstant(
    symbol="c",
    name="Velocidade da luz no vácuo",
    value=299_792_458.0,
    unit="m/s",
    source="SI — valor exato por definição do metro (17ª CGPM, 1983)",
    exact=True,
)

STANDARD_GRAVITY = PhysicalConstant(
    symbol="g0",
    name="Aceleração da gravidade padrão (convencional)",
    value=9.80665,
    unit="m/s^2",
    source="3ª CGPM (1901) — valor convencional exato",
    exact=True,
)

ASTRONOMICAL_UNIT = PhysicalConstant(
    symbol="ua",
    name="Unidade astronômica",
    value=149_597_870_700.0,
    unit="m",
    source="IAU 2012, Resolução B2 — valor exato por definição",
    exact=True,
)

ALL_CONSTANTS: tuple[PhysicalConstant, ...] = (
    GRAVITATIONAL_CONSTANT,
    SPEED_OF_LIGHT,
    STANDARD_GRAVITY,
    ASTRONOMICAL_UNIT,
)

# Atalhos numéricos para uso direto nas fórmulas.
G = GRAVITATIONAL_CONSTANT.value
C = SPEED_OF_LIGHT.value
G0 = STANDARD_GRAVITY.value
AU = ASTRONOMICAL_UNIT.value
