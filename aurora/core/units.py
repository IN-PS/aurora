"""Unidades de medida e conversão para o SI.

Todos os modelos físicos do AURORA trabalham internamente em unidades do SI. Este módulo
converte valores informados pelo usuário em outras unidades (km, ua, min, h, ...) para o SI
na entrada, e do SI para unidades de exibição na saída, verificando a compatibilidade de
dimensões (não é possível informar "400 s" onde se espera um comprimento).

A tabela é deliberadamente pequena e explícita. Uma biblioteca completa de unidades (como
``pint``) pode substituí-la no futuro sem alterar os modelos, que só recebem ``float`` em SI.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from enum import Enum

from ..data.constants import AU
from ..errors import InvalidParameterError, NotFoundError, UnitError


class Dimension(Enum):
    """Grandezas físicas suportadas."""

    DIMENSIONLESS = "adimensional"
    LENGTH = "comprimento"
    TIME = "tempo"
    MASS = "massa"
    SPEED = "velocidade"
    ACCELERATION = "aceleração"
    FORCE = "força"
    GRAVITATIONAL_PARAMETER = "parâmetro gravitacional"
    SPECIFIC_ENERGY = "energia específica"
    ANGLE = "ângulo"


@dataclass(frozen=True)
class Unit:
    """Unidade de medida: ``valor_SI = valor * factor``."""

    symbol: str
    dimension: Dimension
    factor: float
    name: str


_D = Dimension
_UNIT_LIST: tuple[Unit, ...] = (
    Unit("1", _D.DIMENSIONLESS, 1.0, "adimensional"),
    Unit("m", _D.LENGTH, 1.0, "metro"),
    Unit("km", _D.LENGTH, 1e3, "quilômetro"),
    Unit("ua", _D.LENGTH, AU, "unidade astronômica"),
    Unit("s", _D.TIME, 1.0, "segundo"),
    Unit("min", _D.TIME, 60.0, "minuto"),
    Unit("h", _D.TIME, 3600.0, "hora"),
    Unit("d", _D.TIME, 86400.0, "dia (86 400 s)"),
    Unit("kg", _D.MASS, 1.0, "quilograma"),
    Unit("t", _D.MASS, 1e3, "tonelada"),
    Unit("m/s", _D.SPEED, 1.0, "metro por segundo"),
    Unit("km/s", _D.SPEED, 1e3, "quilômetro por segundo"),
    Unit("km/h", _D.SPEED, 1e3 / 3600.0, "quilômetro por hora"),
    Unit("m/s^2", _D.ACCELERATION, 1.0, "metro por segundo ao quadrado"),
    Unit("N", _D.FORCE, 1.0, "newton"),
    Unit("kN", _D.FORCE, 1e3, "quilonewton"),
    Unit("m^3/s^2", _D.GRAVITATIONAL_PARAMETER, 1.0, "metro cúbico por segundo ao quadrado"),
    Unit("km^3/s^2", _D.GRAVITATIONAL_PARAMETER, 1e9, "quilômetro cúbico por segundo ao quadrado"),
    Unit("J/kg", _D.SPECIFIC_ENERGY, 1.0, "joule por quilograma"),
    Unit("MJ/kg", _D.SPECIFIC_ENERGY, 1e6, "megajoule por quilograma"),
    Unit("km^2/s^2", _D.SPECIFIC_ENERGY, 1e6, "quilômetro quadrado por segundo ao quadrado"),
    Unit("rad", _D.ANGLE, 1.0, "radiano"),
    Unit("grau", _D.ANGLE, math.pi / 180.0, "grau"),
)

_UNITS: dict[str, Unit] = {u.symbol: u for u in _UNIT_LIST}

_ALIASES: dict[str, str] = {
    "au": "ua",
    "AU": "ua",
    "UA": "ua",
    "dia": "d",
    "m/s2": "m/s^2",
    "m/s²": "m/s^2",
    "m3/s2": "m^3/s^2",
    "m³/s²": "m^3/s^2",
    "km3/s2": "km^3/s^2",
    "km³/s²": "km^3/s^2",
    "km2/s2": "km^2/s^2",
    "km²/s²": "km^2/s^2",
    "m^2/s^2": "J/kg",
    "m2/s2": "J/kg",
    "m²/s²": "J/kg",
    "°": "grau",
    "deg": "grau",
}

# Unidades de exibição alternativas para cada unidade SI (usadas pela interface).
_DISPLAY_ALTERNATIVES: dict[str, tuple[str, ...]] = {
    "m": ("km", "ua"),
    "s": ("min", "h", "d"),
    "m/s": ("km/s", "km/h"),
    "J/kg": ("MJ/kg",),
    "m^3/s^2": ("km^3/s^2",),
    "rad": ("grau",),
    "N": ("kN",),
}

_SI_BY_DIMENSION: dict[Dimension, Unit] = {u.dimension: u for u in _UNIT_LIST if u.factor == 1.0}


def get_unit(symbol: str) -> Unit:
    """Retorna a unidade pelo símbolo (aceita alguns sinônimos, como ``au`` e ``m/s2``)."""
    key = symbol.strip()
    key = _ALIASES.get(key, key)
    unit = _UNITS.get(key)
    if unit is None:
        known = ", ".join(_UNITS)
        raise NotFoundError(f"Unidade desconhecida: {symbol!r}. Unidades aceitas: {known}.")
    return unit


def si_unit(dimension: Dimension) -> Unit:
    """Unidade SI coerente para a grandeza."""
    return _SI_BY_DIMENSION[dimension]


def list_units() -> tuple[Unit, ...]:
    """Todas as unidades conhecidas."""
    return _UNIT_LIST


def _check_dimension(unit: Unit, expected: Dimension | None) -> None:
    if expected is not None and unit.dimension is not expected:
        raise UnitError(
            f"Unidade '{unit.symbol}' é de {unit.dimension.value}, mas era esperada uma unidade de {expected.value}."
        )


def to_si(value: float, unit: str, expected: Dimension | None = None) -> float:
    """Converte ``value`` expresso em ``unit`` para o SI."""
    u = get_unit(unit)
    _check_dimension(u, expected)
    return value * u.factor


def from_si(value: float, unit: str, expected: Dimension | None = None) -> float:
    """Converte ``value`` (em SI) para ``unit``."""
    u = get_unit(unit)
    _check_dimension(u, expected)
    return value / u.factor


def display_alternatives(si_symbol: str) -> tuple[str, ...]:
    """Unidades alternativas sugeridas para exibir um valor expresso em ``si_symbol``."""
    return _DISPLAY_ALTERNATIVES.get(si_symbol, ())


_QUANTITY_RE = re.compile(
    r"^\s*(?P<num>[+-]?(?:\d+(?:[.,]\d*)?|[.,]\d+)(?:[eE][+-]?\d+)?)\s*(?P<unit>\S.*?)?\s*$"
)


def parse_number(text: str) -> float:
    """Lê um número real finito; aceita vírgula ou ponto como separador decimal (ISO 80000-1)."""
    match = _QUANTITY_RE.match(text)
    if match is None or match.group("unit"):
        raise InvalidParameterError(f"Valor numérico inválido: {text!r}.")
    value = float(match.group("num").replace(",", "."))
    if not math.isfinite(value):
        raise InvalidParameterError(f"Valor numérico fora do intervalo representável: {text!r}.")
    return value


def parse_quantity(text: str, default_unit: str) -> float:
    """Lê uma grandeza como ``"400 km"`` ou ``"400km"`` e retorna o valor em SI.

    Sem unidade explícita, assume-se ``default_unit``. A unidade informada precisa ter a
    mesma dimensão de ``default_unit``. Separador de milhar não é aceito.
    """
    expected = get_unit(default_unit)
    match = _QUANTITY_RE.match(text)
    if match is None:
        raise InvalidParameterError(
            f"Grandeza inválida: {text!r}. Use o formato '<número> [unidade]', por exemplo '400 km'."
        )
    value = float(match.group("num").replace(",", "."))
    unit = get_unit(match.group("unit")) if match.group("unit") else expected
    _check_dimension(unit, expected.dimension)
    result = value * unit.factor
    if not math.isfinite(result):
        raise InvalidParameterError(f"Valor fora do intervalo representável: {text!r}.")
    return result
