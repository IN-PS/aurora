"""Testes da validação de parâmetros e das unidades."""

import math

import numpy as np
import pytest

from aurora.core import units
from aurora.core.units import Dimension
from aurora.core.validation import as_non_negative, as_nonzero, as_positive, as_real, as_vector
from aurora.data.constants import AU
from aurora.errors import InvalidParameterError, NotFoundError, UnitError

# --- validação -----------------------------------------------------------------------


def test_as_real_aceita_inteiros_floats_e_numpy():
    assert as_real("x", 3) == 3.0
    assert as_real("x", 2.5) == 2.5
    assert as_real("x", np.float32(1.5)) == 1.5
    assert isinstance(as_real("x", np.int64(7)), float)


@pytest.mark.parametrize("bad", [True, False, "1.0", None, 1j, math.nan, math.inf, [1.0], 10**400])
def test_as_real_rejeita_entradas_invalidas(bad):
    with pytest.raises(InvalidParameterError):
        as_real("x", bad)


def test_restricoes_de_sinal():
    assert as_positive("x", 1e-300) == 1e-300
    with pytest.raises(InvalidParameterError, match="positivo"):
        as_positive("x", 0.0)
    assert as_non_negative("x", 0.0) == 0.0
    with pytest.raises(InvalidParameterError, match="negativo"):
        as_non_negative("x", -1.0)
    with pytest.raises(InvalidParameterError, match="zero"):
        as_nonzero("x", 0.0)


def test_as_vector():
    v = as_vector("r", [1, 2, 3])
    assert v.dtype == np.float64 and v.tolist() == [1.0, 2.0, 3.0]
    for bad in ([1, 2], [1, 2, math.nan], ["1", "2", "3"], [True, False, True], [1, 2, 3j], "abc"):
        with pytest.raises(InvalidParameterError):
            as_vector("r", bad)


def test_as_vector_retorna_copia():
    original = np.array([1.0, 2.0, 3.0])
    v = as_vector("r", original)
    v[0] = 99.0
    assert original[0] == 1.0


def test_mensagens_em_portugues():
    with pytest.raises(InvalidParameterError, match="O parâmetro 'massa' deve ser um número real"):
        as_real("massa", "dez")


# --- unidades ------------------------------------------------------------------------


def test_conversoes_basicas():
    assert units.to_si(400, "km") == 400_000.0
    assert units.to_si(1, "ua") == AU == 149_597_870_700.0
    assert units.to_si(1, "au") == AU
    assert units.to_si(2, "h") == 7200.0
    assert units.to_si(1, "d") == 86400.0
    assert units.to_si(36, "km/h") == pytest.approx(10.0, rel=1e-15)
    assert units.to_si(398600.4418, "km^3/s^2") == pytest.approx(3.986004418e14, rel=1e-15)
    assert units.to_si(180, "grau") == pytest.approx(math.pi, rel=1e-15)
    assert units.from_si(7200.0, "min") == 120.0


def test_ida_e_volta():
    for unit in ("km", "ua", "min", "h", "d", "t", "km/s", "km/h", "MJ/kg", "km^3/s^2", "grau"):
        assert units.from_si(units.to_si(123.456, unit), unit) == pytest.approx(123.456, rel=1e-15)


def test_unidade_desconhecida():
    with pytest.raises(NotFoundError, match="Unidade desconhecida"):
        units.get_unit("parsec")


def test_dimensao_incompativel():
    with pytest.raises(UnitError, match="tempo"):
        units.to_si(5, "s", expected=Dimension.LENGTH)


@pytest.mark.parametrize(
    "text, default, expected",
    [
        ("400 km", "m", 400_000.0),
        ("400km", "m", 400_000.0),
        ("  6778,137 km ", "m", 6_778_137.0),
        ("1.5e3", "m", 1500.0),
        ("-2,5", "m/s", -2.5),
        (".5 h", "s", 1800.0),
        ("0.5", "1", 0.5),
        ("398600.4418 km³/s²", "m^3/s^2", 3.986004418e14),
    ],
)
def test_parse_quantity(text, default, expected):
    assert units.parse_quantity(text, default) == pytest.approx(expected, rel=1e-15)


@pytest.mark.parametrize("text", ["", "km", "abc", "1.2.3 km", "6.378.137", "nan", "inf m"])
def test_parse_quantity_invalida(text):
    with pytest.raises((InvalidParameterError, NotFoundError)):
        units.parse_quantity(text, "m")


def test_parse_quantity_overflow():
    with pytest.raises(InvalidParameterError):
        units.parse_quantity("1e300 ua", "m")


def test_parse_quantity_dimensao_errada():
    with pytest.raises(UnitError):
        units.parse_quantity("10 s", "m")


def test_parse_number():
    assert units.parse_number("-3,5") == -3.5
    with pytest.raises(InvalidParameterError):
        units.parse_number("3 km")


def test_cada_dimensao_tem_unidade_si():
    for dim in Dimension:
        assert units.si_unit(dim).factor == 1.0
