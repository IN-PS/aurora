"""Testes das manobras impulsivas: equação de Tsiolkovsky e transferência de Hohmann."""

import math

import pytest

from aurora.core.formula import get_formula
from aurora.data.bodies import EARTH
from aurora.data.constants import G0
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.maneuvers import (
    hohmann_delta_v,
    hohmann_time,
    hohmann_transfer,
    propellant_for_delta_v,
    rocket_delta_v,
)

MU = EARTH.mu


def test_tsiolkovsky_razao_de_massa_e():
    # m0/mf = e  =>  Δv = Isp·g0 exatamente.
    assert rocket_delta_v(300.0, math.e * 1000.0, 1000.0) == pytest.approx(300.0 * G0, rel=1e-15)
    assert rocket_delta_v(450.0, 5000.0, 5000.0) == 0.0


def test_propelente_e_a_inversa_de_tsiolkovsky():
    for isp, m0, dv in ((300.0, 5000.0, 2500.0), (450.0, 1e5, 9000.0), (220.0, 12.0, 0.5)):
        mp = propellant_for_delta_v(isp, m0, dv)
        assert 0.0 < mp < m0
        assert rocket_delta_v(isp, m0, m0 - mp) == pytest.approx(dv, rel=1e-12)


def test_propelente_para_dv_pequeno_sem_cancelamento():
    # Para Δv → 0, m_p ≈ m0·Δv/(Isp·g0); 1 − exp(−x) perderia algarismos, expm1 não.
    dv = 1e-9
    assert propellant_for_delta_v(300.0, 1000.0, dv) == pytest.approx(1000.0 * dv / (300.0 * G0), rel=1e-9)
    assert propellant_for_delta_v(300.0, 1000.0, 0.0) == 0.0


@pytest.mark.parametrize(
    ("func", "args", "error"),
    [
        (rocket_delta_v, (300.0, 1000.0, 2000.0), DomainError),  # massa final maior que a inicial
        (rocket_delta_v, (0.0, 1000.0, 500.0), InvalidParameterError),
        (rocket_delta_v, (300.0, 1000.0, 0.0), InvalidParameterError),
        (propellant_for_delta_v, (300.0, 1000.0, -1.0), InvalidParameterError),
        (propellant_for_delta_v, (math.nan, 1000.0, 1.0), InvalidParameterError),
    ],
)
def test_entradas_invalidas(func, args, error):
    with pytest.raises(error):
        func(*args)


def test_hohmann_contra_formulas_fechadas():
    r1, r2 = EARTH.radius + 400e3, EARTH.radius + 35_786e3
    t = hohmann_transfer(MU, r1, r2)
    # Formas fechadas independentes (Curtis, eq. 6.2–6.3; Vallado, alg. 36).
    dv1 = math.sqrt(MU / r1) * (math.sqrt(2 * r2 / (r1 + r2)) - 1)
    dv2 = math.sqrt(MU / r2) * (1 - math.sqrt(2 * r1 / (r1 + r2)))
    assert t.delta_v1 == pytest.approx(dv1, rel=1e-13)
    assert t.delta_v2 == pytest.approx(dv2, rel=1e-13)
    assert t.delta_v_total == pytest.approx(dv1 + dv2, rel=1e-13)
    assert t.transfer_time == pytest.approx(math.pi * math.sqrt(((r1 + r2) / 2) ** 3 / MU), rel=1e-14)
    # Valores de referência deste caso (LEO 400 km → GEO), calculados com μ e R do AURORA.
    assert t.delta_v_total == pytest.approx(3853.96, abs=0.01)
    assert t.transfer_time / 3600 == pytest.approx(5.2913, abs=1e-4)
    # Conservação do momento angular na elipse: r1·v_partida = r2·v_chegada.
    assert r1 * t.speed_transfer_departure == pytest.approx(r2 * t.speed_transfer_arrival, rel=1e-14)


def test_hohmann_descida_tem_queimas_retrogradas_e_mesmo_custo():
    r1, r2 = EARTH.radius + 400e3, EARTH.radius + 2000e3
    up, down = hohmann_transfer(MU, r1, r2), hohmann_transfer(MU, r2, r1)
    assert up.delta_v1 > 0 and up.delta_v2 > 0
    assert down.delta_v1 < 0 and down.delta_v2 < 0
    assert down.delta_v_total == pytest.approx(up.delta_v_total, rel=1e-13)
    assert down.transfer_time == pytest.approx(up.transfer_time, rel=1e-15)


def test_hohmann_entre_orbitas_iguais_nao_custa_nada():
    t = hohmann_transfer(MU, 7e6, 7e6)
    assert t.delta_v_total == pytest.approx(0.0, abs=1e-9)
    assert t.transfer_time == pytest.approx(math.pi * math.sqrt(7e6**3 / MU))


def test_formulas_registradas_com_unidades_si():
    r1, r2 = 7e6, 4.2e7
    assert get_formula("delta_v_hohmann").unit == "m/s"
    assert get_formula("tempo_hohmann").unit == "s"
    assert get_formula("massa_propelente").unit == "kg"
    assert get_formula("delta_v_tsiolkovsky").unit == "m/s"
    assert hohmann_delta_v(MU, r1, r2) == hohmann_transfer(MU, r1, r2).delta_v_total
    assert hohmann_time(MU, r1, r2) == hohmann_transfer(MU, r1, r2).transfer_time
    with pytest.raises(InvalidParameterError):
        hohmann_delta_v(MU, -1.0, r2)
