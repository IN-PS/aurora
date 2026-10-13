"""Testes da gravitação universal de Newton."""

import math

import numpy as np
import pytest

from aurora.data.bodies import EARTH, MOON
from aurora.data.constants import G
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.gravitation import (
    GravityModel,
    NewtonianGravity,
    PointMass,
    gravitational_acceleration,
    gravitational_force,
    gravitational_parameter,
    point_mass_acceleration,
)


def test_forca_formula_direta():
    assert gravitational_force(m1=1000.0, m2=2.0, r=10.0) == pytest.approx(G * 1000.0 * 2.0 / 100.0, rel=1e-15)


def test_forca_lei_do_inverso_do_quadrado_e_simetria():
    f1 = gravitational_force(5e24, 1e3, 7e6)
    assert gravitational_force(5e24, 1e3, 14e6) == pytest.approx(f1 / 4.0, rel=1e-15)
    assert gravitational_force(1e3, 5e24, 7e6) == f1
    assert gravitational_force(0.0, 1.0, 1.0) == 0.0


def test_forca_parametros_invalidos():
    with pytest.raises(InvalidParameterError):
        gravitational_force(1.0, 1.0, 0.0)
    with pytest.raises(InvalidParameterError):
        gravitational_force(-1.0, 1.0, 1.0)


def test_parametro_gravitacional():
    assert gravitational_parameter(M=1.0) == G
    assert gravitational_parameter(M=2e30) == pytest.approx(G * 2e30, rel=1e-15)


def test_gravidade_na_superficie_terrestre_modelo_esferico():
    # μ/R² com o raio equatorial: ≈ 9.798 m/s² (sem efeito centrífugo e sem achatamento;
    # difere do valor convencional g0 = 9.80665 m/s², que não é um valor local calculado).
    g = gravitational_acceleration(EARTH.mu, EARTH.radius)
    assert g == pytest.approx(EARTH.mu / EARTH.radius**2, rel=1e-15)
    assert 9.79 < g < 9.81


def test_forca_igual_massa_vezes_aceleracao():
    m = 1234.0
    M = EARTH.mu / G
    r = 7.0e6
    assert gravitational_force(M, m, r) == pytest.approx(m * gravitational_acceleration(EARTH.mu, r), rel=1e-14)


def test_aceleracao_vetorial_aponta_para_a_fonte():
    a = point_mass_acceleration(EARTH.mu, [7e6, 0, 0])
    assert a[0] < 0 and a[1] == 0 and a[2] == 0
    assert np.linalg.norm(a) == pytest.approx(gravitational_acceleration(EARTH.mu, 7e6), rel=1e-15)

    b = point_mass_acceleration(1.0, [1, 1, 1], source_position=[1, 1, 3])
    assert b == pytest.approx([0.0, 0.0, 0.25])


def test_aceleracao_na_singularidade():
    with pytest.raises(DomainError):
        point_mass_acceleration(1.0, [1, 2, 3], [1, 2, 3])


def test_superposicao_de_duas_massas_iguais_se_anula_no_ponto_medio():
    model = NewtonianGravity([PointMass("A", 5.0, [-1.0, 0, 0]), PointMass("B", 5.0, [1.0, 0, 0])])
    assert model.acceleration(np.zeros(3)) == pytest.approx([0.0, 0.0, 0.0], abs=1e-15)


def test_superposicao_terra_lua_e_soma_das_contribuicoes():
    earth_moon = 384_400e3  # distância média aproximada, apenas para o teste
    earth = PointMass.from_body(EARTH, [0, 0, 0])
    moon = PointMass.from_body(MOON, [earth_moon, 0, 0])
    model = NewtonianGravity([earth, moon])
    p = np.array([1e8, 2e7, -3e6])
    expected = point_mass_acceleration(EARTH.mu, p) + point_mass_acceleration(MOON.mu, p, [earth_moon, 0, 0])
    assert model.acceleration(p) == pytest.approx(expected, rel=1e-15)
    assert isinstance(model, GravityModel)


def test_massa_pontual_e_imutavel():
    pm = PointMass("X", 1.0, [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        pm.position[0] = 5.0
    with pytest.raises(InvalidParameterError):
        PointMass("X", -1.0, [0, 0, 0])
    with pytest.raises(InvalidParameterError):
        NewtonianGravity([])


def test_valores_de_referencia_dos_corpos():
    assert EARTH.mu == 3.986004418e14
    assert EARTH.radius == 6_378_137.0
    assert all(body.sources for body in (EARTH, MOON))
    assert math.isclose(G, 6.67430e-11)
