"""Testes de mecânica orbital: identidades analíticas e casos de referência conhecidos."""

import math

import numpy as np
import pytest

from aurora.data.bodies import EARTH, MOON, SUN
from aurora.data.constants import AU
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.orbits import (
    OrbitType,
    apoapsis_radius,
    circular_orbit_radius,
    circular_orbit_summary,
    circular_orbital_period,
    circular_orbital_speed,
    escape_speed,
    orbit_from_state,
    orbital_period,
    periapsis_radius,
    radius_from_altitude,
    semi_major_axis_from_energy,
    specific_orbital_energy,
    vis_viva_speed,
)

MU = EARTH.mu
R_LEO = EARTH.radius + 400e3


# --- fórmulas da especificação --------------------------------------------------------


def test_velocidade_circular_formula():
    assert circular_orbital_speed(MU, R_LEO) == pytest.approx(math.sqrt(MU / R_LEO), rel=1e-15)


def test_periodo_circular_formula():
    expected = 2 * math.pi * math.sqrt(R_LEO**3 / MU)
    assert circular_orbital_period(MU, R_LEO) == pytest.approx(expected, rel=1e-15)


def test_energia_especifica_formula():
    assert specific_orbital_energy(MU, R_LEO, 7000.0) == pytest.approx(7000.0**2 / 2 - MU / R_LEO, rel=1e-15)


# --- identidades analíticas ------------------------------------------------------------


def test_circunferencia_percorrida_em_um_periodo():
    v = circular_orbital_speed(MU, R_LEO)
    T = circular_orbital_period(MU, R_LEO)
    assert v * T == pytest.approx(2 * math.pi * R_LEO, rel=1e-15)


def test_energia_da_orbita_circular():
    v = circular_orbital_speed(MU, R_LEO)
    assert specific_orbital_energy(MU, R_LEO, v) == pytest.approx(-MU / (2 * R_LEO), rel=1e-14)
    assert semi_major_axis_from_energy(MU, -MU / (2 * R_LEO)) == pytest.approx(R_LEO, rel=1e-15)


def test_vis_viva_com_a_igual_r_reproduz_orbita_circular():
    assert vis_viva_speed(MU, R_LEO, R_LEO) == pytest.approx(circular_orbital_speed(MU, R_LEO), rel=1e-15)


def test_velocidade_de_escape_e_raiz_de_dois_vezes_a_circular():
    assert escape_speed(MU, R_LEO) == pytest.approx(math.sqrt(2) * circular_orbital_speed(MU, R_LEO), rel=1e-15)
    assert specific_orbital_energy(MU, R_LEO, escape_speed(MU, R_LEO)) == pytest.approx(0.0, abs=1e-6)


def test_terceira_lei_de_kepler():
    # T² ∝ a³: dobrar o semieixo multiplica o período por 2^(3/2).
    assert orbital_period(MU, 2 * R_LEO) == pytest.approx(2**1.5 * orbital_period(MU, R_LEO), rel=1e-15)
    assert orbital_period(MU, R_LEO) == circular_orbital_period(MU, R_LEO)


def test_raio_por_periodo_e_inversa_do_periodo():
    for r in (R_LEO, 2.6e7, 4.2e7, 3.8e8):
        assert circular_orbit_radius(MU, circular_orbital_period(MU, r)) == pytest.approx(r, rel=1e-14)


def test_periapse_apoapse_e_conservacao_do_momento_angular():
    a, e = 2.0e7, 0.6
    rp, ra = periapsis_radius(a, e), apoapsis_radius(a, e)
    assert rp + ra == pytest.approx(2 * a, rel=1e-15)
    vp, va = vis_viva_speed(MU, rp, a), vis_viva_speed(MU, ra, a)
    assert rp * vp == pytest.approx(ra * va, rel=1e-14)
    # Energia igual nos dois pontos.
    assert specific_orbital_energy(MU, rp, vp) == pytest.approx(specific_orbital_energy(MU, ra, va), rel=1e-13)


# --- casos de referência conhecidos ------------------------------------------------------


def test_orbita_geoestacionaria():
    # Período = dia sideral (≈ 86 164,09 s). Valor de referência amplamente tabulado: r ≈ 42 164 km.
    r = circular_orbit_radius(MU, 86164.0905)
    assert r / 1e3 == pytest.approx(42164.0, abs=1.0)
    assert circular_orbital_speed(MU, r) == pytest.approx(3074.7, abs=0.1)


def test_orbita_circular_de_1_ua_tem_periodo_de_aproximadamente_um_ano_sideral():
    # Ano sideral ≈ 365,256 d. A diferença (~1e-6 relativa) vem de desprezar a massa da
    # Terra-Lua e de a órbita real não ser exatamente circular de raio 1 ua.
    T = circular_orbital_period(SUN.mu, AU)
    assert T / 86400.0 == pytest.approx(365.256, rel=1e-5)


# --- domínio e validação ------------------------------------------------------------------


def test_vis_viva_fora_do_dominio():
    with pytest.raises(DomainError, match="r > 2a"):
        vis_viva_speed(MU, 3e7, 1e7)
    with pytest.raises(InvalidParameterError):
        vis_viva_speed(MU, 1e7, 0.0)
    # Hipérbole (a < 0) é permitida.
    assert vis_viva_speed(MU, 1e7, -1e7) == pytest.approx(math.sqrt(MU * 3e-7), rel=1e-15)


def test_semieixo_com_energia_nula_e_hiperbolica():
    with pytest.raises(DomainError, match="parabólica"):
        semi_major_axis_from_energy(MU, 0.0)
    assert semi_major_axis_from_energy(MU, 1e6) < 0


def test_excentricidade_invalida_para_elipse():
    with pytest.raises(DomainError):
        periapsis_radius(1e7, 1.0)
    with pytest.raises(InvalidParameterError):
        apoapsis_radius(1e7, -0.1)


@pytest.mark.parametrize("mu, r", [(0.0, 1e7), (-MU, 1e7), (MU, 0.0), (MU, -1.0), (MU, math.nan)])
def test_parametros_orbitais_invalidos(mu, r):
    with pytest.raises(InvalidParameterError):
        circular_orbital_speed(mu, r)


# --- elementos a partir do vetor de estado -------------------------------------------------


def test_estado_circular():
    v = circular_orbital_speed(MU, R_LEO)
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [0, v, 0])
    assert orbit.orbit_type is OrbitType.CIRCULAR
    assert orbit.eccentricity < 1e-14
    assert orbit.semi_major_axis == pytest.approx(R_LEO, rel=1e-14)
    assert orbit.period == pytest.approx(circular_orbital_period(MU, R_LEO), rel=1e-14)
    assert orbit.inclination == 0.0


def test_estado_eliptico_no_periapse():
    rp, ra = 7.0e6, 2.0e7
    a = (rp + ra) / 2
    vp = vis_viva_speed(MU, rp, a)
    orbit = orbit_from_state(MU, [rp, 0, 0], [0, 0, vp])  # plano xz: inclinação de 90°
    assert orbit.orbit_type is OrbitType.ELLIPTIC
    assert orbit.semi_major_axis == pytest.approx(a, rel=1e-13)
    assert orbit.eccentricity == pytest.approx((ra - rp) / (ra + rp), rel=1e-13)
    assert orbit.periapsis_radius == pytest.approx(rp, rel=1e-13)
    assert orbit.apoapsis_radius == pytest.approx(ra, rel=1e-13)
    assert orbit.specific_energy == pytest.approx(-MU / (2 * a), rel=1e-13)
    assert orbit.inclination == pytest.approx(math.pi / 2, rel=1e-15)
    assert orbit.eccentricity_vector == pytest.approx([orbit.eccentricity, 0, 0], abs=1e-15)


def test_estado_retrogrado():
    v = circular_orbital_speed(MU, R_LEO)
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [0, -v, 0])
    assert orbit.inclination == pytest.approx(math.pi, rel=1e-15)


def test_estado_parabolico():
    v = escape_speed(MU, R_LEO)
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [0, v, 0])
    assert orbit.orbit_type is OrbitType.PARABOLIC
    assert orbit.semi_major_axis == math.inf
    assert orbit.periapsis_radius == pytest.approx(R_LEO, rel=1e-14)
    assert orbit.apoapsis_radius == math.inf
    assert orbit.period == math.inf


def test_estado_hiperbolico():
    vp = 1.2 * escape_speed(MU, R_LEO)
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [0, vp, 0])
    assert orbit.orbit_type is OrbitType.HYPERBOLIC
    assert orbit.eccentricity == pytest.approx(R_LEO * vp**2 / MU - 1.0, rel=1e-13)
    assert orbit.semi_major_axis < 0
    assert orbit.periapsis_radius == pytest.approx(R_LEO, rel=1e-13)
    assert orbit.apoapsis_radius == math.inf
    assert orbit.period == math.inf


def test_estado_radial_degenerado():
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [1000.0, 0, 0])
    assert orbit.orbit_type is OrbitType.RADIAL
    assert orbit.eccentricity == pytest.approx(1.0, rel=1e-15)
    assert math.isnan(orbit.inclination)
    assert orbit.periapsis_radius == 0.0
    # Queda livre a partir do repouso: apoapse (2a) é a própria posição inicial.
    rest = orbit_from_state(MU, [R_LEO, 0, 0], [0, 0, 0])
    assert rest.orbit_type is OrbitType.RADIAL
    assert rest.apoapsis_radius == pytest.approx(R_LEO, rel=1e-15)


def test_estado_invalido():
    with pytest.raises(DomainError):
        orbit_from_state(MU, [0, 0, 0], [1, 0, 0])
    with pytest.raises(InvalidParameterError):
        orbit_from_state(MU, [1, 0], [1, 0, 0])


def test_estado_com_numpy_e_listas_da_mesmo_resultado():
    r, v = [7e6, 1e5, -2e5], [100.0, 7400.0, 900.0]
    o1 = orbit_from_state(MU, r, v)
    o2 = orbit_from_state(MU, np.array(r), np.array(v))
    assert o1.eccentricity == o2.eccentricity


# --- resumo de órbita circular ---------------------------------------------------------------


def test_resumo_orbita_circular():
    summary = circular_orbit_summary(EARTH, radius_from_altitude(EARTH, 400e3))
    assert summary.radius == R_LEO
    assert summary.altitude == pytest.approx(400e3, rel=1e-9)
    keys = [r.formula.key for r in summary.results]
    assert keys == [
        "velocidade_orbital_circular",
        "periodo_orbital_circular",
        "energia_especifica",
        "velocidade_escape",
        "aceleracao_gravitacional",
    ]
    by_key = {r.formula.key: r.value for r in summary.results}
    assert by_key["periodo_orbital_circular"] == circular_orbital_period(EARTH.mu, R_LEO)


def test_resumo_rejeita_orbita_abaixo_da_superficie():
    with pytest.raises(DomainError, match="intersectaria"):
        circular_orbit_summary(MOON, radius_from_altitude(MOON, -1.0))


# --- regressões da auditoria técnica (docs/AUDITORIA_ETAPA1.md) ---------------------------


def _state_from_elements(mu, a, e, inc, raan, argp, nu):
    """Vetor de estado a partir de elementos clássicos (fórmulas perifocais; referência independente)."""
    p = a * (1 - e * e)
    r = p / (1 + e * math.cos(nu))
    r_pf = np.array([r * math.cos(nu), r * math.sin(nu), 0.0])
    v_pf = np.array([-math.sin(nu), e + math.cos(nu), 0.0]) * math.sqrt(mu / p)
    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(inc), math.sin(inc)
    cw, sw = math.cos(argp), math.sin(argp)
    Q = np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ]
    )
    return Q @ r_pf, Q @ v_pf


def test_ida_e_volta_elementos_estado_elementos():
    # ORBI-6: orientações, anomalias e excentricidades arbitrárias (inclui e → 1 e hipérboles).
    import random

    rng = random.Random(11)
    for _ in range(500):
        hyper = rng.random() < 0.25
        e = rng.uniform(1.05, 5.0) if hyper else rng.choice((rng.uniform(0.0, 0.99), rng.uniform(0.99, 0.9999)))
        a = (-1.0 if hyper else 1.0) * 10 ** rng.uniform(6.8, 8.5)
        inc = rng.uniform(0.0, math.pi)
        nu_max = 0.9 * math.acos(-1.0 / e) if hyper else math.pi
        r, v = _state_from_elements(MU, a, e, inc, rng.uniform(0, 2 * math.pi), rng.uniform(0, 2 * math.pi),
                                    rng.uniform(-nu_max, nu_max))
        orbit = orbit_from_state(MU, r, v)
        assert orbit.orbit_type is (OrbitType.HYPERBOLIC if hyper else OrbitType.ELLIPTIC)
        assert orbit.semi_major_axis == pytest.approx(a, rel=1e-10)
        assert orbit.eccentricity == pytest.approx(e, abs=1e-13)
        assert orbit.inclination == pytest.approx(inc, abs=1e-13)


def test_trajetoria_quase_radial_ligada_nao_e_parabolica():
    # ORBI-1: com velocidade tangencial de 5 cm/s o estado era rotulado PARABÓLICO, com T = ∞.
    radial = orbit_from_state(MU, [EARTH.radius, 0, 0], [3000.0, 0.0, 0])
    for vt in (0.05, 0.2):
        orbit = orbit_from_state(MU, [EARTH.radius, 0, 0], [3000.0, vt, 0])
        assert orbit.orbit_type is OrbitType.ELLIPTIC
        assert orbit.period == pytest.approx(radial.period, rel=1e-9)
        assert orbit.periapsis_radius + orbit.apoapsis_radius == pytest.approx(2 * orbit.semi_major_axis, rel=1e-14)
    hyper = orbit_from_state(MU, [7e6, 0, 0], [2 * escape_speed(MU, 7e6), 1e-2, 0])
    assert hyper.orbit_type is OrbitType.HYPERBOLIC and hyper.semi_major_axis < 0


def test_escape_radial_tem_semieixo_infinito():
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [escape_speed(MU, R_LEO), 0, 0])
    assert orbit.orbit_type is OrbitType.RADIAL
    assert orbit.semi_major_axis == math.inf and orbit.period == math.inf


def test_inclinacao_pequena_e_recuperada_com_precisao():
    # ORBI-3: com acos, i = 1e-9 rad era reportada como 0.
    v = circular_orbital_speed(MU, R_LEO)
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [0, v * math.cos(1e-9), v * math.sin(1e-9)])
    assert orbit.inclination == pytest.approx(1e-9, rel=1e-6)


def test_tolerancia_da_classificacao_e_limitada():
    # PHYS-5: tol ≥ 1 classificava qualquer trajetória como radial.
    v = circular_orbital_speed(MU, R_LEO)
    for tol in (1e-2, 0.5, 1.0):
        with pytest.raises(InvalidParameterError, match="tol"):
            orbit_from_state(MU, [R_LEO, 0, 0], [0, v, 0], tol=tol)


def test_estado_fora_do_intervalo_do_float64():
    with pytest.raises(DomainError, match="representável"):
        orbit_from_state(MU, [1e200, 0, 0], [1e200, 1e200, 0])


def test_vetores_do_resultado_sao_somente_leitura():
    v = circular_orbital_speed(MU, R_LEO)
    orbit = orbit_from_state(MU, [R_LEO, 0, 0], [0, v, 0])
    with pytest.raises(ValueError):
        orbit.angular_momentum[0] = 1.0


def test_raio_nao_positivo_gera_mensagem_de_altitude_negativa():
    # PHYS-10/CLI--8: antes a mensagem citava o parâmetro interno 'radius'.
    for radius in (0.0, -1.0):
        with pytest.raises(DomainError, match="altitude negativa"):
            circular_orbit_summary(EARTH, radius)


def test_periodo_de_1_ua_concorda_com_o_ano_sideral():
    # Com μ do Sol apenas, a diferença para o ano sideral (365,256363 d) é ≈ 1,5e-6 relativa.
    assert circular_orbital_period(SUN.mu, AU) / 86400.0 == pytest.approx(365.256363, rel=3e-6)
