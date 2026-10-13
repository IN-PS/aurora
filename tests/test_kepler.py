"""Testes da propagação analítica de Kepler (variável universal) e dos elementos orbitais.

As referências são independentes do código testado: a equação de Kepler clássica
(M = E − e·sen E, elipse; M = e·senh H − H, hipérbole) e a equação de Barker (parábola),
resolvidas aqui mesmo por Newton.
"""

import math

import numpy as np
import pytest

from aurora.data.bodies import EARTH
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.kepler import propagate_kepler, state_from_elements, stumpff_c, stumpff_s
from aurora.physics.orbits import orbit_from_state, orbital_period

MU = EARTH.mu


def _closed_c(z):
    if z > 0:
        return (1 - math.cos(math.sqrt(z))) / z
    return (math.cosh(math.sqrt(-z)) - 1) / (-z)


def _closed_s(z):
    if z > 0:
        s = math.sqrt(z)
        return (s - math.sin(s)) / s**3
    s = math.sqrt(-z)
    return (math.sinh(s) - s) / s**3


@pytest.mark.parametrize("z", [-50.0, -1.0, -0.1, 0.1, 1.0, 30.0, 4 * math.pi**2 - 0.5])
def test_stumpff_forma_fechada_fora_da_serie(z):
    assert stumpff_c(z) == pytest.approx(_closed_c(z), rel=1e-14)
    assert stumpff_s(z) == pytest.approx(_closed_s(z), rel=1e-14)


def test_stumpff_serie_proxima_de_zero():
    assert stumpff_c(0.0) == 0.5 and stumpff_s(0.0) == pytest.approx(1 / 6, rel=1e-15)
    # Continuidade na fronteira série/forma fechada (|z| = 0,1).
    for z in (0.1, -0.1):
        below = z * (1 - 1e-12)
        assert stumpff_c(below) == pytest.approx(_closed_c(z), rel=1e-11)
        assert stumpff_s(below) == pytest.approx(_closed_s(z), rel=1e-11)
    # Para z minúsculo a forma fechada perde algarismos; a série não.
    assert stumpff_c(1e-10) == pytest.approx(0.5 - 1e-10 / 24, rel=1e-15)


def _elliptic_reference(a, e, t):
    """Estado no referencial perifocal no instante t após o periapse (equação de Kepler)."""
    n = math.sqrt(MU / a**3)
    mean = math.fmod(n * t, 2 * math.pi)
    ecc = mean if e < 0.8 else math.pi
    for _ in range(60):
        ecc -= (ecc - e * math.sin(ecc) - mean) / (1 - e * math.cos(ecc))
    r = a * (1 - e * math.cos(ecc))
    position = np.array((a * (math.cos(ecc) - e), a * math.sqrt(1 - e * e) * math.sin(ecc), 0.0))
    factor = math.sqrt(MU * a) / r
    velocity = factor * np.array((-math.sin(ecc), math.sqrt(1 - e * e) * math.cos(ecc), 0.0))
    return position, velocity


def _hyperbolic_reference(a, e, t):
    """Idem para a hipérbole (a < 0): M = e·senh H − H."""
    n = math.sqrt(MU / (-a) ** 3)
    mean = n * t
    big_h = math.asinh(mean / e)
    for _ in range(100):
        big_h -= (e * math.sinh(big_h) - big_h - mean) / (e * math.cosh(big_h) - 1)
    r = a * (1 - e * math.cosh(big_h))
    b = -a * math.sqrt(e * e - 1)
    position = np.array((-a * (e - math.cosh(big_h)), b * math.sinh(big_h), 0.0))
    velocity = math.sqrt(-MU * a) / r * np.array((-math.sinh(big_h), math.sqrt(e * e - 1) * math.cosh(big_h), 0.0))
    return position, velocity


@pytest.mark.parametrize("e", [0.0, 0.01, 0.3, 0.7, 0.95])
@pytest.mark.parametrize("fraction", [0.013, 0.25, 0.5, 0.77, 3.4])
def test_elipse_contra_equacao_de_kepler(e, fraction):
    a = 2.5e7
    r0, v0 = state_from_elements(MU, a, e)  # periapse sobre o eixo x
    t = fraction * orbital_period(MU, a)
    r, v = propagate_kepler(MU, r0, v0, t)
    r_ref, v_ref = _elliptic_reference(a, e, t)
    assert np.linalg.norm(r - r_ref) <= 1e-9 * a
    assert np.linalg.norm(v - v_ref) <= 1e-9 * np.linalg.norm(v0)


@pytest.mark.parametrize("e", [1.05, 1.5, 3.0, 10.0])
@pytest.mark.parametrize("t", [60.0, 3600.0, 86400.0, 3e7])
def test_hiperbole_contra_equacao_de_kepler(e, t):
    a = -1.2e7
    r0, v0 = state_from_elements(MU, a, e)
    r, v = propagate_kepler(MU, r0, v0, t)
    r_ref, v_ref = _hyperbolic_reference(a, e, t)
    assert np.linalg.norm(r - r_ref) <= 1e-10 * np.linalg.norm(r_ref)
    assert np.linalg.norm(v - v_ref) <= 1e-10 * np.linalg.norm(v_ref)


def test_parabola_satisfaz_equacao_de_barker():
    rp = 7.0e6
    r0 = np.array((rp, 0.0, 0.0))
    v0 = np.array((0.0, math.sqrt(2 * MU / rp), 0.0))  # velocidade de escape: parábola
    p = 2 * rp
    for t in (100.0, 5000.0, 1e6):
        r, v = propagate_kepler(MU, r0, v0, t)
        nu = math.atan2(r[1], r[0])
        d = math.tan(nu / 2)
        barker = 0.5 * math.sqrt(p**3 / MU) * (d + d**3 / 3)
        assert barker == pytest.approx(t, rel=1e-12)
        # Energia específica nula ao longo da parábola.
        assert abs(0.5 * v @ v - MU / np.linalg.norm(r)) <= 1e-9 * (v0 @ v0)


def test_ida_e_volta_retorna_ao_estado_inicial():
    r0, v0 = state_from_elements(MU, 2.4e7, 0.73, 0.47, 1.0, 4.0, 0.3)
    r1, v1 = propagate_kepler(MU, r0, v0, 12345.6)
    r2, v2 = propagate_kepler(MU, r1, v1, -12345.6)
    assert np.linalg.norm(r2 - r0) <= 1e-6  # m
    assert np.linalg.norm(v2 - v0) <= 1e-9  # m/s


def test_periodo_completo_e_propagacao_longa():
    a = 7.0e6
    r0, v0 = state_from_elements(MU, a, 0.1, 0.9, 0.2, 0.4, 1.0)
    period = orbital_period(MU, a)
    r, v = propagate_kepler(MU, r0, v0, period)
    assert np.linalg.norm(r - r0) <= 1e-6
    # 10 000 períodos: a redução pelo período mantém a anomalia universal limitada.
    r, v = propagate_kepler(MU, r0, v0, 10_000 * period + 0.25 * period)
    r_ref, v_ref = propagate_kepler(MU, r0, v0, 0.25 * period)
    assert np.linalg.norm(r - r_ref) <= 1e-6 * a


def test_conserva_energia_e_momento_angular():
    r0, v0 = state_from_elements(MU, 1.5e7, 0.6, 1.1, 0.3, 2.0, 0.0)
    energy0 = 0.5 * v0 @ v0 - MU / np.linalg.norm(r0)
    h0 = np.cross(r0, v0)
    for t in np.linspace(-50000, 50000, 41):
        r, v = propagate_kepler(MU, r0, v0, float(t))
        assert 0.5 * v @ v - MU / np.linalg.norm(r) == pytest.approx(energy0, rel=1e-11)
        assert np.linalg.norm(np.cross(r, v) - h0) <= 1e-11 * np.linalg.norm(h0)


def test_trajetoria_radial():
    # Lançamento vertical abaixo da velocidade de escape: momento angular nulo.
    r0 = np.array((7.0e6, 0.0, 0.0))
    v0 = np.array((5000.0, 0.0, 0.0))
    energy0 = 0.5 * v0 @ v0 - MU / 7.0e6
    r, v = propagate_kepler(MU, r0, v0, 600.0)
    assert r[1] == 0.0 and r[2] == 0.0 and r[0] > 7.0e6
    assert 0.5 * v @ v - MU / np.linalg.norm(r) == pytest.approx(energy0, rel=1e-10)


def test_intervalo_nulo_devolve_o_estado():
    r0, v0 = np.array((7e6, 0, 0.0)), np.array((0, 7.5e3, 0.0))
    r, v = propagate_kepler(MU, r0, v0, 0.0)
    assert np.array_equal(r, r0) and np.array_equal(v, v0)


@pytest.mark.parametrize(
    ("args", "error"),
    [
        ((MU, (0, 0, 0), (1, 0, 0), 10.0), DomainError),
        ((0.0, (7e6, 0, 0), (0, 7e3, 0), 10.0), InvalidParameterError),
        ((MU, (7e6, 0), (0, 7e3, 0), 10.0), InvalidParameterError),
        ((MU, (7e6, 0, 0), (0, 7e3, 0), math.nan), InvalidParameterError),
        ((MU, (7e6, 0, math.inf), (0, 7e3, 0), 1.0), InvalidParameterError),
    ],
)
def test_propagate_kepler_rejeita_entradas_invalidas(args, error):
    with pytest.raises(error):
        propagate_kepler(*args)


def test_elementos_orbitais_ida_e_volta():
    a, e, inc, raan, argp = 2.66e7, 0.74, math.radians(63.4), math.radians(40), math.radians(270)
    r, v = state_from_elements(MU, a, e, inc, raan, argp, math.radians(30))
    orbit = orbit_from_state(MU, r, v)
    assert orbit.semi_major_axis == pytest.approx(a, rel=1e-12)
    assert orbit.eccentricity == pytest.approx(e, rel=1e-12)
    assert orbit.inclination == pytest.approx(inc, rel=1e-12)
    # No periapse (ν = 0) o raio é a(1 − e) e a velocidade é perpendicular ao raio.
    r, v = state_from_elements(MU, a, e, inc, raan, argp, 0.0)
    assert np.linalg.norm(r) == pytest.approx(a * (1 - e), rel=1e-14)
    assert abs(r @ v) <= 1e-9 * np.linalg.norm(r) * np.linalg.norm(v)
    # Plano xy, Ω = ω = 0: o periapse fica sobre +x.
    r, _ = state_from_elements(MU, 7e6, 0.1)
    assert r[0] == pytest.approx(6.3e6) and r[1] == 0.0 and r[2] == 0.0


@pytest.mark.parametrize(
    ("a", "e", "nu", "error"),
    [
        (7e6, 1.0, 0.0, InvalidParameterError),  # parábola
        (-7e6, 0.5, 0.0, InvalidParameterError),  # a < 0 com elipse
        (7e6, 1.5, 0.0, InvalidParameterError),  # a > 0 com hipérbole
        (-7e6, 2.0, math.radians(150), DomainError),  # fora do ramo: 1 + e·cos ν ≤ 0
        (7e6, -0.1, 0.0, InvalidParameterError),
    ],
)
def test_elementos_orbitais_invalidos(a, e, nu, error):
    with pytest.raises(error):
        state_from_elements(MU, a, e, true_anomaly=nu)


def _random_elliptic_errors(n, seed=1):
    """Maior |Δr|/a em n órbitas elípticas aleatórias (orientação, e, a e Δt arbitrários)."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n):
        a = 10 ** rng.uniform(6.8, 8.5)
        e = rng.uniform(0.0, 0.97)
        inc, raan, argp = rng.uniform(0, math.pi), rng.uniform(0, 2 * math.pi), rng.uniform(0, 2 * math.pi)
        nu0 = rng.uniform(-math.pi, math.pi)
        r0, v0 = state_from_elements(MU, a, e, inc, raan, argp, nu0)
        ecc0 = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(nu0 / 2))
        period = orbital_period(MU, a)
        dt = rng.uniform(-3 * period, 3 * period)
        mean = math.remainder(ecc0 - e * math.sin(ecc0) + 2 * math.pi / period * dt, 2 * math.pi)
        ecc = mean if e < 0.8 else math.pi
        for _ in range(100):
            step = (ecc - e * math.sin(ecc) - mean) / (1 - e * math.cos(ecc))
            ecc -= step
            if abs(step) < 1e-15:
                break
        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ecc / 2), math.sqrt(1 - e) * math.cos(ecc / 2))
        r_ref, _ = state_from_elements(MU, a, e, inc, raan, argp, nu)
        r, _ = propagate_kepler(MU, r0, v0, dt)
        worst = max(worst, float(np.linalg.norm(r - r_ref)) / a)
    return worst


def _random_hyperbolic_errors(n, seed=2):
    """Maior |Δr|/|r| em n hipérboles aleatórias."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n):
        e = rng.uniform(1.05, 4.0)
        a = -(10 ** rng.uniform(6.8, 8.0))
        inc, raan, argp = rng.uniform(0, math.pi), rng.uniform(0, 2 * math.pi), rng.uniform(0, 2 * math.pi)
        nu0 = rng.uniform(-1, 1) * 0.95 * math.acos(-1 / e)
        r0, v0 = state_from_elements(MU, a, e, inc, raan, argp, nu0)
        big_h0 = 2 * math.atanh(math.sqrt((e - 1) / (e + 1)) * math.tan(nu0 / 2))
        dt = rng.uniform(-2e5, 2e5)
        mean = e * math.sinh(big_h0) - big_h0 + math.sqrt(MU / (-a) ** 3) * dt
        big_h = math.asinh(mean / e)
        for _ in range(200):
            step = (e * math.sinh(big_h) - big_h - mean) / (e * math.cosh(big_h) - 1)
            big_h -= step
            if abs(step) < 1e-15 * max(1.0, abs(big_h)):
                break
        nu = 2 * math.atan(math.sqrt((e + 1) / (e - 1)) * math.tanh(big_h / 2))
        r_ref, _ = state_from_elements(MU, a, e, inc, raan, argp, nu)
        r, _ = propagate_kepler(MU, r0, v0, dt)
        worst = max(worst, float(np.linalg.norm(r - r_ref) / np.linalg.norm(r_ref)))
    return worst


def test_varredura_aleatoria_de_elipses_e_hiperboles():
    # 1000 casos de cada tipo, com semente fixa (valores medidos em docs/VALIDACAO_ETAPA2.md).
    assert _random_elliptic_errors(1000) < 1e-10
    assert _random_hyperbolic_errors(1000) < 1e-11
