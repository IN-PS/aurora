"""Testes das condições iniciais idealizadas (fim de queima, órbitas circular e elíptica)."""

import math

import numpy as np
import pytest

from aurora.data.bodies import EARTH, MOON
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.orbits import circular_orbital_speed, orbit_from_state, vis_viva_speed
from aurora.simulation.initial_conditions import burnout_state, circular_orbit_state, elliptic_orbit_state


def test_orbita_circular_inclinada():
    inc, raan = math.radians(51.6), math.radians(30)
    r, v = circular_orbit_state(EARTH, 400e3, inc, raan, math.radians(10))
    radius = EARTH.radius + 400e3
    assert np.linalg.norm(r) == pytest.approx(radius, rel=1e-15)
    assert np.linalg.norm(v) == pytest.approx(circular_orbital_speed(EARTH.mu, radius), rel=1e-15)
    assert abs(r @ v) <= 1e-9 * np.linalg.norm(r) * np.linalg.norm(v)
    h = np.cross(r, v)
    assert math.acos(h[2] / np.linalg.norm(h)) == pytest.approx(inc, rel=1e-13)
    # Linha dos nodos: direção do nodo ascendente = ẑ × ĥ, a Ω do eixo x.
    node = np.cross((0, 0, 1), h)
    assert math.atan2(node[1], node[0]) == pytest.approx(raan, rel=1e-13)
    orbit = orbit_from_state(EARTH.mu, r, v)
    assert orbit.eccentricity < 1e-12


def test_estado_no_fim_da_queima_com_angulo_de_trajetoria():
    gamma = math.radians(5)
    r, v = burnout_state(MOON, 20e3, 1700.0, gamma)
    assert np.linalg.norm(r) == pytest.approx(MOON.radius + 20e3)
    assert np.linalg.norm(v) == pytest.approx(1700.0)
    # sen γ = r̂·v̂ (componente radial da velocidade).
    assert (r @ v) / (np.linalg.norm(r) * np.linalg.norm(v)) == pytest.approx(math.sin(gamma), rel=1e-13)
    # Sem inclinação, movimento progrado no plano xy.
    assert r[2] == 0.0 and v[2] == 0.0 and np.cross(r, v)[2] > 0


def test_orbita_eliptica_comeca_no_periapse():
    r, v = elliptic_orbit_state(EARTH, 200e3, 35_786e3, math.radians(27))
    rp, ra = EARTH.radius + 200e3, EARTH.radius + 35_786e3
    assert np.linalg.norm(r) == pytest.approx(rp, rel=1e-14)
    assert np.linalg.norm(v) == pytest.approx(vis_viva_speed(EARTH.mu, rp, 0.5 * (rp + ra)), rel=1e-14)
    orbit = orbit_from_state(EARTH.mu, r, v)
    assert orbit.apoapsis_radius == pytest.approx(ra, rel=1e-12)
    # Periapse e apoapse iguais: órbita circular.
    r, v = elliptic_orbit_state(EARTH, 500e3, 500e3)
    assert orbit_from_state(EARTH.mu, r, v).eccentricity < 1e-12


@pytest.mark.parametrize(
    ("func", "args", "error"),
    [
        (burnout_state, (EARTH, -1.0, 7000.0), InvalidParameterError),
        (burnout_state, (EARTH, 100e3, -1.0), InvalidParameterError),
        (burnout_state, (EARTH, 100e3, 7000.0, math.radians(91)), DomainError),
        (burnout_state, (EARTH, 100e3, 7000.0, 0.0, math.nan), InvalidParameterError),
        (circular_orbit_state, (EARTH, math.inf), InvalidParameterError),
        (elliptic_orbit_state, (EARTH, 1000e3, 500e3), DomainError),
    ],
)
def test_entradas_invalidas(func, args, error):
    with pytest.raises(error):
        func(*args)
