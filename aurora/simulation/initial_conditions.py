"""Condições iniciais idealizadas para a simulação.

O "lançamento" é representado apenas pelo estado do veículo no fim da queima do lançador
(*burnout*): altitude, velocidade e ângulo de trajetória em um plano orbital dado. Não são
modelados a ascensão atmosférica, o empuxo, o arrasto, a rotação do corpo central nem o
local de lançamento — trata-se de uma condição inicial idealizada, não de um perfil de voo.

Referencial: inercial, centrado no corpo central, com o plano xy como plano de referência
(equador do corpo). Ângulos em radianos.
"""

from __future__ import annotations

import math

import numpy as np

from ..core.validation import as_non_negative, as_real
from ..data.bodies import CelestialBody
from ..errors import DomainError
from ..physics.kepler import state_from_elements
from ..physics.orbits import circular_orbital_speed


def _rotation(inclination: float, raan: float) -> np.ndarray:
    """Matriz que leva o plano xy para o plano orbital (Ω em torno de z, depois i em torno da linha dos nodos)."""
    c_o, s_o = math.cos(raan), math.sin(raan)
    c_i, s_i = math.cos(inclination), math.sin(inclination)
    return np.array(((c_o, -s_o * c_i, s_o * s_i), (s_o, c_o * c_i, -c_o * s_i), (0.0, s_i, c_i)))


def burnout_state(
    body: CelestialBody,
    altitude: float,
    speed: float,
    flight_path_angle: float = 0.0,
    inclination: float = 0.0,
    raan: float = 0.0,
    argument_of_latitude: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Estado no fim da queima: altitude (m), velocidade (m/s) e ângulo de trajetória γ (rad).

    γ = 0 é velocidade horizontal (perpendicular ao raio); γ > 0 aponta para cima. O
    movimento é progrado no plano orbital definido por ``inclination`` e ``raan``;
    ``argument_of_latitude`` posiciona o veículo nesse plano, a partir do nodo ascendente.
    """
    altitude = as_non_negative("altitude", altitude)
    speed = as_non_negative("speed", speed)
    gamma = as_real("flight_path_angle", flight_path_angle)
    if not -math.pi / 2 <= gamma <= math.pi / 2:
        raise DomainError("O ângulo de trajetória deve estar entre −90° e 90°.")
    inc, big_omega, u = (as_real(n, v) for n, v in (("inclination", inclination), ("raan", raan),
                                                       ("argument_of_latitude", argument_of_latitude)))
    r = body.radius + altitude
    radial = np.array((math.cos(u), math.sin(u), 0.0))
    transverse = np.array((-math.sin(u), math.cos(u), 0.0))
    rotation = _rotation(inc, big_omega)
    position = rotation @ (r * radial)
    velocity = rotation @ (speed * (math.sin(gamma) * radial + math.cos(gamma) * transverse))
    return position, velocity


def circular_orbit_state(
    body: CelestialBody, altitude: float, inclination: float = 0.0, raan: float = 0.0, argument_of_latitude: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Estado de uma órbita circular de altitude dada (velocidade circular, γ = 0)."""
    altitude = as_non_negative("altitude", altitude)
    speed = circular_orbital_speed(body.mu, body.radius + altitude)
    return burnout_state(body, altitude, speed, 0.0, inclination, raan, argument_of_latitude)


def elliptic_orbit_state(
    body: CelestialBody,
    periapsis_altitude: float,
    apoapsis_altitude: float,
    inclination: float = 0.0,
    raan: float = 0.0,
    arg_periapsis: float = 0.0,
    true_anomaly: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Estado em uma órbita elíptica definida pelas altitudes de periapse e apoapse (m)."""
    hp = as_non_negative("periapsis_altitude", periapsis_altitude)
    ha = as_non_negative("apoapsis_altitude", apoapsis_altitude)
    if ha < hp:
        raise DomainError("A altitude de apoapse deve ser maior ou igual à de periapse.")
    rp, ra = body.radius + hp, body.radius + ha
    a, e = 0.5 * (rp + ra), (ra - rp) / (ra + rp)
    return state_from_elements(body.mu, a, e, inclination, raan, arg_periapsis, true_anomaly)
