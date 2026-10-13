"""Missões simuladas: transferência de Hohmann com queimas impulsivas.

A missão parte de uma órbita circular equatorial de raio r1, faz a 1ª queima tangencial,
percorre meia elipse de transferência (integração numérica), faz a 2ª queima no raio r2 e
segue na órbita final. As queimas são impulsivas e o propelente é descontado pela equação
de Tsiolkovsky. A órbita alcançada é comparada com a solução analítica: as diferenças
medem o erro numérico da integração (o modelo é o mesmo nos dois cálculos).

Resultados idealizados de simulação; não constituem plano de voo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ..core.validation import as_non_negative, as_positive
from ..data.bodies import CelestialBody
from ..errors import DomainError
from ..physics.gravitation import TwoBodyGravity
from ..physics.maneuvers import HohmannTransfer, hohmann_transfer
from ..physics.orbits import TwoBodyOrbit, orbit_from_state, orbital_period
from .propagator import IntegratorSettings, Termination, Trajectory, concatenate, propagate
from .vehicle import Spacecraft


@dataclass(frozen=True)
class MissionEvent:
    """Evento de uma missão simulada (tempo em s, Δv em m/s com sinal, massas em kg)."""

    time: float
    name: str
    delta_v: float
    propellant_used: float
    mass_after: float


@dataclass(frozen=True)
class HohmannMission:
    """Resultado da missão de Hohmann simulada."""

    body: CelestialBody
    transfer: HohmannTransfer
    trajectory: Trajectory
    segments: tuple[Trajectory, ...]
    events: tuple[MissionEvent, ...]
    final_orbit: TwoBodyOrbit
    spacecraft_initial: Spacecraft
    spacecraft_final: Spacecraft
    arrival_radius_error: float
    """|r| na chegada (antes da 2ª queima) − r2, em m: erro numérico da integração do arco."""
    final_semi_major_axis_error: float
    """a da órbita final − r2, em m."""

    @property
    def propellant_used(self) -> float:
        return self.spacecraft_initial.propellant_mass - self.spacecraft_final.propellant_mass


def _tangential_burn(velocity: np.ndarray, delta_v: float) -> np.ndarray:
    return velocity + delta_v * velocity / math.hypot(*velocity)


def simulate_hohmann(
    body: CelestialBody,
    r1: float,
    r2: float,
    *,
    spacecraft: Spacecraft | None = None,
    settings: IntegratorSettings | None = None,
    final_coast_orbits: float = 1.0,
    n_points: int = 300,
) -> HohmannMission:
    """Simula a transferência de Hohmann entre órbitas circulares de raios ``r1`` e ``r2`` (m)."""
    r1, r2 = as_positive("r1", r1), as_positive("r2", r2)
    for name, radius in (("r1", r1), ("r2", r2)):
        if radius <= body.radius:
            raise DomainError(f"O raio '{name}' ({radius:.6g} m) não está acima da superfície de {body.name}.")
    if r1 == r2:
        raise DomainError("As órbitas inicial e final são iguais: não há transferência a simular.")
    final_coast_orbits = as_non_negative("final_coast_orbits", final_coast_orbits)
    settings = settings if settings is not None else IntegratorSettings()
    craft = spacecraft if spacecraft is not None else Spacecraft(propellant_mass=1e9)
    transfer = hohmann_transfer(body.mu, r1, r2)
    model = TwoBodyGravity(body.mu)

    position = np.array((r1, 0.0, 0.0))
    velocity = np.array((0.0, transfer.speed_initial_circular, 0.0))
    events: list[MissionEvent] = []

    craft_after_1, used_1 = craft.burn(transfer.delta_v1)
    velocity = _tangential_burn(velocity, transfer.delta_v1)
    events.append(
        MissionEvent(0.0, "1ª queima (início da transferência)", transfer.delta_v1, used_1, craft_after_1.mass)
    )

    arc = propagate(
        model, position, velocity, transfer.transfer_time, settings=settings, n_points=n_points,
        collision_radius=body.radius,
    )
    if arc.termination is Termination.COLLISION:
        raise DomainError("A elipse de transferência atinge a superfície do corpo central.")
    arrival_position, arrival_velocity = arc.final_position, arc.final_velocity
    arrival_radius_error = math.hypot(*arrival_position) - r2

    craft_after_2, used_2 = craft_after_1.burn(transfer.delta_v2)
    final_velocity = _tangential_burn(arrival_velocity, transfer.delta_v2)
    events.append(
        MissionEvent(
            transfer.transfer_time, "2ª queima (circularização)", transfer.delta_v2, used_2, craft_after_2.mass
        )
    )
    final_orbit = orbit_from_state(body.mu, arrival_position, final_velocity)

    segments = [arc]
    if final_coast_orbits > 0.0:
        coast = propagate(
            model, arrival_position, final_velocity, final_coast_orbits * orbital_period(body.mu, r2),
            settings=settings, n_points=n_points, t0=transfer.transfer_time, collision_radius=body.radius,
        )
        segments.append(coast)

    return HohmannMission(
        body=body,
        transfer=transfer,
        trajectory=concatenate(segments),
        segments=tuple(segments),
        events=tuple(events),
        final_orbit=final_orbit,
        spacecraft_initial=craft,
        spacecraft_final=craft_after_2,
        arrival_radius_error=arrival_radius_error,
        final_semi_major_axis_error=final_orbit.semi_major_axis - r2,
    )
