"""Testes do veículo virtual (Tsiolkovsky) e da missão de Hohmann simulada."""

import math

import numpy as np
import pytest

from aurora.core.provenance import DataOrigin
from aurora.data.bodies import EARTH, MOON
from aurora.data.constants import G0
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.maneuvers import hohmann_transfer, rocket_delta_v
from aurora.simulation.missions import simulate_hohmann
from aurora.simulation.vehicle import Spacecraft, propellant_required

R_LEO = EARTH.radius + 400e3
R_GEO = EARTH.radius + 35_786e3


def test_veiculo_massas_e_capacidade():
    craft = Spacecraft(dry_mass=1000.0, propellant_mass=1718.28, isp=300.0)
    assert craft.mass == pytest.approx(2718.28)
    assert craft.delta_v_capacity == pytest.approx(rocket_delta_v(300.0, 2718.28, 1000.0), rel=1e-15)
    assert Spacecraft().propellant_mass == 0.0 and Spacecraft().delta_v_capacity == 0.0


def test_queima_desconta_propelente_e_independe_do_sentido():
    craft = Spacecraft(dry_mass=1000.0, propellant_mass=2000.0, isp=320.0)
    after, used = craft.burn(1500.0)
    assert after.propellant_mass == pytest.approx(2000.0 - used)
    assert rocket_delta_v(320.0, craft.mass, after.mass) == pytest.approx(1500.0, rel=1e-12)
    retro, used_retro = craft.burn(-1500.0)
    assert used_retro == used and retro == after
    assert craft.propellant_mass == 2000.0  # imutável


def test_queima_sem_propelente_suficiente():
    craft = Spacecraft(dry_mass=1000.0, propellant_mass=10.0, isp=300.0)
    with pytest.raises(DomainError, match="Propelente insuficiente"):
        craft.burn(500.0)
    # Toda a capacidade pode ser usada (tolerância apenas de arredondamento).
    empty, used = craft.burn(craft.delta_v_capacity)
    assert used == pytest.approx(10.0, rel=1e-12) and empty.propellant_mass == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize(
    "kwargs",
    [{"dry_mass": 0.0}, {"dry_mass": -5.0}, {"propellant_mass": -1.0}, {"isp": 0.0}, {"name": " "}, {"name": 3}],
)
def test_veiculo_invalido(kwargs):
    with pytest.raises(InvalidParameterError):
        Spacecraft(**kwargs)


def test_queima_rejeita_valores_nao_numericos():
    with pytest.raises(InvalidParameterError):
        Spacecraft(propellant_mass=100.0).burn(True)
    with pytest.raises(InvalidParameterError):
        Spacecraft(propellant_mass=100.0).burn(math.nan)


def test_propelente_necessario():
    mp = propellant_required(1000.0, 300.0, 3000.0)
    assert mp == pytest.approx(1000.0 * (math.exp(3000.0 / (300.0 * G0)) - 1.0), rel=1e-14)
    assert Spacecraft(dry_mass=1000.0, propellant_mass=mp, isp=300.0).delta_v_capacity == pytest.approx(3000.0)
    assert propellant_required(1000.0, 300.0, 0.0) == 0.0
    with pytest.raises(InvalidParameterError):
        propellant_required(1000.0, 0.0, 100.0)
    with pytest.raises(DomainError, match="não representável"):
        propellant_required(1000.0, 1e-300, 100.0)


def test_hohmann_simulado_leo_geo():
    craft = Spacecraft(dry_mass=1000.0, propellant_mass=3000.0, isp=320.0)
    mission = simulate_hohmann(EARTH, R_LEO, R_GEO, spacecraft=craft, n_points=200)
    analytic = hohmann_transfer(EARTH.mu, R_LEO, R_GEO)
    assert mission.transfer == analytic
    # O erro numérico da chegada e da órbita final fica bem abaixo de 1 m (medido: ~1 cm).
    assert abs(mission.arrival_radius_error) < 1.0
    assert abs(mission.final_semi_major_axis_error) < 1.0
    assert mission.final_orbit.eccentricity < 1e-7
    # Eventos e consumo de propelente (Tsiolkovsky com Δv total).
    first, second = mission.events
    assert first.time == 0.0 and second.time == analytic.transfer_time
    assert first.delta_v == analytic.delta_v1 and second.delta_v == analytic.delta_v2
    total_mass_ratio = craft.mass / mission.spacecraft_final.mass
    assert 320.0 * G0 * math.log(total_mass_ratio) == pytest.approx(analytic.delta_v_total, rel=1e-12)
    assert mission.propellant_used == pytest.approx(first.propellant_used + second.propellant_used)
    # Trajetória: arco de transferência + uma órbita final; tempos não decrescentes.
    assert len(mission.segments) == 2
    assert mission.trajectory.times.size == 400 and np.all(np.diff(mission.trajectory.times) >= 0)
    assert mission.trajectory.origin is DataOrigin.SIMULATED
    assert np.linalg.norm(mission.segments[0].final_position) == pytest.approx(R_GEO, abs=1.0)


def test_hohmann_veiculo_padrao_usa_exatamente_o_propelente_necessario():
    mission = simulate_hohmann(EARTH, R_LEO, R_GEO, final_coast_orbits=0.0, n_points=50)
    assert len(mission.segments) == 1
    assert mission.spacecraft_initial.dry_mass == 1000.0 and mission.spacecraft_initial.isp == 300.0
    expected = propellant_required(1000.0, 300.0, mission.transfer.delta_v_total)
    assert mission.spacecraft_initial.propellant_mass == expected
    assert mission.spacecraft_final.propellant_mass == pytest.approx(0.0, abs=1e-6)


def test_hohmann_descida_na_lua():
    mission = simulate_hohmann(MOON, MOON.radius + 500e3, MOON.radius + 100e3, n_points=100)
    assert mission.transfer.delta_v1 < 0 and mission.transfer.delta_v2 < 0
    assert abs(mission.arrival_radius_error) < 0.1
    assert mission.final_orbit.semi_major_axis == pytest.approx(MOON.radius + 100e3, abs=0.1)


@pytest.mark.parametrize(
    ("r1", "r2", "kwargs", "match"),
    [
        (R_LEO, R_LEO, {}, "iguais"),
        (EARTH.radius - 1.0, R_GEO, {}, "superfície"),
        (R_LEO, EARTH.radius, {}, "superfície"),
        (R_LEO, R_GEO, {"spacecraft": Spacecraft(propellant_mass=100.0)}, "Propelente insuficiente"),
        (R_LEO, R_GEO, {"final_coast_orbits": -1.0}, "final_coast_orbits"),
    ],
)
def test_hohmann_simulado_entradas_invalidas(r1, r2, kwargs, match):
    with pytest.raises(InvalidParameterError, match=match):
        simulate_hohmann(EARTH, r1, r2, n_points=20, **kwargs)
