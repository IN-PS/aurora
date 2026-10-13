"""Testes do propagador numérico (SciPy) e da validação contra a solução de Kepler."""

import math

import numpy as np
import pytest

from aurora.core.provenance import DataOrigin
from aurora.data.bodies import EARTH, MOON
from aurora.errors import DomainError, InvalidParameterError, SimulationError
from aurora.physics.gravitation import NewtonianGravity, PointMass, TwoBodyGravity
from aurora.physics.orbits import circular_orbital_speed, orbital_period
from aurora.simulation.propagator import IntegratorSettings, Termination, concatenate, propagate
from aurora.simulation.validation import (
    compare_with_kepler,
    conservation_errors,
    kepler_error_series,
    run_validation,
    standard_cases,
)

MU = EARTH.mu
R_LEO = EARTH.radius + 400e3
V_LEO = circular_orbital_speed(MU, R_LEO)
T_LEO = orbital_period(MU, R_LEO)


def leo(duration=T_LEO, **kwargs):
    return propagate(TwoBodyGravity(MU), (R_LEO, 0, 0), (0, V_LEO, 0), duration, **kwargs)


def test_amostras_e_metadados():
    tr = leo(n_points=101)
    assert tr.times.shape == (101,) and tr.positions.shape == (101, 3) and tr.velocities.shape == (101, 3)
    assert tr.times[0] == 0.0 and tr.times[-1] == T_LEO
    assert np.allclose(np.diff(tr.times), T_LEO / 100)
    assert tr.origin is DataOrigin.SIMULATED
    assert tr.termination is Termination.COMPLETED
    assert tr.evaluations > 0 and tr.duration == T_LEO
    assert any("dois corpos" in h.lower() for h in tr.hypotheses)
    assert any("passo adaptativo" in h for h in tr.hypotheses)
    with pytest.raises(ValueError):
        tr.positions[0, 0] = 1.0  # somente leitura


def test_orbita_circular_volta_ao_inicio():
    tr = leo()
    assert np.linalg.norm(tr.final_position - (R_LEO, 0, 0)) < 0.01  # m, após uma órbita
    radii = np.linalg.norm(tr.positions, axis=1)
    assert np.max(np.abs(radii - R_LEO)) < 0.01


def test_erro_contra_kepler_e_conservacao():
    tr = leo(10 * T_LEO, n_points=400)
    comparison = compare_with_kepler(tr, MU)
    assert comparison.max_position_error < 0.05  # m (medido: ~8 mm)
    assert comparison.max_velocity_error < 1e-4  # m/s
    conservation = conservation_errors(tr, MU)
    assert conservation.energy < 1e-8 and conservation.angular_momentum < 1e-8
    dr, dv = kepler_error_series(tr, MU)
    assert dr[0] == 0.0 and dv[0] == 0.0


def test_tolerancia_menor_reduz_o_erro():
    errors = [
        compare_with_kepler(leo(5 * T_LEO, settings=IntegratorSettings(rtol=rtol)), MU).max_position_error
        for rtol in (1e-6, 1e-9, 1e-12)
    ]
    assert errors[0] > 30 * errors[1] > 0 and errors[1] > 3 * errors[2]


def test_rk45_tambem_converge_com_mais_avaliacoes():
    dop = leo(settings=IntegratorSettings(method="DOP853"))
    rk = leo(settings=IntegratorSettings(method="RK45"))
    assert compare_with_kepler(rk, MU).max_position_error < 1.0
    assert rk.evaluations > dop.evaluations


def test_t0_desloca_o_tempo():
    tr = leo(T_LEO / 2, t0=1000.0, n_points=11)
    assert tr.times[0] == 1000.0 and tr.times[-1] == 1000.0 + T_LEO / 2


def test_colisao_com_a_superficie():
    # Velocidade horizontal muito abaixo da circular: a trajetória cruza a superfície.
    tr = propagate(TwoBodyGravity(MU), (R_LEO, 0, 0), (0, 3000.0, 0), 2 * T_LEO, collision_radius=EARTH.radius)
    assert tr.termination is Termination.COLLISION
    assert np.linalg.norm(tr.final_position) == pytest.approx(EARTH.radius, abs=1e-3)
    assert tr.times[-1] < 2 * T_LEO and np.all(np.diff(tr.times) > 0)


def test_posicao_inicial_dentro_do_corpo():
    with pytest.raises(DomainError):
        propagate(TwoBodyGravity(MU), (EARTH.radius, 0, 0), (0, 8000, 0), 100.0, collision_radius=EARTH.radius)


def test_queda_radial_sem_superficie_falha_com_erro_claro():
    # Sem raio de colisão, a queda livre atinge a singularidade r = 0 do modelo.
    with pytest.raises(SimulationError):
        propagate(TwoBodyGravity(MU), (7e6, 0, 0), (0, 0, 0), 3000.0)


def test_limite_de_avaliacoes():
    with pytest.raises(SimulationError, match="Limite de 100 avaliações"):
        leo(100 * T_LEO, settings=IntegratorSettings(max_evaluations=100))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"method": "Euler"},
        {"rtol": 1e-2},
        {"rtol": 1e-15},
        {"rtol": math.nan},
        {"atol_position": 0.0},
        {"atol_velocity": -1.0},
        {"max_step": 0.0},
        {"max_evaluations": 10},
        {"max_evaluations": True},
        {"max_evaluations": 1e6},
    ],
)
def test_configuracao_invalida(kwargs):
    with pytest.raises(InvalidParameterError):
        IntegratorSettings(**kwargs)


def test_configuracao_descrita_em_portugues():
    settings = IntegratorSettings(max_step=60.0)
    assert "DOP853" in settings.describe() and "passo máximo = 60 s" in settings.describe()
    assert settings.to_dict() == {
        "metodo": "DOP853",
        "rtol": 1e-10,
        "atol_posicao_m": 1e-4,
        "atol_velocidade_m_s": 1e-7,
        "passo_maximo_s": 60.0,
        "limite_avaliacoes": 2_000_000,
    }


@pytest.mark.parametrize(
    ("args", "kwargs"),
    [
        ((object(), (R_LEO, 0, 0), (0, V_LEO, 0), 10.0), {}),
        ((TwoBodyGravity(MU), (R_LEO, 0, 0), (0, V_LEO, 0), 0.0), {}),
        ((TwoBodyGravity(MU), (R_LEO, 0, 0), (0, V_LEO, 0), -5.0), {}),
        ((TwoBodyGravity(MU), (R_LEO, 0, 0), (0, V_LEO, 0), 10.0), {"n_points": 1}),
        ((TwoBodyGravity(MU), (R_LEO, 0, 0), (0, V_LEO, 0), 10.0), {"n_points": 10.0}),
        ((TwoBodyGravity(MU), (R_LEO, 0, math.nan), (0, V_LEO, 0), 10.0), {}),
        ((TwoBodyGravity(MU), (R_LEO, 0, 0), (0, V_LEO, 0), 10.0), {"t0": 1e300, "n_points": 3}),
    ],
)
def test_entradas_invalidas(args, kwargs):
    with pytest.raises(InvalidParameterError):
        propagate(*args, **kwargs)


def test_modelo_newtoniano_com_fonte_fixa_equivale_a_dois_corpos():
    # Uma massa pontual fixa na origem com o μ da Terra reproduz o modelo de dois corpos.
    model = NewtonianGravity([PointMass.from_body(EARTH)])
    tr = propagate(model, (R_LEO, 0, 0), (0, V_LEO, 0), T_LEO, n_points=50)
    ref = leo(n_points=50)
    assert np.max(np.linalg.norm(tr.positions - ref.positions, axis=1)) < 1e-3


def test_concatenar_segmentos():
    first = leo(T_LEO / 2, n_points=10)
    second = propagate(
        TwoBodyGravity(MU), first.final_position, first.final_velocity, T_LEO / 2, n_points=10, t0=T_LEO / 2
    )
    joined = concatenate([first, second])
    assert joined.times.size == 20 and joined.evaluations == first.evaluations + second.evaluations
    assert len(set(joined.hypotheses)) == len(joined.hypotheses)
    with pytest.raises(InvalidParameterError):
        concatenate([])


def test_cenarios_padrao_de_validacao():
    results = run_validation(n_points=200)
    assert [r.case.key for r in results] == ["leo", "geo", "gto", "molniya", "lua", "hiperbolica"]
    # Limites com folga de cerca de 10× sobre os valores medidos (docs/VALIDACAO_ETAPA2.md).
    limits = {"leo": 0.1, "geo": 0.2, "gto": 5.0, "molniya": 2.0, "lua": 0.05, "hiperbolica": 0.2}
    for r in results:
        assert r.comparison.max_position_error < limits[r.case.key], r.case.key
        assert r.conservation.energy < 1e-8, r.case.key
        assert r.conservation.angular_momentum < 1e-8, r.case.key
        assert r.samples == 200
    lua = next(r for r in results if r.case.key == "lua")
    assert lua.case.body is MOON


def test_cenario_desconhecido():
    with pytest.raises(InvalidParameterError, match="desconhecido"):
        run_validation(keys=("marte",))
    assert {c.key for c in standard_cases()} >= {"leo", "hiperbolica"}
