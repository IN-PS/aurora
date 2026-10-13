"""Testes de cinemática (movimento com aceleração constante)."""

import numpy as np
import pytest

from aurora.data.constants import G0
from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.kinematics import (
    average_acceleration,
    constant_acceleration_state,
    displacement_constant_acceleration,
    position_constant_acceleration,
    time_to_reach_position,
    torricelli_speed,
    velocity_constant_acceleration,
)


def test_posicao_velocidade_e_deslocamento():
    assert position_constant_acceleration(x0=5.0, v0=2.0, a=4.0, t=3.0) == 5.0 + 6.0 + 18.0
    assert velocity_constant_acceleration(v0=2.0, a=4.0, t=3.0) == 14.0
    assert displacement_constant_acceleration(v0=2.0, a=4.0, t=3.0) == 24.0


def test_movimento_uniforme_e_caso_particular():
    assert position_constant_acceleration(1.0, 3.0, 0.0, 10.0) == 31.0
    assert velocity_constant_acceleration(3.0, 0.0, 10.0) == 3.0


def test_consistencia_entre_formulas():
    v0, a, t = 12.5, -3.2, 7.0
    dx = displacement_constant_acceleration(v0, a, t)
    v = velocity_constant_acceleration(v0, a, t)
    assert average_acceleration(v0, v, t) == pytest.approx(a, rel=1e-14)
    assert torricelli_speed(v0, a, dx) == pytest.approx(abs(v), rel=1e-12)
    # Δx = (v0 + v)/2 · t no MUV
    assert dx == pytest.approx((v0 + v) / 2 * t, rel=1e-14)


def test_torricelli_fora_do_dominio():
    # Lançado a 10 m/s para cima, nunca sobe 100 m.
    with pytest.raises(DomainError, match="nunca é atingido"):
        torricelli_speed(v0=10.0, a=-G0, dx=100.0)


def test_aceleracao_media_com_intervalo_nulo():
    with pytest.raises(InvalidParameterError):
        average_acceleration(0.0, 1.0, 0.0)


def test_estado_vetorial_3d():
    r, v = constant_acceleration_state(r0=[0, 0, 0], v0=[10, 0, 20], a=[0, 0, -G0], t=2.0)
    assert r == pytest.approx([20.0, 0.0, 40.0 - 2.0 * G0])
    assert v == pytest.approx([10.0, 0.0, 20.0 - 2.0 * G0])
    assert isinstance(r, np.ndarray)


def test_estado_vetorial_equivale_ao_escalar_por_componente():
    r0, v0, a, t = [1.0, -2.0, 3.0], [0.5, 4.0, -1.0], [0.1, -0.2, 9.0], 3.7
    r, v = constant_acceleration_state(r0, v0, a, t)
    for i in range(3):
        assert r[i] == position_constant_acceleration(r0[i], v0[i], a[i], t)
        assert v[i] == velocity_constant_acceleration(v0[i], a[i], t)


def test_estado_vetorial_rejeita_dimensao_errada():
    with pytest.raises(InvalidParameterError):
        constant_acceleration_state([0, 0], [0, 0, 0], [0, 0, 0], 1.0)


def test_tempo_de_voo_de_projetil_vertical():
    # Sobe e volta à origem: t = 0 e t = 2·v0/g.
    times = time_to_reach_position(x0=0.0, v0=10.0, a=-G0, x=0.0)
    assert times[0] == 0.0
    assert times[1] == pytest.approx(20.0 / G0, rel=1e-14)


def test_tempo_para_atingir_altura_intermediaria_e_dupla():
    v0 = 20.0
    h = 10.0
    times = time_to_reach_position(0.0, v0, -G0, h)
    assert len(times) == 2
    for t in times:
        assert position_constant_acceleration(0.0, v0, -G0, t) == pytest.approx(h, rel=1e-12)


def test_tempo_posicao_inatingivel_e_somente_futuro():
    assert time_to_reach_position(0.0, 10.0, -G0, 100.0) == ()
    assert time_to_reach_position(0.0, 5.0, 0.0, -1.0) == ()  # ficou no passado
    assert time_to_reach_position(0.0, 5.0, 0.0, 10.0) == (2.0,)
    assert time_to_reach_position(0.0, 0.0, 0.0, 1.0) == ()


def test_tempo_infinitas_solucoes():
    with pytest.raises(DomainError, match="infinitas"):
        time_to_reach_position(1.0, 0.0, 0.0, 1.0)
