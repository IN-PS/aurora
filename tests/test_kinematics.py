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


# --- regressões da auditoria técnica (docs/AUDITORIA_ETAPA1.md) ---------------------------


def test_altura_maxima_e_atingida_uma_vez_no_instante_de_retorno():
    # PHYS-1: x = v0²/(2g) arredondado era declarado inatingível em ~50 % dos casos.
    assert time_to_reach_position(0.0, 20.0, -9.8, 20.408163265306122) == pytest.approx((20.0 / 9.8,), rel=1e-15)


def test_altura_maxima_em_lancamentos_aleatorios():
    import random

    rng = random.Random(5)
    for _ in range(2000):
        v0, g, x0 = rng.uniform(0.1, 500.0), rng.uniform(0.5, 30.0), rng.uniform(-1e4, 1e4)
        apex = x0 + v0 * v0 / (2 * g)
        assert time_to_reach_position(x0, v0, -g, apex) == pytest.approx((v0 / g,), rel=1e-12)
        assert time_to_reach_position(x0, v0, -g, apex + 1e-6 * abs(apex) + 1e-6) == ()  # sem falso positivo
        assert torricelli_speed(v0, -g, v0 * v0 / (2 * g)) == pytest.approx(0.0, abs=1e-6 * v0)


def test_torricelli_com_movel_em_repouso():
    # PHYS-3: com v0 = a = 0 nenhum deslocamento não nulo ocorre (antes retornava |v| = 0).
    with pytest.raises(DomainError, match="repouso"):
        torricelli_speed(v0=0.0, a=0.0, dx=5.0)
    assert torricelli_speed(v0=0.0, a=0.0, dx=0.0) == 0.0


def test_tempo_com_aceleracao_subnormal():
    # NUME-5: a/2 sofria underflow e o erro dizia "a = 0".
    assert time_to_reach_position(0.0, 0.0, 5e-324, 1e-320) == pytest.approx((63.623894882347464,), rel=1e-15)


def test_tempo_nao_representavel_gera_erro():
    with pytest.raises(DomainError, match="representável"):
        time_to_reach_position(-1e308, 1.0, 0.0, 1e308)


def test_estado_vetorial_com_overflow_gera_erro():
    # NUME-6: antes retornava inf/nan silenciosamente.
    with pytest.raises(DomainError):
        constant_acceleration_state([0, 0, 0], [1e200, 0, 0], [0, 0, 0], 1e200)


# --- regressões da revisão das correções ----------------------------------------------------


@pytest.mark.parametrize(
    "args, expected",
    [
        ((1e9, 0.0, -1e-6, 999999999.999999), (1.3810679320049757,)),  # cai 1 µm: atingido
        ((1e9, 0.0, -1e-6, 1000000000.000001), ()),  # alvo acima de corpo em repouso que desce
        ((1e9, -1e-9, -1e-6, 999999999.999999), (1.3800682940436002,)),
        ((10**15, 1, -1, 10**15 + 2), ()),  # 1,5 m acima do ponto de retorno
        ((10**15, 1, -1, 10**15 - 1), (2.732050807568877,)),
        ((6378137.0, 0.0, -9.81, 6378137.0 - 1e-8), (4.570112613295919e-05,)),
        ((6.4e6, 1e-4, 9.8, 6.4e6), (0.0,)),  # já está no alvo em t = 0
    ],
)
def test_tolerancia_do_ponto_de_retorno_nao_depende_da_origem(args, expected):
    # A primeira versão da correção usava uma faixa proporcional a |x0| e errava estes casos.
    assert time_to_reach_position(*args) == pytest.approx(expected, rel=1e-12)


def test_raiz_negativa_fora_do_intervalo_nao_esconde_a_resposta():
    # Raízes ≈ 1 s e ≈ −2e310 s: a segunda é descartada antes da verificação de finitude.
    assert time_to_reach_position(0.0, 1.0, 1e-310, 1.0) == pytest.approx((1.0,), rel=1e-15)
