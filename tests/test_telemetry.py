"""Testes da telemetria simulada (tabela derivada e exportação CSV/JSON)."""

import json

import numpy as np
import pytest

from aurora import __version__
from aurora.core.provenance import DataOrigin
from aurora.data.bodies import EARTH
from aurora.errors import InvalidParameterError
from aurora.physics.gravitation import TwoBodyGravity
from aurora.simulation.initial_conditions import elliptic_orbit_state
from aurora.simulation.propagator import propagate
from aurora.simulation.telemetry import COLUMNS, read_telemetry_csv, telemetry_table, write_telemetry


@pytest.fixture(scope="module")
def trajectory():
    r, v = elliptic_orbit_state(EARTH, 300e3, 2000e3)
    return propagate(TwoBodyGravity(EARTH.mu), r, v, 3600.0, n_points=25, t0=500.0)


def test_tabela_derivada(trajectory):
    table = telemetry_table(trajectory, EARTH)
    assert table.shape == (25, len(COLUMNS))
    names = [c[0] for c in COLUMNS]
    col = {n: table[:, k] for k, n in enumerate(names)}
    assert col["t"][0] == 0.0 and col["t"][-1] == pytest.approx(3600.0)  # tempo desde o início
    assert np.allclose(col["raio"], np.linalg.norm(trajectory.positions, axis=1), rtol=1e-15)
    assert np.allclose(col["altitude"], col["raio"] - EARTH.radius, rtol=1e-15)
    assert col["altitude"][0] == pytest.approx(300e3, abs=1e-6)
    # Energia e momento angular são constantes do modelo (variação = erro numérico).
    assert np.ptp(col["energia_especifica"]) < 1e-8 * abs(col["energia_especifica"][0])
    assert np.ptp(col["momento_angular"]) < 1e-8 * col["momento_angular"][0]


def test_csv_ida_e_volta_exata(trajectory, tmp_path):
    path = write_telemetry(tmp_path / "tel.csv", trajectory, EARTH, {"experimento": "teste"})
    header, data, meta = read_telemetry_csv(path)
    assert header[0] == "t_s" and header[4] == "vx_m_s" and header[-1] == "momento_angular_m2_s"
    assert np.array_equal(data, telemetry_table(trajectory, EARTH))  # repr preserva todos os bits
    assert meta["origem"] == DataOrigin.SIMULATED.value
    assert meta["aviso"] == DataOrigin.SIMULATED.disclaimer
    assert meta["experimento"] == "teste" and meta["versao_aurora"] == __version__
    assert json.loads(meta["integrador"])["metodo"] == "DOP853"
    assert "dois corpos" in meta["hipoteses"].lower()


def test_json_completo_e_valido(trajectory, tmp_path):
    path = write_telemetry(tmp_path / "tel.JSON", trajectory, EARTH)
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["origem"] == "simulacao" and "não representa" in payload["aviso"]
    assert [c["nome"] for c in payload["colunas"]] == [c[0] for c in COLUMNS]
    assert np.array_equal(np.array(payload["dados"]), telemetry_table(trajectory, EARTH))
    assert payload["corpo"]["nome"] == EARTH.name and payload["corpo"]["mu_m3_s2"] == EARTH.mu
    assert payload["encerramento"] == "duração completa"


def test_extensao_invalida(trajectory, tmp_path):
    with pytest.raises(InvalidParameterError, match=".csv ou .json"):
        write_telemetry(tmp_path / "tel.txt", trajectory, EARTH)
    assert not (tmp_path / "tel.txt").exists()
