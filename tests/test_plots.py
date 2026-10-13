"""Testes dos gráficos (Matplotlib e Plotly são opcionais: extra ``[graficos]``)."""

import builtins
import json
import re

import numpy as np
import pytest

from aurora.core.provenance import DataOrigin
from aurora.data.bodies import EARTH
from aurora.errors import DependencyError, InvalidParameterError
from aurora.interface.plots import FigureData, Marker, Segment, figure_for_trajectory, save_figure
from aurora.physics.gravitation import TwoBodyGravity
from aurora.simulation.initial_conditions import circular_orbit_state
from aurora.simulation.propagator import propagate


@pytest.fixture(scope="module")
def data():
    r, v = circular_orbit_state(EARTH, 400e3, 0.9)
    tr = propagate(TwoBodyGravity(EARTH.mu), r, v, 5400.0, n_points=60)
    return figure_for_trajectory(
        "Teste <órbita> & cia",
        EARTH,
        [("trajetória", tr)],
        markers=[Marker("início", 0.0, tuple(tr.positions[0]))],
        reference_orbits=[("referência", EARTH.radius + 400e3)],
        position_error=np.linspace(0, 1e-3, 60),
        energy_error=np.linspace(0, 1e-10, 60),
        summary=[("Linha", "valor <b>"), ("__DATA__", "__PLOTLY__")],
    )


def test_figura_estatica_png_e_svg(data, tmp_path):
    pytest.importorskip("matplotlib")
    png = save_figure(tmp_path / "orbita.png", data)
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    svg = save_figure(tmp_path / "orbita.svg", data).read_text(encoding="utf-8")
    assert "Natureza do resultado" in svg and "simulação numérica" in svg


def test_pagina_html_autocontida(data, tmp_path):
    pytest.importorskip("plotly")
    page = save_figure(tmp_path / "orbita.html", data).read_text(encoding="utf-8")
    assert "<script src=" not in page  # Plotly embutido: funciona sem internet
    assert "Teste &lt;órbita&gt; &amp; cia" in page and "valor &lt;b&gt;" in page
    assert DataOrigin.SIMULATED.disclaimer in page
    # Substituição em uma passagem: o texto "__DATA__" de um resumo não é reinterpretado.
    assert "<th scope=\"row\">__DATA__</th><td>__PLOTLY__</td>" in page
    payload = json.loads(re.search(r"const DATA = (\{.*?\});\n", page).group(1))
    assert payload["segments"][0]["label"] == "trajetória" and len(payload["segments"][0]["x"]) == 60
    assert page.count("prefers-color-scheme: dark") >= 1


def test_formato_nao_suportado(data, tmp_path):
    with pytest.raises(InvalidParameterError, match="não suportado"):
        save_figure(tmp_path / "orbita.jpg", data)
    with pytest.raises(InvalidParameterError, match="sem extensão"):
        save_figure(tmp_path / "orbita", data)


def test_numero_de_segmentos(data):
    tr = data.segments[0].trajectory
    with pytest.raises(InvalidParameterError):
        FigureData("x", EARTH, ())
    with pytest.raises(InvalidParameterError, match="no máximo dois"):
        FigureData("x", EARTH, (Segment("a", tr),) * 3)


def test_biblioteca_ausente_gera_erro_com_instrucoes(data, tmp_path, monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith(("matplotlib", "plotly")):
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(DependencyError, match="pkg install matplotlib"):
        save_figure(tmp_path / "orbita.png", data)
    with pytest.raises(DependencyError, match="plotly"):
        save_figure(tmp_path / "orbita.html", data)
    assert not (tmp_path / "orbita.png").exists() and not (tmp_path / "orbita.html").exists()
