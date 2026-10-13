"""Testes dos comandos de simulação da linha de comando: propagar, hohmann e validar."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from aurora.core.provenance import DataOrigin
from aurora.interface.cli import main
from aurora.simulation.telemetry import read_telemetry_csv

ROOT = Path(__file__).resolve().parents[1]
FOOTER = f"Natureza do resultado: {DataOrigin.SIMULATED.disclaimer}"


def run(capsys, *argv):
    code = main(list(argv))
    out, err = capsys.readouterr()
    return code, out, err


def run_json(capsys, *argv):
    code, out, err = run(capsys, *argv, "--json")
    assert code == 0, err
    return json.loads(out)


# --------------------------------------------------------------------------------------
# propagar
# --------------------------------------------------------------------------------------


def test_propagar_texto(capsys):
    code, out, err = run(capsys, "propagar", "--altitude", "400 km", "--inclinacao", "51.6", "--orbitas", "2")
    assert code == 0 and err == ""
    assert "Simulação numérica de órbita — Terra" in out
    assert "inclinação 51.6°" in out
    assert "Verificação do erro numérico" in out and "contra a solução analítica de Kepler" in out
    assert out.rstrip().endswith(FOOTER)


def test_propagar_json_com_erro_medido(capsys):
    data = run_json(capsys, "propagar", "--altitude", "400 km", "--orbitas", "3", "--pontos", "50")
    assert data["origem"] == "simulacao" and data["aviso"] == DataOrigin.SIMULATED.disclaimer
    assert data["amostras"] == 50 and data["encerramento"] == "duração completa"
    check = data["verificacao_kepler"]
    assert check["maior_erro_posicao_m"] < 0.05
    assert check["maior_variacao_relativa_energia"] < 1e-8
    assert data["orbita_inicial"]["tipo"] == "circular"
    assert data["integrador"]["metodo"] == "DOP853" and data["integrador"]["rtol"] == 1e-10


def test_propagar_elipse_e_vetor_de_estado(capsys):
    data = run_json(capsys, "propagar", "--periapse", "200 km", "--apoapse", "35786 km", "--orbitas", "1")
    assert data["orbita_inicial"]["excentricidade"] == pytest.approx(0.73008, abs=1e-5)
    assert data["verificacao_kepler"]["maior_erro_posicao_m"] < 1.0
    data = run_json(
        capsys, "propagar", "--posicao", "7000 0 0 km", "--velocidade", "0 12 0 km/s", "--duracao", "1 d"
    )
    assert data["orbita_inicial"]["tipo"] == "hiperbolica"
    assert data["condicao_inicial"]["velocidade_m_s"] == [0.0, 12000.0, 0.0]
    assert data["duracao_s"] == 86400.0


def test_propagar_lua_com_outro_metodo(capsys):
    data = run_json(
        capsys, "propagar", "--corpo", "lua", "--altitude", "100 km", "--orbitas", "1", "--metodo", "RK45",
        "--rtol", "1e-9",
    )
    assert data["corpo"] == "lua" and data["integrador"]["metodo"] == "RK45"


def test_propagar_colisao(capsys):
    data = run_json(capsys, "propagar", "--posicao", "7000 0 0 km", "--velocidade", "0 2 0 km/s", "--duracao", "3 h")
    assert data["encerramento"] == "colisão com a superfície do corpo central"
    assert data["estado_final"]["altitude_m"] == pytest.approx(0.0, abs=1e-3)
    assert data["duracao_s"] < 3 * 3600


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["--orbitas", "1"], "exatamente uma condição inicial"),
        (["--altitude", "400 km", "--periapse", "200 km", "--orbitas", "1"], "exatamente uma condição inicial"),
        (["--periapse", "200 km", "--orbitas", "1"], "juntas"),
        (["--posicao", "7000 0 0 km", "--duracao", "1 h"], "juntas"),
        (["--posicao", "7000 0 km", "--velocidade", "0 8 0 km/s", "--duracao", "1 h"], "três componentes"),
        (["--posicao", "7000 0 0 s", "--velocidade", "0 8 0 km/s", "--duracao", "1 h"], "esperava-se"),
        (["--posicao", "7000 0 0 km", "--velocidade", "0 8 0 km/s", "--inclinacao", "10", "--duracao", "1 h"],
         "não se aplica"),
        (["--posicao", "6000 0 0 km", "--velocidade", "0 8 0 km/s", "--duracao", "1 h"], "dentro"),
        (["--posicao", "7000 0 0 km", "--velocidade", "0 12 0 km/s", "--orbitas", "1"], "órbita fechada"),
        (["--altitude", "400 km", "--orbitas", "0"], "positivo"),
        (["--altitude", "-10 km", "--orbitas", "1"], "altitude"),
        (["--altitude", "400 km", "--duracao", "-1 h"], "duration"),
        (["--altitude", "400 km", "--orbitas", "1", "--rtol", "0.5"], "rtol"),
        (["--altitude", "400 km", "--orbitas", "50", "--max-avaliacoes", "500"], "Limite de 500"),
        (["--altitude", "400 km", "--orbitas", "1", "--telemetria", "x.txt"], ".csv ou .json"),
        (["--altitude", "400 km", "--orbitas", "1", "--grafico", "x.jpg"], "não suportado"),
    ],
)
def test_propagar_erros_de_entrada(capsys, argv, message):
    code, out, err = run(capsys, "propagar", *argv)
    assert code == 1 and out == ""
    assert err.startswith("aurora: erro:") and message in err


def test_propagar_erros_de_uso(capsys):
    for argv in (["--altitude", "400 km"], ["--altitude", "400 km", "--orbitas", "1", "--duracao", "1 h"],
                 ["--altitude", "400 km", "--orbitas", "1", "--pontos", "1"],
                 ["--altitude", "400 km", "--orbitas", "1", "--metodo", "Euler"]):
        with pytest.raises(SystemExit) as exc:
            main(["propagar", *argv])
        assert exc.value.code == 2
    capsys.readouterr()


def test_propagar_grava_telemetria(capsys, tmp_path):
    csv_path, json_path = tmp_path / "leo.csv", tmp_path / "leo.json"
    code, out, _ = run(capsys, "propagar", "--altitude", "400 km", "--orbitas", "1", "--pontos", "30",
                       "--telemetria", str(csv_path))
    assert code == 0 and f"telemetria: {csv_path}" in out
    header, rows, meta = read_telemetry_csv(csv_path)
    assert rows.shape == (30, len(header)) and meta["origem"] == "simulacao"
    assert "órbita circular de altitude" in meta["condicao_inicial"]
    data = run_json(capsys, "propagar", "--altitude", "400 km", "--orbitas", "1", "--telemetria", str(json_path))
    assert data["arquivos"] == [{"tipo": "telemetria", "caminho": str(json_path)}]
    assert json.loads(json_path.read_text(encoding="utf-8"))["origem"] == "simulacao"


def test_propagar_grava_graficos(capsys, tmp_path):
    pytest.importorskip("matplotlib")
    pytest.importorskip("plotly")
    png, page = tmp_path / "o.png", tmp_path / "o.html"
    code, out, _ = run(capsys, "propagar", "--altitude", "400 km", "--orbitas", "1", "--pontos", "40",
                       "--grafico", str(png), "--grafico", str(page))
    assert code == 0 and png.stat().st_size > 0 and page.stat().st_size > 0
    assert f"gráfico: {png}" in out and f"gráfico: {page}" in out


# --------------------------------------------------------------------------------------
# hohmann
# --------------------------------------------------------------------------------------


def test_hohmann_analitico(capsys):
    code, out, _ = run(capsys, "hohmann", "--de", "400 km", "--para", "35786 km")
    assert code == 0
    assert "Transferência de Hohmann — Terra (cálculo analítico)" in out
    assert "Δv total = 3853.959216 m/s" in out
    assert out.rstrip().endswith(f"Natureza do resultado: {DataOrigin.ANALYTIC.disclaimer}")
    data = run_json(capsys, "hohmann", "--de", "400 km", "--para", "35786 km")
    assert data["origem"] == DataOrigin.ANALYTIC.value
    assert data["transferencia"]["delta_v_total_m_s"] == pytest.approx(3853.96, abs=0.01)
    assert data["transferencia"]["tempo_transferencia_s"] == pytest.approx(19048.56, abs=0.01)


def test_hohmann_descida_explica_sinal(capsys):
    code, out, _ = run(capsys, "hohmann", "--de", "1000 km", "--para", "300 km")
    assert code == 0 and "queimas retrógradas" in out


def test_hohmann_simulado(capsys, tmp_path):
    csv_path = tmp_path / "missao.csv"
    data = run_json(capsys, "hohmann", "--de", "400 km", "--para", "35786 km", "--simular",
                    "--isp", "320 s", "--pontos", "100", "--telemetria", str(csv_path))
    sim = data["simulacao"]
    assert data["origem"] == "simulacao"
    assert abs(sim["erro_raio_chegada_m"]) < 1.0 and abs(sim["erro_semieixo_final_m"]) < 1.0
    assert sim["orbita_final"]["excentricidade"] < 1e-7
    assert [e["tempo_s"] for e in sim["eventos"]] == [0.0, pytest.approx(19048.56, abs=0.01)]
    # Padrão: exatamente o propelente necessário; a massa final volta à massa seca.
    assert sim["veiculo"]["propelente_consumido_kg"] == pytest.approx(sim["veiculo"]["propelente_inicial_kg"])
    assert sim["eventos"][-1]["massa_depois_kg"] == pytest.approx(1000.0)
    _, rows, meta = read_telemetry_csv(csv_path)
    assert rows.shape[0] == 200 and "1ª queima" in meta["eventos"]


def test_hohmann_simulado_texto(capsys):
    code, out, _ = run(capsys, "hohmann", "--de", "400 km", "--para", "1000 km", "--simular", "--pontos", "20",
                       "--massa-seca", "500 kg", "--propelente", "200 kg")
    assert code == 0
    assert "Missão simulada" in out and "massa seca 500 kg" in out and "Propelente consumido" in out
    assert out.rstrip().endswith(FOOTER)


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["--de", "400 km", "--para", "35786 km", "--simular", "--propelente", "100 kg"], "Propelente insuficiente"),
        (["--de", "400 km", "--para", "35786 km", "--simular", "--isp", "0 s"], "isp"),
        (["--de", "400 km", "--para", "35786 km", "--simular", "--isp", "1e-300 s"], "não representável"),
        (["--de", "400 km", "--para", "35786 km", "--simular", "--massa-seca", "0 kg"], "dry_mass"),
        (["--de", "400 km", "--para", "35786 km", "--telemetria", "x.csv"], "exigem --simular"),
        (["--de", "-400 km", "--para", "35786 km"], "não negativas"),
        (["--de", "400 km", "--para", "400 km", "--simular"], "iguais"),
        (["--de", "400 km", "--para", "35786 s"], "unidade"),
    ],
)
def test_hohmann_erros(capsys, argv, message):
    code, out, err = run(capsys, "hohmann", *argv)
    assert code == 1 and out == ""
    assert message in err


# --------------------------------------------------------------------------------------
# validar
# --------------------------------------------------------------------------------------


def test_validar_json_e_texto(capsys):
    data = run_json(capsys, "validar", "--cenario", "leo", "--cenario", "hiperbolica")
    assert [c["chave"] for c in data["cenarios"]] == ["leo", "hiperbolica"]
    assert data["cenarios"][0]["maior_erro_posicao_m"] < 0.1
    assert data["origem"] == "simulacao"
    code, out, _ = run(capsys, "validar", "--cenario", "lua", "--rtol", "1e-8")
    assert code == 0 and "rtol = 1e-08" in out and "lua:" in out
    assert "mesmo modelo de dois corpos" in out and out.rstrip().endswith(FOOTER)


def test_validar_cenario_desconhecido_e_erro_de_uso(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["validar", "--cenario", "marte"])
    assert exc.value.code == 2
    capsys.readouterr()


# --------------------------------------------------------------------------------------
# Saída canalizada
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("argv", [["validar", "--cenario", "leo"], ["formulas"], ["quadratica", "1", "-3", "2"]])
def test_leitor_que_fecha_o_pipe_nao_gera_rastreamento(argv):
    proc = subprocess.Popen(
        [sys.executable, "-m", "aurora", *argv], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    proc.stdout.close()  # como "aurora ... | head" quando o head já terminou
    err = proc.stderr.read().decode("utf-8", "replace")
    proc.stderr.close()
    assert proc.wait(timeout=60) == 1
    assert "Traceback" not in err and "BrokenPipeError" not in err and "Exception ignored" not in err
