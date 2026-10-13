"""Testes da interface de linha de comando."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from aurora.interface.cli import main

ROOT = Path(__file__).resolve().parents[1]


def run(capsys, *argv):
    code = main(list(argv))
    out, err = capsys.readouterr()
    return code, out, err


def test_quadratica_texto(capsys):
    code, out, _ = run(capsys, "quadratica", "1", "-3", "2")
    assert code == 0
    assert "duas raízes reais distintas" in out
    assert "x₁ = 1" in out and "x₂ = 2" in out


def test_quadratica_aceita_negativos_em_notacao_cientifica_e_virgula(capsys):
    code, out, _ = run(capsys, "quadratica", "1", "-1e8", "1", "--json")
    assert code == 0
    data = json.loads(out)
    assert data["raizes"][0]["real"] == pytest.approx(1e-8, rel=1e-15)
    code, out, _ = run(capsys, "quadratica", "2", "-3,5", "0", "--json")
    assert code == 0
    assert json.loads(out)["coeficientes"]["b"] == -3.5


def test_quadratica_complexa_json(capsys):
    code, out, _ = run(capsys, "quadratica", "1", "2", "5", "--json")
    data = json.loads(out)
    assert data["classificacao"] == "complexas_conjugadas"
    assert data["raizes"] == [{"real": -1.0, "imag": 2.0}, {"real": -1.0, "imag": -2.0}]


def test_quadratica_a_nulo_retorna_erro(capsys):
    code, out, err = run(capsys, "quadratica", "0", "1", "1")
    assert code == 1
    assert "erro" in err and "linear" in err
    assert out == ""


def test_quadratica_entrada_nao_numerica(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["quadratica", "1", "abc", "2"])
    assert exc.value.code == 2
    assert "erro" in capsys.readouterr().err


def test_formulas_lista_e_detalhes(capsys):
    code, out, _ = run(capsys, "formulas")
    assert code == 0 and "velocidade_orbital_circular" in out
    code, out, _ = run(capsys, "formulas", "periodo_orbital_circular")
    assert code == 0
    assert "T = 2*pi*sqrt(r^3/mu)" in out and "Hipóteses" in out


def test_calcular_com_unidades(capsys):
    code, out, _ = run(
        capsys, "calcular", "velocidade_orbital_circular", "mu=398600.4418 km^3/s^2", "r=6778.137km", "--json"
    )
    assert code == 0
    data = json.loads(out)
    assert data["entradas"]["r"]["valor"] == pytest.approx(6778137.0)
    assert data["valor"] == pytest.approx((3.986004418e14 / 6778137.0) ** 0.5, rel=1e-14)
    assert data["origem"] == "calculo_analitico"


def test_calcular_texto_com_conversao_de_saida(capsys):
    code, out, _ = run(capsys, "calcular", "periodo_orbital_circular", "mu=3.986004418e14", "r=42164km", "--em", "h")
    assert code == 0
    assert "Natureza do resultado" in out and "não é medição" in out
    assert " h" in out


@pytest.mark.parametrize(
    "argv, fragment",
    [
        (["calcular", "inexistente"], "desconhecida"),
        (["calcular", "velocidade_escape", "mu=1"], "ausente"),
        (["calcular", "velocidade_escape", "mu=1", "r=2", "q=3"], "não possui"),
        (["calcular", "velocidade_escape", "mu=1", "r=5 s"], "comprimento"),
        (["calcular", "velocidade_escape", "mu=1", "r"], "NOME=VALOR"),
        (["calcular", "velocidade_escape", "mu=1", "r=1", "r=2"], "mais de uma vez"),
        (["calcular", "velocidade_escape", "mu=-1", "r=1"], "positivo"),
        (["calcular", "velocidade_escape", "mu=1", "r=1", "--em", "s"], "tempo"),
    ],
)
def test_calcular_erros(capsys, argv, fragment):
    code, _, err = run(capsys, *argv)
    assert code == 1
    assert fragment in err


def test_orbita_circular(capsys):
    code, out, _ = run(capsys, "orbita-circular", "--corpo", "Terra", "--altitude", "400 km")
    assert code == 0
    assert "Terra" in out and "WGS 84" in out
    assert "Hipóteses" in out and "não é medição" in out


def test_orbita_circular_json_por_raio(capsys):
    code, out, _ = run(capsys, "orbita-circular", "--corpo", "lua", "--raio", "1837.4 km", "--json")
    assert code == 0
    data = json.loads(out)
    assert data["corpo"]["chave"] == "lua"
    assert data["altitude"]["valor"] == pytest.approx(100e3)
    assert len(data["resultados"]) == 5


def test_orbita_circular_abaixo_da_superficie(capsys):
    code, _, err = run(capsys, "orbita-circular", "--altitude", "-10 km")
    assert code == 1
    assert "intersectaria" in err


def test_orbita_circular_corpo_desconhecido(capsys):
    code, _, err = run(capsys, "orbita-circular", "--corpo", "plutao", "--altitude", "10 km")
    assert code == 1
    assert "desconhecido" in err


def test_corpos_e_constantes(capsys):
    code, out, _ = run(capsys, "corpos")
    assert code == 0 and "IAU 2015" in out
    code, out, _ = run(capsys, "constantes")
    assert code == 0 and "CODATA 2018" in out


def test_execucao_como_modulo():
    env = {**os.environ, "PYTHONPATH": str(ROOT), "PYTHONIOENCODING": "utf-8"}
    proc = subprocess.run(
        [sys.executable, "-m", "aurora", "--versao"], cwd=ROOT, env=env, capture_output=True, text=True, check=False
    )
    assert proc.returncode == 0
    assert "AURORA" in proc.stdout
