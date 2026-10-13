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


# --- regressões da auditoria técnica (docs/AUDITORIA_ETAPA1.md) ---------------------------


@pytest.mark.parametrize(
    "argv",
    [
        ["calcular", "velocidade_escape", "--json", "mu=1", "r=1"],  # CLI--1: falhava no Python 3.10–3.12
        ["calcular", "velocidade_escape", "mu=1", "--json", "r=1"],  # falhava em todas as versões
        ["calcular", "--json", "velocidade_escape", "mu=1", "r=1"],
    ],
)
def test_opcoes_em_qualquer_posicao(capsys, argv):
    code, out, _ = run(capsys, *argv)
    assert code == 0
    assert json.loads(out)["valor"] == pytest.approx(2**0.5)


def test_json_com_conversao_de_unidade(capsys):
    argv = ["calcular", "velocidade_escape", "mu=3.986004418e14", "r=6778137", "--em", "km/s", "--json"]
    code, out, _ = run(capsys, *argv)
    assert code == 0
    data = json.loads(out)
    assert data["valor_convertido"]["unidade"] == "km/s"
    assert data["valor_convertido"]["valor"] == pytest.approx(data["valor"] / 1000.0, rel=1e-15)


def test_conversao_de_saida_com_overflow_gera_erro(capsys):
    code, _, err = run(capsys, "calcular", "velocidade_muv", "v0=1e308", "a=0", "t=0", "--em", "km/h")
    assert code == 1 and "excede" in err


def test_quadratica_texto_complexa_e_dupla(capsys):
    _, out, _ = run(capsys, "quadratica", "1", "2", "5")
    assert "x₁ = -1 + 2·i" in out and "x₂ = -1 − 2·i" in out
    _, out, _ = run(capsys, "quadratica", "1", "-4", "4")
    assert "x₁ = x₂ = 2" in out


@pytest.mark.parametrize("coef, root", [(("1", "0.6", "0.09"), -0.3), (("1", "-0,2", "0,01"), 0.1)])
def test_quadratica_decimais_digitados_sao_exatos(capsys, coef, root):
    # ARCH-1: (x + 0,3)² era classificada como "complexas conjugadas".
    code, out, _ = run(capsys, "quadratica", *coef, "--json")
    data = json.loads(out)
    assert data["classificacao"] == "raiz_dupla" and data["sinal_discriminante"] == 0
    assert data["raizes"][0]["real"] == root
    assert data["origem"] == "calculo_analitico"


def test_quadratica_discriminante_fora_do_float64(capsys):
    # ARCH-11: antes exibia "Δ = ∞".
    _, out, _ = run(capsys, "quadratica", "1e200", "-3e200", "2e200")
    assert "excede o intervalo do float64" in out and "= ∞" not in out


def test_ajuda_e_versao_com_nomes_convencionais(capsys):
    for flag in ("--help", "--ajuda"):
        with pytest.raises(SystemExit) as exc:
            main([flag])
        assert exc.value.code == 0
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "AURORA" in capsys.readouterr().out


@pytest.mark.parametrize(
    "argv, fragment",
    [
        ([], "argumentos obrigatórios ausentes"),
        (["quadratica", "1", "2"], "argumentos obrigatórios ausentes"),
        (["quadratica", "1", "2", "3", "4"], "argumentos não reconhecidos"),
        (["inexistente"], "opção inválida"),
        (["orbita-circular"], "informe um dos argumentos"),
    ],
)
def test_erros_de_uso_em_portugues(capsys, argv, fragment):
    with pytest.raises(SystemExit) as exc:
        main(argv)
    assert exc.value.code == 2
    assert fragment in capsys.readouterr().err


def test_mensagens_nao_mostram_o_espaco_de_protecao_de_negativos(capsys):
    code, _, err = run(capsys, "orbita-circular", "--altitude", "-1.2.3 km")
    assert code == 1 and "'-1.2.3 km'" in err


def test_menos_zero_nao_aparece_na_saida(capsys):
    _, out, _ = run(capsys, "calcular", "forca_gravitacional", "m1=-0", "m2=5", "r=1")
    assert "F = 0 N" in out and "-0" not in out


@pytest.mark.parametrize(
    "argv",
    [
        ["quadratica", "1e308", "0", "1e308"],
        ["quadratica", "1e308", "1e308", "-1e308"],
        ["quadratica", "4.9e-324", "4.9e-324", "0"],
        ["quadratica", "1e-300", "1e300", "1"],
        ["calcular", "aceleracao_gravitacional", "mu=3.986e14", "r=1e-170"],
        ["calcular", "forca_gravitacional", "m1=1", "m2=1", "r=5e-324"],
        ["calcular", "velocidade_torricelli", "v0=1e200", "a=1", "dx=1"],
        ["calcular", "periodo_orbital_eliptico", "mu=1", "a=1e300"],
        ["orbita-circular", "--raio", "1e308 ua"],
        ["orbita-circular", "--altitude", "-7000 km"],
    ],
)
def test_entradas_extremas_nunca_vazam_traceback(capsys, argv):
    # Antes da auditoria, OverflowError e ZeroDivisionError escapavam como traceback.
    code, _, err = run(capsys, *argv)
    assert code in (0, 1)
    assert "Traceback" not in err


def test_json_preserva_caracteres_com_saida_nao_utf8():
    # CLI--4: com saída ASCII, μ/☉/ç viravam '?' silenciosamente no JSON exportado.
    env = {**os.environ, "PYTHONPATH": str(ROOT), "PYTHONIOENCODING": "ascii"}
    proc = subprocess.run(
        [sys.executable, "-m", "aurora", "orbita-circular", "--corpo", "sol", "--raio", "1 ua", "--json"],
        cwd=ROOT, env=env, capture_output=True, check=False,
    )
    assert proc.returncode == 0
    data = json.loads(proc.stdout.decode("ascii"))
    assert "☉" in data["corpo"]["fontes"][0] and "Resolução" in data["corpo"]["fontes"][0]
