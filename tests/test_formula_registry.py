"""Testes do registro de fórmulas (metadados, unidades e avaliação)."""

import inspect

import pytest

import aurora.physics  # noqa: F401  (registra as fórmulas)
from aurora.core import units
from aurora.core.formula import get_formula, list_formulas, register_formula
from aurora.core.provenance import DataOrigin
from aurora.errors import DomainError, InvalidParameterError, NotFoundError

FORMULAS_DA_ESPECIFICACAO = {
    "forca_gravitacional": "F = G*m1*m2/r^2",
    "parametro_gravitacional": "mu = G*M",
    "velocidade_orbital_circular": "v = sqrt(mu/r)",
    "periodo_orbital_circular": "T = 2*pi*sqrt(r^3/mu)",
    "energia_especifica": "epsilon = v^2/2 - mu/r",
}


def test_formulas_da_especificacao_estao_registradas():
    for key, expression in FORMULAS_DA_ESPECIFICACAO.items():
        assert get_formula(key).expression == expression


def test_metadados_completos_e_coerentes():
    formulas = list_formulas()
    assert len(formulas) >= 15
    for f in formulas:
        assert f.title and f.expression and f.hypotheses, f.key
        assert units.get_unit(f.unit).factor == 1.0, f.key
        signature = list(inspect.signature(f.function).parameters)
        assert signature == [p.name for p in f.parameters], f.key
        for p in f.parameters:
            assert p.description, (f.key, p.name)
            assert units.get_unit(p.unit).factor == 1.0, (f.key, p.name)


def test_avaliacao_retorna_resultado_analitico_com_metadados():
    result = get_formula("velocidade_escape").evaluate(mu=2.0, r=1.0)
    assert result.value == 2.0
    assert result.unit == "m/s"
    assert result.origin is DataOrigin.ANALYTIC
    data = result.to_dict()
    assert data["origem"] == "calculo_analitico"
    assert data["entradas"]["r"] == {"valor": 1.0, "unidade": "m", "descricao": "Distância ao centro do corpo central"}
    assert "não é medição" in data["aviso"]
    assert data["hipoteses"]


def test_parametros_ausentes_ou_desconhecidos():
    f = get_formula("velocidade_escape")
    with pytest.raises(InvalidParameterError, match="ausente"):
        f.evaluate(mu=1.0)
    with pytest.raises(InvalidParameterError, match="não reconhecido"):
        f.evaluate(mu=1.0, r=1.0, x=3.0)
    with pytest.raises(InvalidParameterError, match="não possui"):
        f.parameter("x")


def test_formula_desconhecida():
    with pytest.raises(NotFoundError, match="desconhecida"):
        get_formula("dobra_espacial")


def test_resultado_nao_finito_vira_erro_de_dominio():
    with pytest.raises(DomainError, match="não é representável"):
        get_formula("forca_gravitacional").evaluate(m1=1e300, m2=1e300, r=1e-10)


def test_registro_rejeita_assinatura_incoerente():
    with pytest.raises(ValueError, match="diferem da assinatura"):

        @register_formula(
            key="_teste_assinatura",
            title="t",
            expression="y = x",
            unit="m",
            parameters={"z": ("z", "m")},
            hypotheses=("h",),
        )
        def _f(x):
            return x


def test_registro_rejeita_unidade_fora_do_si():
    with pytest.raises(ValueError, match="SI"):

        @register_formula(
            key="_teste_unidade",
            title="t",
            expression="y = x",
            unit="km",
            parameters={"x": ("x", "m")},
            hypotheses=("h",),
        )
        def _f(x):
            return x


# --- regressões da auditoria técnica (docs/AUDITORIA_ETAPA1.md) ---------------------------


def test_rede_de_seguranca_converte_divisao_por_zero_em_erro_de_dominio():
    from aurora.core import formula as formula_module

    @register_formula(
        key="_teste_divisao",
        title="teste",
        expression="y = 1/x",
        unit="m",
        parameters={"x": ("x", "m")},
        hypotheses=("h",),
    )
    def _f(x):
        return 1.0 / (x * 0.0)

    try:
        with pytest.raises(DomainError, match="não é representável"):
            get_formula("_teste_divisao").evaluate(x=1.0)
    finally:
        formula_module._REGISTRY.pop("_teste_divisao", None)


def test_recarregar_modulo_nao_duplica_nem_quebra_o_registro():
    # ARCH-4: reload (ex.: autoreload em notebooks) levantava ValueError. Executado em outro
    # processo porque reload recria as classes do módulo e afetaria os demais testes.
    import subprocess
    import sys
    from pathlib import Path

    code = (
        "import importlib, aurora.physics, aurora.physics.orbits as o\n"
        "from aurora.core.formula import list_formulas\n"
        "n = len(list_formulas()); importlib.reload(o); assert len(list_formulas()) == n; print(n)"
    )
    root = Path(__file__).resolve().parents[1]
    proc = subprocess.run([sys.executable, "-c", code], cwd=root, capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "18"


def test_mensagens_de_busca_de_formula():
    with pytest.raises(NotFoundError, match="Fórmulas disponíveis: .*velocidade_escape"):
        get_formula("inexistente")
    with pytest.raises(InvalidParameterError):
        get_formula(None)


def test_formulas_orbitais_trazem_referencias_e_hipoteses_coerentes():
    for key in ("velocidade_orbital_circular", "periodo_orbital_eliptico", "velocidade_vis_viva"):
        assert get_formula(key).references
    # ARCH-14: periapse/apoapse não dependem de μ.
    for key in ("raio_periapse", "raio_apoapse"):
        assert not any("μ = G·M" in h for h in get_formula(key).hypotheses)
