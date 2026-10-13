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
