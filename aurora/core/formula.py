"""Registro de fórmulas com unidades, hipóteses e referências.

Os modelos físicos são funções Python comuns que recebem e retornam ``float`` em unidades do
SI. O decorador :func:`register_formula` associa a cada função os metadados que a interface
precisa exibir — expressão, unidade do resultado, descrição e unidade de cada parâmetro e
hipóteses do modelo — sem acoplar o modelo à interface.

Ao registrar, conferimos que os parâmetros declarados coincidem com a assinatura da função e
que todas as unidades estão no SI; incoerências são detectadas na importação do módulo.

O registro é preenchido quando ``aurora.physics`` é importado. Referências bibliográficas
gerais (livros-texto) acompanham as fórmulas orbitais; a bibliografia completa está em
``docs/MODELOS_FISICOS.md``.
"""

from __future__ import annotations

import functools
import inspect
import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

from ..errors import DomainError, InvalidParameterError, NotFoundError
from .provenance import DataOrigin
from .units import get_unit


@dataclass(frozen=True)
class Parameter:
    """Parâmetro de entrada de uma fórmula (sempre em unidade SI)."""

    name: str
    description: str
    unit: str


@dataclass(frozen=True)
class Formula:
    """Fórmula registrada: função pura + metadados para exibição e validação."""

    key: str
    title: str
    expression: str
    unit: str
    parameters: tuple[Parameter, ...]
    hypotheses: tuple[str, ...]
    function: Callable[..., float] = field(repr=False)
    references: tuple[str, ...] = ()

    @property
    def symbol(self) -> str:
        """Símbolo do resultado (lado esquerdo da expressão)."""
        return self.expression.split("=", 1)[0].strip()

    def parameter(self, name: str) -> Parameter:
        for p in self.parameters:
            if p.name == name:
                return p
        names = ", ".join(p.name for p in self.parameters)
        raise InvalidParameterError(
            f"A fórmula '{self.key}' não possui o parâmetro '{name}'. Parâmetros: {names}."
        )

    def evaluate(self, **inputs: float) -> FormulaResult:
        """Avalia a fórmula com entradas em SI e retorna o resultado com metadados."""
        expected = [p.name for p in self.parameters]
        missing = [name for name in expected if name not in inputs]
        extra = [name for name in inputs if name not in expected]
        if extra:
            raise InvalidParameterError(
                f"Parâmetro(s) não reconhecido(s) para '{self.key}': {', '.join(extra)}. "
                f"Esperados: {', '.join(expected)}."
            )
        if missing:
            raise InvalidParameterError(
                f"Parâmetro(s) ausente(s) para '{self.key}': {', '.join(missing)}."
            )
        value = self.function(**inputs)
        return FormulaResult(formula=self, value=value, inputs={name: float(inputs[name]) for name in expected})


@dataclass(frozen=True)
class FormulaResult:
    """Resultado de uma fórmula, com entradas, unidade, hipóteses e procedência."""

    formula: Formula
    value: float
    inputs: Mapping[str, float]
    origin: DataOrigin = DataOrigin.ANALYTIC

    @property
    def unit(self) -> str:
        return self.formula.unit

    def to_dict(self) -> dict[str, Any]:
        """Representação serializável (chaves em português) para exportação."""
        f = self.formula
        return {
            "formula": f.key,
            "titulo": f.title,
            "expressao": f.expression,
            "simbolo": f.symbol,
            "valor": self.value,
            "unidade": f.unit,
            "entradas": {
                p.name: {"valor": self.inputs[p.name], "unidade": p.unit, "descricao": p.description}
                for p in f.parameters
            },
            "hipoteses": list(f.hypotheses),
            "referencias": list(f.references),
            "origem": self.origin.value,
            "aviso": self.origin.disclaimer,
        }


_REGISTRY: dict[str, Formula] = {}


def register_formula(
    *,
    key: str,
    title: str,
    expression: str,
    unit: str,
    parameters: Mapping[str, tuple[str, str]],
    hypotheses: tuple[str, ...],
    references: tuple[str, ...] = (),
) -> Callable[[Callable[..., float]], Callable[..., float]]:
    """Registra uma função como fórmula consultável pela interface.

    ``parameters`` mapeia o nome de cada argumento para ``(descrição, unidade SI)``. A
    função decorada passa a verificar que o resultado é finito (detecta overflow).
    """

    def decorator(func: Callable[..., float]) -> Callable[..., float]:
        previous = _REGISTRY.get(key)
        if previous is not None:
            original = getattr(previous.function, "__wrapped__", previous.function)
            same_origin = (original.__module__, original.__qualname__) == (func.__module__, func.__qualname__)
            if not same_origin:  # recarregar o mesmo módulo (ex.: autoreload em notebooks) é permitido
                raise ValueError(f"Fórmula já registrada: {key!r}")
        signature_names = list(inspect.signature(func).parameters)
        if signature_names != list(parameters):
            raise ValueError(
                f"Parâmetros declarados para {key!r} ({list(parameters)}) diferem da assinatura {signature_names}"
            )
        for symbol in (unit, *(u for _, u in parameters.values())):
            if get_unit(symbol).factor != 1.0:
                raise ValueError(f"A fórmula {key!r} deve usar unidades SI; recebido {symbol!r}")

        message = (
            f"O resultado de '{title}' não é representável em ponto flutuante de 64 bits "
            "(valores de entrada extremos)."
        )

        @functools.wraps(func)
        def checked(*args: Any, **kwargs: Any) -> float:
            try:
                result = func(*args, **kwargs)
            except (ZeroDivisionError, OverflowError):  # rede de segurança para over/underflow
                raise DomainError(message) from None
            if not math.isfinite(result):
                raise DomainError(message)
            return result

        _REGISTRY[key] = Formula(
            key=key,
            title=title,
            expression=expression,
            unit=unit,
            parameters=tuple(Parameter(name, desc, u) for name, (desc, u) in parameters.items()),
            hypotheses=tuple(hypotheses),
            function=checked,
            references=tuple(references),
        )
        return checked

    return decorator


def get_formula(key: str) -> Formula:
    """Retorna a fórmula registrada com a chave dada.

    As fórmulas são registradas na importação de ``aurora.physics``.
    """
    if not isinstance(key, str):
        raise InvalidParameterError(f"A chave da fórmula deve ser um texto; recebido {key!r}.")
    formula = _REGISTRY.get(key.strip())
    if formula is None:
        if not _REGISTRY:
            raise NotFoundError("Nenhuma fórmula registrada: importe o pacote aurora.physics.")
        raise NotFoundError(
            f"Fórmula desconhecida: {key.strip()!r}. Fórmulas disponíveis: {', '.join(sorted(_REGISTRY))}."
        )
    return formula


def list_formulas() -> tuple[Formula, ...]:
    """Todas as fórmulas registradas, na ordem de registro."""
    return tuple(_REGISTRY.values())
