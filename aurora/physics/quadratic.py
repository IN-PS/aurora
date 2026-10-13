"""Calculadora de equações do segundo grau ``a·x² + b·x + c = 0``.

Fórmulas (com ``a ≠ 0``)::

    Δ = b² − 4·a·c
    x = (−b ± √Δ) / (2·a)

Cuidados numéricos (ponto flutuante IEEE 754 de 64 bits):

* **Discriminante exato.** ``b²`` e ``4·a·c`` podem ser quase iguais e a subtração em ponto
  flutuante perderia todos os algarismos significativos. Como todo ``float`` é um número
  racional exato, Δ é calculado em aritmética racional (:class:`fractions.Fraction`) e só
  então arredondado. Assim, a classificação (duas raízes reais, raiz dupla ou complexas) é
  exata para os coeficientes fornecidos.
* **Cancelamento catastrófico.** Quando ``|b| ≫ |4ac|``, ``−b + √Δ`` (ou ``−b − √Δ``)
  subtrai números quase iguais. Usamos a forma estável ``q = −(b + sinal(b)·√Δ)/2``,
  ``x₁ = q/a`` e ``x₂ = c/q`` (Press et al., *Numerical Recipes*, §5.6).
* **Overflow.** ``√Δ`` é obtida reescalando Δ por potências de 4, de modo que coeficientes
  grandes (ex.: 10²⁰⁰) não produzem ``inf`` em resultados intermediários.

Limitação inerente (não algorítmica): perto de uma raiz dupla o problema é mal condicionado;
uma perturbação relativa ε nos coeficientes desloca as raízes em ordem de √ε.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction

from ..core.validation import as_real
from ..errors import DomainError


class RootKind(Enum):
    """Natureza das raízes, determinada pelo sinal do discriminante."""

    TWO_REAL = "duas raízes reais distintas (Δ > 0)"
    DOUBLE = "uma raiz real dupla (Δ = 0)"
    COMPLEX = "duas raízes complexas conjugadas (Δ < 0)"


@dataclass(frozen=True)
class QuadraticSolution:
    """Solução de ``a·x² + b·x + c = 0``.

    ``roots`` contém sempre duas raízes, contadas com multiplicidade: reais em ordem
    crescente, ou complexas conjugadas com a parte imaginária positiva primeiro.
    ``discriminant`` é Δ arredondado para ``float`` (``±inf`` se exceder o intervalo
    representável; a classificação em ``kind`` continua exata).
    """

    a: float
    b: float
    c: float
    discriminant: float
    kind: RootKind
    roots: tuple[float, float] | tuple[complex, complex]

    @property
    def real_roots(self) -> tuple[float, ...]:
        """Raízes reais distintas: duas, uma (raiz dupla) ou nenhuma."""
        if self.kind is RootKind.TWO_REAL:
            return self.roots  # type: ignore[return-value]
        if self.kind is RootKind.DOUBLE:
            return (self.roots[0],)  # type: ignore[return-value]
        return ()


def _exact_discriminant(a: float, b: float, c: float) -> Fraction:
    fa, fb, fc = Fraction(a), Fraction(b), Fraction(c)
    return fb * fb - 4 * fa * fc


def _fraction_to_float(value: Fraction) -> float:
    try:
        return float(value)
    except OverflowError:
        return math.inf if value > 0 else -math.inf


def _sqrt_fraction(value: Fraction) -> float:
    """√value (value > 0) sem overflow/underflow intermediário; erro de ~1 ulp."""
    e = (value.numerator.bit_length() - value.denominator.bit_length()) // 2
    scaled = value / Fraction(4) ** e if e >= 0 else value * Fraction(4) ** (-e)
    return math.ldexp(math.sqrt(float(scaled)), e)


def discriminant(a: float, b: float, c: float) -> float:
    """Δ = b² − 4·a·c, calculado exatamente e arredondado uma única vez para ``float``."""
    return _fraction_to_float(_exact_discriminant(as_real("a", a), as_real("b", b), as_real("c", c)))


def solve_quadratic(a: float, b: float, c: float) -> QuadraticSolution:
    """Resolve ``a·x² + b·x + c = 0`` para coeficientes reais finitos com ``a ≠ 0``.

    Levanta :class:`~aurora.errors.DomainError` se ``a = 0`` (equação linear) ou se alguma
    raiz exceder o intervalo representável em ponto flutuante de 64 bits.
    """
    a = as_real("a", a)
    b = as_real("b", b)
    c = as_real("c", c)
    if a == 0.0:
        raise DomainError(
            "O coeficiente 'a' deve ser diferente de zero; com a = 0 a equação é linear "
            "(b·x + c = 0), não do segundo grau."
        )

    delta = _exact_discriminant(a, b, c)
    half_b = 0.5 * b

    if delta > 0:
        kind = RootKind.TWO_REAL
        half_sqrt = 0.5 * _sqrt_fraction(delta)
        q = -(half_b + math.copysign(half_sqrt, b))
        if q == 0.0:
            raise DomainError("Coeficientes pequenos demais para o cálculo em ponto flutuante de 64 bits.")
        x1, x2 = sorted((q / a, c / q))
        roots: tuple[float, float] | tuple[complex, complex] = (x1 + 0.0, x2 + 0.0)
        parts = roots
    elif delta == 0:
        kind = RootKind.DOUBLE
        x = -half_b / a + 0.0  # "+ 0.0" normaliza −0.0 para 0.0
        roots = (x, x)
        parts = roots
    else:
        kind = RootKind.COMPLEX
        real = -half_b / a + 0.0
        imag = 0.5 * _sqrt_fraction(-delta) / abs(a)
        roots = (complex(real, imag), complex(real, -imag))
        parts = (real, imag)

    if not all(math.isfinite(p) for p in parts):
        raise DomainError("Raiz fora do intervalo representável em ponto flutuante de 64 bits.")

    return QuadraticSolution(
        a=a, b=b, c=c, discriminant=_fraction_to_float(delta), kind=kind, roots=roots
    )
