"""Calculadora de equações do segundo grau ``a·x² + b·x + c = 0``.

Fórmulas (com ``a ≠ 0``)::

    Δ = b² − 4·a·c
    x = (−b ± √Δ) / (2·a)

Cuidados numéricos (ponto flutuante IEEE 754 de 64 bits):

* **Discriminante exato.** ``b²`` e ``4·a·c`` podem ser quase iguais e a subtração em ponto
  flutuante perderia todos os algarismos significativos. Δ é calculado em aritmética
  racional (:class:`fractions.Fraction`) e só então arredondado. Coeficientes ``int`` ou
  ``Fraction`` são usados com seu valor exato — a interface converte os decimais digitados
  (``0.1``, ``0,6``) em ``Fraction`` —, e ``float`` com seu valor binário exato. Assim, a
  classificação (duas raízes reais, raiz dupla ou complexas) é exata para os coeficientes
  fornecidos: ``x² + 0,6x + 0,09`` digitado na interface tem raiz dupla.
* **Cancelamento catastrófico.** Quando ``|b| ≫ |4ac|``, ``−b + √Δ`` (ou ``−b − √Δ``)
  subtrai números quase iguais. Usamos a forma estável ``q = −(b + sinal(b)·√Δ)/2``,
  ``x₁ = q/a`` e ``x₂ = c/q`` (Press et al., *Numerical Recipes*, §5.6).
* **Overflow e underflow.** ``√Δ`` é aproximada (erro de ~1 ulp) reescalando Δ por potências
  de 4, e ``q``, ``q/a`` e ``c/q`` são calculados em aritmética racional e arredondados uma
  única vez. Não há overflow nem underflow intermediário: coeficientes subnormais ou
  próximos de 1,8 × 10³⁰⁸ dão raízes corretas (erro de ~2 ulp), e só é levantado erro
  quando a própria raiz não cabe em ``float64``.

Limitação inerente (não algorítmica): perto de uma raiz dupla o problema é mal condicionado;
uma perturbação relativa ε nos coeficientes desloca as raízes em ordem de √ε.
"""

from __future__ import annotations

import math
import numbers
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


def _exact(value: object, as_float: float) -> Fraction:
    """Valor exato do coeficiente: inteiros e frações são mantidos; ``float`` é exato como binário."""
    return Fraction(value) if isinstance(value, numbers.Rational) else Fraction(as_float)


def _fraction_to_float(value: Fraction) -> float:
    try:
        return float(value)
    except OverflowError:
        return math.inf if value > 0 else -math.inf


def _sqrt_fraction(value: Fraction) -> Fraction:
    """Aproximação racional de √value (value > 0), com erro relativo de ~1 ulp de float64.

    Δ é reescalado por potências de 4 antes da raiz, e o fator 2^e é reaplicado em
    aritmética racional: não há overflow nem underflow intermediário.
    """
    e = (value.numerator.bit_length() - value.denominator.bit_length()) // 2
    scaled = value / Fraction(4) ** e if e >= 0 else value * Fraction(4) ** (-e)
    return Fraction(math.sqrt(float(scaled))) * Fraction(2) ** e


def _solve_rational(fa: Fraction, fb: Fraction, fc: Fraction) -> tuple[Fraction, RootKind, tuple]:
    """Núcleo do solver para coeficientes racionais exatos (``fa ≠ 0``): retorna (Δ exato, tipo, raízes)."""
    delta = fb * fb - 4 * fa * fc
    if delta > 0:
        kind = RootKind.TWO_REAL
        s = _sqrt_fraction(delta)
        q = -(fb + s) / 2 if fb >= 0 else -(fb - s) / 2  # |q| ≥ √Δ/2 > 0: sem cancelamento
        x1, x2 = sorted((_fraction_to_float(q / fa), _fraction_to_float(fc / q)))
        roots: tuple = (x1, x2)
        parts: tuple[float, ...] = roots
    elif delta == 0:
        kind = RootKind.DOUBLE
        x = _fraction_to_float(-fb / (2 * fa))
        roots = (x, x)
        parts = roots
    else:
        kind = RootKind.COMPLEX
        real = _fraction_to_float(-fb / (2 * fa))
        imag = _fraction_to_float(_sqrt_fraction(-delta) / (2 * abs(fa)))
        roots = (complex(real, imag), complex(real, -imag))
        parts = (real, imag)
    if not all(math.isfinite(p) for p in parts):
        raise DomainError("Raiz fora do intervalo representável em ponto flutuante de 64 bits.")
    return delta, kind, roots


def discriminant(a: float, b: float, c: float) -> float:
    """Δ = b² − 4·a·c, calculado exatamente e arredondado uma única vez para ``float``."""
    fa, fb, fc = (_exact(v, as_real(n, v)) for n, v in (("a", a), ("b", b), ("c", c)))
    return _fraction_to_float(fb * fb - 4 * fa * fc)


def solve_quadratic(a: float, b: float, c: float) -> QuadraticSolution:
    """Resolve ``a·x² + b·x + c = 0`` para coeficientes reais finitos com ``a ≠ 0``.

    Coeficientes ``int`` ou :class:`fractions.Fraction` (por exemplo, decimais digitados na
    interface) são usados com seu valor exato; ``float`` é usado com seu valor binário exato.
    Levanta :class:`~aurora.errors.DomainError` se ``a = 0`` (equação linear) ou se alguma
    raiz exceder o intervalo representável em ponto flutuante de 64 bits.
    """
    fa = _exact(a, af := as_real("a", a))
    fb = _exact(b, bf := as_real("b", b))
    fc = _exact(c, cf := as_real("c", c))
    if fa == 0:
        raise DomainError(
            "O coeficiente 'a' deve ser diferente de zero; com a = 0 a equação é linear "
            "(b·x + c = 0), não do segundo grau."
        )
    delta, kind, roots = _solve_rational(fa, fb, fc)
    return QuadraticSolution(a=af, b=bf, c=cf, discriminant=_fraction_to_float(delta), kind=kind, roots=roots)
