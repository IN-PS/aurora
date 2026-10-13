"""Validação de parâmetros numéricos.

Todas as funções retornam o valor convertido para ``float`` (ou ``numpy.ndarray`` de
``float64``) e levantam :class:`~aurora.errors.InvalidParameterError` com mensagem em
português quando a entrada é inválida. Valores lógicos (``True``/``False``), textos, números
complexos, NaN e infinitos são sempre rejeitados.
"""

from __future__ import annotations

import math
import numbers

import numpy as np

from ..errors import InvalidParameterError


def as_real(name: str, value: object) -> float:
    """Exige um número real finito e o converte para ``float``."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise InvalidParameterError(
            f"O parâmetro '{name}' deve ser um número real; recebido {value!r} ({type(value).__name__})."
        )
    try:
        x = float(value)
    except OverflowError:
        raise InvalidParameterError(
            f"O parâmetro '{name}' excede o intervalo representável em ponto flutuante de 64 bits."
        ) from None
    if not math.isfinite(x):
        raise InvalidParameterError(f"O parâmetro '{name}' deve ser finito; recebido {x!r}.")
    return x


def as_positive(name: str, value: object) -> float:
    """Exige um número real finito estritamente positivo."""
    x = as_real(name, value)
    if not x > 0.0:
        raise InvalidParameterError(f"O parâmetro '{name}' deve ser estritamente positivo; recebido {x!r}.")
    return x


def as_non_negative(name: str, value: object) -> float:
    """Exige um número real finito maior ou igual a zero."""
    x = as_real(name, value)
    if x < 0.0:
        raise InvalidParameterError(f"O parâmetro '{name}' não pode ser negativo; recebido {x!r}.")
    return x


def as_nonzero(name: str, value: object) -> float:
    """Exige um número real finito diferente de zero."""
    x = as_real(name, value)
    if x == 0.0:
        raise InvalidParameterError(f"O parâmetro '{name}' deve ser diferente de zero.")
    return x


def as_vector(name: str, value: object, size: int = 3) -> np.ndarray:
    """Exige um vetor de ``size`` componentes reais finitas; retorna cópia em ``float64``."""
    try:
        raw = np.asarray(value)
    except (TypeError, ValueError):
        raw = None
    if raw is None or raw.dtype.kind not in "iuf":
        raise InvalidParameterError(f"O parâmetro '{name}' deve ser um vetor de números reais; recebido {value!r}.")
    if raw.shape != (size,):
        raise InvalidParameterError(
            f"O parâmetro '{name}' deve ter exatamente {size} componentes; recebido formato {raw.shape}."
        )
    vec = raw.astype(np.float64, copy=True)
    if not np.all(np.isfinite(vec)):
        raise InvalidParameterError(f"Todas as componentes de '{name}' devem ser finitas; recebido {vec!r}.")
    return vec
