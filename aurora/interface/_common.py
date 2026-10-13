"""Formatação comum da interface de linha de comando (texto e JSON)."""

from __future__ import annotations

import json
import math
from collections.abc import Sequence
from typing import Any

from ..core import units
from ..core.provenance import DataOrigin
from ..errors import AuroraError


def _fmt(x: float, digits: int = 10) -> str:
    if math.isnan(x):
        return "indefinido"
    if math.isinf(x):
        return "∞" if x > 0 else "−∞"
    return f"{x + 0.0:.{digits}g}"  # "+ 0.0" evita exibir "-0"


def _fmt_quantity(value: float, unit: str) -> str:
    """Valor na unidade SI seguido de unidades alternativas, ex.: '7668.56 m/s (7.66856 km/s)'."""
    text = f"{_fmt(value)} {unit}"
    if math.isfinite(value):
        alts = []
        for alt in units.display_alternatives(unit):
            try:
                converted = units.from_si(value, alt)
            except AuroraError:  # conversão apenas ilustrativa: omitida se exceder o float64
                continue
            if 0.01 <= abs(converted) < 1e9:  # omite conversões pouco legíveis (ex.: 3e-11 ua)
                alts.append(f"{_fmt(converted, 7)} {alt}")
        if alts:
            text += f"  ({'; '.join(alts)})"
    return text


def _json_safe(obj: Any) -> Any:
    """Prepara para JSON estrito: não finitos viram os textos "inf", "-inf" ou "nan"; −0.0 vira 0.0."""
    if isinstance(obj, float):
        if not math.isfinite(obj):
            return "nan" if math.isnan(obj) else ("inf" if obj > 0 else "-inf")
        return obj + 0.0
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _print_json(data: Any) -> None:
    # ensure_ascii: caracteres como μ, ☉ e ç viram escapes \\uXXXX e não são corrompidos quando a
    # saída padrão não é UTF-8 (ex.: redirecionamento no Windows).
    print(json.dumps(_json_safe(data), ensure_ascii=True, indent=2, allow_nan=False))


def _print_list(title: str, items: Sequence[str]) -> None:
    print(title)
    indent = " " * (len(title) - len(title.lstrip()) + 2)
    for item in items:
        print(f"{indent}- {item}")


def _print_footer(origin: DataOrigin) -> None:
    print(f"Natureza do resultado: {origin.disclaimer}")
