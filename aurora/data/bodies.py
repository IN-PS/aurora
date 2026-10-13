"""Corpos celestes de referência para os cálculos analíticos.

Usamos o parâmetro gravitacional ``mu = G*M`` publicado diretamente, e não ``G`` vezes uma
massa: ``mu`` é determinado por rastreamento de órbitas com precisão muito superior à de
``G`` ou ``M`` isoladamente.

O conjunto é propositalmente pequeno. A etapa 4 (exploração espacial) prevê a importação de
efemérides e constantes de fontes confiáveis em vez de ampliar esta tabela manualmente.
"""

from __future__ import annotations

import unicodedata
from dataclasses import dataclass

from ..errors import NotFoundError


@dataclass(frozen=True)
class CelestialBody:
    """Corpo celeste modelado como esfera com simetria esférica de massa."""

    key: str
    name: str
    mu: float
    """Parâmetro gravitacional G*M, em m^3/s^2."""
    radius: float
    """Raio de referência, em m (ver ``radius_description``)."""
    radius_description: str
    sources: tuple[str, ...]


SUN = CelestialBody(
    key="sol",
    name="Sol",
    mu=1.3271244e20,
    radius=6.957e8,
    radius_description="raio nominal",
    sources=("IAU 2015, Resolução B3: parâmetro gravitacional nominal (GM)☉ e raio nominal R☉",),
)

EARTH = CelestialBody(
    key="terra",
    name="Terra",
    mu=3.986004418e14,
    radius=6_378_137.0,
    radius_description="raio equatorial",
    sources=("WGS 84 (NGA TR8350.2): GM da Terra e semieixo maior do elipsoide",),
)

MOON = CelestialBody(
    key="lua",
    name="Lua",
    mu=4.9028e12,
    radius=1_737_400.0,
    radius_description="raio médio",
    sources=(
        "GM: efemérides JPL DE430 (Folkner et al., 2014), arredondado a 5 algarismos significativos",
        "Raio médio: IAU WGCCRE 2009 (Archinal et al., 2011)",
    ),
)

_BODIES: dict[str, CelestialBody] = {body.key: body for body in (SUN, EARTH, MOON)}


def _normalize(text: str) -> str:
    decomposed = unicodedata.normalize("NFKD", text.strip().casefold())
    return "".join(ch for ch in decomposed if not unicodedata.combining(ch))


def get_body(key: str) -> CelestialBody:
    """Retorna o corpo pela chave (``"terra"``, ``"lua"``, ``"sol"``), sem distinção de caixa ou acento."""
    body = _BODIES.get(_normalize(key))
    if body is None:
        known = ", ".join(sorted(_BODIES))
        raise NotFoundError(f"Corpo celeste desconhecido: {key!r}. Disponíveis: {known}.")
    return body


def list_bodies() -> tuple[CelestialBody, ...]:
    """Lista os corpos disponíveis."""
    return tuple(_BODIES.values())
