"""Hierarquia de exceções do AURORA.

Módulo sem dependências, utilizável por todas as camadas. Todas as mensagens são escritas
em português brasileiro, pois são exibidas diretamente ao usuário pela interface.
"""

from __future__ import annotations


class AuroraError(Exception):
    """Base de todas as exceções do AURORA."""


class InvalidParameterError(AuroraError, ValueError):
    """Parâmetro com tipo ou valor inválido (não numérico, NaN, infinito, fora de faixa)."""


class DomainError(InvalidParameterError):
    """Parâmetros válidos isoladamente, mas fora do domínio matemático ou físico do modelo."""


class UnitError(AuroraError, ValueError):
    """Unidade desconhecida ou incompatível com a grandeza esperada."""


class NotFoundError(AuroraError, LookupError):
    """Item (fórmula, corpo celeste, unidade) não encontrado."""

    def __str__(self) -> str:  # LookupError herdaria a representação entre aspas de KeyError
        return str(self.args[0]) if self.args else ""


class SimulationError(AuroraError):
    """Falha da simulação numérica (integrador não convergiu, estado não representável…)."""


class DependencyError(AuroraError):
    """Dependência opcional ausente (por exemplo, a biblioteca de gráficos)."""
