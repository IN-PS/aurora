"""Procedência dos dados: distingue cálculo analítico, simulação e medição real.

Todo resultado exposto ao usuário carrega uma :class:`DataOrigin`. O pacote ``aurora`` (o
simulador) só produz ``ANALYTIC`` e ``SIMULATED``. ``MEASURED`` está definido aqui apenas
para que um futuro software de hardware, desenvolvido em pacote separado, possa usar o mesmo
contrato de dados; nenhum código do simulador deve criar dados com essa origem.
"""

from __future__ import annotations

from enum import Enum


class DataOrigin(Enum):
    """Origem de um dado ou resultado."""

    ANALYTIC = "calculo_analitico"
    SIMULATED = "simulacao"
    MEASURED = "medicao"

    @property
    def label(self) -> str:
        """Descrição curta em português."""
        return _LABELS[self]

    @property
    def disclaimer(self) -> str:
        """Aviso a ser exibido junto com o resultado."""
        return _DISCLAIMERS[self]


_LABELS = {
    DataOrigin.ANALYTIC: "Cálculo analítico",
    DataOrigin.SIMULATED: "Simulação numérica",
    DataOrigin.MEASURED: "Medição",
}

_DISCLAIMERS = {
    DataOrigin.ANALYTIC: "cálculo analítico com modelo idealizado — não é medição nem dado observacional.",
    DataOrigin.SIMULATED: (
        "simulação numérica com modelo idealizado — não representa detecção ou medição de objeto real."
    ),
    DataOrigin.MEASURED: (
        "dado obtido por instrumento real (reservado ao software de hardware; nunca produzido pelo simulador)."
    ),
}
