"""Veículo espacial virtual: massa seca, propelente e impulso específico.

O consumo de propelente em cada manobra impulsiva segue a equação de Tsiolkovsky
(:mod:`aurora.physics.maneuvers`). O veículo é imutável: uma queima devolve um novo veículo
com menos propelente.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from ..core.validation import as_non_negative, as_positive, as_real
from ..errors import DomainError, InvalidParameterError
from ..physics.maneuvers import propellant_for_delta_v, rocket_delta_v


@dataclass(frozen=True)
class Spacecraft:
    """Veículo virtual (massas em kg, impulso específico em s)."""

    name: str = "Veículo virtual"
    dry_mass: float = 1000.0
    propellant_mass: float = 0.0
    isp: float = 300.0

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise InvalidParameterError("O veículo precisa de um nome.")
        object.__setattr__(self, "dry_mass", as_positive("dry_mass", self.dry_mass))
        object.__setattr__(self, "propellant_mass", as_non_negative("propellant_mass", self.propellant_mass))
        object.__setattr__(self, "isp", as_positive("isp", self.isp))

    @property
    def mass(self) -> float:
        """Massa total atual (kg)."""
        return self.dry_mass + self.propellant_mass

    @property
    def delta_v_capacity(self) -> float:
        """Δv máximo (m/s) obtido queimando todo o propelente restante."""
        return rocket_delta_v(self.isp, self.mass, self.dry_mass)

    def burn(self, delta_v: float) -> tuple[Spacecraft, float]:
        """Executa uma queima de módulo ``delta_v`` (m/s); retorna o novo veículo e o propelente gasto (kg)."""
        delta_v = abs(as_real("delta_v", delta_v))  # o sentido da queima não altera o consumo
        needed = propellant_for_delta_v(self.isp, self.mass, delta_v)
        if needed > self.propellant_mass * (1.0 + 1e-12):
            raise DomainError(
                f"Propelente insuficiente para Δv = {delta_v:.6g} m/s: necessário {needed:.6g} kg, "
                f"disponível {self.propellant_mass:.6g} kg (capacidade total ≈ {self.delta_v_capacity:.6g} m/s)."
            )
        used = min(needed, self.propellant_mass)
        return replace(self, propellant_mass=self.propellant_mass - used), used
