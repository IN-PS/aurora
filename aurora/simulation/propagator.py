"""Propagação numérica de trajetórias (integração das equações de movimento).

Equações de movimento de uma partícula de teste, em referencial inercial::

    dr⃗/dt = v⃗
    dv⃗/dt = a⃗(r⃗, v⃗, t)        (dado por um GravityModel)

A integração usa :func:`scipy.integrate.solve_ivp` com o método de Dormand–Prince de ordem
8(5,3) (``DOP853``, padrão) ou 5(4) (``RK45``). O passo é ajustado automaticamente para que o
erro **local** estimado de cada componente fique abaixo de ``atol + rtol·|y|``, com
tolerâncias absolutas separadas para posição (m) e velocidade (m/s). O erro **global**
acumula ao longo da simulação; ele é medido, para o modelo de dois corpos, comparando com a
solução analítica de Kepler (:mod:`aurora.simulation.validation`).

Os resultados são **simulações numéricas de modelos idealizados** (``DataOrigin.SIMULATED``);
não representam medições nem trajetórias reais certificadas.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from ..core.provenance import DataOrigin
from ..core.validation import as_positive, as_real, as_vector
from ..errors import AuroraError, DomainError, InvalidParameterError, SimulationError
from ..physics.gravitation import GravityModel

METHODS = ("DOP853", "RK45")
_MAX_POINTS = 1_000_000


@dataclass(frozen=True)
class IntegratorSettings:
    """Configuração do integrador numérico (tolerâncias em SI)."""

    method: str = "DOP853"
    rtol: float = 1e-10
    atol_position: float = 1e-4
    """Tolerância absoluta por componente de posição, em m."""
    atol_velocity: float = 1e-7
    """Tolerância absoluta por componente de velocidade, em m/s."""
    max_step: float = math.inf
    """Maior passo permitido, em s (``inf`` = sem limite)."""
    max_evaluations: int = 2_000_000
    """Limite de avaliações da aceleração; protege contra simulações longas demais (da ordem de
    meio minuto num computador de mesa com o valor padrão)."""

    def __post_init__(self) -> None:
        if self.method not in METHODS:
            raise InvalidParameterError(
                f"Método de integração desconhecido: {self.method!r}. Use: {', '.join(METHODS)}."
            )
        rtol = as_positive("rtol", self.rtol)
        if not 1e-13 <= rtol <= 1e-3:
            raise InvalidParameterError("A tolerância relativa 'rtol' deve estar entre 1e-13 e 1e-3.")
        as_positive("atol_position", self.atol_position)
        as_positive("atol_velocity", self.atol_velocity)
        if not (self.max_step == math.inf or as_positive("max_step", self.max_step) > 0):
            raise InvalidParameterError("O passo máximo deve ser positivo.")
        limit = self.max_evaluations
        if isinstance(limit, bool) or not isinstance(limit, int) or limit < 100:
            raise InvalidParameterError("O limite de avaliações deve ser um inteiro maior ou igual a 100.")

    def describe(self) -> str:
        return (
            f"{self.method}, rtol = {self.rtol:g}, atol = {self.atol_position:g} m / {self.atol_velocity:g} m/s"
            + ("" if self.max_step == math.inf else f", passo máximo = {self.max_step:g} s")
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "metodo": self.method,
            "rtol": self.rtol,
            "atol_posicao_m": self.atol_position,
            "atol_velocidade_m_s": self.atol_velocity,
            "passo_maximo_s": None if self.max_step == math.inf else self.max_step,
            "limite_avaliacoes": self.max_evaluations,
        }


class Termination(Enum):
    """Motivo do fim da integração."""

    COMPLETED = "duração completa"
    COLLISION = "colisão com a superfície do corpo central"


def _readonly(array: np.ndarray) -> np.ndarray:
    array = np.ascontiguousarray(array, dtype=np.float64)
    array.flags.writeable = False
    return array


@dataclass(frozen=True, eq=False)
class Trajectory:
    """Trajetória simulada: amostras de tempo (s), posição (m) e velocidade (m/s).

    Posições e velocidades são relativas ao centro do corpo central (modelo de dois corpos)
    ou à origem do referencial do modelo. Os arrays são somente leitura.
    """

    times: np.ndarray
    positions: np.ndarray
    velocities: np.ndarray
    model: str
    settings: IntegratorSettings
    hypotheses: tuple[str, ...]
    evaluations: int
    termination: Termination
    origin: DataOrigin = field(default=DataOrigin.SIMULATED)

    def __post_init__(self) -> None:
        object.__setattr__(self, "times", _readonly(self.times))
        object.__setattr__(self, "positions", _readonly(self.positions))
        object.__setattr__(self, "velocities", _readonly(self.velocities))

    @property
    def duration(self) -> float:
        return float(self.times[-1] - self.times[0])

    @property
    def final_position(self) -> np.ndarray:
        return self.positions[-1]

    @property
    def final_velocity(self) -> np.ndarray:
        return self.velocities[-1]


class _EvaluationLimitReached(Exception):
    pass


_NUMERICAL_HYPOTHESIS = (
    "Integração numérica com passo adaptativo: o erro local é limitado pelas tolerâncias do "
    "integrador; o erro global acumula com o tempo de simulação (ver a verificação contra Kepler)."
)


def propagate(
    model: GravityModel,
    position: object,
    velocity: object,
    duration: float,
    *,
    settings: IntegratorSettings | None = None,
    n_points: int = 500,
    t0: float = 0.0,
    collision_radius: float | None = None,
) -> Trajectory:
    """Integra a trajetória por ``duration`` segundos a partir do estado inicial dado.

    Retorna ``n_points`` amostras igualmente espaçadas no tempo (incluindo os instantes
    inicial e final). Se ``collision_radius`` (m) for informado, a integração termina quando
    a distância à origem cai a esse valor (colisão com a superfície), e o instante da colisão
    é acrescentado como última amostra.
    """
    from scipy.integrate import solve_ivp  # importado só aqui: a CLI parte rápido sem o SciPy

    if not isinstance(model, GravityModel):
        raise InvalidParameterError("O modelo deve implementar acceleration(position, velocity, time).")
    settings = settings if settings is not None else IntegratorSettings()
    r0 = as_vector("position", position)
    v0 = as_vector("velocity", velocity)
    duration = as_positive("duration", duration)
    t0 = as_real("t0", t0)
    if isinstance(n_points, bool) or not isinstance(n_points, int) or not 2 <= n_points <= _MAX_POINTS:
        raise InvalidParameterError(f"O número de amostras deve ser um inteiro entre 2 e {_MAX_POINTS}.")
    t_end = t0 + duration
    if not math.isfinite(t_end) or t_end == t0:
        raise InvalidParameterError("Intervalo de tempo não representável.")

    events = None
    if collision_radius is not None:
        radius = as_positive("collision_radius", collision_radius)
        if math.hypot(*r0) <= radius:
            raise DomainError("A posição inicial está dentro (ou na superfície) do corpo central.")

        def collision(t: float, y: np.ndarray) -> float:
            return math.hypot(y[0], y[1], y[2]) - radius

        collision.terminal = True  # type: ignore[attr-defined]
        collision.direction = -1  # type: ignore[attr-defined]
        events = [collision]

    evaluations = 0

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        nonlocal evaluations
        evaluations += 1
        if evaluations > settings.max_evaluations:
            raise _EvaluationLimitReached
        acceleration = model.acceleration(y[:3], y[3:], t)
        return np.concatenate((y[3:], acceleration))

    t_eval = np.linspace(t0, t_end, n_points)
    t_eval[-1] = t_end
    atol = np.array([settings.atol_position] * 3 + [settings.atol_velocity] * 3)
    try:
        with np.errstate(over="raise", invalid="raise"):
            solution = solve_ivp(
                rhs,
                (t0, t_end),
                np.concatenate((r0, v0)),
                method=settings.method,
                t_eval=t_eval,
                rtol=settings.rtol,
                atol=atol,
                max_step=settings.max_step,
                events=events,
            )
    except AuroraError as exc:
        raise SimulationError(f"A simulação foi interrompida: {exc}") from None
    except FloatingPointError:
        raise SimulationError("A simulação produziu valores fora do intervalo representável.") from None
    except _EvaluationLimitReached:
        raise SimulationError(
            f"Limite de {settings.max_evaluations} avaliações da função atingido antes do fim da simulação: "
            "reduza a duração, aumente a tolerância ou eleve o limite."
        ) from None
    if solution.status == -1:
        raise SimulationError(f"O integrador numérico falhou: {solution.message}")

    times = solution.t
    states = solution.y.T
    termination = Termination.COMPLETED
    if solution.status == 1 and events is not None and len(solution.t_events[0]) > 0:
        termination = Termination.COLLISION
        t_hit = float(solution.t_events[0][0])
        if times.size == 0 or t_hit > times[-1]:
            times = np.append(times, t_hit)
            states = np.vstack((states, solution.y_events[0][0]))
    if not np.all(np.isfinite(states)):
        raise SimulationError("A simulação produziu valores não finitos.")

    describe = getattr(model, "describe", None)
    model_text = describe() if callable(describe) else type(model).__name__
    hypotheses = (*getattr(model, "hypotheses", ()), _NUMERICAL_HYPOTHESIS)
    return Trajectory(
        times=times,
        positions=states[:, :3],
        velocities=states[:, 3:],
        model=model_text,
        settings=settings,
        hypotheses=tuple(hypotheses),
        evaluations=int(solution.nfev),
        termination=termination,
    )


def concatenate(trajectories: list[Trajectory] | tuple[Trajectory, ...]) -> Trajectory:
    """Une segmentos consecutivos (ex.: arcos de uma missão) em uma só trajetória.

    Amostras com o mesmo instante na emenda são mantidas (antes e depois de uma queima
    impulsiva a velocidade muda, então ambas as amostras são informativas).
    """
    if not trajectories:
        raise InvalidParameterError("É necessário ao menos um segmento de trajetória.")
    first = trajectories[0]
    hypotheses: list[str] = []
    for segment in trajectories:
        hypotheses.extend(h for h in segment.hypotheses if h not in hypotheses)
    return Trajectory(
        times=np.concatenate([s.times for s in trajectories]),
        positions=np.vstack([s.positions for s in trajectories]),
        velocities=np.vstack([s.velocities for s in trajectories]),
        model=first.model,
        settings=first.settings,
        hypotheses=tuple(hypotheses),
        evaluations=sum(s.evaluations for s in trajectories),
        termination=trajectories[-1].termination,
    )
