"""Propagação analítica do problema de dois corpos (formulação por variável universal).

Dado o estado (r⃗₀, v⃗₀) de um veículo em torno de um corpo de parâmetro gravitacional μ, a
posição e a velocidade após um intervalo Δt são obtidas sem integração numérica, pelos
coeficientes de Lagrange (Curtis, *Orbital Mechanics for Engineering Students*, alg. 3.3 e
3.4; Vallado, *Fundamentals of Astrodynamics and Applications*, alg. 8)::

    α = 2/r₀ − v₀²/μ                       (inverso do semieixo maior; vale para toda cônica)
    F(χ) = (r₀·v_r0/√μ)·χ²·C(z) + (1 − α·r₀)·χ³·S(z) + r₀·χ − √μ·Δt = 0,   z = α·χ²
    f = 1 − χ²·C(z)/r₀             g = Δt − χ³·S(z)/√μ
    ḟ = √μ·(α·χ³·S(z) − χ)/(r·r₀)  ġ = 1 − χ²·C(z)/r
    r⃗ = f·r⃗₀ + g·v⃗₀               v⃗ = ḟ·r⃗₀ + ġ·v⃗₀

C(z) e S(z) são as funções de Stumpff. A mesma formulação cobre órbitas elípticas,
parabólicas e hiperbólicas, inclusive trajetórias radiais.

Cuidados numéricos:

* para |z| < 0,1, C e S são calculadas por série de Taylor (a forma fechada
  ``(1 − cos √z)/z`` perderia algarismos por cancelamento);
* F(χ) é estritamente crescente (F′(χ) = r > 0), por isso a raiz é encontrada por Newton
  protegido por um intervalo que sempre contém a solução (com bisseção de reserva);
* em órbitas elípticas o intervalo Δt é reduzido pelo número inteiro de períodos antes da
  solução, pois o movimento é periódico; isso mantém χ limitado em propagações longas.

Este é um cálculo **analítico** do modelo idealizado de dois corpos; serve de referência
para validar o integrador numérico (:mod:`aurora.simulation`).
"""

from __future__ import annotations

import math
import sys

import numpy as np

from ..core.validation import as_non_negative, as_positive, as_real, as_vector
from ..errors import DomainError, InvalidParameterError

_SERIES_LIMIT = 0.1
_SERIES_TERMS = 12
_EPS = sys.float_info.epsilon


def _stumpff_series(z: float, start: int) -> float:
    """Σ (−z)^k / (start + 2k)! — série de Taylor de C (start = 2) e S (start = 3)."""
    term = 1.0 / math.factorial(start)
    total = term
    for n in range(1, _SERIES_TERMS):
        term *= -z / ((start + 2 * n - 1) * (start + 2 * n))
        total += term
    return total


def stumpff_c(z: float) -> float:
    """Função de Stumpff C(z) = (1 − cos √z)/z (z > 0), (cosh √−z − 1)/(−z) (z < 0), 1/2 (z = 0)."""
    if abs(z) < _SERIES_LIMIT:
        return _stumpff_series(z, 2)
    if z > 0.0:
        return (1.0 - math.cos(math.sqrt(z))) / z
    return (math.cosh(math.sqrt(-z)) - 1.0) / (-z)


def stumpff_s(z: float) -> float:
    """Função de Stumpff S(z) = (√z − sin √z)/√z³ (z > 0), (sinh √−z − √−z)/√−z³ (z < 0), 1/6 (z = 0)."""
    if abs(z) < _SERIES_LIMIT:
        return _stumpff_series(z, 3)
    if z > 0.0:
        s = math.sqrt(z)
        return (s - math.sin(s)) / (s * s * s)
    s = math.sqrt(-z)
    return (math.sinh(s) - s) / (s * s * s)


def _solve_universal_anomaly(r0: float, vr0: float, alpha: float, sqrt_mu: float, dt: float) -> float:
    """Raiz χ de F(χ) = 0 (ver docstring do módulo); χ tem o mesmo sinal de Δt."""
    a1 = r0 * vr0 / sqrt_mu
    a2 = 1.0 - alpha * r0
    target = sqrt_mu * dt

    def residual_and_noise(chi: float) -> tuple[float, float]:
        """F(χ) e uma estimativa do erro de arredondamento ao avaliá-la."""
        z = alpha * chi * chi
        try:
            terms = (a1 * chi * chi * stumpff_c(z), a2 * chi**3 * stumpff_s(z), r0 * chi, -target)
        except OverflowError:  # cosh/sinh de argumento enorme: muito além da raiz
            return math.copysign(math.inf, chi), 0.0
        return sum(terms), 8.0 * _EPS * sum(abs(t) for t in terms)

    def residual(chi: float) -> float:
        return residual_and_noise(chi)[0]

    def derivative(chi: float) -> float:  # F′(χ) = raio no instante correspondente (> 0)
        z = alpha * chi * chi
        return a1 * chi * (1.0 - z * stumpff_s(z)) + a2 * chi * chi * stumpff_c(z) + r0

    # Intervalo inicial [lo, hi] com F(lo) ≤ 0 ≤ F(hi), expandido por duplicação.
    guess = target / r0  # primeira ordem: F ≈ r₀·χ − √μ·Δt
    if dt > 0.0:
        lo, hi = 0.0, max(guess, _EPS)
        while residual(hi) < 0.0:
            lo, hi = hi, 2.0 * hi
            if not math.isfinite(hi):
                raise DomainError("Intervalo de propagação fora do alcance do cálculo em ponto flutuante.")
    else:
        lo, hi = min(guess, -_EPS), 0.0
        while residual(lo) > 0.0:
            lo, hi = 2.0 * lo, lo
            if not math.isfinite(lo):
                raise DomainError("Intervalo de propagação fora do alcance do cálculo em ponto flutuante.")

    # Newton protegido (estratégia "rtsafe", Press et al., Numerical Recipes, §9.4): usa bisseção
    # quando o passo de Newton sairia do intervalo ou não o reduziria rápido o bastante — o que
    # acontece longe da raiz em hipérboles, onde F cresce exponencialmente com χ.
    chi = 0.5 * (lo + hi)
    step = step_old = hi - lo
    for _ in range(300):
        value, noise = residual_and_noise(chi)
        if abs(value) <= noise:  # resíduo já no nível do erro de arredondamento de F
            return chi
        if value < 0.0:
            lo = chi
        else:
            hi = chi
        try:
            slope = derivative(chi)
        except OverflowError:
            slope = math.inf
        newton_ok = (
            math.isfinite(slope)
            and slope > 0.0
            and lo < chi - value / slope < hi
            and abs(2.0 * value) <= abs(step_old * slope)
        )
        step_old = step
        if newton_ok:
            step = value / slope
            candidate = chi - step
        else:
            step = 0.5 * (hi - lo)
            candidate = lo + step
        if abs(candidate - chi) <= 4.0 * _EPS * max(abs(candidate), 1.0) or hi - lo <= 4.0 * _EPS * max(
            abs(lo), abs(hi)
        ):
            return candidate
        chi = candidate
    raise DomainError("A equação de Kepler universal não convergiu.")


def propagate_kepler(mu: float, position: object, velocity: object, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Estado (m, m/s) após ``dt`` segundos no problema de dois corpos (``dt`` pode ser negativo).

    ``position`` e ``velocity`` são relativos ao centro do corpo central, em referencial
    inercial. Retorna ``(r, v)`` como vetores ``float64``.
    """
    mu = as_positive("mu", mu)
    r0_vec = as_vector("position", position)
    v0_vec = as_vector("velocity", velocity)
    dt = as_real("dt", dt)
    r0 = math.hypot(*r0_vec)
    if r0 == 0.0:
        raise DomainError("A posição inicial não pode coincidir com o centro do corpo central.")
    if dt == 0.0:
        return r0_vec, v0_vec

    sqrt_mu = math.sqrt(mu)
    v0_sq = float(np.dot(v0_vec, v0_vec))
    vr0 = float(np.dot(r0_vec, v0_vec)) / r0
    alpha = 2.0 / r0 - v0_sq / mu
    if alpha > 0.0:  # elipse: o movimento é periódico
        period = 2.0 * math.pi / (sqrt_mu * alpha * math.sqrt(alpha))
        dt = math.fmod(dt, period)
        if dt == 0.0:
            return r0_vec, v0_vec

    chi = _solve_universal_anomaly(r0, vr0, alpha, sqrt_mu, dt)
    z = alpha * chi * chi
    try:
        c, s = stumpff_c(z), stumpff_s(z)
        f = 1.0 - chi * chi / r0 * c
        g = dt - chi**3 / sqrt_mu * s
        r_vec = f * r0_vec + g * v0_vec
        r = math.hypot(*r_vec)
        if r == 0.0:
            raise DomainError("A trajetória passa pelo centro do corpo central (colisão).")
        f_dot = sqrt_mu / (r * r0) * (alpha * chi**3 * s - chi)
        g_dot = 1.0 - chi * chi / r * c
    except OverflowError:
        raise DomainError("Estado propagado fora do intervalo representável em ponto flutuante.") from None
    v_vec = f_dot * r0_vec + g_dot * v0_vec
    if not (np.all(np.isfinite(r_vec)) and np.all(np.isfinite(v_vec))):
        raise DomainError("Estado propagado fora do intervalo representável em ponto flutuante.")
    return r_vec, v_vec


def state_from_elements(
    mu: float,
    a: float,
    e: float,
    inclination: float = 0.0,
    raan: float = 0.0,
    arg_periapsis: float = 0.0,
    true_anomaly: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Vetor de estado (m, m/s) a partir dos elementos orbitais clássicos.

    ``a`` em m (positivo para elipse, negativo para hipérbole), ``e`` ≥ 0 e diferente de 1
    (a parábola não é definida pelo semieixo maior), ângulos em radianos: inclinação,
    longitude do nodo ascendente (``raan``), argumento do periapse e anomalia verdadeira.
    Usa as fórmulas no referencial perifocal e a rotação 3-1-3 (Ω, i, ω).
    """
    mu = as_positive("mu", mu)
    a = as_real("a", a)
    e = as_non_negative("e", e)
    angles = [as_real(n, v) for n, v in (("inclination", inclination), ("raan", raan),
                                          ("arg_periapsis", arg_periapsis), ("true_anomaly", true_anomaly))]
    inc, big_omega, small_omega, nu = angles
    if e == 1.0:
        raise InvalidParameterError("Órbita parabólica (e = 1) não pode ser definida pelo semieixo maior.")
    p = a * (1.0 - e * e)
    if not p > 0.0:
        raise InvalidParameterError(
            "Combinação inválida: use a > 0 com 0 ≤ e < 1 (elipse) ou a < 0 com e > 1 (hipérbole)."
        )
    denom = 1.0 + e * math.cos(nu)
    if denom <= 0.0:
        raise DomainError("Anomalia verdadeira fora do ramo da hipérbole (1 + e·cos ν ≤ 0).")
    r = p / denom
    r_pf = np.array((r * math.cos(nu), r * math.sin(nu), 0.0))
    v_pf = math.sqrt(mu / p) * np.array((-math.sin(nu), e + math.cos(nu), 0.0))
    c_o, s_o = math.cos(big_omega), math.sin(big_omega)
    c_i, s_i = math.cos(inc), math.sin(inc)
    c_w, s_w = math.cos(small_omega), math.sin(small_omega)
    rotation = np.array(
        (
            (c_o * c_w - s_o * s_w * c_i, -c_o * s_w - s_o * c_w * c_i, s_o * s_i),
            (s_o * c_w + c_o * s_w * c_i, -s_o * s_w + c_o * c_w * c_i, -c_o * s_i),
            (s_w * s_i, c_w * s_i, c_i),
        )
    )
    return rotation @ r_pf, rotation @ v_pf
