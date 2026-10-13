"""Testes da calculadora de equações do segundo grau."""

import math
import random
from fractions import Fraction

import pytest

from aurora.errors import DomainError, InvalidParameterError
from aurora.physics.quadratic import RootKind, discriminant, solve_quadratic


def test_duas_raizes_reais_distintas():
    sol = solve_quadratic(1, -3, 2)
    assert sol.kind is RootKind.TWO_REAL
    assert sol.discriminant == 1.0
    assert sol.roots == (1.0, 2.0)
    assert sol.real_roots == (1.0, 2.0)


def test_raizes_reais_em_ordem_crescente_com_a_negativo():
    sol = solve_quadratic(-2, 0, 8)  # -2x² + 8 = 0 -> x = ±2
    assert sol.roots == (-2.0, 2.0)


def test_raiz_dupla():
    sol = solve_quadratic(1, -4, 4)
    assert sol.kind is RootKind.DOUBLE
    assert sol.discriminant == 0.0
    assert sol.roots == (2.0, 2.0)
    assert sol.real_roots == (2.0,)


def test_raiz_dupla_em_zero_nao_e_menos_zero():
    sol = solve_quadratic(3, 0, 0)
    assert sol.kind is RootKind.DOUBLE
    assert math.copysign(1.0, sol.roots[0]) == 1.0


def test_raizes_complexas_conjugadas():
    sol = solve_quadratic(1, 2, 5)
    assert sol.kind is RootKind.COMPLEX
    assert sol.discriminant == -16.0
    assert sol.roots == (complex(-1, 2), complex(-1, -2))
    assert sol.real_roots == ()


def test_raizes_satisfazem_a_equacao():
    for a, b, c in [(1, -3, 2), (2, 1, -7), (1, 2, 5), (-0.5, 3.25, 1.5)]:
        sol = solve_quadratic(a, b, c)
        for x in sol.roots:
            assert abs(a * x * x + b * x + c) <= 1e-12 * max(abs(a * x * x), abs(b * x), abs(c), 1.0)


def test_relacoes_de_girard_vieta_em_coeficientes_aleatorios():
    rng = random.Random(20261013)
    for _ in range(2000):
        a = rng.uniform(-1e3, 1e3) or 1.0
        b = rng.uniform(-1e3, 1e3)
        c = rng.uniform(-1e3, 1e3)
        x1, x2 = solve_quadratic(a, b, c).roots
        scale = abs(b / a) + abs(complex(x1)) + abs(complex(x2))
        assert abs((x1 + x2) - (-b / a)) <= 1e-12 * scale
        assert abs(x1 * x2 - c / a) <= 1e-12 * (abs(c / a) + abs(complex(x1)) * abs(complex(x2)))


def test_cancelamento_catastrofico_evitado():
    # A fórmula ingênua (-b + sqrt(Δ))/(2a) daria ≈ -7.45e-9 (erro de ~25 %) para a raiz pequena.
    sol = solve_quadratic(1.0, 1e8, 1.0)
    small, large = sol.roots[1], sol.roots[0]
    assert small == pytest.approx(-1e-8, rel=1e-15)
    assert large == pytest.approx(-1e8, rel=1e-15)


def test_discriminante_exato_exemplo_de_kahan():
    # Exemplo clássico de W. Kahan: b'² − a·c = 1.890625 exatamente, mas b'² e a·c excedem
    # 2^53 e a subtração direta em float64 perde todos os algarismos corretos.
    a, b, c = 94906265.625, -2 * 94906267.0, 94906268.375
    assert b * b - 4 * a * c == 0.0  # conta ingênua: falsa raiz dupla
    assert discriminant(a, b, c) == 4 * 1.890625
    sol = solve_quadratic(a, b, c)
    assert sol.kind is RootKind.TWO_REAL
    assert sol.roots[0] == 1.0
    assert sol.roots[1] == pytest.approx(1.0 + 5.5 / 189812531.25, rel=1e-15)


def test_raizes_proximas_nao_sao_confundidas_com_raiz_dupla():
    # (x − (r−1))(x − (r+1)) = x² − 2r·x + (r² − 1), com coeficientes exatos em float64.
    # A conta ingênua b·b − 4·a·c dá 0 (falsa raiz dupla); o valor exato é Δ = 4.
    r = 101343649.0
    a, b, c = 1.0, -2 * r, r * r - 1.0
    assert b * b - 4 * a * c == 0.0
    sol = solve_quadratic(a, b, c)
    assert sol.kind is RootKind.TWO_REAL
    assert sol.discriminant == 4.0
    assert sol.roots == (r - 1, r + 1)


def test_raiz_dupla_exata():
    r = 2.0**25 + 1
    sol = solve_quadratic(3.0, -6.0 * r, 3.0 * r * r)
    assert sol.kind is RootKind.DOUBLE
    assert sol.roots == (r, r)


def test_coeficientes_muito_grandes_nao_causam_overflow():
    sol = solve_quadratic(1e200, -3e200, 2e200)
    assert sol.roots == pytest.approx((1.0, 2.0), rel=1e-15)
    assert sol.discriminant == math.inf  # Δ = 1e400 não cabe em float64; a classificação continua exata


def test_coeficientes_muito_pequenos_nao_causam_underflow():
    sol = solve_quadratic(1e-200, -3e-200, 2e-200)
    assert sol.roots == pytest.approx((1.0, 2.0), rel=1e-15)


def test_raiz_que_nao_cabe_em_float64_gera_erro():
    with pytest.raises(DomainError):
        solve_quadratic(1e-300, 1e300, 1.0)


def test_a_igual_a_zero_e_rejeitado():
    with pytest.raises(DomainError, match="linear"):
        solve_quadratic(0, 2, 1)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf, "1", None, True, 1 + 2j])
def test_coeficientes_invalidos(bad):
    with pytest.raises(InvalidParameterError):
        solve_quadratic(1, bad, 1)


# --- regressões da auditoria técnica (docs/AUDITORIA_ETAPA1.md) ---------------------------


@pytest.mark.parametrize(
    "coef, expected",
    [
        ((1e308, 0.0, 1e308), (1j, -1j)),  # NUME-1: antes, OverflowError
        ((1e308, 0.0, -1e308), (-1.0, 1.0)),
        ((-1.7e308, 1.7e308, 1.7e308), (-0.6180339887498948, 1.618033988749895)),
        ((5e-324, 5e-324, -5e-324), (-1.618033988749895, 0.6180339887498948)),  # NUME-2: antes, (-1, 1)
        ((5e-324, 5e-324, 0.0), (-1.0, 0.0)),  # NUME-2: antes, erro espúrio
    ],
)
def test_coeficientes_extremos_dao_raizes_corretas(coef, expected):
    sol = solve_quadratic(*coef)
    for got, exp in zip(sol.roots, expected, strict=True):
        assert got == pytest.approx(exp, rel=1e-15, abs=0.0)


def test_parte_real_de_raiz_complexa_com_b_subnormal():
    sol = solve_quadratic(1e-300, 5e-324, 1.0)
    assert sol.roots[0].real == pytest.approx(-2.4703282292062327e-24, rel=1e-15)
    assert sol.roots[0].imag == pytest.approx(1e150, rel=1e-15)


@pytest.mark.parametrize(
    "coef, root",
    [
        ((Fraction(1), Fraction("0.6"), Fraction("0.09")), -0.3),  # (x + 0,3)²
        ((1, Fraction("-0.2"), Fraction("0.01")), 0.1),  # (x − 0,1)²
        ((Fraction("0.01"), Fraction("0.2"), 1), -10.0),
    ],
)
def test_coeficientes_racionais_exatos_dao_raiz_dupla(coef, root):
    # ARCH-1: decimais lidos exatamente (como faz a interface) preservam a raiz dupla.
    sol = solve_quadratic(*coef)
    assert sol.kind is RootKind.DOUBLE
    assert sol.roots == (root, root)


def test_floats_binarios_sao_tratados_pelo_valor_exato_do_float():
    # 0.6 e 0.09 em float64 não são os decimais exatos; para esses binários Δ < 0 de fato.
    sol = solve_quadratic(1.0, 0.6, 0.09)
    assert sol.kind is RootKind.COMPLEX
    assert sol.discriminant == float(Fraction(0.6) ** 2 - 4 * Fraction(0.09))


def _reference_roots(a, b, c):
    """Raízes em precisão estendida (60 dígitos) pela forma estável, para comparação."""
    from decimal import Decimal, localcontext

    with localcontext() as ctx:
        ctx.prec = 60
        a, b, c = Decimal(a), Decimal(b), Decimal(c)
        d = b * b - 4 * a * c
        if d > 0:
            q = -(b + d.sqrt()) / 2 if b >= 0 else -(b - d.sqrt()) / 2
            return sorted((q / a, c / q))
        if d == 0:
            return [-b / (2 * a)] * 2
        return [-b / (2 * a), (-d).sqrt() / (2 * abs(a))]


def test_precisao_de_poucos_ulp_em_varias_escalas():
    # Normal, perto de 1,8e308 e subnormal: erro de no máximo 4 ulp em relação à referência.
    rng = random.Random(99)
    for lo, hi in ((-3, 3), (305, 307), (-322, -310)):
        for _ in range(300):
            a, b, c = (rng.choice((-1, 1)) * rng.uniform(1, 9.9) * 10.0 ** rng.randint(lo, hi) for _ in range(3))
            sol = solve_quadratic(a, b, c)
            got = [sol.roots[0].real, sol.roots[0].imag] if sol.kind is RootKind.COMPLEX else list(sol.roots)
            for x, ref in zip(got, _reference_roots(a, b, c), strict=True):
                assert abs(x - float(ref)) <= 4 * math.ulp(float(ref)), (a, b, c)


# --- regressões da revisão das correções ----------------------------------------------------


def test_inteiros_do_numpy_sao_tratados_com_precisao_arbitraria():
    # Fraction(np.int64) mantinha o inteiro de 64 bits: AttributeError ou estouro silencioso.
    np = pytest.importorskip("numpy")
    assert solve_quadratic(np.int64(1), np.int64(-5), np.int64(6)).roots == (2.0, 3.0)
    assert solve_quadratic(np.int64(2**31), np.int64(0), np.int64(-(2**31))).roots == (-1.0, 1.0)
    assert solve_quadratic(np.uint64(1), np.uint64(5), np.uint64(6)).roots == (-3.0, -2.0)
    assert discriminant(np.int64(1), np.int64(4_000_000_000), np.int64(1)) == 1.6e19


def test_raizes_nao_retornam_menos_zero():
    small = solve_quadratic(1.0, 1e10, 5e-324).roots[1]
    assert small == 0.0 and math.copysign(1.0, small) == 1.0
    z = solve_quadratic(1e10, 5e-324, 1e10).roots[0]
    assert math.copysign(1.0, z.real) == 1.0
