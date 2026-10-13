"""Interface de linha de comando do AURORA (painel de cálculos da etapa 1).

Exemplos::

    python -m aurora quadratica 1 -3 2
    python -m aurora formulas
    python -m aurora formulas velocidade_orbital_circular
    python -m aurora calcular velocidade_orbital_circular mu=3.986004418e14 "r=6778.137 km"
    python -m aurora orbita-circular --corpo terra --altitude "400 km"
    python -m aurora corpos
    python -m aurora constantes

Toda saída numérica usa ponto como separador decimal (a entrada aceita ponto ou vírgula).
A opção ``--json`` exporta o resultado com entradas, unidades, hipóteses e procedência.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from collections.abc import Sequence
from fractions import Fraction
from typing import Any

from .. import (  # noqa: F401  (importar ``physics`` registra as fórmulas)
    __version__,
    physics,
)
from ..core import units
from ..core.formula import Formula, FormulaResult, get_formula, list_formulas
from ..core.provenance import DataOrigin
from ..data.bodies import get_body, list_bodies
from ..data.constants import ALL_CONSTANTS
from ..errors import AuroraError, InvalidParameterError
from ..physics.orbits import (
    CircularOrbitSummary,
    circular_orbit_summary,
    radius_from_altitude,
)
from ..physics.quadratic import QuadraticSolution, RootKind, solve_quadratic
from ._common import _fmt, _fmt_quantity, _print_footer, _print_json, _print_list
from .cli_simulation import add_simulation_commands

PROG = "aurora"


# --------------------------------------------------------------------------------------
# Formatação
# --------------------------------------------------------------------------------------


def _fmt_complex(z: complex) -> str:
    sign = "+" if z.imag >= 0 else "−"
    return f"{_fmt(z.real)} {sign} {_fmt(abs(z.imag))}·i"



# --------------------------------------------------------------------------------------
# Comandos
# --------------------------------------------------------------------------------------


_ROOT_KIND_CODES = {
    RootKind.TWO_REAL: "duas_reais_distintas",
    RootKind.DOUBLE: "raiz_dupla",
    RootKind.COMPLEX: "complexas_conjugadas",
}


_DISCRIMINANT_SIGN = {RootKind.TWO_REAL: 1, RootKind.DOUBLE: 0, RootKind.COMPLEX: -1}


def _fmt_discriminant(sol: QuadraticSolution) -> str:
    """Relação para Δ arredondado, sem sugerir 'infinito' ou 'zero' quando Δ só excede o float64."""
    d = sol.discriminant
    if math.isinf(d):
        bound = _fmt(sys.float_info.max)
        return f"{'> ' if d > 0 else '< −'}{bound} (excede o intervalo do float64; o sinal é exato)"
    if d == 0.0 and sol.kind is not RootKind.DOUBLE:
        relation = "> 0" if sol.kind is RootKind.TWO_REAL else "< 0"
        return f"{relation}, com módulo menor que o menor float64 (≈ 4.9e-324)"
    return f"= {_fmt(d)}"


def _cmd_quadratic(args: argparse.Namespace) -> int:
    sol: QuadraticSolution = solve_quadratic(args.a, args.b, args.c)
    if args.json:
        roots = [
            {"real": r.real, "imag": r.imag} if isinstance(r, complex) else {"real": r, "imag": 0.0}
            for r in sol.roots
        ]
        _print_json(
            {
                "equacao": "a*x^2 + b*x + c = 0",
                "coeficientes": {"a": sol.a, "b": sol.b, "c": sol.c},
                "discriminante": sol.discriminant,
                "sinal_discriminante": _DISCRIMINANT_SIGN[sol.kind],
                "classificacao": _ROOT_KIND_CODES[sol.kind],
                "descricao": sol.kind.value,
                "raizes": roots,
                "origem": DataOrigin.ANALYTIC.value,
            }
        )
        return 0

    print(f"Equação: ({_fmt(sol.a)})·x² + ({_fmt(sol.b)})·x + ({_fmt(sol.c)}) = 0")
    print(f"Discriminante: Δ = b² − 4·a·c {_fmt_discriminant(sol)}")
    print(f"Classificação: {sol.kind.value}")
    if sol.kind is RootKind.DOUBLE:
        print(f"  x₁ = x₂ = {_fmt(sol.roots[0])}")  # type: ignore[arg-type]
    else:
        for label, root in zip(("x₁", "x₂"), sol.roots, strict=True):
            text = _fmt_complex(root) if isinstance(root, complex) else _fmt(root)
            print(f"  {label} = {text}")
    print(
        "Método: Δ calculado em aritmética racional exata a partir dos coeficientes digitados "
        "(decimais lidos sem arredondamento); raízes pela forma numericamente estável "
        "q = −(b + sinal(b)·√Δ)/2, x = q/a e x = c/q, arredondadas uma única vez."
    )
    return 0


def _cmd_formulas(args: argparse.Namespace) -> int:
    if args.chave is None:
        formulas = list_formulas()
        width = max(len(f.key) for f in formulas)
        print("Fórmulas disponíveis (use 'formulas <chave>' para detalhes):")
        for f in formulas:
            print(f"  {f.key:<{width}}  {f.expression}  [{f.unit}]  — {f.title}")
        return 0
    _print_formula_details(get_formula(args.chave))
    return 0


def _print_formula_details(f: Formula) -> None:
    print(f"{f.title}  ({f.key})")
    print(f"  Expressão: {f.expression}")
    print(f"  Unidade do resultado: {f.unit}")
    print("  Parâmetros:")
    for p in f.parameters:
        print(f"    {p.name} [{p.unit}]: {p.description}")
    _print_list("  Hipóteses:", f.hypotheses)
    if f.references:
        _print_list("  Referências:", f.references)


def _parse_assignments(formula: Formula, items: Sequence[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for item in items:
        name, sep, raw = item.partition("=")
        name = name.strip()
        if not sep or not name:
            raise InvalidParameterError(f"Entrada inválida {item.strip()!r}: use o formato NOME=VALOR, ex.: r=6778km.")
        if name in values:
            raise InvalidParameterError(f"Parâmetro '{name}' informado mais de uma vez.")
        param = formula.parameter(name)
        values[name] = units.parse_quantity(raw, param.unit)
    return values


def _print_result(result: FormulaResult, out_unit: str | None = None) -> None:
    f = result.formula
    print(f.title)
    print(f"  Expressão: {f.expression}")
    print("  Entradas (convertidas para o SI):")
    width = max(len(p.name) for p in f.parameters)
    for p in f.parameters:
        print(f"    {p.name:<{width}} = {_fmt(result.inputs[p.name])} {p.unit}  ({p.description})")
    print(f"  Resultado: {f.symbol} = {_fmt_quantity(result.value, f.unit)}")
    if out_unit is not None:
        dim = units.get_unit(f.unit).dimension
        converted = units.from_si(result.value, out_unit, expected=dim)
        print(f"             {f.symbol} = {_fmt(converted)} {units.get_unit(out_unit).symbol}")
    _print_list("  Hipóteses:", f.hypotheses)
    _print_footer(result.origin)


def _cmd_calculate(args: argparse.Namespace) -> int:
    formula = get_formula(args.formula)
    inputs = _parse_assignments(formula, args.valores)
    result = formula.evaluate(**inputs)
    if args.em is not None:
        args.em = args.em.strip()
        # valida a unidade e a conversão (dimensão e overflow) antes de qualquer saída
        units.from_si(result.value, args.em, expected=units.get_unit(formula.unit).dimension)
    if args.json:
        data = result.to_dict()
        if args.em is not None:
            unit = units.get_unit(args.em)
            data["valor_convertido"] = {"valor": units.from_si(result.value, unit.symbol), "unidade": unit.symbol}
        _print_json(data)
        return 0
    _print_result(result, args.em)
    return 0


def _cmd_circular_orbit(args: argparse.Namespace) -> int:
    body = get_body(args.corpo)
    if args.raio is not None:
        radius = units.parse_quantity(args.raio, "m")
    else:
        radius = radius_from_altitude(body, units.parse_quantity(args.altitude, "m"))
    summary: CircularOrbitSummary = circular_orbit_summary(body, radius)

    if args.json:
        _print_json(
            {
                "corpo": {
                    "chave": body.key,
                    "nome": body.name,
                    "mu": {"valor": body.mu, "unidade": "m^3/s^2"},
                    "raio_referencia": {"valor": body.radius, "unidade": "m", "tipo": body.radius_description},
                    "fontes": list(body.sources),
                },
                "altitude": {"valor": summary.altitude, "unidade": "m"},
                "raio_orbital": {"valor": summary.radius, "unidade": "m"},
                "resultados": [r.to_dict() for r in summary.results],
                "origem": DataOrigin.ANALYTIC.value,
                "aviso": DataOrigin.ANALYTIC.disclaimer,
            }
        )
        return 0

    print(f"Órbita circular em torno de: {body.name}")
    print(f"  μ = {_fmt(body.mu)} m^3/s^2; {body.radius_description} R = {_fmt_quantity(body.radius, 'm')}")
    _print_list("  Fontes dos dados do corpo:", body.sources)
    print(f"  Altitude h = {_fmt_quantity(summary.altitude, 'm')}")
    print(f"  Raio orbital r = R + h = {_fmt_quantity(summary.radius, 'm')}")
    print()
    width = max(len(r.formula.title) for r in summary.results)
    for r in summary.results:
        print(f"  {r.formula.title:<{width}}  {r.formula.symbol} = {_fmt_quantity(r.value, r.unit)}")
    print()
    hypotheses: list[str] = []
    for r in summary.results:
        hypotheses.extend(h for h in r.formula.hypotheses if h not in hypotheses)
    _print_list("Hipóteses:", hypotheses)
    _print_footer(DataOrigin.ANALYTIC)
    return 0


def _cmd_bodies(args: argparse.Namespace) -> int:
    print("Corpos celestes de referência (modelados como esferas):")
    for body in list_bodies():
        print(f"  {body.key}: {body.name}")
        print(f"    μ = {_fmt(body.mu)} m^3/s^2")
        print(f"    {body.radius_description} = {_fmt_quantity(body.radius, 'm')}")
        for source in body.sources:
            print(f"    Fonte: {source}")
    return 0


def _cmd_constants(args: argparse.Namespace) -> int:
    print("Constantes físicas utilizadas:")
    for c in ALL_CONSTANTS:
        status = "exata" if c.exact else f"incerteza relativa ≈ {c.relative_uncertainty:.1e}"
        print(f"  {c.symbol} = {_fmt(c.value, 12)} {c.unit}  — {c.name} ({status})")
        print(f"    Fonte: {c.source}")
    return 0


# --------------------------------------------------------------------------------------
# Analisador de argumentos
# --------------------------------------------------------------------------------------


class _Formatter(argparse.RawDescriptionHelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):  # type: ignore[no-untyped-def]
        super().add_usage(usage, actions, groups, prefix if prefix is not None else "uso: ")


# Tradução das mensagens geradas pelo próprio argparse (em inglês). Mensagens sem
# correspondência são exibidas sem alteração.
_ARGPARSE_MESSAGES = (
    (re.compile(r"^the following arguments are required: (.+)$"), r"argumentos obrigatórios ausentes: \1"),
    (re.compile(r"^unrecognized arguments: (.+)$"), r"argumentos não reconhecidos: \1"),
    (
        re.compile(r"^argument (.+?): invalid choice: (.+?) \(choose from (.+)\)$"),
        r"argumento \1: opção inválida: \2 (escolha entre \3)",
    ),
    (re.compile(r"^argument (.+?): expected one argument$"), r"argumento \1: esperava um valor"),
    (re.compile(r"^argument (.+?): not allowed with argument (.+)$"), r"argumento \1: não pode ser usado com \2"),
    (re.compile(r"^one of the arguments (.+) is required$"), r"informe um dos argumentos: \1"),
    (re.compile(r"^ambiguous option: (.+?) could match (.+)$"), r"opção ambígua: \1 pode corresponder a \2"),
    (re.compile(r"^argument (.+?): (.+)$"), r"argumento \1: \2"),
)


def _translate_argparse(message: str) -> str:
    for pattern, replacement in _ARGPARSE_MESSAGES:
        if pattern.match(message):
            return pattern.sub(replacement, message)
    return message


class _Parser(argparse.ArgumentParser):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("formatter_class", _Formatter)
        kwargs.setdefault("add_help", False)
        super().__init__(*args, **kwargs)
        self._positionals.title = "argumentos"
        self._optionals.title = "opções"
        self.add_argument("-h", "--ajuda", "--help", action="help", help="mostra esta ajuda e sai")

    def error(self, message: str) -> None:  # type: ignore[override]
        self.print_usage(sys.stderr)
        self.exit(2, f"{self.prog}: erro: {_translate_argparse(message)}\n")


def _number(text: str) -> Fraction:
    """Coeficiente numérico com o valor decimal exato digitado (``0.1`` é exatamente 1/10).

    Valores que não são números normais em ``float64`` usam o valor do ``float`` (o mesmo que é
    exibido), de modo que a equação mostrada é sempre a equação resolvida; valores não nulos
    abaixo do menor ``float64`` são rejeitados. Isso também limita o custo da aritmética exata.
    """
    try:
        value = units.parse_number(text)  # valida o formato e o intervalo representável
    except AuroraError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from None
    text = text.strip().replace(",", ".")
    if value == 0.0:
        if any(ch in "123456789" for ch in re.split("[eE]", text)[0]):
            raise argparse.ArgumentTypeError(
                f"Valor {text!r} é menor, em módulo, que o menor número representável em float64 (≈ 4,9e-324)."
            )
        return Fraction(0)
    if abs(value) < sys.float_info.min:  # subnormal: usa o valor exibido
        return Fraction(value)
    try:
        return Fraction(text)
    except ValueError:  # ex.: mais dígitos do que o limite de conversão de inteiros do Python
        return Fraction(value)


def build_parser() -> argparse.ArgumentParser:
    parser = _Parser(
        prog=PROG,
        description=(
            "AURORA — cálculos físicos e simulação orbital (etapas 1 e 2).\n"
            "Todos os resultados são calculados com modelos idealizados; nenhum é medição."
        ),
    )
    parser.add_argument(
        "--versao", "--version", action="version", version=f"AURORA {__version__}", help="mostra a versão e sai"
    )
    sub = parser.add_subparsers(dest="comando", title="comandos", metavar="COMANDO", parser_class=_Parser)
    sub.required = True

    p = sub.add_parser("quadratica", help="resolve a·x² + b·x + c = 0", description="Resolve a·x² + b·x + c = 0.")
    p.add_argument("a", type=_number, help="coeficiente de x² (não nulo)")
    p.add_argument("b", type=_number, help="coeficiente de x")
    p.add_argument("c", type=_number, help="termo independente")
    p.add_argument("--json", action="store_true", help="saída em JSON")
    p.set_defaults(handler=_cmd_quadratic)

    p = sub.add_parser("formulas", help="lista as fórmulas ou mostra detalhes de uma delas")
    p.add_argument("chave", nargs="?", help="chave da fórmula (opcional)")
    p.set_defaults(handler=_cmd_formulas)

    p = sub.add_parser(
        "calcular",
        help="avalia uma fórmula registrada",
        description=(
            "Avalia uma fórmula. Cada parâmetro é informado como NOME=VALOR, com unidade opcional\n"
            "(sem unidade, assume-se a unidade SI do parâmetro). Ex.:\n"
            '  aurora calcular velocidade_orbital_circular mu=398600.4418km^3/s^2 "r=6778.137 km"'
        ),
    )
    p.add_argument("formula", help="chave da fórmula (veja o comando 'formulas')")
    p.add_argument("valores", nargs="*", metavar="NOME=VALOR", help="parâmetros, ex.: r=6778.137km")
    p.add_argument("--em", metavar="UNIDADE", help="exibe também o resultado nesta unidade")
    p.add_argument("--json", action="store_true", help="saída em JSON")
    p.set_defaults(handler=_cmd_calculate)

    p = sub.add_parser(
        "orbita-circular",
        help="grandezas de uma órbita circular em torno de um corpo",
        description="Velocidade, período, energia, velocidade de escape e gravidade de uma órbita circular.",
    )
    p.add_argument("--corpo", default="terra", help="terra, lua ou sol (padrão: terra)")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--altitude", metavar="GRANDEZA", help='altitude acima do raio de referência, ex.: "400 km"')
    group.add_argument("--raio", metavar="GRANDEZA", help='raio orbital a partir do centro, ex.: "42164 km"')
    p.add_argument("--json", action="store_true", help="saída em JSON")
    p.set_defaults(handler=_cmd_circular_orbit)

    p = sub.add_parser("corpos", help="lista os corpos celestes de referência e suas fontes")
    p.set_defaults(handler=_cmd_bodies)

    p = sub.add_parser("constantes", help="lista as constantes físicas e suas fontes")
    p.set_defaults(handler=_cmd_constants)

    add_simulation_commands(sub)

    return parser


_NEGATIVE_NUMBER = re.compile(r"^-(?:\d|[.,]\d)")


def _protect_negative_numbers(argv: Sequence[str]) -> list[str]:
    """Evita que valores como ``-1e8`` sejam interpretados como opções pelo argparse."""
    return [" " + arg if _NEGATIVE_NUMBER.match(arg) else arg for arg in argv]


def _configure_output() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")  # type: ignore[union-attr]
        except (AttributeError, ValueError):
            pass


def _silence_closed_stdout() -> None:
    """Redireciona a saída padrão para /dev/null depois que o leitor fechou o pipe.

    Sem isso, o Python tentaria descarregar o buffer ao encerrar e imprimiria outro erro
    (receita da documentação do módulo ``signal``, seção "Note on SIGPIPE").
    """
    try:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
    except (AttributeError, OSError, ValueError):  # saída capturada (testes) ou sem descritor
        pass


def main(argv: Sequence[str] | None = None) -> int:
    """Ponto de entrada. Retorna 0 em caso de sucesso e 1 em erro de cálculo ou de entrada.

    Como no ``argparse``, erros de uso (código 2) e as opções ``--ajuda``/``--versao``
    (código 0) encerram a execução levantando ``SystemExit``. Se a saída for canalizada para
    um programa que encerra antes de ler tudo (ex.: ``aurora ... | head``), a execução termina
    em silêncio com código 1, sem rastreamento de pilha.
    """
    try:
        code = _main(argv)
        sys.stdout.flush()  # erros de escrita aparecem aqui, e não no encerramento do interpretador
        return code
    except BrokenPipeError:
        _silence_closed_stdout()
        return 1


def _main(argv: Sequence[str] | None) -> int:
    _configure_output()
    parser = build_parser()
    args, extras = parser.parse_known_args(_protect_negative_numbers(sys.argv[1:] if argv is None else argv))
    if extras and getattr(args, "comando", None) == "calcular":
        # Antes do Python 3.13, uma opção entre a chave da fórmula e os valores (ex.:
        # "calcular velocidade_escape --json mu=1 r=1") deixa os valores como "não reconhecidos".
        args.valores = [*args.valores, *(e for e in extras if not e.startswith("-"))]
        extras = [e for e in extras if e.startswith("-") and e != "--"]
    if extras:
        parser.error(f"unrecognized arguments: {' '.join(e.strip() for e in extras)}")
    try:
        return args.handler(args)
    except AuroraError as exc:
        print(f"{PROG}: erro: {exc}", file=sys.stderr)
        return 1
