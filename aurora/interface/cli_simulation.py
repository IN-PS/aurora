"""Comandos de simulação orbital da interface de linha de comando (etapa 2).

* ``propagar``: integra numericamente uma órbita no modelo de dois corpos e mede o erro
  numérico contra a solução analítica de Kepler; exporta telemetria e gráficos.
* ``hohmann``: transferência de Hohmann analítica e, com ``--simular``, executada
  numericamente com consumo de propelente.
* ``validar``: executa os cenários padrão de validação do propagador.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import numpy as np

from ..core import units
from ..core.provenance import DataOrigin
from ..data.bodies import CelestialBody, get_body
from ..errors import AuroraError, DomainError, InvalidParameterError
from ..physics.gravitation import TwoBodyGravity
from ..physics.maneuvers import HohmannTransfer, hohmann_transfer
from ..physics.orbits import TwoBodyOrbit, orbit_from_state
from ..simulation.initial_conditions import circular_orbit_state, elliptic_orbit_state
from ..simulation.missions import HohmannMission, simulate_hohmann
from ..simulation.propagator import METHODS, IntegratorSettings, Trajectory, propagate
from ..simulation.telemetry import write_telemetry
from ..simulation.validation import (
    angular_momentum_series,
    kepler_error_series,
    run_validation,
    specific_energy_series,
    standard_cases,
)
from ..simulation.vehicle import Spacecraft
from ._common import _fmt, _fmt_quantity, _print_footer, _print_json, _print_list
from .plots import SUPPORTED_FORMATS, Marker, figure_for_trajectory, save_figure

_KM = 1e-3


# --------------------------------------------------------------------------------------
# Tipos de argumento
# --------------------------------------------------------------------------------------


def _positive_int(text: str) -> int:
    try:
        value = int(text.strip())
    except ValueError:
        raise argparse.ArgumentTypeError(f"Número inteiro inválido: {text.strip()!r}.") from None
    if value < 2:
        raise argparse.ArgumentTypeError("O número de amostras deve ser pelo menos 2.")
    return value


def _real(text: str) -> float:
    try:
        return units.parse_number(text)
    except AuroraError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from None


def _parse_vector(text: str, default_unit: str, name: str) -> np.ndarray:
    """Lê "x y z [unidade]" (componentes separados por espaço) e retorna o vetor em SI."""
    tokens = text.split()
    expected = units.get_unit(default_unit)
    unit = expected
    if len(tokens) == 4:
        unit = units.get_unit(tokens[3])
        if unit.dimension is not expected.dimension:
            raise InvalidParameterError(
                f"Unidade '{unit.symbol}' em {name} é de {unit.dimension.value}; "
                f"esperava-se {expected.dimension.value}."
            )
        tokens = tokens[:3]
    if len(tokens) != 3:
        raise InvalidParameterError(
            f"{name} deve ter três componentes separados por espaço, com unidade opcional; ex.: \"7000 0 0 km\"."
        )
    return np.array([units.to_si(units.parse_number(t), unit.symbol) for t in tokens])


def _check_output_paths(telemetry: str | None, graphs: list[str] | None) -> None:
    if telemetry is not None and Path(telemetry).suffix.lower() not in (".csv", ".json"):
        raise InvalidParameterError("O arquivo de telemetria deve terminar em .csv ou .json.")
    for graph in graphs or ():
        if Path(graph).suffix.lower() not in SUPPORTED_FORMATS:
            raise InvalidParameterError(
                f"Formato de gráfico não suportado em {graph!r}. Use: {', '.join(SUPPORTED_FORMATS)}."
            )


def _settings(args: argparse.Namespace) -> IntegratorSettings:
    return IntegratorSettings(method=args.metodo, rtol=args.rtol)


# --------------------------------------------------------------------------------------
# Descrição de órbitas
# --------------------------------------------------------------------------------------


def _orbit_text(orbit: TwoBodyOrbit) -> str:
    parts = [orbit.orbit_type.value]
    if math.isfinite(orbit.semi_major_axis):
        parts.append(f"a = {_fmt(orbit.semi_major_axis * _KM, 8)} km")
    parts.append(f"e = {_fmt(orbit.eccentricity, 6)}")
    if math.isfinite(orbit.inclination):
        parts.append(f"i = {_fmt(math.degrees(orbit.inclination), 6)}°")
    if math.isfinite(orbit.period):
        parts.append(f"período = {_fmt_quantity(orbit.period, 's')}")
    return "; ".join(parts)


def _orbit_dict(orbit: TwoBodyOrbit) -> dict[str, Any]:
    return {
        "tipo": orbit.orbit_type.name.lower(),
        "descricao": orbit.orbit_type.value,
        "semieixo_maior_m": orbit.semi_major_axis,
        "excentricidade": orbit.eccentricity,
        "inclinacao_graus": math.degrees(orbit.inclination) if math.isfinite(orbit.inclination) else math.nan,
        "periodo_s": orbit.period,
        "energia_especifica_J_kg": orbit.specific_energy,
        "raio_periapse_m": orbit.periapsis_radius,
        "raio_apoapse_m": orbit.apoapsis_radius,
    }


def _diagnostics(trajectory: Trajectory, mu: float) -> dict[str, Any]:
    dr, dv = kepler_error_series(trajectory, mu)
    energy = specific_energy_series(trajectory, mu)
    h = angular_momentum_series(trajectory)
    rel = lambda s: float(np.max(np.abs(s - s[0])) / abs(s[0])) if s[0] != 0 else float(np.max(np.abs(s - s[0])))  # noqa: E731
    return {
        "dr": dr,
        "energy_rel": np.abs(energy - energy[0]) / abs(energy[0]) if energy[0] != 0 else np.abs(energy - energy[0]),
        "max_dr": float(dr.max()),
        "final_dr": float(dr[-1]),
        "max_dv": float(dv.max()),
        "energy": rel(energy),
        "h": rel(h),
    }


def _write_outputs(
    args: argparse.Namespace, body: CelestialBody, trajectory: Trajectory, figure_kwargs: dict[str, Any],
    telemetry_extra: dict[str, Any] | None = None,
) -> list[tuple[str, str]]:
    written: list[tuple[str, str]] = []
    if args.telemetria:
        written.append(("telemetria", str(write_telemetry(args.telemetria, trajectory, body, telemetry_extra))))
    for graph in args.grafico or ():
        written.append(("gráfico", str(save_figure(graph, figure_for_trajectory(body=body, **figure_kwargs)))))
    return written


# --------------------------------------------------------------------------------------
# propagar
# --------------------------------------------------------------------------------------


def _initial_state(args: argparse.Namespace, body: CelestialBody) -> tuple[np.ndarray, np.ndarray, str]:
    inclination = units.parse_quantity(args.inclinacao, "grau") if args.inclinacao is not None else 0.0
    modes = [args.altitude is not None, args.periapse is not None or args.apoapse is not None,
             args.posicao is not None or args.velocidade is not None]
    if sum(modes) != 1:
        raise InvalidParameterError(
            "Informe exatamente uma condição inicial: --altitude; ou --periapse e --apoapse; "
            "ou --posicao e --velocidade."
        )
    incl_text = f"inclinação {_fmt(math.degrees(inclination), 6)}°"
    if args.altitude is not None:
        altitude = units.parse_quantity(args.altitude, "m")
        r, v = circular_orbit_state(body, altitude, inclination)
        return r, v, f"órbita circular de altitude {_fmt_quantity(altitude, 'm')}, {incl_text}"
    if args.periapse is not None or args.apoapse is not None:
        if args.periapse is None or args.apoapse is None:
            raise InvalidParameterError("Informe --periapse e --apoapse juntas (altitudes).")
        hp = units.parse_quantity(args.periapse, "m")
        ha = units.parse_quantity(args.apoapse, "m")
        r, v = elliptic_orbit_state(body, hp, ha, inclination)
        return r, v, (
            f"elipse com periapse a {_fmt(hp * _KM, 8)} km e apoapse a {_fmt(ha * _KM, 8)} km de altitude, "
            f"{incl_text}, partindo do periapse"
        )
    if args.posicao is None or args.velocidade is None:
        raise InvalidParameterError("Informe --posicao e --velocidade juntas.")
    if args.inclinacao is not None:
        raise InvalidParameterError("--inclinacao não se aplica a --posicao/--velocidade (o vetor já define o plano).")
    r = _parse_vector(args.posicao, "m", "--posicao")
    v = _parse_vector(args.velocidade, "m/s", "--velocidade")
    return r, v, "vetor de estado informado (referencial inercial centrado no corpo)"


def _cmd_propagate(args: argparse.Namespace) -> int:
    body = get_body(args.corpo)
    _check_output_paths(args.telemetria, args.grafico)
    r0, v0, description = _initial_state(args, body)
    if math.hypot(*r0) <= body.radius:
        raise DomainError(f"A posição inicial está dentro (ou na superfície) de {body.name}.")
    initial_orbit = orbit_from_state(body.mu, r0, v0)
    if args.orbitas is not None:
        if args.orbitas <= 0:
            raise InvalidParameterError("--orbitas deve ser positivo.")
        if not math.isfinite(initial_orbit.period):
            raise DomainError("--orbitas exige uma órbita fechada (elíptica); para esta trajetória use --duracao.")
        duration = args.orbitas * initial_orbit.period
    else:
        duration = units.parse_quantity(args.duracao, "s")
    settings = _settings(args)
    trajectory = propagate(
        TwoBodyGravity(body.mu), r0, v0, duration, settings=settings, n_points=args.pontos,
        collision_radius=body.radius,
    )
    diag = _diagnostics(trajectory, body.mu)
    final_r = math.hypot(*trajectory.final_position)
    final_v = math.hypot(*trajectory.final_velocity)
    summary = [
        ("Corpo central", body.name),
        ("Condição inicial", description),
        ("Órbita inicial", _orbit_text(initial_orbit)),
        ("Integrador", settings.describe()),
        ("Duração simulada", _fmt_quantity(trajectory.duration, "s")),
        ("Encerramento", trajectory.termination.value),
        ("Maior erro de posição contra Kepler", f"{diag['max_dr']:.3g} m"),
        ("Maior variação relativa da energia", f"{diag['energy']:.3g}"),
    ]
    written = _write_outputs(
        args, body, trajectory,
        {
            "title": f"Simulação orbital — {body.name}",
            "segments": [("trajetória", trajectory)],
            "position_error": diag["dr"],
            "energy_error": diag["energy_rel"],
            "summary": summary,
        },
        {"condicao_inicial": description, "orbita_inicial": _orbit_dict(initial_orbit)},
    )

    if args.json:
        _print_json(
            {
                "corpo": body.key,
                "condicao_inicial": {"descricao": description, "posicao_m": r0.tolist(), "velocidade_m_s": v0.tolist()},
                "orbita_inicial": _orbit_dict(initial_orbit),
                "modelo": trajectory.model,
                "integrador": settings.to_dict(),
                "avaliacoes_da_funcao": trajectory.evaluations,
                "amostras": int(trajectory.times.size),
                "duracao_s": trajectory.duration,
                "encerramento": trajectory.termination.value,
                "estado_final": {
                    "posicao_m": trajectory.final_position.tolist(),
                    "velocidade_m_s": trajectory.final_velocity.tolist(),
                    "altitude_m": final_r - body.radius,
                },
                "verificacao_kepler": {
                    "maior_erro_posicao_m": diag["max_dr"],
                    "erro_posicao_final_m": diag["final_dr"],
                    "maior_erro_velocidade_m_s": diag["max_dv"],
                    "maior_variacao_relativa_energia": diag["energy"],
                    "maior_variacao_relativa_momento_angular": diag["h"],
                },
                "arquivos": [{"tipo": k, "caminho": v} for k, v in written],
                "hipoteses": list(trajectory.hypotheses),
                "origem": trajectory.origin.value,
                "aviso": trajectory.origin.disclaimer,
            }
        )
        return 0

    print(f"Simulação numérica de órbita — {body.name}")
    print(f"  Condição inicial (idealizada): {description}")
    print(f"  Órbita inicial: {_orbit_text(initial_orbit)}")
    print(f"  Modelo: {trajectory.model}")
    print(f"  Integrador: {settings.describe()}; {trajectory.evaluations} avaliações da função")
    print(f"  Duração simulada: {_fmt_quantity(trajectory.duration, 's')}; {trajectory.times.size} amostras")
    print(f"  Encerramento: {trajectory.termination.value}")
    print(f"  Estado final: altitude {_fmt_quantity(final_r - body.radius, 'm')}; "
          f"velocidade {_fmt_quantity(final_v, 'm/s')}")
    print()
    print("Verificação do erro numérico (contra a solução analítica de Kepler, mesmo modelo):")
    print(f"  maior erro de posição: {diag['max_dr']:.3g} m (no fim: {diag['final_dr']:.3g} m); "
          f"maior erro de velocidade: {diag['max_dv']:.3g} m/s")
    print(f"  maior variação relativa da energia: {diag['energy']:.3g}; do momento angular: {diag['h']:.3g}")
    if written:
        print()
        _print_list("Arquivos gerados:", [f"{kind}: {path}" for kind, path in written])
    print()
    _print_list("Hipóteses:", trajectory.hypotheses)
    _print_footer(DataOrigin.SIMULATED)
    return 0


# --------------------------------------------------------------------------------------
# hohmann
# --------------------------------------------------------------------------------------


def _print_transfer(body: CelestialBody, h1: float, h2: float) -> None:
    t = hohmann_transfer(body.mu, body.radius + h1, body.radius + h2)
    print(f"Transferência de Hohmann — {body.name} (cálculo analítico)")
    print(f"  Órbita inicial: circular, altitude {_fmt_quantity(h1, 'm')} (r1 = {_fmt(t.r1 * _KM, 10)} km)")
    print(f"  Órbita final:   circular, altitude {_fmt_quantity(h2, 'm')} (r2 = {_fmt(t.r2 * _KM, 10)} km)")
    print(f"  Elipse de transferência: a = {_fmt(t.transfer_semi_major_axis * _KM, 10)} km")
    print(f"  1ª queima: Δv1 = {_fmt_quantity(t.delta_v1, 'm/s')}  ({_fmt(t.speed_initial_circular, 8)} → "
          f"{_fmt(t.speed_transfer_departure, 8)} m/s)")
    print(f"  2ª queima: Δv2 = {_fmt_quantity(t.delta_v2, 'm/s')}  ({_fmt(t.speed_transfer_arrival, 8)} → "
          f"{_fmt(t.speed_final_circular, 8)} m/s)")
    print(f"  Δv total = {_fmt_quantity(t.delta_v_total, 'm/s')}")
    print(f"  Tempo de transferência = {_fmt_quantity(t.transfer_time, 's')}")
    if t.delta_v1 < 0:
        print("  (Δv negativos indicam queimas retrógradas: redução de órbita.)")


def _transfer_dict(t: HohmannTransfer) -> dict[str, Any]:
    return {
        "r1_m": t.r1, "r2_m": t.r2, "semieixo_transferencia_m": t.transfer_semi_major_axis,
        "delta_v1_m_s": t.delta_v1, "delta_v2_m_s": t.delta_v2, "delta_v_total_m_s": t.delta_v_total,
        "tempo_transferencia_s": t.transfer_time,
        "velocidades_m_s": {
            "circular_inicial": t.speed_initial_circular, "partida_transferencia": t.speed_transfer_departure,
            "chegada_transferencia": t.speed_transfer_arrival, "circular_final": t.speed_final_circular,
        },
    }


def _cmd_hohmann(args: argparse.Namespace) -> int:
    body = get_body(args.corpo)
    h1 = units.parse_quantity(args.de, "m")
    h2 = units.parse_quantity(args.para, "m")
    if h1 < 0 or h2 < 0:
        raise DomainError("As altitudes devem ser não negativas.")
    if not args.simular and (args.telemetria or args.grafico):
        raise InvalidParameterError("--telemetria e --grafico exigem --simular.")
    _check_output_paths(args.telemetria, args.grafico)
    transfer = hohmann_transfer(body.mu, body.radius + h1, body.radius + h2)

    if not args.simular:
        if args.json:
            _print_json({"corpo": body.key, "transferencia": _transfer_dict(transfer),
                         "origem": DataOrigin.ANALYTIC.value, "aviso": DataOrigin.ANALYTIC.disclaimer})
            return 0
        _print_transfer(body, h1, h2)
        print()
        _print_list("Hipóteses:", (
            "Órbitas inicial e final circulares e coplanares; queimas impulsivas e tangenciais.",
            "Problema de dois corpos; perturbações desprezadas; mecânica newtoniana.",
        ))
        _print_footer(DataOrigin.ANALYTIC)
        return 0

    dry = units.parse_quantity(args.massa_seca, "kg")
    isp = units.parse_quantity(args.isp, "s")
    if args.propelente is not None:
        propellant = units.parse_quantity(args.propelente, "kg")
    else:  # exatamente o necessário (Tsiolkovsky), para a simulação não falhar por falta de propelente
        from ..data.constants import G0

        propellant = dry * math.expm1(transfer.delta_v_total / (isp * G0)) * (1.0 + 1e-12)
    craft = Spacecraft(name="Veículo virtual", dry_mass=dry, propellant_mass=propellant, isp=isp)
    mission: HohmannMission = simulate_hohmann(
        body, transfer.r1, transfer.r2, spacecraft=craft, settings=_settings(args), n_points=args.pontos
    )
    events_meta = [
        {"tempo_s": e.time, "evento": e.name, "delta_v_m_s": e.delta_v, "propelente_kg": e.propellant_used,
         "massa_depois_kg": e.mass_after}
        for e in mission.events
    ]
    markers = [Marker(e.name.split(" (")[0], e.time, tuple(p)) for e, p in zip(
        mission.events, (mission.segments[0].positions[0], mission.segments[0].positions[-1]), strict=True)]
    written = _write_outputs(
        args, body, mission.trajectory,
        {
            "title": f"Transferência de Hohmann — {body.name}",
            "segments": [("transferência", mission.segments[0])]
            + ([("órbita final", mission.segments[1])] if len(mission.segments) > 1 else []),
            "markers": markers,
            "reference_orbits": [("órbita inicial (referência)", transfer.r1)],
            "summary": [
                ("Δv total (analítico)", f"{transfer.delta_v_total:.6g} m/s"),
                ("Tempo de transferência", _fmt_quantity(transfer.transfer_time, "s")),
                ("Erro do raio na chegada (numérico)", f"{mission.arrival_radius_error:.3g} m"),
                ("Propelente consumido", f"{mission.propellant_used:.6g} kg"),
            ],
        },
        {"eventos": events_meta, "transferencia_analitica": _transfer_dict(transfer)},
    )

    if args.json:
        _print_json(
            {
                "corpo": body.key,
                "transferencia": _transfer_dict(transfer),
                "simulacao": {
                    "eventos": events_meta,
                    "erro_raio_chegada_m": mission.arrival_radius_error,
                    "orbita_final": _orbit_dict(mission.final_orbit),
                    "erro_semieixo_final_m": mission.final_semi_major_axis_error,
                    "veiculo": {"massa_seca_kg": dry, "propelente_inicial_kg": craft.propellant_mass,
                                "propelente_consumido_kg": mission.propellant_used, "isp_s": isp},
                    "integrador": mission.trajectory.settings.to_dict(),
                },
                "arquivos": [{"tipo": k, "caminho": v} for k, v in written],
                "hipoteses": list(mission.trajectory.hypotheses),
                "origem": DataOrigin.SIMULATED.value,
                "aviso": DataOrigin.SIMULATED.disclaimer,
            }
        )
        return 0

    _print_transfer(body, h1, h2)
    print()
    print("Missão simulada (integração numérica, queimas impulsivas):")
    print(f"  Veículo: massa seca {_fmt(dry, 8)} kg; propelente {_fmt(craft.propellant_mass, 8)} kg; "
          f"Isp {_fmt(isp, 6)} s")
    for e in mission.events:
        print(f"  t = {_fmt_quantity(e.time, 's')}: {e.name}; Δv = {_fmt(e.delta_v, 8)} m/s; "
              f"propelente gasto {_fmt(e.propellant_used, 8)} kg; massa {_fmt(e.mass_after, 8)} kg")
    print(f"  Propelente consumido: {_fmt(mission.propellant_used, 8)} kg "
          f"(restam {_fmt(mission.spacecraft_final.propellant_mass, 6)} kg)")
    print(f"  Órbita final simulada: {_orbit_text(mission.final_orbit)}")
    print(f"  Erro numérico do raio na chegada: {mission.arrival_radius_error:.3g} m; "
          f"do semieixo final: {mission.final_semi_major_axis_error:.3g} m")
    if written:
        print()
        _print_list("Arquivos gerados:", [f"{kind}: {path}" for kind, path in written])
    print()
    _print_list("Hipóteses:", mission.trajectory.hypotheses)
    _print_footer(DataOrigin.SIMULATED)
    return 0


# --------------------------------------------------------------------------------------
# validar
# --------------------------------------------------------------------------------------


def _cmd_validate(args: argparse.Namespace) -> int:
    settings = _settings(args)
    keys = tuple(args.cenario) if args.cenario else None
    results = run_validation(settings, keys)
    if args.json:
        _print_json(
            {
                "integrador": settings.to_dict(),
                "cenarios": [
                    {
                        "chave": r.case.key,
                        "descricao": r.case.description,
                        "duracao_s": r.case.duration,
                        "maior_erro_posicao_m": r.comparison.max_position_error,
                        "maior_erro_velocidade_m_s": r.comparison.max_velocity_error,
                        "variacao_relativa_energia": r.conservation.energy,
                        "variacao_relativa_momento_angular": r.conservation.angular_momentum,
                        "avaliacoes_da_funcao": r.evaluations,
                    }
                    for r in results
                ],
                "origem": DataOrigin.SIMULATED.value,
                "aviso": DataOrigin.SIMULATED.disclaimer,
            }
        )
        return 0
    print("Validação do propagador numérico contra a solução analítica de Kepler (modelo de dois corpos)")
    print(f"  Integrador: {settings.describe()}")
    print()
    header = (
        f"  {'cenário':<12} {'duração':>10} {'máx |Δr|':>11} {'máx |Δv|':>13}"
        f" {'|Δε/ε|':>9} {'|Δh/h|':>9} {'avaliações':>11}"
    )
    print(header)
    for r in results:
        hours = r.case.duration / 3600.0
        print(
            f"  {r.case.key:<12} {hours:>8.2f} h {r.comparison.max_position_error:>9.2e} m "
            f"{r.comparison.max_velocity_error:>9.2e} m/s {r.conservation.energy:>9.1e} "
            f"{r.conservation.angular_momentum:>9.1e} {r.evaluations:>11d}"
        )
    print()
    _print_list("Cenários:", [f"{r.case.key}: {r.case.description}" for r in results])
    print()
    print("O erro medido é apenas o da integração numérica: a referência usa o mesmo modelo de dois corpos.")
    _print_footer(DataOrigin.SIMULATED)
    return 0


# --------------------------------------------------------------------------------------
# Registro dos subcomandos
# --------------------------------------------------------------------------------------


def _add_integrator_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("--rtol", type=_real, default=1e-10, help="tolerância relativa do integrador (padrão: 1e-10)")
    p.add_argument("--metodo", choices=METHODS, default="DOP853", help="método de integração (padrão: DOP853)")


def _add_output_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("--pontos", type=_positive_int, default=500, help="amostras por trecho da trajetória (padrão: 500)")
    p.add_argument("--telemetria", metavar="ARQUIVO", help="grava a telemetria simulada (.csv ou .json)")
    p.add_argument("--grafico", metavar="ARQUIVO", action="append",
                   help="grava gráfico .png, .svg, .pdf ou .html interativo (pode repetir)")
    p.add_argument("--json", action="store_true", help="saída em JSON")


def add_simulation_commands(sub: Any) -> None:
    """Acrescenta os comandos ``propagar``, ``hohmann`` e ``validar`` ao analisador."""
    p = sub.add_parser(
        "propagar",
        help="simula numericamente uma órbita (modelo de dois corpos)",
        description=(
            "Integra numericamente uma órbita no modelo de dois corpos (DOP853, passo adaptativo) e mede\n"
            "o erro numérico contra a solução analítica de Kepler. Exemplos:\n"
            '  aurora propagar --altitude "400 km" --inclinacao 51.6 --orbitas 3 --grafico orbita.png\n'
            '  aurora propagar --periapse "200 km" --apoapse "35786 km" --orbitas 2 --telemetria gto.csv\n'
            '  aurora propagar --posicao "7000 0 0 km" --velocidade "0 8.5 0 km/s" --duracao "5 h"'
        ),
    )
    p.add_argument("--corpo", default="terra", help="terra, lua ou sol (padrão: terra)")
    initial = p.add_argument_group("condição inicial (escolha uma)")
    initial.add_argument("--altitude", metavar="GRANDEZA", help='órbita circular nessa altitude, ex.: "400 km"')
    initial.add_argument("--periapse", metavar="GRANDEZA", help="altitude do periapse (use com --apoapse)")
    initial.add_argument("--apoapse", metavar="GRANDEZA", help="altitude do apoapse (use com --periapse)")
    initial.add_argument("--posicao", metavar='"X Y Z [unid]"', help='posição inicial, ex.: "7000 0 0 km"')
    initial.add_argument("--velocidade", metavar='"VX VY VZ [unid]"', help='velocidade inicial, ex.: "0 7.5 0 km/s"')
    p.add_argument("--inclinacao", metavar="ÂNGULO",
                   help="inclinação do plano orbital (graus por padrão; aceita '1.2 rad')")
    duration = p.add_mutually_exclusive_group(required=True)
    duration.add_argument("--duracao", metavar="GRANDEZA", help='tempo simulado, ex.: "6 h" ou "2 d"')
    duration.add_argument("--orbitas", type=_real, metavar="N", help="número de períodos orbitais (órbitas fechadas)")
    _add_integrator_options(p)
    _add_output_options(p)
    p.set_defaults(handler=_cmd_propagate)

    p = sub.add_parser(
        "hohmann",
        help="transferência de Hohmann entre órbitas circulares (analítica e simulada)",
        description=(
            "Calcula a transferência de Hohmann entre órbitas circulares coplanares e, com --simular,\n"
            "executa-a numericamente com consumo de propelente (Tsiolkovsky). Exemplo:\n"
            '  aurora hohmann --de "400 km" --para "35786 km" --simular --grafico hohmann.html'
        ),
    )
    p.add_argument("--corpo", default="terra", help="terra, lua ou sol (padrão: terra)")
    p.add_argument("--de", required=True, metavar="ALTITUDE", help='altitude da órbita inicial, ex.: "400 km"')
    p.add_argument("--para", required=True, metavar="ALTITUDE", help='altitude da órbita final, ex.: "35786 km"')
    p.add_argument("--simular", action="store_true", help="executa a missão por integração numérica")
    p.add_argument("--massa-seca", default="1000 kg", metavar="GRANDEZA", help="massa sem propelente (padrão: 1000 kg)")
    p.add_argument("--propelente", metavar="GRANDEZA",
                   help="massa de propelente (padrão: exatamente a necessária para a missão)")
    p.add_argument("--isp", default="300 s", metavar="GRANDEZA", help="impulso específico (padrão: 300 s)")
    _add_integrator_options(p)
    _add_output_options(p)
    p.set_defaults(handler=_cmd_hohmann)

    p = sub.add_parser(
        "validar",
        help="valida o propagador numérico contra a solução analítica de Kepler",
        description="Executa os cenários padrão de validação e mostra os erros numéricos medidos.",
    )
    p.add_argument("--cenario", action="append", choices=[c.key for c in standard_cases()],
                   help="executa só este cenário (pode repetir)")
    _add_integrator_options(p)
    p.add_argument("--json", action="store_true", help="saída em JSON")
    p.set_defaults(handler=_cmd_validate)

