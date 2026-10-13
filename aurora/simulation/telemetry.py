"""Telemetria simulada: tabela de grandezas derivadas da trajetória e exportação CSV/JSON.

Todo arquivo exportado declara a origem ``simulacao`` e o aviso de que os dados vêm de uma
simulação numérica de modelo idealizado — nunca de medição de um veículo real. O cabeçalho
registra também o modelo, o integrador, as tolerâncias e as hipóteses, para que o
experimento possa ser reproduzido.
"""

from __future__ import annotations

import csv
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from .. import __version__
from ..data.bodies import CelestialBody
from ..errors import InvalidParameterError
from .propagator import Trajectory

COLUMNS: tuple[tuple[str, str, str], ...] = (
    ("t", "s", "Tempo desde o início da simulação"),
    ("x", "m", "Posição x (referencial inercial centrado no corpo central)"),
    ("y", "m", "Posição y"),
    ("z", "m", "Posição z"),
    ("vx", "m/s", "Velocidade x"),
    ("vy", "m/s", "Velocidade y"),
    ("vz", "m/s", "Velocidade z"),
    ("raio", "m", "Distância ao centro do corpo central"),
    ("altitude", "m", "Altitude acima do raio de referência do corpo"),
    ("velocidade", "m/s", "Módulo da velocidade"),
    ("energia_especifica", "J/kg", "Energia orbital específica v²/2 − μ/r"),
    ("momento_angular", "m^2/s", "Módulo do momento angular específico |r × v|"),
)


def telemetry_table(trajectory: Trajectory, body: CelestialBody) -> np.ndarray:
    """Matriz (amostras × colunas de :data:`COLUMNS`) com as grandezas em SI."""
    r = np.linalg.norm(trajectory.positions, axis=1)
    speed = np.linalg.norm(trajectory.velocities, axis=1)
    energy = 0.5 * speed**2 - body.mu / r
    h = np.linalg.norm(np.cross(trajectory.positions, trajectory.velocities), axis=1)
    return np.column_stack(
        (trajectory.times - trajectory.times[0], trajectory.positions, trajectory.velocities, r, r - body.radius, speed,
         energy, h)
    )


def telemetry_metadata(
    trajectory: Trajectory, body: CelestialBody, extra: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """Metadados que acompanham a telemetria exportada (chaves em português)."""
    metadata: dict[str, Any] = {
        "titulo": "AURORA — telemetria simulada",
        "origem": trajectory.origin.value,
        "aviso": trajectory.origin.disclaimer,
        "versao_aurora": __version__,
        "corpo": {
            "nome": body.name,
            "mu_m3_s2": body.mu,
            "raio_referencia_m": body.radius,
            "tipo_raio": body.radius_description,
            "fontes": list(body.sources),
        },
        "modelo": trajectory.model,
        "integrador": trajectory.settings.to_dict(),
        "avaliacoes_da_funcao": trajectory.evaluations,
        "encerramento": trajectory.termination.value,
        "hipoteses": list(trajectory.hypotheses),
    }
    if extra:
        metadata.update(extra)
    return metadata


def _json_number(x: float) -> float | str:
    return x + 0.0 if math.isfinite(x) else ("inf" if x > 0 else "-inf" if x < 0 else "nan")


def write_telemetry(
    path: str | Path, trajectory: Trajectory, body: CelestialBody, extra: Mapping[str, Any] | None = None
) -> Path:
    """Grava a telemetria em ``.csv`` (cabeçalho comentado com ``#``) ou ``.json``; retorna o caminho."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in (".csv", ".json"):
        raise InvalidParameterError("O arquivo de telemetria deve terminar em .csv ou .json.")
    table = telemetry_table(trajectory, body)
    metadata = telemetry_metadata(trajectory, body, extra)
    if suffix == ".json":
        payload = {
            **metadata,
            "colunas": [{"nome": n, "unidade": u, "descricao": d} for n, u, d in COLUMNS],
            "dados": [[_json_number(float(x)) for x in row] for row in table],
        }
        path.write_text(json.dumps(payload, ensure_ascii=True, indent=1, allow_nan=False) + "\n", encoding="utf-8")
        return path
    with path.open("w", newline="", encoding="utf-8") as handle:
        for key, value in metadata.items():
            text = json.dumps(value, ensure_ascii=False) if not isinstance(value, str) else value
            handle.write(f"# {key}: {text}\n")
        writer = csv.writer(handle)
        writer.writerow([f"{n}_{u.replace('/', '_').replace('^', '')}" for n, u, _ in COLUMNS])
        for row in table:
            writer.writerow([repr(float(x)) for x in row])
    return path


def read_telemetry_csv(path: str | Path) -> tuple[list[str], np.ndarray, dict[str, str]]:
    """Lê um CSV gravado por :func:`write_telemetry`: (nomes das colunas, dados, metadados brutos)."""
    metadata: dict[str, str] = {}
    rows: list[Sequence[str]] = []
    header: list[str] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# "):
                key, _, value = line[2:].partition(": ")
                metadata[key] = value.rstrip("\n")
            elif not header:
                header = next(csv.reader([line]))
            elif line.strip():
                rows.append(next(csv.reader([line])))
    return header, np.array(rows, dtype=float), metadata
