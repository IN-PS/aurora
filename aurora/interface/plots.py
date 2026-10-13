"""Gráficos de trajetórias simuladas.

* Arquivos ``.png``, ``.svg`` ou ``.pdf``: figura estática (Matplotlib, API orientada a
  objetos, sem depender de interface gráfica) com projeção da órbita no plano xy, vista 3D,
  altitude e velocidade em função do tempo e, quando disponível, o erro numérico medido
  contra a solução de Kepler e a variação relativa da energia.
* Arquivos ``.html``: página interativa autocontida (Plotly embutido, funciona sem
  internet, inclusive no navegador de um celular), com órbita 3D, altitude e velocidade,
  tooltips, tabela-resumo e modo escuro.

Ambas as bibliotecas são opcionais (extra ``[graficos]``). Cores: paleta categórica de
referência validada para daltonismo (slots 1 e 2) e tinta neutra; uma série por eixo
(sem eixo y duplo).
"""

from __future__ import annotations

import html
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..core.provenance import DataOrigin
from ..data.bodies import CelestialBody
from ..errors import DependencyError, InvalidParameterError
from ..simulation.propagator import Trajectory

# Paleta de referência (categórica, slots 1–2) e tinta/superfície, modo claro e escuro.
LIGHT = {
    "surface": "#fcfcfb",
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "series": ("#2a78d6", "#eb6834"),
    "body": "#e1e0d9",
}
DARK = {
    "surface": "#1a1a19",
    "ink": "#f0efec",
    "ink2": "#c3c2b7",
    "muted": "#898781",
    "grid": "#2c2c2a",
    "axis": "#383835",
    "series": ("#3987e5", "#d95926"),
    "body": "#383835",
}
STATIC_FORMATS = (".png", ".svg", ".pdf")
SUPPORTED_FORMATS = (*STATIC_FORMATS, ".html")


@dataclass(frozen=True)
class Segment:
    """Parte da trajetória com rótulo próprio (ex.: arco de transferência e órbita final)."""

    label: str
    trajectory: Trajectory


@dataclass(frozen=True)
class Marker:
    """Evento marcado no gráfico (posição em m, tempo em s)."""

    label: str
    time: float
    position: tuple[float, float, float]


@dataclass(frozen=True)
class FigureData:
    """Tudo o que os gráficos precisam, em SI."""

    title: str
    body: CelestialBody
    segments: tuple[Segment, ...]
    markers: tuple[Marker, ...] = ()
    reference_orbits: tuple[tuple[str, float], ...] = ()
    """Círculos de referência no plano xy: (rótulo, raio em m)."""
    position_error: np.ndarray | None = None
    """|Δr| contra Kepler em cada amostra (m), na ordem das amostras concatenadas."""
    energy_error: np.ndarray | None = None
    """|Δε/ε₀| em cada amostra."""
    summary: tuple[tuple[str, str], ...] = field(default=())
    """Linhas (rótulo, valor) da tabela-resumo do HTML."""

    def __post_init__(self) -> None:
        if not self.segments:
            raise InvalidParameterError("É necessário ao menos um segmento de trajetória para o gráfico.")
        if len(self.segments) > 2:
            raise InvalidParameterError("O gráfico admite no máximo dois segmentos com cores distintas.")


def _concat_times(data: FigureData) -> np.ndarray:
    return np.concatenate([s.trajectory.times for s in data.segments])


def save_figure(path: str | Path, data: FigureData) -> Path:
    """Grava o gráfico no formato indicado pela extensão (.png, .svg, .pdf ou .html)."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_FORMATS:
        raise InvalidParameterError(
            f"Formato de gráfico não suportado: {suffix or '(sem extensão)'}. Use: {', '.join(SUPPORTED_FORMATS)}."
        )
    if suffix == ".html":
        _save_html(path, data)
    else:
        _save_static(path, data)
    return path


# --------------------------------------------------------------------------------------
# Figura estática (Matplotlib)
# --------------------------------------------------------------------------------------


def _import_matplotlib():  # type: ignore[no-untyped-def]
    try:
        from matplotlib.figure import Figure
    except ImportError:
        raise DependencyError(
            "Gráficos estáticos exigem o Matplotlib: pip install \"aurora-simulador[graficos]\" "
            "(no Termux: pkg install matplotlib)."
        ) from None
    return Figure


def _style_axes(ax, c: dict) -> None:  # type: ignore[no-untyped-def]
    ax.set_facecolor(c["surface"])
    ax.grid(True, color=c["grid"], linewidth=0.75, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(c["axis"])
        ax.spines[side].set_linewidth(0.75)
    ax.tick_params(colors=c["muted"], labelcolor=c["ink2"], labelsize=8)
    ax.xaxis.label.set_color(c["ink2"])
    ax.yaxis.label.set_color(c["ink2"])
    ax.title.set_color(c["ink"])


def _save_static(path: Path, data: FigureData) -> None:
    Figure = _import_matplotlib()  # noqa: N806
    c = LIGHT
    has_errors = data.position_error is not None or data.energy_error is not None
    rows = 3 if has_errors else 2
    fig = Figure(figsize=(11, 4.2 * rows), dpi=150, facecolor=c["surface"], layout="constrained")
    km, hours = 1e-3, 1.0 / 3600.0
    t_all = _concat_times(data) * hours
    multi = len(data.segments) > 1

    # (a) Projeção no plano da órbita inicial (forma verdadeira da cônica)
    ax = fig.add_subplot(rows, 2, 1)
    _style_axes(ax, c)
    basis, plane_title, labels = _projection_basis(data.segments[0].trajectory)
    ax.add_patch(_circle(data.body.radius * km, c["body"]))
    ax.annotate(data.body.name, (0, 0), ha="center", va="center", color=c["ink2"], fontsize=8)
    for label, radius in data.reference_orbits:
        theta = np.linspace(0, 2 * math.pi, 361)
        ax.plot(radius * km * np.cos(theta), radius * km * np.sin(theta), color=c["axis"], linewidth=0.75, label=label)
    for k, seg in enumerate(data.segments):
        p = seg.trajectory.positions @ basis.T * km
        ax.plot(p[:, 0], p[:, 1], color=c["series"][k], linewidth=1.5, solid_capstyle="round", label=seg.label)
    _plot_markers_2d(ax, data, c, km, basis)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(plane_title, fontsize=10, loc="left")
    if multi or data.reference_orbits:
        _legend(ax, c)

    # (b) Vista 3D
    ax3 = fig.add_subplot(rows, 2, 2, projection="3d")
    ax3.set_facecolor(c["surface"])
    u, v = np.mgrid[0 : 2 * math.pi : 40j, 0 : math.pi : 20j]
    rb = data.body.radius * km
    ax3.plot_surface(rb * np.cos(u) * np.sin(v), rb * np.sin(u) * np.sin(v), rb * np.cos(v), color=c["body"],
                     linewidth=0, alpha=0.6, shade=False)
    extent = rb
    for k, seg in enumerate(data.segments):
        p = seg.trajectory.positions * km
        ax3.plot(p[:, 0], p[:, 1], p[:, 2], color=c["series"][k], linewidth=1.5, label=seg.label)
        extent = max(extent, float(np.max(np.abs(p))))
    for axis in (ax3.set_xlim, ax3.set_ylim, ax3.set_zlim):
        axis(-extent, extent)
    ax3.set_box_aspect((1, 1, 1))
    for pane_axis in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
        pane_axis.set_pane_color(c["surface"])
        pane_axis._axinfo["grid"]["color"] = c["grid"]  # noqa: SLF001 (estilo da grade 3D)
        pane_axis._axinfo["grid"]["linewidth"] = 0.6  # noqa: SLF001
    for setter, name in ((ax3.set_xlabel, "x (km)"), (ax3.set_ylabel, "y (km)"), (ax3.set_zlabel, "z (km)")):
        setter(name, color=c["ink2"], fontsize=8)
    ax3.tick_params(colors=c["muted"], labelsize=7)
    ax3.set_title("Vista 3D", fontsize=10, loc="left", color=c["ink"])

    # (c) Altitude e (d) velocidade — um gráfico por grandeza (sem eixo duplo)
    for col, (title, ylabel, values_of) in enumerate(
        (
            ("Altitude", "altitude (km)", lambda tr: (np.linalg.norm(tr.positions, axis=1) - data.body.radius) * km),
            ("Velocidade", "velocidade (km/s)", lambda tr: np.linalg.norm(tr.velocities, axis=1) * km),
        )
    ):
        ax = fig.add_subplot(rows, 2, 3 + col)
        _style_axes(ax, c)
        for k, seg in enumerate(data.segments):
            ax.plot(seg.trajectory.times * hours, values_of(seg.trajectory), color=c["series"][k], linewidth=1.5,
                    label=seg.label)
        for marker in data.markers:
            ax.axvline(marker.time * hours, color=c["axis"], linewidth=0.75)
        ax.set_xlabel("tempo (h)")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10, loc="left")
        if multi:
            _legend(ax, c)

    # (e) e (f) erros numéricos
    if has_errors:
        for col, (title, ylabel, series) in enumerate(
            (
                ("Erro de posição contra Kepler", "|Δr| (m)", data.position_error),
                ("Variação relativa da energia", "|Δε/ε₀|", data.energy_error),
            )
        ):
            ax = fig.add_subplot(rows, 2, 5 + col)
            _style_axes(ax, c)
            if series is None:
                ax.text(0.5, 0.5, "não disponível", transform=ax.transAxes, ha="center", color=c["muted"])
            else:
                values = np.asarray(series, dtype=float)
                positive = values > 0
                ax.plot(t_all[positive], values[positive], color=c["series"][0], linewidth=1.5)
                if np.any(positive):
                    ax.set_yscale("log")
            ax.set_xlabel("tempo (h)")
            ax.set_ylabel(ylabel)
            ax.set_title(title, fontsize=10, loc="left")

    origin = data.segments[0].trajectory.origin
    fig.suptitle(data.title, color=c["ink"], fontsize=13, x=0.01, ha="left")
    fig.supxlabel(f"Natureza do resultado: {origin.disclaimer}", color=c["ink2"], fontsize=8, x=0.01, ha="left")
    fig.savefig(path, facecolor=c["surface"])


def _circle(radius: float, color: str):  # type: ignore[no-untyped-def]
    from matplotlib.patches import Circle

    return Circle((0, 0), radius, facecolor=color, edgecolor="none", zorder=0)


def _projection_basis(trajectory: Trajectory) -> tuple[np.ndarray, str, tuple[str, str]]:
    """Base ortonormal do plano da órbita inicial (p̂ ao longo de r₀, q̂ = ŵ × p̂); xy se h ≈ 0."""
    r0, v0 = trajectory.positions[0], trajectory.velocities[0]
    h = np.cross(r0, v0)
    h_norm = float(np.linalg.norm(h))
    if h_norm <= 1e-12 * float(np.linalg.norm(r0)) * float(np.linalg.norm(v0)) or h_norm == 0.0:
        return np.eye(3)[:2], "Projeção no plano xy", ("x (km)", "y (km)")
    p_hat = r0 / np.linalg.norm(r0)
    w_hat = h / h_norm
    q_hat = np.cross(w_hat, p_hat)
    return np.vstack((p_hat, q_hat)), "Plano da órbita inicial", ("p (km), na direção da posição inicial", "q (km)")


def _plot_markers_2d(ax, data: FigureData, c: dict, km: float, basis: np.ndarray) -> None:  # type: ignore[no-untyped-def]
    first = basis @ data.segments[0].trajectory.positions[0] * km
    last = basis @ data.segments[-1].trajectory.positions[-1] * km
    extent = max(float(np.max(np.abs(seg.trajectory.positions))) for seg in data.segments) * km
    candidates = [("início", first)] + [(m.label, basis @ np.asarray(m.position) * km) for m in data.markers]
    candidates.append(("fim", last))
    points: list[tuple[str, np.ndarray]] = []
    for label, p in candidates:  # pontos coincidentes (< 1 % da escala) recebem um só rótulo
        for k, (other_label, q) in enumerate(points):
            if float(np.linalg.norm(p - q)) < 0.01 * extent:
                points[k] = (f"{other_label} / {label}", q)
                break
        else:
            points.append((label, p))
    for label, p in points:
        ax.plot(p[0], p[1], marker="o", markersize=6, color=c["ink"], markeredgecolor=c["surface"], markeredgewidth=1.5,
                linestyle="none", zorder=5)
        ax.annotate(label, (p[0], p[1]), textcoords="offset points", xytext=(6, 6), color=c["ink2"], fontsize=8)


def _legend(ax, c: dict) -> None:  # type: ignore[no-untyped-def]
    legend = ax.legend(frameon=False, fontsize=8, loc="best")
    for text in legend.get_texts():
        text.set_color(c["ink2"])


# --------------------------------------------------------------------------------------
# Página HTML interativa (Plotly)
# --------------------------------------------------------------------------------------


def _plotly_js() -> str:
    try:
        from plotly.offline import get_plotlyjs
    except ImportError:
        raise DependencyError(
            "Gráficos interativos exigem o Plotly: pip install \"aurora-simulador[graficos]\" (ou pip install plotly)."
        ) from None
    return get_plotlyjs()


def _decimate(n: int, limit: int = 4000) -> slice:
    return slice(None, None, max(1, math.ceil(n / limit)))


def _html_payload(data: FigureData) -> dict:
    km, hours = 1e-3, 1.0 / 3600.0
    segments = []
    for seg in data.segments:
        tr = seg.trajectory
        step = _decimate(tr.times.size)
        r = np.linalg.norm(tr.positions, axis=1)
        segments.append(
            {
                "label": seg.label,
                "t_h": (tr.times[step] * hours).round(6).tolist(),
                "x": (tr.positions[step, 0] * km).round(3).tolist(),
                "y": (tr.positions[step, 1] * km).round(3).tolist(),
                "z": (tr.positions[step, 2] * km).round(3).tolist(),
                "alt_km": ((r[step] - data.body.radius) * km).round(3).tolist(),
                "v_kms": (np.linalg.norm(tr.velocities, axis=1)[step] * km).round(6).tolist(),
            }
        )
    return {
        "body": {"name": data.body.name, "radius_km": data.body.radius * km},
        "segments": segments,
        "markers": [
            {"label": m.label, "t_h": m.time * hours, "x": m.position[0] * km, "y": m.position[1] * km,
             "z": m.position[2] * km}
            for m in data.markers
        ],
        "palette": {"light": LIGHT, "dark": DARK},
    }


_HTML_TEMPLATE = """<!doctype html>
<html lang="pt-BR">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --border: rgba(11,11,11,0.10);
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --surface: #1a1a19; --page: #0d0d0d; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --border: rgba(255,255,255,0.10);
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --page: #0d0d0d; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --border: rgba(255,255,255,0.10);
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1100px; margin: 0 auto; padding: 16px; }
h1 { font-size: 1.25rem; margin: 0 0 4px; }
.aviso { color: var(--ink2); font-size: 0.85rem; margin: 0 0 16px; }
.grade { display: grid; gap: 16px; grid-template-columns: 1fr 1fr; }
.cartao { background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 8px; min-width: 0; }
.cartao.largo { grid-column: 1 / -1; }
.cartao h2 { font-size: 0.95rem; margin: 4px 8px; color: var(--ink); }
.grafico { width: 100%; height: 340px; }
.grafico.alto { height: 520px; }
table { border-collapse: collapse; width: 100%; font-size: 0.85rem; }
th, td { text-align: left; padding: 6px 8px; border-top: 1px solid var(--grid); }
td { font-variant-numeric: tabular-nums; color: var(--ink2); }
@media (max-width: 720px) { .grade { grid-template-columns: 1fr; } .grafico.alto { height: 380px; } }
</style>
</head>
<body>
<main>
<h1>__TITLE__</h1>
<p class="aviso">__DISCLAIMER__</p>
<div class="grade">
  <section class="cartao largo"><h2>Trajetória em 3D (km)</h2><div id="orbita" class="grafico alto"></div></section>
  <section class="cartao"><h2>Altitude (km)</h2><div id="altitude" class="grafico"></div></section>
  <section class="cartao"><h2>Velocidade (km/s)</h2><div id="velocidade" class="grafico"></div></section>
  <section class="cartao largo"><h2>Resumo</h2><table>__SUMMARY__</table></section>
</div>
</main>
<script>__PLOTLY__</script>
<script>
const DATA = __DATA__;
function palette() {
  const forced = document.documentElement.getAttribute("data-theme");
  const dark = forced ? forced === "dark" : window.matchMedia("(prefers-color-scheme: dark)").matches;
  return dark ? DATA.palette.dark : DATA.palette.light;
}
function axis(c, title) {
  return { title: { text: title, font: { color: c.ink2, size: 12 } }, gridcolor: c.grid, linecolor: c.axis,
           zerolinecolor: c.axis, tickfont: { color: c.ink2, size: 11 } };
}
function draw() {
  const c = palette();
  const cfg = { responsive: true, displaylogo: false };
  const base = { paper_bgcolor: c.surface, plot_bgcolor: c.surface, font: { color: c.ink2,
                 family: "system-ui, -apple-system, Segoe UI, sans-serif" }, margin: { l: 56, r: 16, t: 8, b: 48 },
                 showlegend: DATA.segments.length > 1, legend: { orientation: "h", y: 1.12, font: { color: c.ink2 } },
                 hoverlabel: { bgcolor: c.surface, bordercolor: c.axis, font: { color: c.ink } } };
  // Esfera do corpo central
  const n = 32, R = DATA.body.radius_km, sx = [], sy = [], sz = [];
  for (let i = 0; i <= n; i++) {
    const th = Math.PI * i / n, rx = [], ry = [], rz = [];
    for (let j = 0; j <= 2 * n; j++) {
      const ph = Math.PI * j / n;
      rx.push(R * Math.sin(th) * Math.cos(ph)); ry.push(R * Math.sin(th) * Math.sin(ph)); rz.push(R * Math.cos(th));
    }
    sx.push(rx); sy.push(ry); sz.push(rz);
  }
  const traces3d = [{ type: "surface", x: sx, y: sy, z: sz, showscale: false, opacity: 0.6, hoverinfo: "skip",
                      colorscale: [[0, c.body], [1, c.body]], name: DATA.body.name, showlegend: false }];
  DATA.segments.forEach((s, k) => traces3d.push({
    type: "scatter3d", mode: "lines", x: s.x, y: s.y, z: s.z, name: s.label,
    line: { color: c.series[k], width: 4 }, customdata: s.t_h.map((t, i) => [t, s.alt_km[i], s.v_kms[i]]),
    hovertemplate: "t = %{customdata[0]:.3f} h<br>altitude = %{customdata[1]:,.1f} km<br>" +
                   "velocidade = %{customdata[2]:.4f} km/s<extra>" + s.label + "</extra>" }));
  if (DATA.markers.length) traces3d.push({ type: "scatter3d", mode: "markers+text", name: "eventos", showlegend: false,
    x: DATA.markers.map(m => m.x), y: DATA.markers.map(m => m.y), z: DATA.markers.map(m => m.z),
    text: DATA.markers.map(m => m.label), textfont: { color: c.ink2 }, marker: { size: 5, color: c.ink },
    hovertemplate: "%{text}<extra></extra>" });
  Plotly.react("orbita", traces3d, Object.assign({}, base, { scene: {
    aspectmode: "data", bgcolor: c.surface,
    xaxis: axis(c, "x (km)"), yaxis: axis(c, "y (km)"), zaxis: axis(c, "z (km)") },
    margin: { l: 0, r: 0, t: 8, b: 0 } }), cfg);
  const series2d = [["altitude", "alt_km", "altitude (km)", ",.1f"],
                    ["velocidade", "v_kms", "velocidade (km/s)", ".4f"]];
  for (const [id, key, label, fmt] of series2d) {
    const traces = DATA.segments.map((s, k) => ({ type: "scatter", mode: "lines", x: s.t_h, y: s[key], name: s.label,
      line: { color: c.series[k], width: 2 },
      hovertemplate: "t = %{x:.3f} h<br>" + label + " = %{y:" + fmt + "}<extra>" + s.label + "</extra>" }));
    const shapes = DATA.markers.map(m => ({ type: "line", xref: "x", yref: "paper", x0: m.t_h, x1: m.t_h, y0: 0, y1: 1,
                                             line: { color: c.axis, width: 1 } }));
    Plotly.react(id, traces, Object.assign({}, base, { hovermode: "x unified", shapes,
      xaxis: axis(c, "tempo (h)"), yaxis: axis(c, label) }), cfg);
  }
}
draw();
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
new MutationObserver(draw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
</script>
</body>
</html>
"""


def _save_html(path: Path, data: FigureData) -> None:
    plotly_js = _plotly_js()
    origin: DataOrigin = data.segments[0].trajectory.origin
    rows = "".join(
        f"<tr><th scope=\"row\">{html.escape(label)}</th><td>{html.escape(value)}</td></tr>"
        for label, value in data.summary
    )
    page = (
        _HTML_TEMPLATE.replace("__TITLE__", html.escape(data.title))
        .replace("__DISCLAIMER__", html.escape(f"Natureza do resultado: {origin.disclaimer}"))
        .replace("__SUMMARY__", rows or "<tr><td>—</td></tr>")
        .replace("__DATA__", json.dumps(_html_payload(data), ensure_ascii=True).replace("</", "<\\/"))
        .replace("__PLOTLY__", plotly_js.replace("</script", "<\\/script"))
    )
    path.write_text(page, encoding="utf-8")


def figure_for_trajectory(
    title: str,
    body: CelestialBody,
    segments: Sequence[tuple[str, Trajectory]],
    *,
    markers: Sequence[Marker] = (),
    reference_orbits: Sequence[tuple[str, float]] = (),
    position_error: np.ndarray | None = None,
    energy_error: np.ndarray | None = None,
    summary: Sequence[tuple[str, str]] = (),
) -> FigureData:
    """Atalho para montar :class:`FigureData`."""
    return FigureData(
        title=title,
        body=body,
        segments=tuple(Segment(label, tr) for label, tr in segments),
        markers=tuple(markers),
        reference_orbits=tuple(reference_orbits),
        position_error=position_error,
        energy_error=energy_error,
        summary=tuple(summary),
    )
