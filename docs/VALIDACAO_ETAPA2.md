# Validação da etapa 2 — erros numéricos medidos

Este documento registra **medições** do erro numérico da simulação orbital, feitas com o
próprio código do repositório (versão 0.2.0) em um computador Linux x86-64. Os valores podem
variar ligeiramente com a versão do SciPy e com o processador, mas a ordem de grandeza não.
Cada tabela indica o comando ou o teste que a reproduz.

> Todas as medições comparam a simulação com a **solução analítica do mesmo modelo** (dois
> corpos, massas pontuais, sem perturbações). Elas medem a qualidade da integração numérica;
> **não** medem a distância entre o modelo e um satélite real, que é muito maior: arrasto,
> achatamento da Terra, Lua e Sol afastam uma órbita baixa real da órbita kepleriana em
> quilômetros ou mais ao longo de um único dia.

## 1. Referência analítica: propagação de Kepler

A propagação por variável universal (`aurora.physics.kepler.propagate_kepler`) é comparada
com soluções independentes, escritas nos próprios testes (`tests/test_kepler.py`):

| Verificação | Casos | Maior erro medido | Limite exigido no teste |
|---|---|---|---|
| Elipses aleatórias (a de 6 300 a 316 000 km, 0 ≤ e < 0,97, orientação qualquer, Δt de até ±3 períodos) contra a equação de Kepler M = E − e·sen E | 1 000 (semente fixa) | \|Δr\|/a = 3,5 × 10⁻¹³ | 10⁻¹⁰ |
| Hipérboles aleatórias (1,05 ≤ e ≤ 4, Δt de até ±2 × 10⁵ s) contra M = e·senh H − H | 1 000 (semente fixa) | \|Δr\|/\|r\| = 1,8 × 10⁻¹³ | 10⁻¹¹ |
| Elipses com e de 0 a 0,95 em cinco instantes (até 3,4 períodos) | 25 | — | 10⁻⁹ · a |
| Hipérboles com e de 1,05 a 10 até 3 × 10⁷ s | 16 | — | 10⁻¹⁰ relativo |
| Parábola contra a equação de Barker | 3 instantes (até 10⁶ s) | — | 10⁻¹² relativo no tempo |
| Ida e volta: propagar 12 345,6 s e voltar | 1 | — | 10⁻⁶ m e 10⁻⁹ m/s |
| 10 000 períodos + ¼ de período contra ¼ de período | 1 | — | 10⁻⁶ · a |
| Conservação da energia e do momento angular (41 instantes) | 41 | — | 10⁻¹¹ relativo |

Erros de 2 a 4 × 10⁻¹³ correspondem a cerca de 800 a 1 600 vezes o épsilon de máquina
(2,2 × 10⁻¹⁶), o esperado para dezenas de operações em ponto flutuante sobre números da
ordem de a. Como essa solução é exata para o modelo e muito mais precisa que a integração
numérica (seção 2: erros relativos da ordem de 10⁻¹⁰ a 10⁻⁸), ela serve de referência para medir o
erro do integrador.

## 2. Integrador numérico: cenários padrão

Reprodução: `aurora validar` (padrão: DOP853, `rtol = 1e-10`, `atol = 1e-4 m / 1e-7 m/s`,
400 amostras por cenário). Saída obtida:

| Cenário | Duração | Maior \|Δr\| | Maior \|Δv\| | \|Δε/ε\| | \|Δh/h\| | Avaliações |
|---|---|---|---|---|---|---|
| `leo`: circular 400 km, 10 órbitas | 15,43 h | 7,61 × 10⁻³ m | 8,37 × 10⁻⁶ m/s | 2,8 × 10⁻¹⁰ | 1,4 × 10⁻¹⁰ | 3 983 |
| `geo`: geoestacionária, 2 órbitas | 47,87 h | 1,90 × 10⁻² m | 1,32 × 10⁻⁶ m/s | 2,4 × 10⁻¹⁰ | 1,2 × 10⁻¹⁰ | 878 |
| `gto`: 200 × 35 786 km, i = 27°, 5 órbitas | 52,59 h | 5,07 × 10⁻¹ m | 4,58 × 10⁻⁴ m/s | 1,4 × 10⁻⁹ | 3,0 × 10⁻¹⁰ | 4 901 |
| `molniya`: 500 × 39 800 km, i = 63,4°, 3 órbitas | 35,83 h | 1,40 × 10⁻¹ m | 1,17 × 10⁻⁴ m/s | 2,1 × 10⁻⁹ | 4,9 × 10⁻¹⁰ | 3 056 |
| `lua`: circular 100 km na Lua, 5 órbitas | 9,82 h | 2,64 × 10⁻³ m | 2,27 × 10⁻⁶ m/s | 3,9 × 10⁻¹⁰ | 1,9 × 10⁻¹⁰ | 1 943 |
| `hiperbolica`: 1,2·v_escape a 400 km, 2 dias | 48,00 h | 1,83 × 10⁻² m | 3,20 × 10⁻⁷ m/s | 5,0 × 10⁻¹¹ | 1,1 × 10⁻¹⁰ | 611 |

Δε/ε e Δh/h são as maiores variações relativas da energia específica e do momento angular
específico, que o modelo conserva exatamente. Os seis cenários levam cerca de 0,4 s no total.

### Efeito da tolerância e do método

Reprodução: `aurora validar --rtol 1e-8`, `aurora validar --rtol 1e-12` e
`aurora validar --metodo RK45`. Maior \|Δr\| (avaliações da função entre parênteses):

| Cenário | DOP853, rtol 10⁻⁸ | DOP853, rtol 10⁻¹⁰ (padrão) | DOP853, rtol 10⁻¹² | RK45, rtol 10⁻¹⁰ |
|---|---|---|---|---|
| `leo` | 0,613 m (2 303) | 7,61 mm (3 983) | 1,54 mm (4 973) | 0,936 m (9 878) |
| `geo` | 1,95 m (518) | 19,0 mm (878) | 0,62 mm (1 328) | 0,159 m (2 156) |
| `gto` | 149 m (3 185) | 0,507 m (4 901) | 0,171 m (5 495) | 2,67 m (7 424) |
| `molniya` | 38,3 m (1 913) | 0,140 m (3 056) | 26,7 mm (3 410) | 1,34 m (4 508) |
| `lua` | 0,113 m (1 178) | 2,64 mm (1 943) | 1,07 mm (2 183) | 0,103 m (4 358) |
| `hiperbolica` | 0,601 m (431) | 18,3 mm (611) | 2,10 mm (755) | 0,197 m (788) |

Observações:

- De 10⁻⁸ para 10⁻¹⁰ o erro cai de 33 a 294 vezes; de 10⁻¹⁰ para 10⁻¹² cai bem menos (de
  2,5 a 31 vezes), sinal de que o arredondamento de ponto flutuante começa a pesar. Por isso
  tolerâncias abaixo de 10⁻¹³ são recusadas.
- Com a mesma tolerância, o RK45 (ordem 5) precisou de 1,3 a 2,5 vezes mais avaliações e
  errou de 5 a 123 vezes mais que o DOP853 (ordem 8). O DOP853 é o padrão.
- O teste `tests/test_propagator.py::test_cenarios_padrao_de_validacao` repete estes cenários
  (com 200 amostras) e exige erros abaixo de limites com folga de 10 a 20 vezes sobre os
  valores medidos: `leo` 0,1 m; `geo` 0,2 m; `gto` 5 m; `molniya` 2 m; `lua` 0,05 m;
  `hiperbolica` 0,2 m; energia e momento angular abaixo de 10⁻⁸.

## 3. Crescimento do erro com o tempo simulado

O erro **global** não é limitado pela tolerância: ele acumula. Medição com o script abaixo
(DOP853, tolerâncias padrão, órbitas equatoriais partindo do periapse, 200 amostras):

| Órbita | 1 órbita | 10 órbitas | 100 órbitas | 1 000 órbitas |
|---|---|---|---|---|
| LEO 400 km (período 1,54 h) — maior \|Δr\| | 2,32 mm | 7,57 mm | 1,49 m | 169 m |
| LEO 400 km — \|Δε/ε\| | 2,9 × 10⁻¹⁰ | 2,7 × 10⁻¹⁰ | 6,3 × 10⁻¹⁰ | 5,4 × 10⁻⁹ |
| GTO 200 × 35 786 km (período 10,5 h) — maior \|Δr\| | 0,147 m | 13,3 m | 1,49 km | 150 km |
| GTO — \|Δε/ε\| | 1,3 × 10⁻⁹ | 5,3 × 10⁻⁹ | 5,2 × 10⁻⁸ | 5,2 × 10⁻⁷ |

Leitura dos números:

- A energia deriva lentamente: na GTO, a variação cresce 10 vezes a cada década de órbitas
  a partir de 10 órbitas; na LEO, 8,6 vezes de 100 para 1 000 órbitas. Um pequeno erro de
  energia é um pequeno erro no período, e o atraso ou adiantamento ao longo da órbita se
  acumula a cada volta. Por isso, a partir de 10 órbitas, o erro de posição cresce
  aproximadamente com o **quadrado** do número de órbitas (de 100 a 200 vezes por década
  nestas medições).
- Quase todo o erro é **ao longo da trajetória** (fase). Na LEO após 1 000 órbitas, a
  componente radial do erro final foi de 3,8 cm, contra cerca de 170 m na direção do
  movimento; a componente normal ao plano foi nula nestas órbitas equatoriais.
- O controle de passo atua componente a componente (x, y, z), por isso o erro depende também
  da orientação da órbita nos eixos. Para a GTO após 5 órbitas, medimos 2,84 m no plano
  equatorial e 0,507 m com i = 27° e ω = 270° (cenário `gto` da seção 2).
- Para simulações longas, reduza `rtol` (por exemplo, `--rtol 1e-12`) e confira o erro
  informado. No modelo de dois corpos a solução de Kepler é exata e preferível; o integrador
  existe para os modelos com perturbações das próximas etapas, que não têm solução fechada.

Script usado (executar a partir da raiz do repositório):

```python
from aurora.data.bodies import EARTH
from aurora.physics.gravitation import TwoBodyGravity
from aurora.physics.orbits import orbital_period
from aurora.simulation.initial_conditions import circular_orbit_state, elliptic_orbit_state
from aurora.simulation.propagator import propagate
from aurora.simulation.validation import compare_with_kepler, conservation_errors

model = TwoBodyGravity(EARTH.mu)
casos = {
    "LEO 400 km": (circular_orbit_state(EARTH, 400e3), EARTH.radius + 400e3),
    "GTO 200 × 35 786 km": (elliptic_orbit_state(EARTH, 200e3, 35_786e3), EARTH.radius + (200e3 + 35_786e3) / 2),
}
for nome, ((r, v), a) in casos.items():
    for n in (1, 10, 100, 1000):
        tr = propagate(model, r, v, n * orbital_period(EARTH.mu, a), n_points=200)
        c, e = compare_with_kepler(tr, EARTH.mu), conservation_errors(tr, EARTH.mu)
        print(f"{nome} | {n} | {c.max_position_error:.2e} m | {e.energy:.1e} | {tr.evaluations}")
```

Tempo de execução medido: 1 000 órbitas da GTO (810 650 avaliações) levaram 9,2 s; o custo
foi de cerca de 11 a 15 µs por avaliação da função neste computador. O limite padrão de
2 000 000 de avaliações por trecho (`--max-avaliacoes`) interrompeu uma simulação de 10⁶
órbitas LEO após 30,4 s, com mensagem de erro em vez de travar.

## 4. Missão de Hohmann simulada

Reprodução: `aurora hohmann --de "400 km" --para "35786 km" --simular --isp "320 s"`.

| Grandeza | Valor |
|---|---|
| Δv1 (analítico) | 2 397,4725 m/s |
| Δv2 (analítico) | 1 456,4867 m/s |
| Δv total (analítico) | 3 853,9592 m/s |
| Tempo de transferência (analítico) | 19 048,56 s (5,2913 h) |
| Erro numérico do raio na chegada (antes da 2ª queima) | 0,0119 m |
| Erro do semieixo maior da órbita final | 0,0122 m |
| Excentricidade da órbita final simulada (ideal: 0) | 7,0 × 10⁻¹⁰ |
| Propelente (massa seca 1 000 kg, Isp 320 s; padrão: exatamente o necessário) | 2 414,76 kg |
| Propelente por queima | 1 824,14 kg e 590,62 kg |

O Δv total confere com as formas fechadas Δv1 = √(μ/r1)·(√(2r2/(r1 + r2)) − 1) e
Δv2 = √(μ/r2)·(1 − √(2r1/(r1 + r2))) a 10⁻¹³ (`tests/test_maneuvers.py`). O consumo confere
com Tsiolkovsky aplicado ao Δv total: Isp·g0·ln(m0/mf) = Δv1 + Δv2 com erro relativo abaixo
de 10⁻¹² (`tests/test_vehicle_missions.py`).

Estes números **não** são um plano de voo: as queimas são instantâneas (uma queima real de
2,4 km/s dura minutos e tem perdas gravitacionais), o lançador e a fase de subida não são
modelados, e a órbita final não sofre perturbações.

## 5. Ambientes em que a suíte foi executada

Medido em Linux x86-64 (`pytest`, 428 testes):

| Ambiente | Resultado |
|---|---|
| Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10.9, Plotly 7.1.0 | 428 aprovados |
| Python 3.11.17, NumPy 2.4.6, SciPy 1.17.1, Matplotlib 3.11.2, Plotly 7.1.0 | 428 aprovados |
| Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2, Plotly 7.1.0 | 428 aprovados |
| Python 3.13.16, NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2, Plotly 7.1.0 | 428 aprovados |
| Python 3.14.6, NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2, Plotly 7.1.0 | 428 aprovados |
| Versões mínimas: Python 3.11.17, NumPy 1.24.4, SciPy 1.10.1, Matplotlib 3.6.3, Plotly 5.0.0 | 428 aprovados (avisos de depreciação internos do Matplotlib 3.6) |
| Python 3.12.3, NumPy 2.5.4, SciPy 1.18.1, **sem** Matplotlib e Plotly | 425 aprovados, 3 pulados (gráficos) |

A página HTML interativa foi aberta no Chromium em tela de celular (390 × 844 px, modo claro)
e de computador (1280 × 900 px, modo escuro): três gráficos desenhados, sem rolagem
horizontal e sem erros de JavaScript. **Não** houve execução em aparelho Android.
