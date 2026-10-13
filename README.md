# AURORA v1.0

*Em desenvolvimento — etapas 1 e 2 concluídas (versão do pacote: 0.2.0).*

**Plataforma educacional de simulação de radar, satélites e exploração espacial**
Instituto Regional de Pesquisas Científicas, Tecnologia e Engenharia do Sul de Minas Gerais

> **Aviso.** Todos os resultados do AURORA vêm de cálculos analíticos ou simulações com
> modelos idealizados. Nenhum resultado é medição de objeto real, o simulador não detecta
> objetos reais e nenhuma trajetória calculada constitui plano de voo certificado. O pacote
> `aurora` não contém código para controlar hardware (ver [Arquitetura](docs/ARQUITETURA.md)).

## Estado do projeto

| Etapa | Conteúdo | Situação |
|---|---|---|
| 1. Fundamentos matemáticos | Estrutura do projeto, equação do 2º grau, motor de fórmulas físicas com unidades, testes | **Concluída** (pacote `0.1.0`) |
| 2. Simulação orbital | Dois corpos (analítico e numérico com controle de passo), validação contra Kepler, gráficos 2D/3D, veículo virtual, Hohmann, telemetria simulada | **Concluída** (pacote `0.2.0`) |
| 3. Radar virtual | Varredura simulada, alvos virtuais, distância e velocidade radial com incerteza | Planejada |
| 4. Exploração espacial | Sistema Solar, sondas interplanetárias, registro de missões, efemérides | Planejada |
| 5. Preparação para hardware | Requisitos, arquitetura separada de sensores/telemetria, plano de testes | Planejada |

## O que a etapa 2 oferece

- **Propagação analítica de dois corpos** (formulação por variável universal, coeficientes
  de Lagrange): elipses, parábolas, hipérboles e trajetórias radiais, para qualquer
  intervalo de tempo, positivo ou negativo. É a referência exata do modelo, usada para medir
  o erro do integrador.
- **Propagação numérica com controle de passo** (`scipy.integrate.solve_ivp`, Dormand–Prince
  DOP853 ou RK45) com tolerâncias em SI (padrão: `rtol = 1e-10`, `atol = 0,1 mm` e
  `0,1 µm/s`), detecção de colisão com a superfície e limite de avaliações. O erro global de
  cada simulação é **medido** contra a solução analítica e exibido junto com o resultado.
- **Validação contra casos analíticos**: seis cenários padrão (órbita baixa, geoestacionária,
  elipses de transferência e tipo Molniya, órbita lunar, hipérbole). Com as tolerâncias
  padrão, o maior erro de posição medido vai de 2,6 mm (órbita lunar, 5 voltas) a 0,51 m
  (elipse 200 × 35 786 km, 5 voltas). Ver [Validação da etapa 2](docs/VALIDACAO_ETAPA2.md).
- **Gráficos**: figura estática (`.png`, `.svg`, `.pdf`) com projeção no plano da órbita,
  vista 3D, altitude, velocidade, erro contra Kepler e variação da energia; página HTML
  interativa autocontida (órbita 3D, altitude, velocidade, tabela-resumo, modo escuro), que
  funciona sem internet, inclusive no navegador de um celular.
- **Veículo virtual** com massa seca, propelente e impulso específico; consumo pela equação
  de Tsiolkovsky.
- **Condições iniciais idealizadas**: estado no fim da queima do lançador (altitude,
  velocidade, ângulo de trajetória, plano orbital), órbitas circulares e elípticas.
- **Transferência de Hohmann** analítica e simulada (duas queimas impulsivas), com o erro
  numérico da chegada e o propelente consumido.
- **Telemetria simulada** em CSV ou JSON, com modelo, integrador, tolerâncias, hipóteses e a
  marca de origem `simulacao` em todo arquivo.
- **4 fórmulas novas** no motor de fórmulas (22 no total): Δv de Tsiolkovsky, massa de
  propelente, Δv e tempo de Hohmann.

Nada disso é plano de voo: as manobras são impulsivas e o modelo despreza arrasto,
achatamento da Terra, terceiros corpos e pressão de radiação (lista completa em
[Modelos físicos](docs/MODELOS_FISICOS.md)).

## O que a etapa 1 oferece

- **Calculadora de equações do 2º grau** com discriminante exato (aritmética racional) e
  fórmula numericamente estável: os decimais digitados são lidos sem arredondamento
  (`x² + 0,6x + 0,09` tem raiz dupla), raízes próximas não são confundidas com raiz dupla, e
  coeficientes de 10⁻³²⁴ a 10³⁰⁸ dão raízes com erro de no máximo ≈ 2 ulp.
- **Motor de fórmulas físicas**: 18 fórmulas registradas na etapa 1, em unidades do SI, cada uma com
  expressão, unidade do resultado, unidade e descrição de cada parâmetro e lista de
  hipóteses — acessíveis pela interface (`aurora formulas` / `aurora calcular`):
  - cinemática do movimento uniformemente variado em 1D (posição, velocidade,
    deslocamento, aceleração média, Torricelli);
  - gravitação universal de Newton (força, μ = G·M, aceleração);
  - órbitas circulares e elípticas (velocidade circular, período, raio por período, energia
    específica, vis-viva, velocidade de escape, semieixo maior, periapse/apoapse).
- **API Python adicional** (somente como biblioteca, sem comando na interface): estado 3D
  no MUV, instantes de passagem por uma posição, superposição gravitacional de várias
  massas pontuais e classificação da cônica a partir do vetor de estado.
- **Unidades**: entrada em `km`, `ua`, `h`, `km/s`, `km^3/s^2` etc., com conversão para o SI e
  verificação de dimensão (não aceita `400 s` onde se espera um comprimento).
- **Validação de entradas** com mensagens em português (NaN, infinito, texto, valores fora
  do domínio, órbita abaixo da superfície, valores fora do intervalo do `float64`…).
- **Procedência explícita**: todo resultado físico (fórmulas e órbitas) é marcado como
  *cálculo analítico* e exibido com aviso de que não é medição.
- **Interface de linha de comando** em português, com exportação em JSON.
- **Testes automatizados** (pytest), incluindo identidades analíticas, casos de referência
  conhecidos (órbita geoestacionária, período de uma órbita de 1 ua), conversão ida e volta
  elementos → estado → elementos e entradas extremas.

## Instalação

Requisitos: Python 3.10 ou superior, NumPy e SciPy. Os gráficos usam Matplotlib e Plotly,
que são opcionais (extra `graficos`): sem eles, tudo funciona, exceto `--grafico`.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[graficos,dev]"   # AURORA + Matplotlib/Plotly + pytest
```

Sem instalar, os comandos também funcionam a partir da raiz do repositório com
`python -m aurora ...` (desde que as dependências estejam disponíveis).

## Uso pela linha de comando

```text
aurora --ajuda                     # lista os comandos
aurora quadratica A B C            # resolve A·x² + B·x + C = 0
aurora formulas [CHAVE]            # lista as fórmulas ou mostra detalhes/hipóteses de uma delas
aurora calcular CHAVE NOME=VALOR…  # avalia uma fórmula (unidades opcionais em cada valor)
aurora orbita-circular --corpo terra --altitude "400 km"
aurora corpos                      # corpos de referência (Sol, Terra, Lua) e fontes dos dados
aurora constantes                  # constantes físicas e fontes

# etapa 2 — simulação orbital
aurora propagar --altitude "400 km" --inclinacao 51.6 --orbitas 3 --grafico orbita.html
aurora propagar --periapse "200 km" --apoapse "35786 km" --orbitas 2 --telemetria gto.csv
aurora propagar --posicao "7000 0 0 km" --velocidade "0 12 0 km/s" --duracao "1 d"
aurora hohmann --de "400 km" --para "35786 km"                     # cálculo analítico
aurora hohmann --de "400 km" --para "35786 km" --simular --isp "320 s" --grafico hohmann.png
aurora validar                     # cenários de validação contra a solução de Kepler
```

Em `propagar`, a condição inicial é uma de três: `--altitude` (órbita circular),
`--periapse` com `--apoapse` (elipse, partindo do periapse) ou `--posicao` com
`--velocidade` (vetor de estado no referencial inercial centrado no corpo). A duração é
`--duracao` ou `--orbitas`. Opções do integrador: `--rtol`, `--metodo {DOP853,RK45}` e
`--max-avaliacoes`. Saídas: `--pontos`, `--telemetria arquivo.csv|.json`,
`--grafico arquivo.png|.svg|.pdf|.html` (pode repetir) e `--json`.

Exemplos (saídas reais, abreviadas):

```text
$ aurora quadratica 1 -3 2
Equação: (1)·x² + (-3)·x + (2) = 0
Discriminante: Δ = b² − 4·a·c = 1
Classificação: duas raízes reais distintas (Δ > 0)
  x₁ = 1
  x₂ = 2

$ aurora orbita-circular --corpo terra --altitude "400 km"
Órbita circular em torno de: Terra
  μ = 3.986004418e+14 m^3/s^2; raio equatorial R = 6378137 m  (6378.137 km)
  ...
  Velocidade orbital circular                 v = 7668.558175 m/s  (7.668558 km/s; 27606.81 km/h)
  Período orbital circular                    T = 5553.624271 s  (92.5604 min; 1.542673 h; 0.06427806 d)
  Energia orbital específica (dois corpos)    epsilon = -29403392.24 J/kg  (-29.40339 MJ/kg)
  ...
Hipóteses:
  - Órbita perfeitamente circular (excentricidade nula).
  - Problema de dois corpos: massa do veículo desprezível frente à do corpo central, ...
  ...
Natureza do resultado: cálculo analítico com modelo idealizado — não é medição nem dado observacional.

$ aurora calcular raio_orbita_circular mu=3.986004418e14 "T=86164.0905 s" --em km
  Resultado: r = 42164169.62 m  (42164.17 km)

$ aurora propagar --altitude "400 km" --inclinacao 51.6 --orbitas 3
Simulação numérica de órbita — Terra
  Condição inicial (idealizada): órbita circular de altitude 400000 m  (400 km), inclinação 51.6°
  ...
  Integrador: DOP853, rtol = 1e-10, atol = 0.0001 m / 1e-07 m/s; 1256 avaliações da função
  Duração simulada: 16660.87281 s  (277.6812 min; 4.62802 h; 0.1928342 d); 500 amostras
  ...
Verificação do erro numérico (contra a solução analítica de Kepler, mesmo modelo):
  maior erro de posição: 0.00463 m (no fim: 0.0045 m); maior erro de velocidade: 5.11e-06 m/s
  maior variação relativa da energia: 3.01e-10; do momento angular: 1.5e-10
  ...
Natureza do resultado: simulação numérica com modelo idealizado — não representa detecção ou medição de objeto real.

$ aurora hohmann --de "400 km" --para "35786 km" --simular --isp "320 s"
Transferência de Hohmann — Terra (cálculo analítico)
  ...
  Δv total = 3853.959216 m/s  (3.853959 km/s; 13874.25 km/h)
  Tempo de transferência = 19048.56251 s  (317.476 min; 5.291267 h; 0.2204695 d)

Missão simulada (integração numérica, queimas impulsivas):
  Veículo: massa seca 1000 kg; propelente 2414.7617 kg; Isp 320 s
  ...
  Erro numérico do raio na chegada: 0.0119 m; do semieixo final: 0.0122 m
```

O "erro numérico" é a diferença entre a integração e a solução analítica **do mesmo
modelo**: mede a qualidade do integrador, não a distância até um satélite real.

Observações:

- Valores sem unidade são interpretados na unidade SI do parâmetro (veja `aurora formulas CHAVE`).
- A entrada aceita ponto ou vírgula decimal (`9,81` ou `9.81`); a saída usa ponto, como é
  usual em software científico (a ISO 80000-1 admite ambos).
- **Não use separador de milhar.** Um único `.` ou `,` é sempre o separador decimal:
  `"1.000 km"` vale **1 km**, não 1000 km. Números com mais de um separador
  (`"6.778.137 km"`) são rejeitados com mensagem explicativa. As entradas são sempre
  ecoadas convertidas para o SI, para conferência.
- Números negativos funcionam diretamente (`aurora quadratica 1 -1e8 1`).
- As opções podem vir em qualquer posição (`aurora calcular velocidade_escape --json mu=1 r=1`)
  e podem ser abreviadas quando não houver ambiguidade (`--alt` para `--altitude`).
- Em `quadratica`, coeficientes não nulos menores que o menor `float64` (≈ 4,9 × 10⁻³²⁴) são
  rejeitados, para que a equação exibida seja sempre a resolvida.
- Saídas longas podem ser canalizadas (`aurora validar | head`) sem mensagens de erro.
- `--json` exporta o resultado com entradas, unidades, hipóteses e procedência. O texto é
  escrito em ASCII com escapes `\uXXXX` (μ, ç, ☉…), que qualquer leitor de JSON decodifica,
  para não ser corrompido em terminais que não usam UTF-8. Valores não finitos, como um Δ
  que excede o `float64`, aparecem como os textos `"inf"`, `"-inf"` ou `"nan"`.

## Uso no Android

O AURORA **não tem aplicativo Android**. No celular ele roda em um terminal Linux, como o
[Termux](https://termux.dev), e os gráficos interativos são páginas HTML que abrem no
navegador do aparelho:

1. Instale o Termux pelo F-Droid ou pelas versões publicadas no GitHub do projeto (os
   mantenedores recomendam essas fontes; a versão da Google Play é um ramo experimental,
   com funcionalidades ausentes).
2. No Termux, instale as bibliotecas científicas **pré-compiladas** do repositório do Termux
   (via `pip`, NumPy, SciPy e Matplotlib exigiriam compilação no aparelho); o Plotly é
   Python puro e vem do `pip`:

   ```bash
   pkg update
   pkg install git python python-numpy python-scipy matplotlib
   pip install plotly                  # opcional: só para os gráficos .html
   git clone https://github.com/fishing-ainpm/Radar-system-.git
   cd Radar-system-
   python -m aurora --ajuda
   python -m aurora validar
   python -m aurora hohmann --de "400 km" --para "35786 km" --simular --grafico hohmann.html
   ```

3. Para ver o gráfico: `termux-open hohmann.html` abre a página no navegador padrão. Outra
   opção é dar acesso ao armazenamento (`termux-setup-storage`) e gravar direto na pasta de
   downloads: `--grafico ~/storage/downloads/hohmann.html`. A página tem o Plotly embutido
   (≈ 5 MB) e funciona sem internet; é possível girar e ampliar a órbita 3D com os dedos e
   tocar nas curvas para ler os valores.
4. Opcional: `pip install pytest` e `python -m pytest` para executar os testes;
   `pip install -e .` para ter o comando `aurora`.

**Situação da verificação:** estes passos **ainda não foram executados em um aparelho
Android**. Foi verificado em Linux x86-64: a suíte de testes passa com Python 3.10 a 3.14
(o Termux distribui Python 3.14, NumPy 2.4, SciPy 1.18 e Matplotlib 3.11 em outubro de 2026);
a página HTML foi aberta no Chromium em tela de celular (390 × 844, modo claro) e de
computador (modo escuro), sem rolagem horizontal nem erros de JavaScript.

## Uso como biblioteca

```python
from aurora.physics.quadratic import solve_quadratic
from aurora.physics.orbits import circular_orbital_speed, orbit_from_state
from aurora.core.formula import get_formula
from aurora.data.bodies import EARTH

solve_quadratic(1, -3, 2).roots                 # (1.0, 2.0)
circular_orbital_speed(EARTH.mu, 6_778_137.0)   # 7668.55... (m/s)

resultado = get_formula("periodo_orbital_circular").evaluate(mu=EARTH.mu, r=6_778_137.0)
resultado.value, resultado.unit                  # (5553.62..., 's')
resultado.formula.hypotheses                     # hipóteses do modelo
resultado.to_dict()                              # exportação com procedência

orbita = orbit_from_state(EARTH.mu, [7.0e6, 0, 0], [0, 8.5e3, 0])
orbita.orbit_type, orbita.eccentricity, orbita.period
```

Simulação orbital (etapa 2):

```python
import math
from aurora.data.bodies import EARTH
from aurora.physics.gravitation import TwoBodyGravity
from aurora.physics.kepler import propagate_kepler
from aurora.simulation.initial_conditions import circular_orbit_state
from aurora.simulation.propagator import IntegratorSettings, propagate
from aurora.simulation.validation import compare_with_kepler
from aurora.simulation.missions import simulate_hohmann
from aurora.simulation.telemetry import write_telemetry

r0, v0 = circular_orbit_state(EARTH, 400e3, inclination=math.radians(51.6))
traj = propagate(TwoBodyGravity(EARTH.mu), r0, v0, 6 * 3600.0,
                 settings=IntegratorSettings(rtol=1e-10), collision_radius=EARTH.radius)
traj.times, traj.positions, traj.velocities      # amostras (s, m, m/s), somente leitura
compare_with_kepler(traj, EARTH.mu)              # erro numérico medido contra a solução exata
propagate_kepler(EARTH.mu, r0, v0, 3600.0)       # estado analítico após 1 h
write_telemetry("leo.csv", traj, EARTH)          # telemetria simulada com metadados

missao = simulate_hohmann(EARTH, EARTH.radius + 400e3, EARTH.radius + 35_786e3)
missao.transfer.delta_v_total, missao.arrival_radius_error, missao.propellant_used
```

Todas as funções recebem e retornam valores no SI; erros de entrada levantam exceções da
hierarquia `aurora.errors.AuroraError` com mensagens em português (`SimulationError` para
falhas da integração, `DependencyError` quando falta uma biblioteca opcional). As fórmulas
são registradas ao importar `aurora.physics` (ou qualquer um de seus módulos). Coeficientes
`int` ou `fractions.Fraction` passados a `solve_quadratic` são usados com valor exato;
`float` é usado com seu valor binário exato (`0.1` em `float` não é exatamente 1/10).

## Testes

```bash
pytest
```

O fluxo `.github/workflows/testes.yml` executa a suíte no GitHub Actions com Python 3.10 a
3.14, com as versões mínimas declaradas (NumPy 1.24, SciPy 1.10, Matplotlib 3.6,
Plotly 5.0) e sem as bibliotecas de gráficos (os testes de gráficos são então pulados).

## Estrutura do repositório

```text
aurora/
  errors.py            exceções (mensagens em pt-BR)
  data/                dados de referência com fontes: constantes, corpos celestes
  core/                validação, unidades, registro de fórmulas, procedência dos dados
  physics/             modelos: quadratic, kinematics, gravitation, orbits, kepler, maneuvers
  simulation/          propagador numérico, validação, condições iniciais, veículo,
                       telemetria e missões (etapa 2)
  interface/           linha de comando (cli.py, cli_simulation.py) e gráficos (plots.py)
tests/                 testes automatizados (pytest)
docs/
  ARQUITETURA.md       arquitetura, decisões tecnológicas, separação simulador × hardware, roteiro
  MODELOS_FISICOS.md   fórmulas, unidades, hipóteses, precisão numérica, validação e limitações
  VALIDACAO_ETAPA2.md  erros numéricos medidos do propagador e da missão de Hohmann
  AUDITORIA_ETAPA1.md  relatório da auditoria técnica da etapa 1 e correções aplicadas
.github/workflows/     integração contínua (pytest em Python 3.10–3.14)
```

## Documentação

- [Arquitetura e roteiro](docs/ARQUITETURA.md)
- [Modelos físicos, hipóteses e validação](docs/MODELOS_FISICOS.md)
- [Validação da etapa 2: erros numéricos medidos](docs/VALIDACAO_ETAPA2.md)
- [Auditoria técnica da etapa 1](docs/AUDITORIA_ETAPA1.md)
