# Arquitetura do AURORA

Este documento descreve a organização do código, as decisões tecnológicas (com justificativa),
a separação entre o simulador e um futuro software de hardware, e o roteiro das próximas
etapas.

## 1. Princípios

1. **Modelos físicos separados de interface e de dados.** Fórmulas são funções puras em SI,
   sem `print` nem leitura de arquivos; a interface apenas lê entradas e exibe resultados; os
   dados de referência ficam em módulos próprios, sempre com a fonte.
2. **Unidades do SI internamente.** Conversões (km, ua, h…) acontecem somente nas bordas
   (entrada e exibição), com verificação de dimensão.
3. **Hipóteses explícitas.** Toda fórmula registrada declara suas hipóteses, exibidas junto
   com o resultado.
4. **Procedência explícita.** Todo resultado carrega sua origem (cálculo analítico, simulação
   ou medição), e o simulador nunca produz dados marcados como medição.
5. **Validação na entrada.** Parâmetros inválidos ou fora do domínio geram erros claros, em
   português, em vez de `NaN` silenciosos.
6. **Primeiro uma versão executável e testada.** Cada etapa só começa depois que a anterior
   está validada por testes; não há módulos vazios ou funcionalidades simuladas de fachada.

## 2. Camadas

```text
┌──────────────────────────────────────────────────────────────┐
│ aurora.interface   linha de comando (etapa 1); painel web     │
│                    (etapas 2–4, a definir — ver seção 5)       │
├──────────────────────────────────────────────────────────────┤
│ aurora.physics     modelos: quadratic, kinematics,            │
│                    gravitation, orbits (+ hypotheses)          │
├──────────────────────────────────────────────────────────────┤
│ aurora.core        validation, units, formula (registro),     │
│                    provenance (DataOrigin)                     │
├──────────────────────────────────────────────────────────────┤
│ aurora.data        constants, bodies (valores com fonte)      │
├──────────────────────────────────────────────────────────────┤
│ aurora.errors      hierarquia de exceções (sem dependências)  │
└──────────────────────────────────────────────────────────────┘
```

Regra de dependência: cada camada importa apenas as camadas **abaixo** dela. `physics` não
conhece a interface; `data` não conhece os modelos. Isso permite trocar ou acrescentar
interfaces (CLI, painel web, notebooks) sem tocar nos modelos.

## 3. Componentes da etapa 1

| Módulo | Responsabilidade |
|---|---|
| `errors.py` | `AuroraError` e subclasses: `InvalidParameterError`, `DomainError`, `UnitError`, `NotFoundError`. |
| `data/constants.py` | G (CODATA 2018), c, g0, ua — com fonte, exatidão e incerteza relativa. |
| `data/bodies.py` | Sol, Terra, Lua: μ e raio de referência, com fontes (IAU 2015, WGS 84, JPL DE430, IAU WGCCRE). |
| `core/validation.py` | Conversão e validação de escalares e vetores (rejeita NaN, ±∞, `bool`, texto, complexos). |
| `core/units.py` | Tabela de unidades por dimensão, conversão SI ↔ unidade, leitura de `"400 km"`. |
| `core/formula.py` | Decorador `register_formula`: associa a cada função sua expressão, unidades, parâmetros e hipóteses; confere a assinatura na importação; `evaluate()` devolve `FormulaResult`. |
| `core/provenance.py` | `DataOrigin` (`ANALYTIC`, `SIMULATED`, `MEASURED`) com rótulos e avisos. |
| `physics/quadratic.py` | Equação do 2º grau com discriminante exato e fórmula estável. |
| `physics/kinematics.py` | MUV em 1D (fórmulas registradas) e 3D (vetorial); instantes de passagem por uma posição. |
| `physics/gravitation.py` | Lei de Newton, μ = G·M, aceleração vetorial; protocolo `GravityModel` e `NewtonianGravity` (superposição de massas pontuais). |
| `physics/orbits.py` | Fórmulas de órbitas circulares/elípticas; `orbit_from_state` (cônica a partir de r⃗, v⃗); resumo de órbita circular. |
| `interface/cli.py` | Comandos `quadratica`, `formulas`, `calcular`, `orbita-circular`, `corpos`, `constantes`; saída em texto ou JSON. |

### Pontos de extensão já previstos

- **Novas fórmulas** (ex.: equação do radar na etapa 3) entram com `@register_formula` e
  aparecem automaticamente na interface, com unidades e hipóteses.
- **Modelos gravitacionais** implementam `acceleration(position, velocity, time)`
  (`GravityModel`). O integrador da etapa 2 dependerá apenas desse protocolo, de modo que
  modelos mais completos — fontes móveis lidas de efemérides, harmônicos zonais (J2),
  correções relativísticas pós-newtonianas (que dependem da velocidade) — possam ser
  acrescentados sem alterá-lo.
- **Procedência** (`DataOrigin`) será usada pela telemetria simulada (etapa 2), pelos alvos
  e ecos do radar virtual (etapa 3) e pelo contrato de dados com o hardware (etapa 5).
- **Efemérides** (etapa 4) substituirão posições fixas por fontes dependentes do tempo, por
  meio de uma camada de importação em `aurora.data`.

### O que funciona e o que ainda é protótipo (etapa 1)

| Componente | Situação |
|---|---|
| Calculadora do 2º grau, 18 fórmulas registradas, unidades, validação, CLI | **Funcional e testado** |
| `orbit_from_state`, `constant_acceleration_state`, `time_to_reach_position`, `NewtonianGravity` | **Funcional e testado**, mas só como API Python (sem comando na CLI) |
| `GravityModel` | **Gancho de arquitetura**: protocolo sem consumidor até existir o integrador da etapa 2. Por ser `runtime_checkable`, `isinstance` só verifica a existência do método `acceleration`, não sua assinatura |
| `NewtonianGravity` | Fontes **estáticas**; `velocity` é validada mas não usada e `time` é ignorado (reservados a modelos futuros) |
| `DataOrigin.SIMULATED` / `MEASURED` | **Reservados**: nada no pacote produz esses valores ainda. A regra "o simulador nunca cria `MEASURED`" é convenção documentada, não imposta por código |
| Painel gráfico (dashboard), radar virtual, propagação numérica, efemérides | **Não implementados** (etapas 2–4) |
| Software de hardware | **Fora deste pacote** por decisão de arquitetura (seção 4) |

O registro de fórmulas é preenchido como efeito da importação de `aurora.physics`; quem usa
apenas `aurora.core.formula` deve importar `aurora.physics` antes. Recarregar um módulo de
física (por exemplo, com *autoreload* em notebooks) é aceito sem duplicar registros.

## 4. Separação entre simulador e hardware

O pacote `aurora` é **exclusivamente um simulador**. Requisitos de arquitetura:

1. **Pacotes distintos.** O software de controle de hardware (etapa 5; nome provisório
   `aurora_hw`) será um pacote separado, com distribuição, testes e revisão próprios. O
   simulador não importa nada dele, não contém drivers e não abre portas seriais, sockets
   de equipamentos ou barramentos.
2. **Dependência em um só sentido.** O software de hardware poderá reutilizar contratos de
   dados do simulador (por exemplo, `DataOrigin` e formatos de telemetria), nunca o inverso.
3. **Dados rotulados.** Dados de instrumentos reais são `MEASURED`; dados do simulador são
   `ANALYTIC` ou `SIMULATED`. A interface deve exibir a origem de cada série, e nenhum código
   do simulador pode criar dados `MEASURED`.
4. **Sem comandos a partir do simulador.** Resultados do simulador não são enviados como
   comandos a equipamentos. Qualquer integração futura passará por um plano de testes em
   bancada e ambiente simulado (etapa 5), com análise de segurança.

## 5. Decisões tecnológicas

| Decisão | Escolha | Justificativa |
|---|---|---|
| Linguagem | Python ≥ 3.10 | Preferência do projeto; ecossistema científico maduro. |
| Cálculo numérico | `float64` (IEEE 754) via `math` e NumPy | Precisão de ~15–16 algarismos significativos, suficiente para os modelos idealizados; NumPy para vetores 3D. Aritmética racional (`fractions`) apenas onde necessária (discriminante). |
| SciPy | **Ainda não incluído** | Não é necessário na etapa 1. Será introduzido na etapa 2 para integração numérica com controle de passo (`scipy.integrate.solve_ivp`, métodos RK45/DOP853), com comparação contra soluções analíticas. |
| Unidades | Tabela própria (`core/units.py`) | Cobre as grandezas usadas, sem dependência extra, e é totalmente testada. Os modelos recebem apenas `float` em SI, então uma biblioteca como `pint` pode substituí-la no futuro sem alterar a física. |
| Interface (etapa 1) | Linha de comando com `argparse` (biblioteca padrão) | O repositório estava vazio, sem tecnologia de interface a seguir. A CLI não adiciona dependências, é testável e permite automatizar e registrar experimentos (saída JSON). |
| Gráficos (etapa 2+) | Matplotlib (figuras estáticas de validação) e Plotly (visualização interativa 2D/3D) | Preferências do projeto. Serão introduzidos quando houver trajetórias propagadas a exibir. |
| Painel (etapas 2–4) | **Proposta**: aplicação web local em Python com Plotly Dash | Mantém uma única linguagem, reaproveita as figuras Plotly e roda no navegador sem servidor externo. Alternativa: Streamlit (mais simples, menos controle de layout). A decisão será confirmada antes de adicionar a dependência. |
| Testes | pytest | Preferência do projeto. |
| Documentação | Markdown em pt-BR | Preferência do projeto. |

### Convenções de código

- Identificadores (módulos, funções, classes) em inglês, como os nomes de módulos da
  especificação (`physics_engine`, `radar_simulator`…); docstrings, mensagens de erro,
  chaves de fórmulas e toda a interface em português brasileiro. As mensagens de erro de uso
  geradas pelo `argparse` (em inglês na biblioteca padrão) são traduzidas por uma tabela de
  padrões; uma mensagem sem padrão correspondente é exibida como está.
- Integração contínua (`.github/workflows/testes.yml`): a suíte roda em Python 3.10–3.14 e
  com a versão mínima declarada do NumPy, porque o comportamento do `argparse` muda entre
  versões do Python (um defeito real da CLI só aparecia no 3.10–3.12).
- Parâmetros das fórmulas usam os símbolos físicos (`mu`, `r`, `v0`, `T`…).
- Grandezas ilimitadas (apoapse de hipérbole, período de trajetória aberta) são `math.inf`;
  grandezas indefinidas (inclinação de trajetória radial) são `math.nan`.

## 6. Roteiro

As etapas seguintes só serão iniciadas após validação da anterior.

- **Etapa 2 — Simulação orbital.** Integrador numérico com controle de passo e tolerâncias
  configuráveis sobre `GravityModel`; validação contra órbitas analíticas (conservação de
  energia e momento angular, período, retorno à posição inicial); gráficos de posição,
  velocidade e altitude; telemetria simulada marcada como `SIMULATED`; transferência de
  Hohmann como caso analítico de referência.
- **Etapa 3 — Radar virtual.** Antena, feixe, alcance e varredura; alvos estáticos e móveis;
  distância (tempo de ida e volta, c·Δt/2) e velocidade radial (efeito Doppler) com ruído,
  resolução e incerteza; visualização de incertezas. Documentação explícita de que o
  simulador não detecta objetos reais.
- **Etapa 4 — Exploração espacial.** Visualizador do Sistema Solar com escalas ajustáveis;
  missões virtuais; registro de missões; importação opcional de efemérides de fontes
  confiáveis (por exemplo, arquivos SPICE/JPL), sem embutir valores não verificados.
- **Etapa 5 — Preparação para hardware.** Documento de requisitos de um satélite
  experimental; arquitetura do pacote separado de sensores e telemetria; plano de testes em
  bancada e ambiente simulado; análise de segurança, orçamento e viabilidade.
