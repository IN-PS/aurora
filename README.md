# AURORA v1.0

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
| 2. Simulação orbital | Propagação numérica de dois corpos com controle de passo, gráficos, validação | Planejada |
| 3. Radar virtual | Varredura simulada, alvos virtuais, distância e velocidade radial com incerteza | Planejada |
| 4. Exploração espacial | Sistema Solar, sondas interplanetárias, registro de missões, efemérides | Planejada |
| 5. Preparação para hardware | Requisitos, arquitetura separada de sensores/telemetria, plano de testes | Planejada |

## O que a etapa 1 oferece

- **Calculadora de equações do 2º grau** com discriminante exato (aritmética racional) e
  fórmula numericamente estável: não confunde raízes próximas com raiz dupla e não perde
  precisão por cancelamento.
- **Motor de fórmulas físicas** (18 fórmulas) em unidades do SI, cada uma com expressão,
  unidade do resultado, unidade e descrição de cada parâmetro e lista de hipóteses:
  - cinemática do movimento uniformemente variado (posição, velocidade, deslocamento,
    aceleração média, Torricelli), em 1D e em 3D (vetorial);
  - gravitação universal de Newton (força, μ = G·M, aceleração, superposição de várias massas);
  - órbitas circulares e elípticas (velocidade circular, período, raio por período, energia
    específica, vis-viva, velocidade de escape, semieixo maior, periapse/apoapse) e
    classificação da cônica a partir do vetor de estado.
- **Unidades**: entrada em `km`, `ua`, `h`, `km/s`, `km^3/s^2` etc., com conversão para o SI e
  verificação de dimensão (não aceita `400 s` onde se espera um comprimento).
- **Validação de entradas** com mensagens em português (NaN, infinito, texto, valores fora
  do domínio, órbita abaixo da superfície…).
- **Procedência explícita**: todo resultado é marcado como *cálculo analítico* e exibido com
  aviso de que não é medição.
- **Interface de linha de comando** em português, com exportação em JSON.
- **143 testes automatizados** (pytest), incluindo identidades analíticas e casos de
  referência conhecidos (órbita geoestacionária, período de uma órbita de 1 ua).

## Instalação

Requisitos: Python 3.10 ou superior. Única dependência de execução: NumPy.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[dev]"            # instala o AURORA e o pytest
```

Sem instalar, os comandos também funcionam a partir da raiz do repositório com
`python -m aurora ...` (desde que o NumPy esteja disponível).

## Uso pela linha de comando

```text
aurora --ajuda                     # lista os comandos
aurora quadratica A B C            # resolve A·x² + B·x + C = 0
aurora formulas [CHAVE]            # lista as fórmulas ou mostra detalhes/hipóteses de uma delas
aurora calcular CHAVE NOME=VALOR…  # avalia uma fórmula (unidades opcionais em cada valor)
aurora orbita-circular --corpo terra --altitude "400 km"
aurora corpos                      # corpos de referência (Sol, Terra, Lua) e fontes dos dados
aurora constantes                  # constantes físicas e fontes
```

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
```

Observações:

- Valores sem unidade são interpretados na unidade SI do parâmetro (veja `aurora formulas CHAVE`).
- A entrada aceita ponto ou vírgula decimal (`9,81` ou `9.81`); a saída usa ponto, como é
  usual em software científico (a ISO 80000-1 admite ambos). Não use separador de milhar.
- Números negativos funcionam diretamente (`aurora quadratica 1 -1e8 1`).
- `--json` exporta o resultado com entradas, unidades, hipóteses e procedência.

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

Todas as funções recebem e retornam valores no SI; erros de entrada levantam exceções da
hierarquia `aurora.errors.AuroraError` com mensagens em português.

## Testes

```bash
pytest
```

## Estrutura do repositório

```text
aurora/
  errors.py            exceções (mensagens em pt-BR)
  data/                dados de referência com fontes: constantes, corpos celestes
  core/                validação, unidades, registro de fórmulas, procedência dos dados
  physics/             modelos: quadratic, kinematics, gravitation, orbits, hypotheses
  interface/cli.py     linha de comando (pt-BR)
tests/                 testes automatizados (pytest)
docs/
  ARQUITETURA.md       arquitetura, decisões tecnológicas, separação simulador × hardware, roteiro
  MODELOS_FISICOS.md   fórmulas, unidades, hipóteses, precisão numérica, validação e limitações
```

## Documentação

- [Arquitetura e roteiro](docs/ARQUITETURA.md)
- [Modelos físicos, hipóteses e validação](docs/MODELOS_FISICOS.md)
