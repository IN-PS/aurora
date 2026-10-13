# AURORA v1.0

*Em desenvolvimento — etapa 1 concluída (versão do pacote: 0.1.0).*

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
  fórmula numericamente estável: os decimais digitados são lidos sem arredondamento
  (`x² + 0,6x + 0,09` tem raiz dupla), raízes próximas não são confundidas com raiz dupla, e
  coeficientes de 10⁻³²⁴ a 10³⁰⁸ dão raízes com erro de poucos ulp.
- **Motor de fórmulas físicas**: 18 fórmulas registradas, em unidades do SI, cada uma com
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
- **220 testes automatizados** (pytest), incluindo identidades analíticas, casos de
  referência conhecidos (órbita geoestacionária, período de uma órbita de 1 ua), conversão
  ida e volta elementos → estado → elementos e entradas extremas. Integração contínua no
  GitHub Actions com Python 3.10 a 3.14.

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
  usual em software científico (a ISO 80000-1 admite ambos).
- **Não use separador de milhar.** Um único `.` ou `,` é sempre o separador decimal:
  `"1.000 km"` vale **1 km**, não 1000 km. Números com mais de um separador
  (`"6.778.137 km"`) são rejeitados com mensagem explicativa. As entradas são sempre
  ecoadas convertidas para o SI, para conferência.
- Números negativos funcionam diretamente (`aurora quadratica 1 -1e8 1`).
- As opções podem vir em qualquer posição (`aurora calcular F --json mu=1 r=1`).
- `--json` exporta o resultado com entradas, unidades, hipóteses e procedência. O texto é
  escrito em ASCII com escapes `\uXXXX` (μ, ç, ☉…), que qualquer leitor de JSON decodifica,
  para não ser corrompido em terminais que não usam UTF-8. Valores não finitos, como um Δ
  que excede o `float64`, aparecem como os textos `"inf"`, `"-inf"` ou `"nan"`.

## Uso no Android

A etapa 1 **não tem interface gráfica nem aplicativo Android**: a interface é a linha de
comando, com saída em texto ou JSON. No Android ela pode ser usada em um terminal Linux,
como o [Termux](https://termux.dev):

1. Instale o Termux pelo F-Droid ou pelas versões publicadas no GitHub do projeto (os
   mantenedores não recomendam a versão da Google Play, que está desatualizada).
2. No Termux:

   ```bash
   pkg update
   pkg install git python python-numpy   # NumPy pré-compilado do repositório do Termux
   git clone https://github.com/fishing-ainpm/Radar-system-.git
   cd Radar-system-
   python -m aurora --ajuda
   python -m aurora orbita-circular --corpo terra --altitude "400 km"
   python -m aurora quadratica 1 -3 2
   ```

3. Opcional: `pip install pytest` e `python -m pytest` para executar os testes;
   `pip install -e .` para ter o comando `aurora`. Prefira o pacote `python-numpy` a
   `pip install numpy`, que no Termux exige compilação.
4. Para abrir um resultado em outro aplicativo, exporte em JSON
   (`python -m aurora ... --json > resultado.json`) e use `termux-setup-storage` para dar
   acesso ao armazenamento compartilhado do aparelho.

**Situação da verificação:** estes passos ainda não foram executados em um aparelho Android.
Foi verificado, em Linux x86-64, que todos os testes passam com Python 3.14.6 e NumPy 2.5.4
(o repositório do Termux distribui Python 3.14 e NumPy 2.4 em outubro de 2026) e com Python
3.10 a 3.13. Um painel gráfico (aplicação web local, que poderia ser aberta no navegador do
celular) está apenas proposto para as etapas 2–4.

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
hierarquia `aurora.errors.AuroraError` com mensagens em português. As fórmulas são
registradas ao importar `aurora.physics` (ou qualquer um de seus módulos). Coeficientes
`int` ou `fractions.Fraction` passados a `solve_quadratic` são usados com valor exato;
`float` é usado com seu valor binário exato (`0.1` em `float` não é exatamente 1/10).

## Testes

```bash
pytest
```

O fluxo `.github/workflows/testes.yml` executa a suíte no GitHub Actions com Python 3.10,
3.11, 3.12, 3.13 e 3.14, e com a versão mínima declarada do NumPy (1.24).

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
  AUDITORIA_ETAPA1.md  relatório da auditoria técnica da etapa 1 e correções aplicadas
.github/workflows/     integração contínua (pytest em Python 3.10–3.14)
```

## Documentação

- [Arquitetura e roteiro](docs/ARQUITETURA.md)
- [Modelos físicos, hipóteses e validação](docs/MODELOS_FISICOS.md)
- [Auditoria técnica da etapa 1](docs/AUDITORIA_ETAPA1.md)
