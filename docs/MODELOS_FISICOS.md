# Modelos físicos e matemáticos — etapas 1 e 2

Este documento lista cada modelo implementado, com unidades, domínio de validade, hipóteses,
cuidados numéricos e a forma como foi validado. As mesmas hipóteses são exibidas pela
interface (`aurora formulas <chave>` e na saída de `propagar` e `hohmann`).

> Todos os resultados são **cálculos analíticos ou simulações numéricas de modelos
> idealizados**. Não são medições e não devem ser usados como planos de voo.

## 1. Convenções

- Unidades do SI em todos os cálculos: m, s, kg, m/s, m/s², N, m³/s², J/kg (= m²/s²), rad.
- Aritmética de ponto flutuante IEEE 754 de 64 bits (`float64`): cerca de 15–16 algarismos
  significativos. A precisão **exibida** não equivale à exatidão física: as incertezas das
  constantes e as simplificações do modelo dominam o erro.
- Raios orbitais são medidos a partir do **centro** do corpo central; altitude é medida a
  partir do raio de referência do corpo (modelado como esfera). Para a Terra o raio de
  referência é o **equatorial** (WGS 84): fora do equador, a "altitude de 400 km" difere da
  altitude geodésica (acima do elipsoide) em até ≈ 21 km nos polos.
- Semieixo maior: `a > 0` para elipses, `a < 0` para hipérboles, `a = ∞` para a parábola.

## 2. Constantes e dados de referência

| Grandeza | Valor (SI) | Fonte | Observação |
|---|---|---|---|
| G | 6,67430 × 10⁻¹¹ m³ kg⁻¹ s⁻² | CODATA 2018 | incerteza relativa ≈ 2,2 × 10⁻⁵ (não propagada) |
| c | 299 792 458 m/s | SI (definição do metro, 1983) | exata |
| g0 | 9,80665 m/s² | 3ª CGPM (1901) | valor convencional exato; **não** é a gravidade local |
| ua | 149 597 870 700 m | IAU 2012, Resolução B2 | exata |
| μ do Sol | 1,3271244 × 10²⁰ m³/s² | IAU 2015, Resolução B3 (valor nominal) | |
| R do Sol | 6,957 × 10⁸ m | IAU 2015, Resolução B3 (raio nominal) | |
| μ da Terra | 3,986004418 × 10¹⁴ m³/s² | WGS 84 | inclui a atmosfera |
| R da Terra | 6 378 137 m | WGS 84 (semieixo maior do elipsoide) | raio **equatorial** |
| μ da Lua | 4,9028 × 10¹² m³/s² | JPL DE430 (Folkner et al., 2014) | arredondado a 5 algarismos |
| R da Lua | 1 737 400 m | IAU WGCCRE 2009 (Archinal et al., 2011) | raio **médio** |

Usamos μ = G·M diretamente, e não G vezes a massa: μ é determinado por rastreamento orbital
com precisão muito superior à de G e de M isoladamente.

## 3. Equação do segundo grau (`physics/quadratic.py`)

`a·x² + b·x + c = 0`, com `Δ = b² − 4ac` e `x = (−b ± √Δ)/(2a)`.

- **Domínio:** coeficientes reais finitos e `a ≠ 0` (com `a = 0` a equação é linear e o
  programa informa isso).
- **Classificação:** Δ > 0 duas raízes reais distintas; Δ = 0 raiz dupla; Δ < 0 raízes
  complexas conjugadas.
- **Cuidados numéricos:**
  1. Δ é calculado **exatamente** em aritmética racional e arredondado uma única vez. A
     interface lê os decimais digitados sem arredondamento (como frações exatas: `0,1` é
     exatamente 1/10), de modo que a classificação é exata para a equação que o usuário vê:
     `x² + 0,6x + 0,09 = (x + 0,3)²` tem raiz dupla. Na API, coeficientes `int` e
     `Fraction` também são exatos; um `float` é tratado pelo seu valor binário exato (em
     `float64`, `0.6` e `0.09` não são esses decimais, e para esses binários Δ < 0 de fato).
     Exemplo de teste: para `x² − 202687298·x + (101343649² − 1)`, a conta ingênua
     `b·b − 4·a·c` dá 0 (falsa raiz dupla); o valor correto é Δ = 4.
  2. As raízes usam a forma estável `q = −(b + sinal(b)·√Δ)/2`, `x₁ = q/a`, `x₂ = c/q`, que
     evita o cancelamento catastrófico. Exemplo: para `x² + 10⁸x + 1`, a fórmula ingênua
     erra a raiz pequena em ~25 %; a forma estável acerta até o último algarismo.
  3. `√Δ` é aproximada (erro de ~1 ulp) reescalando Δ por potências de 4; `q`, `q/a` e `c/q`
     são calculados em aritmética racional e arredondados uma única vez. Não há overflow nem
     underflow intermediário: coeficientes subnormais (≈ 10⁻³²⁴) ou próximos de 1,8 × 10³⁰⁸
     dão raízes com erro de no máximo ≈ 2 ulp (máximo medido: 1,8 ulp contra referência de
     60 dígitos; os testes usam margem de 4 ulp).
     Só é levantado erro quando a própria raiz não cabe em `float64`. Quando Δ excede o
     intervalo do `float64`, a interface informa "Δ > 1,8 × 10³⁰⁸" — o sinal continua exato.
- **Limitação inerente:** perto de uma raiz dupla o problema é mal condicionado — uma
  perturbação relativa ε nos coeficientes desloca as raízes em ordem de √ε. Isso é
  propriedade do problema, não do algoritmo.

## 4. Cinemática — movimento uniformemente variado (`physics/kinematics.py`)

| Chave | Expressão | Unidade |
|---|---|---|
| `posicao_muv` | x = x0 + v0·t + a·t²/2 | m |
| `velocidade_muv` | v = v0 + a·t | m/s |
| `deslocamento_muv` | Δx = v0·t + a·t²/2 | m |
| `aceleracao_media` | a_m = (v − v0)/Δt, com Δt ≠ 0 | m/s² |
| `velocidade_torricelli` | \|v\| = √(v0² + 2·a·Δx), exige v0² + 2aΔx ≥ 0 | m/s |

Também disponíveis na API (sem comando na interface): `constant_acceleration_state` (as
mesmas relações em 3D, vetoriais) e `time_to_reach_position`, que resolve
`a/2·t² + v0·t + (x0 − x) = 0` em aritmética racional para obter os instantes t ≥ 0 em que o
móvel passa por uma posição.

**Ponto de retorno.** A altura máxima `x = x0 + v0²/(2|a|)`, quando escrita em `float64`, fica
ligeiramente acima ou abaixo do valor exato. Por isso, quando o instante de retorno
`t = −v0/a` não está no passado e `x ≠ x0`, se `|Δ| ≤ 2|a|·(ulp(x0) + ulp(x)) + 4ε·v0²`
(ε = 2,2 × 10⁻¹⁶) a posição é considerada atingida uma única vez, nesse instante. A faixa
corresponde ao erro de representação das próprias entradas: cerca de 1 ulp de `x0` e de `x`
em posição. Fora dela vale a solução exata. A equação de Torricelli trata da mesma forma um
radicando negativo apenas pelo arredondamento (até 8ε·max(v0², |2aΔx|); resultado `|v| = 0`).

**Torricelli:** a relação vale para qualquer instante, passado ou futuro, em que o
deslocamento ocorre; com `v0 = a = 0` nenhum deslocamento não nulo ocorre e é levantado erro.

**Hipóteses:** aceleração constante no intervalo; movimento retilíneo (nas fórmulas 1D);
referencial inercial; mecânica clássica (v ≪ c).

## 5. Gravitação universal (`physics/gravitation.py`)

| Chave | Expressão | Unidade |
|---|---|---|
| `forca_gravitacional` | F = G·m1·m2/r² (r > 0, massas ≥ 0) | N |
| `parametro_gravitacional` | μ = G·M | m³/s² |
| `aceleracao_gravitacional` | g = μ/r² | m/s² |

API vetorial: `point_mass_acceleration` (a⃗ = −μ·r⃗/\|r⃗\|³) e `NewtonianGravity`, que soma as
acelerações de várias massas pontuais fixas (planetas e luas), permitindo estudar a
influência de mais de um corpo.

**Hipóteses:** massas pontuais ou esferas com simetria esférica de massa, com o ponto fora
das esferas (teorema das cascas; para dentro de um corpo a fórmula não vale); mecânica
newtoniana, sem correções relativísticas; fontes estáticas; `g = μ/r²` não inclui o efeito
centrífugo da rotação do corpo (por isso μ/R² na Terra ≈ 9,798 m/s² difere de g0).

**Domínio:** massas e μ ≥ 0 (resultado nulo para massa nula), r > 0. As operações são
intercaladas (`μ/r/r` e `G·(m1/r)·(m2/r)`), o que evita under/overflow espúrio de `r²` e de
`G·m1·m2`; resultados fora do intervalo do `float64` — maiores que ≈ 1,8 × 10³⁰⁸ ou, com
massas não nulas, menores que o menor `float64` — geram erro em vez de `inf`, `NaN` ou 0.

**Extensão prevista:** qualquer modelo com o método `acceleration(position, velocity, time)`
segue o protocolo `GravityModel`. Correções relativísticas (dependentes da velocidade) e
fontes móveis vindas de efemérides poderão ser incluídas sem alterar o integrador da etapa 2.

## 6. Órbitas — problema de dois corpos (`physics/orbits.py`)

| Chave | Expressão | Unidade | Domínio |
|---|---|---|---|
| `velocidade_orbital_circular` | v = √(μ/r) | m/s | μ, r > 0 |
| `periodo_orbital_circular` | T = 2π·√(r³/μ) | s | μ, r > 0 |
| `periodo_orbital_eliptico` | T = 2π·√(a³/μ) | s | μ, a > 0 |
| `raio_orbita_circular` | r = (μ·T²/(4π²))^(1/3) | m | μ, T > 0 |
| `energia_especifica` | ε = v²/2 − μ/r | J/kg | μ, r > 0; v ≥ 0 |
| `velocidade_vis_viva` | v = √(μ·(2/r − 1/a)) | m/s | μ, r > 0; a ≠ 0; r ≤ 2a se a > 0 |
| `velocidade_escape` | v_esc = √(2μ/r) | m/s | μ, r > 0 |
| `semieixo_maior_por_energia` | a = −μ/(2ε) | m | μ > 0; ε ≠ 0 |
| `raio_periapse` | r_p = a·(1 − e) | m | a > 0, 0 ≤ e < 1 |
| `raio_apoapse` | r_a = a·(1 + e) | m | a > 0, 0 ≤ e < 1 |

Nas fórmulas orbitais exige-se μ > 0 (corpo central atrator); nas de gravitação, μ ≥ 0.

`orbit_from_state(mu, r⃗, v⃗)` (somente API) calcula, a partir do vetor de estado: energia ε,
momento angular h⃗ = r⃗ × v⃗, vetor excentricidade e⃗ = ((v² − μ/r)·r⃗ − (r⃗·v⃗)·v⃗)/μ, semieixo
maior a = −μ/(2ε), semilatus rectum p = h²/μ, periapse p/(1+e), apoapse a·(1+e), período e
inclinação i = atan2(√(h_x² + h_y²), h_z). Classificação, com tolerância adimensional `tol`
(padrão 10⁻¹⁰; valores acima de 10⁻³ são rejeitados):

- **radial** (degenerada) se |h| ≤ tol·|r|·|v| — com energia desprezível, a = ∞;
- **circular** se e < tol;
- **parabólica** se |e − 1| ≤ tol **e** |ε| ≤ tol·μ/r;
- caso contrário, **elíptica** se ε < 0 e **hiperbólica** se ε > 0.

O teste de energia é necessário porque trajetórias quase radiais têm e ≈ 1 qualquer que seja
a energia; sem ele, um lançamento quase vertical e ligado seria rotulado parabólico, com
semieixo e período infinitos.

**Hipóteses:** problema de dois corpos (massa do veículo desprezível; corpo central em
repouso em referencial inercial); perturbações desprezadas (arrasto atmosférico,
achatamento J2, pressão de radiação solar, terceiros corpos); mecânica newtoniana; μ sem
incerteza propagada.

**Massa do segundo corpo.** As fórmulas usam o μ do corpo central, supondo a massa do
veículo desprezível. Quando a massa do segundo corpo não é desprezível para a precisão
desejada — Lua–Terra (efeito de ≈ 0,6 % no período) ou Terra–Sol (≈ 1,5 × 10⁻⁶) —, a órbita
**relativa** usa μ = G·(M₁ + M₂), que pode ser informado diretamente:
`aurora calcular periodo_orbital_eliptico mu=403503.2418km^3/s^2 a=384400km --em d`
(μ_Terra + μ_Lua), que dá 27,28 d. Com o μ da Terra apenas, o resultado é ≈ 0,6 % maior
(27,45 d). Nenhum dos dois reproduz exatamente o mês sideral (27,32 d): o valor de `a` e a
perturbação do Sol, desprezada pelo modelo, também contam.

**Cuidados numéricos:** os períodos são calculados como `2π·a·√(a/μ)` para não formar `a³`
(evita overflow); para órbitas quase circulares a excentricidade calculada tem erro absoluto
da ordem de 10⁻¹⁵ (alguns ulp)·(v²r/μ), por isso a classificação usa tolerância; a apoapse
usa a·(1+e), bem condicionada mesmo com e → 1 (p/(1−e) perderia precisão); as normas usam
`math.hypot`, que não sofre overflow/underflow ao elevar ao quadrado.

## 7. Propagação analítica de dois corpos (`physics/kepler.py`)

`propagate_kepler(mu, r⃗₀, v⃗₀, Δt)` devolve o estado após Δt (positivo ou negativo) pela
formulação de **variável universal**, válida para elipses, parábolas, hipérboles e
trajetórias radiais (Curtis, alg. 3.3 e 3.4; Vallado, alg. 8):

```text
α = 2/r₀ − v₀²/μ                          (inverso do semieixo maior; zero na parábola)
F(χ) = (r₀·v_r0/√μ)·χ²·C(z) + (1 − α·r₀)·χ³·S(z) + r₀·χ − √μ·Δt = 0,   z = α·χ²
f = 1 − χ²·C(z)/r₀        g = Δt − χ³·S(z)/√μ
ḟ = √μ·(α·χ³·S(z) − χ)/(r·r₀)        ġ = 1 − χ²·C(z)/r
r⃗ = f·r⃗₀ + g·v⃗₀          v⃗ = ḟ·r⃗₀ + ġ·v⃗₀
```

C(z) e S(z) são as funções de Stumpff. **Cuidados numéricos:**

- Para |z| < 0,1, C e S são calculadas pela série de Taylor (12 termos); a forma fechada
  `(1 − cos √z)/z` perderia algarismos por cancelamento perto de z = 0.
- F′(χ) = r > 0, então F é crescente e a raiz é única. Ela é cercada por um intervalo
  [lo, hi] com mudança de sinal (expandido por duplicação) e refinada por Newton protegido
  (estratégia *rtsafe*, Press et al., §9.4): quando o passo de Newton sairia do intervalo ou
  não o reduziria o bastante — o que acontece em hipérboles, onde F cresce exponencialmente
  —, usa-se bisseção.
- O critério de parada compara |F| com uma estimativa do erro de arredondamento ao avaliá-la
  (8ε·Σ\|termos\|), em vez de uma tolerância fixa que poderia ser inatingível.
- Em órbitas elípticas, Δt é reduzido pelo número inteiro de períodos (`fmod`), o que mantém
  χ limitado mesmo em propagações de milhares de períodos.
- Estados que não cabem em `float64` (hipérboles propagadas por tempos enormes) geram
  `DomainError`, nunca `inf` ou `NaN`.

`state_from_elements(mu, a, e, i, Ω, ω, ν)` converte elementos clássicos em vetor de estado
(referencial perifocal e rotação 3-1-3). Exige a > 0 com 0 ≤ e < 1 ou a < 0 com e > 1; a
parábola (e = 1) não é definida pelo semieixo maior e é recusada, assim como uma anomalia
verdadeira fora do ramo da hipérbole (1 + e·cos ν ≤ 0).

**Validação** (detalhes e números em [VALIDACAO_ETAPA2.md](VALIDACAO_ETAPA2.md), seção 1):
contra a equação de Kepler clássica em 1 000 elipses aleatórias (maior erro medido:
3,5 × 10⁻¹³·a), contra a equação hiperbólica em 1 000 hipérboles (1,8 × 10⁻¹³ relativo),
contra a equação de Barker na parábola, ida e volta, 10 000 períodos e conservação da energia
e do momento angular.

## 8. Manobras impulsivas (`physics/maneuvers.py`)

| Chave | Expressão | Unidade | Domínio |
|---|---|---|---|
| `delta_v_tsiolkovsky` | Δv = Isp·g0·ln(m0/mf) | m/s | Isp, m0, mf > 0; mf ≤ m0 |
| `massa_propelente` | m_p = m0·(1 − exp(−Δv/(Isp·g0))) | kg | Isp, m0 > 0; Δv ≥ 0 |
| `delta_v_hohmann` | Δv = \|Δv1\| + \|Δv2\| | m/s | μ, r1, r2 > 0 |
| `tempo_hohmann` | t = π·√(((r1 + r2)/2)³/μ) | s | μ, r1, r2 > 0 |

Transferência de Hohmann entre órbitas circulares coplanares de raios r1 e r2 (a partir do
centro): elipse de transferência com a_t = (r1 + r2)/2;
Δv1 = √(μ(2/r1 − 1/a_t)) − √(μ/r1) e Δv2 = √(μ/r2) − √(μ(2/r2 − 1/a_t)) (vis-viva).
Valores negativos indicam queimas **retrógradas** (descida para uma órbita menor); o custo
total usa os módulos. `hohmann_transfer` devolve também as velocidades antes e depois de cada
queima.

**Cuidados numéricos:** `massa_propelente` usa `expm1`, preciso também para Δv muito pequeno
(1 − exp(−x) perderia todos os algarismos para x ≈ 10⁻¹²); o período usa `a·√(a/μ)`.

**Hipóteses:** queimas **impulsivas** (variação instantânea de velocidade, sem perdas
gravitacionais nem de arrasto); impulso específico constante; g0 = 9,80665 m/s² é o valor
convencional usado na definição de Isp, não a gravidade local; sem forças externas durante a
queima; problema de dois corpos. Não são planos de voo: uma queima real de quilômetros por
segundo dura minutos e tem perdas que o modelo ignora.

## 9. Propagação numérica (`simulation/propagator.py`)

Equações de movimento de uma partícula de teste em referencial inercial:

```text
dr⃗/dt = v⃗          dv⃗/dt = a⃗(r⃗, v⃗, t)
```

A aceleração vem de um modelo que segue o protocolo `GravityModel`. Na etapa 2 o modelo é
`TwoBodyGravity`: a⃗ = −μ·r⃗/r³, com o corpo central fixo na origem (calculado como
`(μ/r/r)/r` para evitar overflow de r³). `NewtonianGravity` (fontes fixas) também pode ser
integrado.

**Integrador.** `scipy.integrate.solve_ivp` com o método de Dormand–Prince de ordem 8
(`DOP853`, padrão; Hairer, Nørsett e Wanner) ou de ordem 5 (`RK45`; Dormand e Prince, 1980).
O passo é escolhido automaticamente para que o erro **local** estimado de cada componente do
estado fique abaixo de `atol_i + rtol·|y_i|`. Padrões:

| Parâmetro | Padrão | Justificativa |
|---|---|---|
| `rtol` | 10⁻¹⁰ | Seis ordens acima do épsilon de máquina: deixa margem para o arredondamento acumulado. Aceito entre 10⁻¹³ e 10⁻³. |
| `atol` de posição | 10⁻⁴ m | Abaixo de `rtol·r` em qualquer órbita da Terra ou da Lua (≈ 7 × 10⁻⁴ m em LEO): a tolerância relativa domina e `atol` só atua quando uma componente passa por zero. |
| `atol` de velocidade | 10⁻⁷ m/s | Idem (`rtol·v` ≈ 7,7 × 10⁻⁷ m/s em LEO). |
| passo máximo | sem limite | O controle de erro já reduz o passo perto do periapse. |
| limite de avaliações | 2 000 000 por trecho | Cerca de 30 s num computador de mesa; uma simulação mais longa termina com erro claro (`SimulationError`) em vez de travar. Ajustável (`--max-avaliacoes`). |

As amostras de saída (`--pontos`, igualmente espaçadas no tempo) vêm da interpolação de alta
ordem do próprio método (saída densa), sem forçar passos pequenos.

**Tipos de erro numérico:**

1. **Erro de truncamento local** — controlado pelas tolerâncias a cada passo.
2. **Erro global** — o acumulado dos erros locais ao longo da simulação. **Não** é limitado
   pela tolerância: em órbitas fechadas, um pequeno erro de energia vira um pequeno erro de
   período, e o atraso de fase se acumula a cada volta. Medimos crescimento aproximadamente
   quadrático com o número de órbitas, quase todo ao longo da trajetória (LEO, tolerâncias
   padrão: 7,6 mm em 10 órbitas, 1,5 m em 100, 169 m em 1 000).
3. **Arredondamento** — `float64` tem ε ≈ 2,2 × 10⁻¹⁶; com tolerâncias muito pequenas, o
   arredondamento acumulado passa a dominar e reduzir `rtol` deixa de ajudar (medido: de
   10⁻¹⁰ para 10⁻¹² o erro cai apenas 2,5 a 31 vezes). Por isso `rtol < 10⁻¹³` é recusado.
4. **Dependência da orientação** — o controle é feito por componente cartesiana, então a
   mesma órbita, girada, tem erro um pouco diferente (medido: fator de até ≈ 6).

**Como o erro é medido.** Para o modelo de dois corpos, `aurora.simulation.validation`
compara cada amostra com a solução analítica (seção 7) a partir do mesmo estado inicial e
calcula as variações relativas da energia e do momento angular. A CLI **sempre** exibe esse
erro medido junto com a simulação; `aurora validar` executa seis cenários de referência. Os
números estão em [VALIDACAO_ETAPA2.md](VALIDACAO_ETAPA2.md).

**Colisão.** Com `collision_radius` (a CLI usa o raio de referência do corpo), um evento
terminal interrompe a integração quando r cai a esse valor; o instante é localizado por busca
de raiz na saída densa (testado: |r| − R ≤ 1 mm) e acrescentado como última amostra. Sem
raio de colisão, uma queda até r = 0 termina com `SimulationError` (singularidade do modelo).

**Hipótese numérica**, exibida com cada simulação: "Integração numérica com passo
adaptativo: o erro local é limitado pelas tolerâncias do integrador; o erro global acumula
com o tempo de simulação."

## 10. Condições iniciais, veículo, missões e telemetria (`simulation/`)

**Referencial.** Inercial, com origem no centro do corpo central e o plano xy no equador do
corpo. Não há data/época associada: os resultados não estão ligados a nenhum instante real, e
o eixo x não aponta para nenhuma direção astronômica específica.

**Condições iniciais idealizadas** (`initial_conditions.py`). O "lançamento" é representado
apenas pelo estado no fim da queima do lançador (*burnout*): altitude, velocidade inercial,
ângulo de trajetória γ (0 = horizontal; positivo para cima), inclinação, longitude do nodo
ascendente e argumento de latitude. Não são modelados subida atmosférica, empuxo, arrasto,
rotação do corpo (a velocidade do local de lançamento não é somada) nem o local de
lançamento. Atalhos: `circular_orbit_state` (velocidade circular, γ = 0) e
`elliptic_orbit_state` (altitudes de periapse e apoapse; parte do periapse).

**Veículo virtual** (`vehicle.py`). Massa seca, massa de propelente e impulso específico;
imutável. Cada queima consome m_p = m·(1 − exp(−|Δv|/(Isp·g0))) (Tsiolkovsky), e uma queima
sem propelente suficiente gera erro. `propellant_required(m_seca, Isp, Δv)` dá o propelente
para somar Δv; para várias queimas, Δv é a soma dos módulos, porque as razões de massa se
multiplicam.

**Missão de Hohmann** (`missions.py`). Parte de (r1, 0, 0) em órbita circular equatorial,
aplica a 1ª queima na direção da velocidade, integra numericamente meia elipse de
transferência, aplica a 2ª queima no raio alcançado e segue na órbita final (uma volta, por
padrão). A diferença entre o raio alcançado e r2, e entre o semieixo da órbita final e r2,
mede o erro numérico. Sem veículo informado, usa 1 000 kg secos, Isp = 300 s e exatamente o
propelente necessário.

**Telemetria simulada** (`telemetry.py`). Tabela com t, posição, velocidade, raio, altitude,
módulo da velocidade, energia específica e momento angular, em SI, exportada em CSV
(metadados em linhas `# chave: valor` e números com todos os algarismos de `float64`) ou JSON.
Todo arquivo registra `origem: simulacao`, o aviso de que não é medição, o corpo e suas
fontes, o modelo, o integrador com tolerâncias, o número de avaliações e as hipóteses — o
suficiente para reproduzir o experimento.

## 11. Validação

Os testes em `tests/` verificam quatro tipos de propriedade:

1. **Fórmulas diretas** — o resultado coincide com a expressão da especificação.
2. **Identidades analíticas independentes** — por exemplo: v·T = 2πr na órbita circular;
   ε = −μ/(2r) na órbita circular; vis-viva com a = r reproduz a velocidade circular;
   v_esc = √2·v_c; T² ∝ a³; r_p + r_a = 2a; conservação do momento angular
   (r_p·v_p = r_a·v_a) e da energia entre periapse e apoapse; raio por período é a inversa
   do período; a solução vetorial 3D coincide com a escalar por componente; relações de
   Girard–Vieta em 2 000 equações aleatórias; ida e volta elementos → estado → elementos
   (fórmulas perifocais independentes) em 500 órbitas com orientação, anomalia e
   excentricidade arbitrárias, inclusive e → 1 e hipérboles.
3. **Casos de referência conhecidos** — raio geoestacionário ≈ 42 164 km para período igual
   ao dia sideral; período de órbita circular de 1 ua em torno do Sol ≈ 365,256 dias
   (concordância de 1,5 × 10⁻⁶ com o ano sideral de 365,256363 d, explicada quase toda por
   se desprezar a massa Terra–Lua); exemplo de Kahan para o discriminante.
4. **Precisão numérica** — raízes do 2º grau com erro ≤ 4 ulp contra referência de 60
   dígitos em escalas normal, próxima de 10³⁰⁸ e subnormal; propagação de Kepler contra as
   equações de Kepler clássica, hiperbólica e de Barker; erro do integrador numérico contra
   a solução analítica, sua redução com a tolerância e a conservação de energia e momento
   angular; consumo de propelente contra Tsiolkovsky; Hohmann contra formas fechadas.

Além disso, os testes cobrem entradas inválidas (NaN, ±∞, texto, valores lógicos, unidades
incompatíveis, separador de milhar), domínios (a = 0, r > 2a, ε = 0, e ≥ 1, órbita abaixo da
superfície), entradas extremas que não podem gerar traceback na interface e as regressões
dos problemas encontrados na auditoria ([AUDITORIA_ETAPA1.md](AUDITORIA_ETAPA1.md)).

## 12. Limitações conhecidas

- Incertezas das constantes não são propagadas aos resultados.
- Órbitas usam o μ do corpo central (veículo de massa desprezível); ver "Massa do segundo
  corpo" (seção 6).
- Corpos modelados como esferas; não há achatamento, rotação, atmosfera ou marés.
- Fontes gravitacionais em `NewtonianGravity` são fixas (sem movimento próprio).
- A simulação usa apenas o modelo de dois corpos: sem achatamento (J2), arrasto, terceiros
  corpos (Lua e Sol sobre satélites da Terra), pressão de radiação nem efeitos relativísticos.
  Para esse modelo a solução de Kepler é exata; o integrador numérico existe para receber
  modelos com perturbações nas próximas etapas, que exigirão nova validação.
- O erro global do integrador cresce com o tempo simulado (aproximadamente com o quadrado do
  número de órbitas); em simulações longas, reduza `rtol` e confira o erro exibido.
- Manobras são impulsivas; não há queimas de duração finita, controle de atitude, estágios
  de lançador nem subida atmosférica.
- O referencial não tem época: não é possível comparar uma simulação com a posição de um
  satélite real em uma data.
- A tabela de corpos tem apenas Sol, Terra e Lua; outros corpos virão de fontes confiáveis
  importadas na etapa 4.

## 13. Referências

- Tiesinga, E. et al. *CODATA recommended values of the fundamental physical constants: 2018*.
  Rev. Mod. Phys. 93, 025010 (2021).
- IAU 2012, Resolução B2 — redefinição da unidade astronômica.
- Prša, A. et al. *Nominal values for selected solar and planetary quantities: IAU 2015
  Resolution B3*. Astron. J. 152, 41 (2016).
- National Geospatial-Intelligence Agency. *World Geodetic System 1984* (TR8350.2).
- Folkner, W. M. et al. *The Planetary and Lunar Ephemerides DE430 and DE431*. IPN Progress
  Report 42-196 (2014).
- Archinal, B. A. et al. *Report of the IAU Working Group on Cartographic Coordinates and
  Rotational Elements: 2009*. Celest. Mech. Dyn. Astron. 109, 101–135 (2011).
- Curtis, H. D. *Orbital Mechanics for Engineering Students*. Elsevier.
- Vallado, D. A. *Fundamentals of Astrodynamics and Applications*. Microcosm Press.
- Press, W. H. et al. *Numerical Recipes*, seções 5.6 (equações quadráticas e cúbicas) e 9.4
  (Newton protegido por bisseção, *rtsafe*). Cambridge University Press.
- Dormand, J. R.; Prince, P. J. *A family of embedded Runge-Kutta formulae*. J. Comput. Appl.
  Math. 6, 19–26 (1980).
- Hairer, E.; Nørsett, S. P.; Wanner, G. *Solving Ordinary Differential Equations I: Nonstiff
  Problems*, 2ª ed. Springer (1993) — método DOP853.
- Virtanen, P. et al. *SciPy 1.0: fundamental algorithms for scientific computing in
  Python*. Nature Methods 17, 261–272 (2020).
