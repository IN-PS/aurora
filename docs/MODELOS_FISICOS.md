# Modelos físicos e matemáticos — etapa 1

Este documento lista cada modelo implementado, com unidades, domínio de validade, hipóteses,
cuidados numéricos e a forma como foi validado. As mesmas hipóteses são exibidas pela
interface (`aurora formulas <chave>`).

> Todos os resultados são **cálculos analíticos com modelos idealizados**. Não são medições
> e não devem ser usados como planos de voo.

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
     dão raízes com erro de no máximo ~2 ulp (verificado contra referência de 60 dígitos).
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
ligeiramente acima ou abaixo do valor exato. Por isso, se `|Δ| ≤ 8ε·max(v0², 2|a|·max(|x0|, |x|))`
(ε = 2,2 × 10⁻¹⁶), a posição é considerada atingida uma única vez, no instante `t = −v0/a`; a
equação de Torricelli trata da mesma forma um radicando negativo apenas pelo arredondamento
(resultado `|v| = 0`). Fora dessa faixa, de largura relativa ~10⁻¹⁵, nada muda.

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

**Domínio:** massas e μ ≥ 0 (resultado nulo para massa nula), r > 0. As divisões são feitas
em sequência (`μ/r/r`), o que evita underflow de `r²` para distâncias minúsculas; resultados
fora do intervalo do `float64` geram erro em vez de `inf`/`NaN`.

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

**Massas comparáveis.** As fórmulas usam o μ do corpo central, supondo a massa do veículo
desprezível. Para a órbita **relativa** de dois corpos de massas comparáveis (Lua em torno
da Terra, Terra em torno do Sol), use μ = G·(M₁ + M₂), informando-o diretamente:
`aurora calcular periodo_orbital_eliptico mu=403503.2418km^3/s^2 a=384400km --em d`
(μ_Terra + μ_Lua), que dá 27,28 d. Com o μ da Terra apenas, o resultado é ≈ 0,6 % maior
(27,45 d). Nenhum dos dois reproduz exatamente o mês sideral (27,32 d): o valor de `a` e a
perturbação do Sol, desprezada pelo modelo, também contam.

**Cuidados numéricos:** os períodos são calculados como `2π·a·√(a/μ)` para não formar `a³`
(evita overflow); para órbitas quase circulares a excentricidade calculada tem erro absoluto
da ordem de 10⁻¹⁵ (alguns ulp)·(v²r/μ), por isso a classificação usa tolerância; a apoapse
usa a·(1+e), bem condicionada mesmo com e → 1 (p/(1−e) perderia precisão); as normas usam
`math.hypot`, que não sofre overflow/underflow ao elevar ao quadrado.

## 7. Validação

Os testes em `tests/` verificam três tipos de propriedade:

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
   dígitos em escalas normal, próxima de 10³⁰⁸ e subnormal.

Além disso, os testes cobrem entradas inválidas (NaN, ±∞, texto, valores lógicos, unidades
incompatíveis, separador de milhar), domínios (a = 0, r > 2a, ε = 0, e ≥ 1, órbita abaixo da
superfície), entradas extremas que não podem gerar traceback na interface e as regressões
dos problemas encontrados na auditoria ([AUDITORIA_ETAPA1.md](AUDITORIA_ETAPA1.md)).

## 8. Limitações conhecidas

- Incertezas das constantes não são propagadas aos resultados.
- Órbitas usam o μ do corpo central (veículo de massa desprezível); ver "Massas comparáveis".
- Corpos modelados como esferas; não há achatamento, rotação, atmosfera ou marés.
- Fontes gravitacionais em `NewtonianGravity` são fixas (sem movimento próprio).
- Não há ainda integração numérica de trajetórias (etapa 2) — somente relações analíticas.
- A tabela de corpos tem apenas Sol, Terra e Lua; outros corpos virão de fontes confiáveis
  importadas na etapa 4.

## 9. Referências

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
- Press, W. H. et al. *Numerical Recipes*, seção 5.6 (equações quadráticas e cúbicas).
  Cambridge University Press.
