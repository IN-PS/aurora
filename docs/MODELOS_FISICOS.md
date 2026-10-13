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
  partir do raio de referência do corpo (modelado como esfera).
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
  1. Δ é calculado **exatamente** em aritmética racional (todo `float` é um racional exato)
     e arredondado uma única vez. A classificação é, portanto, exata para os coeficientes
     fornecidos. Exemplo de teste: para `x² − 202687298·x + (101343649² − 1)`, a conta
     ingênua `b·b − 4·a·c` dá 0 (falsa raiz dupla); o valor correto é Δ = 4.
  2. As raízes usam a forma estável `q = −(b + sinal(b)·√Δ)/2`, `x₁ = q/a`, `x₂ = c/q`, que
     evita o cancelamento catastrófico. Exemplo: para `x² + 10⁸x + 1`, a fórmula ingênua
     erra a raiz pequena em ~25 %; a forma estável acerta até o último algarismo.
  3. `√Δ` é calculada reescalando Δ por potências de 4, evitando overflow intermediário
     com coeficientes como 10²⁰⁰. Se uma raiz não couber em `float64`, é levantado erro.
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

Também disponíveis na API: `constant_acceleration_state` (as mesmas relações em 3D,
vetoriais) e `time_to_reach_position`, que usa a calculadora do 2º grau para obter os
instantes t ≥ 0 em que o móvel passa por uma posição.

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

**Hipóteses:** massas pontuais ou esferas com simetria esférica de massa; mecânica
newtoniana, sem correções relativísticas; fontes estáticas; `g = μ/r²` não inclui o efeito
centrífugo da rotação do corpo (por isso μ/R² na Terra ≈ 9,798 m/s² difere de g0).

**Extensão prevista:** qualquer modelo com o método `acceleration(position, velocity, time)`
segue o protocolo `GravityModel`. Correções relativísticas (dependentes da velocidade) e
fontes móveis vindas de efemérides poderão ser incluídas sem alterar o integrador da etapa 2.

## 6. Órbitas — problema de dois corpos (`physics/orbits.py`)

| Chave | Expressão | Unidade | Domínio |
|---|---|---|---|
| `velocidade_orbital_circular` | v = √(μ/r) | m/s | μ, r > 0 |
| `periodo_orbital_circular` | T = 2π·√(r³/μ) | s | μ, r > 0 |
| `periodo_orbital_eliptico` | T = 2π·√(a³/μ) | s | a > 0 |
| `raio_orbita_circular` | r = (μ·T²/(4π²))^(1/3) | m | T > 0 |
| `energia_especifica` | ε = v²/2 − μ/r | J/kg | v ≥ 0 |
| `velocidade_vis_viva` | v = √(μ·(2/r − 1/a)) | m/s | a ≠ 0; r ≤ 2a se a > 0 |
| `velocidade_escape` | v_esc = √(2μ/r) | m/s | |
| `semieixo_maior_por_energia` | a = −μ/(2ε) | m | ε ≠ 0 |
| `raio_periapse` | r_p = a·(1 − e) | m | a > 0, 0 ≤ e < 1 |
| `raio_apoapse` | r_a = a·(1 + e) | m | a > 0, 0 ≤ e < 1 |

`orbit_from_state(mu, r⃗, v⃗)` calcula, a partir do vetor de estado: energia ε, momento
angular h⃗ = r⃗ × v⃗, vetor excentricidade e⃗ = ((v² − μ/r)·r⃗ − (r⃗·v⃗)·v⃗)/μ, semieixo maior,
semilatus rectum p = h²/μ, periapse p/(1+e), apoapse p/(1−e), período e inclinação, e
classifica a trajetória em circular, elíptica, parabólica, hiperbólica ou radial
(degenerada), com tolerância configurável (padrão 10⁻¹⁰ na excentricidade).

**Hipóteses:** problema de dois corpos (massa do veículo desprezível; corpo central em
repouso em referencial inercial); perturbações desprezadas (arrasto atmosférico,
achatamento J2, pressão de radiação solar, terceiros corpos); mecânica newtoniana; μ sem
incerteza propagada.

**Cuidados numéricos:** os períodos são calculados como `2π·a·√(a/μ)` para não formar `a³`
(evita overflow); para órbitas quase circulares a excentricidade calculada tem erro absoluto
da ordem de 10⁻¹⁶·(v²r/μ), por isso a classificação usa tolerância.

## 7. Validação

Os testes em `tests/` verificam três tipos de propriedade:

1. **Fórmulas diretas** — o resultado coincide com a expressão da especificação.
2. **Identidades analíticas independentes** — por exemplo: v·T = 2πr na órbita circular;
   ε = −μ/(2r) na órbita circular; vis-viva com a = r reproduz a velocidade circular;
   v_esc = √2·v_c; T² ∝ a³; r_p + r_a = 2a; conservação do momento angular
   (r_p·v_p = r_a·v_a) e da energia entre periapse e apoapse; raio por período é a inversa
   do período; a solução vetorial 3D coincide com a escalar por componente; relações de
   Girard–Vieta em 2 000 equações aleatórias.
3. **Casos de referência conhecidos** — raio geoestacionário ≈ 42 164 km para período igual
   ao dia sideral; período de órbita circular de 1 ua em torno do Sol ≈ 365,256 dias
   (concordância de ~10⁻⁶, coerente com as simplificações do modelo); exemplo de Kahan para
   o discriminante.

Além disso, os testes cobrem entradas inválidas (NaN, ±∞, texto, valores lógicos, unidades
incompatíveis) e domínios (a = 0, r > 2a, ε = 0, e ≥ 1, órbita abaixo da superfície).

## 8. Limitações conhecidas

- Incertezas das constantes não são propagadas aos resultados.
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
