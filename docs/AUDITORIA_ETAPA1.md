# Auditoria técnica — AURORA v1.0, etapa 1

- **Data:** 13/10/2026
- **Código auditado:** commit `06219b9` (pacote 0.1.0, etapa 1 — fundamentos matemáticos)
- **Natureza:** auditoria de software e de modelos físicos. Todos os números abaixo foram
  obtidos executando o código; nenhum é medição de fenômeno real.

## 1. Escopo e método

Cinco auditores independentes examinaram o projeto, cada um com um enfoque: (1) fórmulas
físicas, unidades SI, constantes e condições iniciais; (2) comparação dos cálculos orbitais com
resultados analíticos conhecidos e com uma implementação independente; (3) robustez numérica
e validação de entradas; (4) interface de linha de comando, códigos de saída, codificação,
empacotamento e compatibilidade de versões do Python; (5) arquitetura real, exatidão da
documentação e inventário de protótipos.

Cada achado precisava vir acompanhado do comando executado e da saída observada. Em seguida,
um **verificador adversarial** por enfoque tentou refutar cada achado, reproduzindo-o de forma
independente e reavaliando a gravidade. Resultado: **62 achados, nenhum refutado** (8
confirmados apenas em parte e alguns com gravidade reduzida). Achados repetidos entre
enfoques foram agrupados em **43 problemas distintos** (seção 5).

Também foram verificados: compilação do pacote (wheel), instalação limpa em ambiente virtual,
execução com Python 3.10, 3.11, 3.12, 3.13 e 3.14, e com NumPy 1.24.4 (mínimo declarado) até
2.5.4.

## 2. Arquitetura real

O grafo de importações (40 arestas internas) confirmou a arquitetura documentada, sem ciclos
e sem dependências "para cima":

```text
aurora.errors  ←  aurora.data  ←  aurora.core  ←  aurora.physics  ←  aurora.interface (CLI)
```

- `aurora.errors`: exceções; `aurora.data`: constantes e corpos com fontes;
  `aurora.core`: validação, unidades, registro de fórmulas, procedência;
  `aurora.physics`: modelos (funções puras em SI); `aurora.interface.cli`: única interface.
- O registro de fórmulas é preenchido como efeito da importação de `aurora.physics`.
- Não há código de hardware, GUI, servidor web ou aplicativo Android no repositório.

## 3. Testes

| Momento | Resultado |
|---|---|
| Antes das correções (commit `06219b9`) | **143 passaram, 0 falharam** em Python 3.10.20, 3.11.17, 3.12.3, 3.13.16 e 3.14.6; também com NumPy 1.24.4 |
| Depois das correções | **239 passaram, 0 falharam** (143 originais + 96 novos); números por versão na seção 10 |
| Testes novos executados contra o código original | **69 dos 96 falham** — confirmam que reproduzem os defeitos corrigidos |

Os 27 testes novos que já passavam no código original cobrem comportamentos que estavam
corretos e passam a ser protegidos (inclusive contra as regressões encontradas na revisão
da seção 9), ou defeitos que só aparecem no Python 3.10–3.12 (ordem das opções da CLI),
agora exercitados pela integração contínua.

A suíte original passava inteira apesar dos defeitos porque testava apenas casos centrais
(coeficientes até 10²⁰⁰, estados orbitais alinhados aos eixos, opções no fim da linha de
comando, saída sempre UTF-8).

## 4. Verificação física e comparação com resultados analíticos

As 18 fórmulas registradas foram conferidas por análise dimensional (todas coerentes em SI)
e comparadas com expressões escritas à mão em 2 000 entradas aleatórias cada (erro relativo
máximo 9,5 × 10⁻¹⁶). As constantes e os dados dos corpos coincidem com as fontes citadas
(CODATA 2018, IAU 2012/2015, WGS 84, JPL DE430, IAU WGCCRE 2009).

Comparações com valores conhecidos, recalculadas com o código final (referência de alta
precisão em `decimal` com 50 dígitos):

| Caso | AURORA | Valor de referência | Erro numérico relativo |
|---|---|---|---|
| GEO: raio para T = dia sideral (86 164,0905 s) | 42 164,1696 km | ≈ 42 164 km | 1,2 × 10⁻¹⁵ |
| GPS: semieixo para T = meio dia sideral | 26 561,762 km | ≈ 26 562 km | 9,3 × 10⁻¹⁶ |
| LEO 400 km: velocidade circular | 7 668,558 m/s | ≈ 7,67 km/s | 6,4 × 10⁻¹⁷ |
| LEO 400 km: período | 92,5604 min | ≈ 92,6 min | 1,7 × 10⁻¹⁶ |
| Escape na superfície da Terra (raio equatorial) | 11,1799 km/s | 11,18 km/s | 4,6 × 10⁻¹⁷ |
| Escape na superfície da Lua | 2,3757 km/s | 2,38 km/s | 6,2 × 10⁻¹⁸ |
| Escape do Sol a 1 ua | 42,122 km/s | ≈ 42,1 km/s | 4,5 × 10⁻¹⁷ |
| Velocidade circular a 1 ua | 29,785 km/s | 29,78 km/s | — |
| Período circular a 1 ua | 365,256898 d | 365,256363 d (ano sideral) | desvio físico 1,5 × 10⁻⁶ (massa Terra–Lua desprezada) |
| Lua, a = 384 400 km, só μ_Terra | 27,4519 d | 27,32 d (mês sideral) | +0,48 % (aproximação de dois corpos) |
| Lua, a = 384 400 km, μ_Terra + μ_Lua | 27,2846 d | 27,32 d | −0,14 % (a e perturbação solar) |
| Hohmann 6 678 km → 42 164 km: Δv total | 3,8926 km/s | ≈ 3,9 km/s | — |
| Hohmann: tempo de transferência | 5,2750 h | ≈ 5,3 h | — |
| g = μ/R² no equador (sem rotação) | 9,79829 m/s² | ≈ 9,798 m/s² | — |

Outras verificações: 3ª lei de Kepler em 10 000 pares de órbitas (desvio ≤ 8,9 × 10⁻¹⁶);
conservação de energia e momento angular entre periapse e apoapse; ida e volta
elementos → estado → elementos em 5 000 órbitas com orientação arbitrária, incluindo e → 1 e
hipérboles (nenhuma classificação errada; erro em a ≤ 3,4 × 10⁻¹², em e ≤ 4,4 × 10⁻¹⁵, em i ≤
1,1 × 10⁻¹⁵ rad); invariância da classificação a mudanças de escala de 10⁻⁶⁰ a 10⁶⁰.

**Conclusão:** a física e as unidades estavam corretas. Os defeitos estavam no tratamento
numérico de casos-limite, na classificação de trajetórias quase radiais, na interface e na
documentação.

## 5. Problemas encontrados e correções

Gravidade após a verificação adversarial. Os identificadores remetem aos relatórios dos
auditores (PHYS = física/unidades, ORBI = órbitas, NUME = numérico, CLI = interface,
ARCH = arquitetura/documentação).

### 5.1 Defeitos e robustez

| # | Gravidade | Achados | Problema | Situação |
|---|---|---|---|---|
| 1 | Alta | ARCH-1 | O Δ "exato" era calculado sobre os `float` binários, não sobre os decimais digitados: `quadratica 1 0.6 0.09` = (x + 0,3)² era classificada como "complexas conjugadas". | **Corrigido.** A CLI lê os decimais como frações exatas; `solve_quadratic` usa o valor exato de `int`/`Fraction`. |
| 2 | Média | NUME-1, ARCH-2, CLI-2 | Coeficientes ≈ 10³⁰⁸ com raízes representáveis (ex.: 10³⁰⁸x² − 10³⁰⁸ = 0) causavam `OverflowError` e traceback. | **Corrigido.** q, q/a e c/q calculados em aritmética racional e arredondados uma vez. |
| 3 | Média | NUME-2 | Coeficientes subnormais davam raízes erradas ((−1, 1) em vez de (−1,618; 0,618)) ou erro espúrio. | **Corrigido** pela mesma mudança; erro máximo medido de 1,8 ulp em todas as escalas testadas. |
| 4 | Média | PHYS-1 | `time_to_reach_position` declarava a altura máxima inatingível em ~50 % dos casos (arredondamento de v0²/2g). | **Corrigido.** Tolerância igual ao erro de representação das entradas (≈ 1 ulp de x0 e x), aplicada só se o retorno não está no passado: 0/20 000 falhas e 0/20 000 falsos positivos com x0 até ±10⁶ m. |
| 5 | Média | PHYS-3 | Torricelli retornava \|v\| = 0 para móvel em repouso (v0 = a = 0) com Δx ≠ 0. | **Corrigido** (erro de domínio). |
| 6 | Média | ORBI-1 | Trajetórias quase radiais (ex.: lançamento vertical com 5 cm/s horizontais) eram rotuladas "parabólicas" com a = T = ∞, mesmo fortemente ligadas. | **Corrigido.** Parabólica exige também energia desprezível; senão, o sinal de ε decide. |
| 7 | Média | NUME-3, PHYS-2, ARCH-3, CLI-3 | Gravitação: r² sofria underflow (ZeroDivisionError, traceback) ou overflow (resultado 0); a versão vetorial retornava inf/NaN. | **Corrigido.** Operações intercaladas (μ/r/r, G·(m1/r)·(m2/r)), erro claro para resultado fora do `float64` (inclusive underflow com massas não nulas) e rede de segurança no registro de fórmulas. |
| 8 | Média | CLI-1 | `aurora calcular velocidade_escape --json mu=1 r=1` falhava no Python 3.10–3.12; opção entre dois valores falhava em todas as versões. | **Corrigido** (`parse_known_args`); verificado no 3.10. |
| 9 | Média | CLI-4 | `--json` trocava μ, ☉, ç… por "?" quando a saída não era UTF-8 (ex.: redirecionamento no Windows). | **Corrigido** (JSON em ASCII com escapes `\uXXXX`). |
| 10 | Média | NUME-4, ARCH-12 | "1.000 km" era lido silenciosamente como 1 km; "6.778.137 km" gerava "unidade desconhecida '.137 km'". | **Corrigido em parte.** Números com mais de um separador são rejeitados com mensagem clara; o separador único continua decimal (decisão de projeto, documentada no README; as entradas são ecoadas em SI). |
| 11 | Baixa | ORBI-2 | Apoapse calculada como p/(1−e) perdia precisão com e → 1 e quebrava r_p + r_a = 2a. | **Corrigido** (a·(1+e)). |
| 12 | Baixa | ORBI-3 | Inclinação por `acos` perdia ~10⁻⁸ rad perto de 0 e π (i = 10⁻⁹ rad saía 0). | **Corrigido** (`atan2`). |
| 13 | Baixa | ORBI-4, PHYS-4, NUME-5, NUME-6 | Funções fora do registro retornavam inf/NaN silenciosamente; mensagens citavam coeficientes internos ('c', "a = 0" para a subnormal). | **Corrigido** (verificações de finitude; Δx calculado exatamente). |
| 14 | Baixa | PHYS-5 | `orbit_from_state` aceitava qualquer `tol`; tol ≥ 1 rotulava tudo como radial. | **Corrigido** (0 < tol ≤ 10⁻³). |
| 15 | Baixa | PHYS-6 | `PointMass`: `==` devolvia um array e `hash()` falhava; vetores de `TwoBodyOrbit` eram graváveis. | **Corrigido.** |
| 16 | Baixa | PHYS-10, CLI-8 | Erro de `--raio`/altitude citava o parâmetro interno "radius". | **Corrigido.** |
| 17 | Baixa | NUME-7 | Leitura de grandezas com tempo quadrático em textos longos (5 s com 40 000 caracteres). | **Corrigido** (0,1 ms). |
| 18 | Baixa | NUME-8 | Saídas exibiam "-0 N". | **Corrigido.** |
| 19 | Baixa | NUME-9 | Conversões de unidade não validavam o valor; tipos errados geravam exceções Python genéricas. | **Corrigido** (escalares e arrays NumPy validados). |
| 20 | Baixa | CLI-5, ARCH-7 | Mensagens de erro de uso do `argparse` em inglês. | **Corrigido** (tabela de tradução, com retorno à mensagem original se não houver padrão). |
| 21 | Baixa | CLI-6 | `--help` e `--version` não eram reconhecidos. | **Corrigido** (aliases). Consequência: as abreviações `--ver`/`--vers` passaram a ser ambíguas (mensagem traduzida); as demais abreviações continuam aceitas. |
| 22 | Baixa | CLI-7 | O espaço inserido para proteger números negativos aparecia nas mensagens. | **Corrigido.** |
| 23 | Baixa | CLI-9 | `--em` exibia ∞ para resultado finito; sentinelas "inf"/"nan" do JSON não documentadas. | **Corrigido** (erro claro) e documentado. |
| 24 | Baixa | ARCH-11 | CLI exibia "Δ = ∞" para Δ finito que só excede o `float64`. | **Corrigido** ("Δ > 1,8 × 10³⁰⁸ …"). |
| 25 | Baixa | ARCH-4 | Recarregar um módulo de física (autoreload) levantava `ValueError`; mensagem de fórmula desconhecida voltada à CLI. | **Corrigido**; dependência de `import aurora.physics` documentada. |
| 26 | Info | NUME-10 | Tipos numéricos aceitos inconsistentes (`Decimal` rejeitado, `Fraction` aceito). | **Documentado.** |
| 27 | Info | NUME-11 | Torricelli com \|v0\| > 1,3 × 10¹⁵⁴ m/s gera erro espúrio de overflow. | **Não corrigido** — fora de qualquer regime físico; registrado aqui. |

### 5.2 Documentação e hipóteses

| # | Gravidade | Achados | Problema | Situação |
|---|---|---|---|---|
| 28 | Baixa | PHYS-7, ARCH-14 | Periapse/apoapse exibiam hipótese sobre μ (não usam μ); faltava a condição do teorema das cascas. | **Corrigido.** (A hipótese de referencial inercial na cinemática foi mantida.) |
| 29 | Baixa | PHYS-8, CLI-13, ARCH-10 | Campo de referências das fórmulas sempre vazio. | **Corrigido em parte:** as 10 fórmulas orbitais citam livros-texto (Curtis; Vallado); as demais não têm referência específica. |
| 30 | Baixa | PHYS-9 | Domínio de μ (≥ 0 na gravitação, > 0 nas órbitas) incompleto na documentação. | **Documentado.** |
| 31 | Baixa | ARCH-5 | Docs afirmavam que "todo resultado" carrega procedência — falso para a quadrática. | **Corrigido** (texto ajustado; JSON da quadrática com `origem`). |
| 32 | Baixa | ARCH-6 | README misturava fórmulas da CLI com funções só da API. | **Corrigido.** |
| 33 | Baixa | ARCH-8 | Docstrings descreviam as camadas incorretamente. | **Corrigido.** |
| 34 | Baixa | CLI-10 | Docstring de `main()` prometia retornar 2 em erro de uso (na verdade, `SystemExit`). | **Corrigido.** |
| 35 | Baixa | CLI-11, ARCH-9 | "AURORA v1.0" × pacote 0.1.0; versão definida em dois lugares. | **Corrigido** (título explicita a etapa; versão única em `aurora.__version__`). |
| 36 | Baixa | ARCH-15 | Comentário de teste citava o ramo errado da fórmula ingênua. | **Corrigido.** |
| 37 | Info | PHYS-12 | Erro de arredondamento da excentricidade subestimado em ~10× na documentação. | **Corrigido.** |
| 38 | Info | PHYS-11 | Altitude medida a partir do raio equatorial (até ≈ 21 km de diferença nos polos). | **Documentado.** |
| 39 | Info | ORBI-5 | μ do corpo central sem massa reduzida; falta de orientação para massas comparáveis. | **Documentado** (como usar μ = G(M₁+M₂)). |
| 40 | Info | ARCH-17 | Arredondamento decimal → binário não documentado. | **Corrigido** (item 1) e documentado. |
| 41 | Info | ARCH-16 | Ganchos de arquitetura sem consumidor (ver seção 7). | **Documentado**; `NewtonianGravity` passou a validar fontes e velocidade. |
| 42 | Info | CLI-14, ARCH-18 | Não há interface gráfica nem instruções para Android. | **Documentado** (README, seção "Uso no Android"); não verificado em aparelho. |

### 5.3 Testes

| # | Gravidade | Achados | Problema | Situação |
|---|---|---|---|---|
| 43 | Baixa | ORBI-6, NUME-12, CLI-12, ARCH-13 | Lacunas: estados orbitais genéricos, entradas extremas, posição das opções, codificação, versões antigas do Python. | **Corrigido:** 96 testes novos e integração contínua (Python 3.10–3.14, NumPy 1.24). Não foi adicionado teste que compare as tabelas Markdown com o registro (considerado frágil). |

## 6. Aproximações do modelo (todas documentadas; as que se aplicam à CLI são exibidas por ela)

- Problema de dois corpos com o corpo central em repouso e veículo de massa desprezível.
- Sem perturbações: arrasto, achatamento (J2), pressão de radiação, terceiros corpos.
- Corpos esféricos; altitude a partir do raio de referência (equatorial, para a Terra).
- g = μ/r² sem o efeito centrífugo da rotação.
- Mecânica clássica (sem correções relativísticas).
- Incertezas de G e μ não propagadas; μ da Lua arredondado a 5 algarismos.
- Fontes gravitacionais fixas em `NewtonianGravity` (somente API).
- Classificação de trajetórias com tolerância, padrão 10⁻¹⁰ (somente API).

## 7. Funcionalidades que ainda são protótipos ou não existem

- **Ganchos sem consumidor:** `GravityModel` (aguarda o integrador da etapa 2); parâmetros
  `velocity`/`time` de `NewtonianGravity`; `DataOrigin.SIMULATED` e `MEASURED` (nada os
  produz; "o simulador nunca cria `MEASURED`" é convenção, não imposição de código).
- **Somente API, sem comando na CLI:** `orbit_from_state`, `constant_acceleration_state`,
  `time_to_reach_position`, `NewtonianGravity`.
- **Não implementados:** propagação numérica de trajetórias, gráficos, painel gráfico,
  radar virtual, efemérides, missões interplanetárias, qualquer software de hardware
  (etapas 2 a 5).

## 8. Limitações remanescentes

- Separador decimal único: "1.000 km" continua significando 1 km (ver item 10).
- Torricelli com velocidades acima de ~10¹⁵⁴ m/s (item 27).
- Mensagens do `argparse` sem padrão na tabela de tradução continuam em inglês.
- Não houve execução em aparelho Android.

## 9. Revisão adversarial das correções

Antes da entrega, as próprias correções (diferença `06219b9..69f83dc`) passaram por uma
segunda revisão independente, também com verificação adversarial (3 revisores e 3
verificadores). Ela encontrou **16 problemas nas correções**, todos confirmados ao menos em
parte e todos tratados, com testes de regressão:

| Gravidade | Problema introduzido ou deixado pelas correções | Tratamento |
|---|---|---|
| Alta | Inteiros do NumPy em `solve_quadratic`: `Fraction(np.int64)` mantinha aritmética de 64 bits (AttributeError ou raízes erradas por estouro silencioso). | Numerador e denominador convertidos para `int` do Python. |
| Baixa | A tolerância do ponto de retorno crescia com \|x0\| (dependia da origem do eixo): falsos "atingido" e falsos "nunca atingido". | Faixa reduzida ao erro de representação das entradas e aplicada só se o retorno não está no passado e x ≠ x0. |
| Baixa | Decimais exatos sem limite: `quadratica 1e-100000000 1 1` levava minutos; coeficientes abaixo do `float64` eram resolvidos mas exibidos como 0. | Valores não normais usam o `float` exibido; não nulos abaixo do menor `float64` são rejeitados. |
| Baixa | Raízes −0.0 voltaram a aparecer. | Normalização restaurada. |
| Baixa | Erro quando só a raiz negativa (descartada) excedia o `float64`. | Raízes negativas descartadas antes da verificação. |
| Baixa | `G·m1·m2` ainda sofria underflow (0 N silencioso). | Operações intercaladas e erro em underflow. |
| Baixa | `to_si`/`from_si` deixaram de aceitar arrays NumPy. | Arrays aceitos e validados. |
| Baixa | Abreviações de opções (`--alt`, `--js`) desativadas sem aviso. | Reativadas. |
| Baixa | Sobras da linha de comando: `--` ou uma opção inválida faziam todos os NOME=VALOR serem rejeitados no Python 3.10–3.12. | Só as sobras com "-" são rejeitadas; `--` é ignorado. |
| Baixa | O espaço de proteção de negativos ainda aparecia em uma mensagem de `calcular`. | Corrigido. |
| Baixa | Documentação: limites de ulp inconsistentes; exemplo do README com fórmula inexistente; motivo errado para evitar o Termux da Google Play; "massas comparáveis" para Terra–Sol; título da seção 6 deste relatório. | Textos corrigidos. |

## 10. Como reproduzir

```bash
pip install -e ".[dev]"
python -m pytest -q          # suíte completa
```

Resultado da suíte final (239 testes) por interpretador, em Linux x86-64:

| Interpretador | NumPy | Resultado |
|---|---|---|
| Python 3.10.20 | 2.2.6 | 239 passaram, 0 falharam |
| Python 3.10.20 | 1.24.4 | 239 passaram, 0 falharam |
| Python 3.11.17 | 2.4.6 | 239 passaram, 0 falharam |
| Python 3.11.17 | 1.24.4 | 239 passaram, 0 falharam |
| Python 3.12.3 | 2.5.4 | 239 passaram, 0 falharam |
| Python 3.13.16 | 2.5.3 | 239 passaram, 0 falharam |
| Python 3.14.6 | 2.5.4 | 239 passaram, 0 falharam |

A integração contínua (`.github/workflows/testes.yml`) repete a suíte a cada envio, com
Python 3.10 a 3.14 e com NumPy 1.24.
