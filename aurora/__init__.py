"""AURORA — plataforma educacional de simulação de mecânica orbital, radar e exploração espacial.

Instituto Regional de Pesquisas Científicas, Tecnologia e Engenharia do Sul de Minas Gerais.

Organização em camadas (cada camada só depende das camadas acima dela na lista):

- ``aurora.errors``: hierarquia de exceções (sem dependências);
- ``aurora.data``: constantes físicas e dados de corpos celestes, com fontes citadas;
- ``aurora.core``: validação de parâmetros, unidades, registro de fórmulas e procedência
  dos dados (analítico / simulado / medido);
- ``aurora.physics``: modelos físicos e matemáticos (funções puras, unidades do SI);
- ``aurora.simulation``: simulação numérica (propagação de trajetórias, validação contra
  soluções analíticas, veículo virtual, telemetria simulada e missões idealizadas);
- ``aurora.interface``: interfaces com o usuário (linha de comando e gráficos).

Todo resultado produzido por este pacote é calculado ou simulado a partir de modelos
idealizados. Nenhum resultado constitui medição de objetos ou fenômenos reais, e o pacote
não contém código para controlar hardware.
"""

__version__ = "0.2.0"
