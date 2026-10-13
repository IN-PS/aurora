"""AURORA — plataforma educacional de simulação de mecânica orbital, radar e exploração espacial.

Instituto Regional de Pesquisas Científicas, Tecnologia e Engenharia do Sul de Minas Gerais.

Organização em camadas (cada camada só depende das camadas acima dela na lista):

- ``aurora.data``: constantes físicas e dados de corpos celestes, com fontes citadas;
- ``aurora.core``: erros, validação de parâmetros, unidades, registro de fórmulas e
  procedência dos dados (analítico / simulado / medido);
- ``aurora.physics``: modelos físicos e matemáticos (funções puras, unidades do SI);
- ``aurora.interface``: interfaces com o usuário (por ora, linha de comando).

Todo resultado produzido por este pacote é calculado ou simulado a partir de modelos
idealizados. Nenhum resultado constitui medição de objetos ou fenômenos reais, e o pacote
não contém código para controlar hardware.
"""

__version__ = "0.1.0"
