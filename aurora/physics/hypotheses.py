"""Hipóteses de modelagem compartilhadas entre as fórmulas (texto exibido ao usuário)."""

CLASSICAL = "Mecânica clássica newtoniana: velocidades muito menores que c, sem correções relativísticas."
INERTIAL_FRAME = "Grandezas medidas em um referencial inercial."
CONSTANT_ACCELERATION = "Aceleração constante durante todo o intervalo considerado."
ONE_DIMENSIONAL = "Movimento retilíneo (uma dimensão); o sinal indica o sentido no eixo."

POINT_MASSES = (
    "Corpos tratados como massas pontuais ou esferas com distribuição de massa esfericamente simétrica."
)
G_VALUE = (
    "G = 6,67430 × 10⁻¹¹ m³ kg⁻¹ s⁻² (CODATA 2018); a incerteza relativa de G (≈ 2,2 × 10⁻⁵) não é propagada."
)
MU_GIVEN = (
    "μ = G·M informado diretamente (conhecido com precisão muito maior que G e M separados); "
    "sua incerteza não é propagada."
)
TWO_BODY = (
    "Problema de dois corpos: massa do veículo desprezível frente à do corpo central, "
    "que permanece em repouso em um referencial inercial."
)
NO_PERTURBATIONS = (
    "Perturbações desprezadas: arrasto atmosférico, achatamento do corpo central (J2), "
    "pressão de radiação solar e atração de terceiros corpos."
)
RADIUS_FROM_CENTER = "Raios medidos a partir do centro de massa do corpo central (não a partir da superfície)."
POTENTIAL_ZERO_AT_INFINITY = "Energia potencial gravitacional nula no infinito."
