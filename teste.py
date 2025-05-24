# Vou apresentar o código modificado do MATLAB convertido para Python, adaptado ao modelo descrito por Gianini.
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parâmetros do secador rotativo direto (com base na tese de Gianini)
L = 1.47        # Comprimento do secador (m)
De = 0.16       # Diâmetro do secador (m)
N = 4           # Rotação do secador (rpm)
S = np.tan(1.57 * np.pi / 180)  # Inclinação do secador (rad)
dp = 0.002      # Diâmetro médio das partículas (m)
Xs_ini = 0.22   # Umidade inicial (base seca)
Ta = 55         # Temperatura do ar (°C)
ur = 0.19       # Umidade relativa do ar
ms_dot = 0.00333 # Vazão mássica do sólido (kg/s)

# Propriedades físicas
rho_s = 700     # Massa específica do sólido seco (kg/m³)
rho_ap = 500    # Massa específica aparente (kg/m³)
cps = 1.67      # Calor específico do sólido (kJ/kg°C)
cpw = 4.18      # Calor específico da água (kJ/kg°C)

# Parâmetros auxiliares
epsilon = 1 - rho_ap / rho_s
psi = 10 ** (-1.82 * (1 - epsilon))
tr = 0.19 * L / (N**0.19 * De * S) # Tempo de residência (s)
vs = L / tr                         # Velocidade axial do sólido (m/s)

# Condições iniciais
y0 = [Xs_ini, 25]  # Umidade inicial e temperatura inicial do sólido (°C)

# Equações diferenciais baseadas em Gianini
def drying_ode(t, y):
    Xs, Ts = y

    # Equações da tese Gianini
    
    K = (-4.7e-3 * Ta + 0.77) * Xs ** 2 + (2.2e-3 * Ta - 0.25) * Xs + 2.7e-3 * np.exp(71.81 / Ta)
    Xe = 0.834 / (1 + 0.036 * (Ta + 273.15) * np.log(1 / ur))
    R = K * (Xs - Xe)

    # Transferência de calor
    as_ = 1000
    hc = 0.05
    lambda_ = 2500 - cpw * Ts

    # Balanço de massa e energia
    dXsdt = -R
    denominator = cps + 0.22 * cpw
    Q_dot = as_ * hc * (Ta - Ts)
    dTsdt = (1 / denominator) * (Q_dot / rho_s - lambda_ * R)

    return [dXsdt, dTsdt]

# Tempo da simulação
t_span = (0, tr)
t_eval = np.linspace(0, tr, 100)

# Resolver as EDOs
sol = solve_ivp(drying_ode, t_span, y0, t_eval=t_eval)

# Plotar resultados
fig, ax = plt.subplots(1, 2, figsize=(14, 5))

# Umidade vs tempo
ax[0].plot(sol.t, sol.y[0], 'b-', linewidth=2)
ax[0].set_xlabel('Tempo (s)')
ax[0].set_ylabel('Umidade (base seca)')
ax[0].set_title('Umidade do Farelo com o Tempo')
ax[0].grid(True)

# Umidade vs comprimento do secador
z = vs * sol.t
ax[1].plot(z, sol.y[0], 'r-', linewidth=2)
ax[1].set_xlabel('Comprimento do Secador (m)')
ax[1].set_ylabel('Umidade (base seca)')
ax[1].set_title('Umidade ao Longo do Secador')
ax[1].grid(True)

plt.tight_layout()
plt.show()

# Exibir resultados finais
Xs_final = sol.y[0, -1]
Ts_final = sol.y[1, -1]

(tr, vs, Xs_final, Ts_final)
