import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Função de simulação
def simulate_soybean_drying(Ta, Xs_ini, Ts_ini):
    # Parâmetros do secador e propriedades físicas
    L = 1.47
    De = 0.16
    N = 4
    S = np.tan(1.57 * np.pi / 180)
    dp = 0.002
    ur = 0.19
    ms_dot = 0.00333
    rho_s = 700
    rho_ap = 500
    cps = 1.67
    cpw = 4.18
    epsilon = 1 - rho_ap / rho_s
    tr = 0.19 * L / (N**0.19 * De * S)
    vs = L / tr
    def K_adjusted(Ta, Xs):
        # Modelo simplificado com dependência linear em Ta e Xs
        return 0.01 * (Ta / 40) * (Xs + 0.1)
    def K_complete(Ta, u, dp, rho_a=1.2, mu=1.8e-5, D_AB=2.5e-5, rho_s=700, c_p=1.67):
         # Número de Reynolds
        Re = rho_a * u * dp / mu
        # Número de Schmidt
        Sc = mu / (rho_a * D_AB)
        # Número de Sherwood
        Sh = 2 + (0.6 * Re**0.5 * Sc**(1/3))
        # Correlação empírica clássica: Sherwood
        C, m = 2, 0.6  # valores típicos para partículas em fluxo de ar
        Sh = C * (Re ** m) * (Sc ** (1/3))
        # Coeficiente convectivo de massa
        k_m = Sh * D_AB / dp  # [m/s]
        # Coeficiente global de transferência de massa volumétrico
        K = k_m / (rho_s * c_p)  # [1/s]
        return K
    
    def drying_ode(t, y):
        Xs, Ts = y
        #K = (-4.7e-3 * Ta + 0.77) * Xs**2 + (2.2e-3 * Ta - 0.25) * Xs + 2.7e-3 * np.exp(71.81 / Ta)
        K= K_adjusted(Ta, Xs)  # Chama a função ajustada
        Xe = 0.834 / (1 + 0.036 * (Ta + 273.15) * np.log(1 / ur))
        R = K * (Xs - Xe)
        as_ = 1000
        hc = 0.05
        lambda_ = 2500 - cpw * Ts
        dXsdt = -R
        denominator = cps + 0.22 * cpw
        Q_dot = as_ * hc * (Ta - Ts)
        dTsdt = (1 / denominator) * (Q_dot / rho_s - lambda_ * R)
        return [dXsdt, dTsdt]
    
    y0 = [Xs_ini, Ts_ini]
    t_span = (0, tr)
    t_eval = np.linspace(0, tr, 100)
    sol = solve_ivp(drying_ode, t_span, y0, t_eval=t_eval)
    z = vs * sol.t
    return z, sol.y[0], sol.y[1], tr, vs

# Condições de simulação
conditions = [
    {'Ta': 40, 'Xs_ini': 0.22, 'Ts_ini': 25},
    {'Ta': 50, 'Xs_ini': 0.22, 'Ts_ini': 25},
    {'Ta': 60, 'Xs_ini': 0.22, 'Ts_ini': 25},
    {'Ta': 40, 'Xs_ini': 0.18, 'Ts_ini': 25},
    {'Ta': 50, 'Xs_ini': 0.18, 'Ts_ini': 25},
    {'Ta': 60, 'Xs_ini': 0.18, 'Ts_ini': 25},
    {'Ta': 40, 'Xs_ini': 0.14, 'Ts_ini': 25},
    {'Ta': 50, 'Xs_ini': 0.14, 'Ts_ini': 25},
    {'Ta': 60, 'Xs_ini': 0.14, 'Ts_ini': 25}
]  

# Plot dos perfis de umidade
plt.figure(figsize=(10, 6))

for cond in conditions:
    z, Xs, Ts, tr, vs = simulate_soybean_drying(cond['Ta'], cond['Xs_ini'], cond['Ts_ini'])
    label = f"Ta={cond['Ta']}°C, Xs_ini={cond['Xs_ini']}, Ts_ini={cond['Ts_ini']}°C"
    plt.plot(z, Xs, label=label)

plt.xlabel("Comprimento do Secador (m)")
plt.ylabel("Umidade (b.s.)")
plt.title("Perfis de Umidade ao Longo do Secador")
plt.legend()
plt.grid(True)
plt.show()
