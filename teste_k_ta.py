# 1. Plotar K vs Ta para diferentes valores de Xs
import numpy as np
import matplotlib.pyplot as plt

Ta_range = np.linspace(20, 200, 100)  # Faixa de temperatura do ar
Xs_values = [0.14, 0.18, 0.22]  # Diferentes umidades iniciais

plt.figure(figsize=(10, 6))

def K_adjusted(Ta, Xs):
    # Modelo simplificado com dependência linear em Ta e Xs
    return 0.01 * (Ta / 40) * (Xs + 0.1)

for Xs in Xs_values:
    ##K_values = (-4.7e-3 * Ta_range + 0.77) * Xs**2 + (2.2e-3 * Ta_range - 0.25) * Xs + 2.7e-3 * np.exp(71.81 / Ta_range)
    K_values = K_adjusted(Ta_range, Xs)  # Chama a função ajustada
    plt.plot(Ta_range, K_values, label=f"Xs = {Xs}")


plt.xlabel('Temperatura do Ar (°C)')
plt.ylabel('Coeficiente de Transferência de Massa (K)')
plt.title('Variação de K com Temperatura do Ar para Diferentes Xs')
plt.grid(True)
plt.legend()
plt.show()
