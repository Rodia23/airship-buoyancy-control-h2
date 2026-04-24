import numpy as np
import matplotlib.pyplot as plt
from t2_simulacion import _Params  # Importamos tu lógica real

# 1. Obtener parámetros oficiales de tu simulación
p = _Params()
L_v = 2257.0      # kJ/kg (Calor latente del agua)
t_objetivo = 10.0 # Umbral de 10 minutos definido en el paper

# 2. Cálculo dinámico de masa de agua (basado en tu Eq. 14)
# Usamos las densidades calculadas por tu propio modelo
rho_aire = p.rho_aire_std * np.exp(-p.h0 / p.H_escala)
rho_h2 = p.P_atmos / (p.R_H2 * p.T1)
m_h2_quemada = p.m_carga * (rho_h2 / (rho_aire - rho_h2))
m_h2o = m_h2_quemada * 9.0  # Relación estequiométrica

# Cálculo del umbral de potencia
Q_min_kW = (m_h2o * L_v) / (t_objetivo * 60)

# 3. Generación de datos para la curva
Q_range = np.logspace(1.5, 4, 500) # De ~30kW a 10,000kW
t_range = (m_h2o * L_v) / (Q_range * 60)

# 4. Gráfica estética
plt.figure(figsize=(9, 6.5), dpi=150)
ax = plt.gca()

# Fondo y rejilla sutil
ax.set_facecolor('#fcfcfc')
plt.grid(True, which="both", ls="-", alpha=0.1, color='gray')

# Dibujar la curva y zonas
plt.plot(Q_range, t_range, label='Condensation curve', color='#2c7bb6', lw=3, zorder=3)
plt.axhspan(0, t_objetivo, color='green', alpha=0.08, label='Operational zone ($t < 10$ min)')

# Líneas de umbral naranja
plt.axvline(Q_min_kW, color='darkorange', ls='--', lw=1.5, alpha=0.6)
plt.axhline(t_objetivo, color='darkorange', ls='--', lw=1.5, alpha=0.6)
plt.scatter([Q_min_kW], [t_objetivo], color='darkorange', s=80, edgecolors='white', zorder=5)

# MEJORA: Etiqueta visible con offset
ax.annotate(f'Critical threshold:\n$Q_{{min}} \\approx {Q_min_kW:.0f}$ kW\n$t = 10$ min',
            xy=(Q_min_kW, t_objetivo), xytext=(Q_min_kW * 1.6, t_objetivo * 4),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='darkorange'),
            fontsize=11, fontweight='bold', color='darkorange',
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="darkorange", alpha=0.9))

# Formato de ejes y "Offset" inferior
plt.xscale('log')
plt.ylim(-8, 500) # Offset negativo para que el eje X respire
plt.xlim(30, 10000)

plt.xlabel('Condensation power $Q_{cond}$ [kW]', fontsize=12, labelpad=12)
plt.ylabel('Condensation time [min]', fontsize=12, labelpad=12)
plt.title(f'Water vapor condensation time\n($m_{{H_2O}} = {m_h2o:.1f}$ kg for a $1000$ kg payload)', 
          fontsize=14, pad=20, fontweight='bold')

plt.legend(frameon=True, loc='upper right', facecolor='white', framealpha=1)
plt.tight_layout()

plt.savefig("grafica_condensacion_final.png")
plt.show()

print(f"Synchronized water mass: {m_h2o:.2f} kg")
print(f"Calculated minimum power: {Q_min_kW:.2f} kW")