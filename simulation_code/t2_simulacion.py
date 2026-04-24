# =============================================================================
# t2_simulacion.py (v12.0 — Versión Final con Escala Log-X y Superposición)
# =============================================================================

import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import parametros as _p_default

def correr(n_entregas=1, toneladas_por_drop=1, graficar=True):
    """
    Simula una misión de 'n' entregas. 
    Para n=1, usa escala logarítmica en tiempo y superpone masa/flujo.
    """
    # 1. Ajuste de carga y tiempos
    carga_total_kg = n_entregas * toneladas_por_drop * 1000.0
    p = _Params({'m_carga': carga_total_kg})
    
    # Si es 1 entrega, drop al segundo 1.0 para log-x. Si no, cada 30s.
    intervalo = 1.0 if n_entregas == 1 else 30.0
    t_final = 250.0 if n_entregas == 1 else (n_entregas + 3) * intervalo
    
    # 2. Cálculos iniciales de flotabilidad [cite: 12, 13]
    rho_h0 = p.rho_aire_std * np.exp(-p.h0 / p.H_escala)
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini = (p.m_estructura + p.m_carga) / (rho_h0 / rho_H2_ini - 1)

    ultimo_drop_anunciado = -1

    print("\n" + "="*60)
    print(f"MISIÓN CONFIGURADA: {n_entregas} entregas de {toneladas_por_drop} t")
    print(f"Inercia Inicial: {p.m_estructura + p.m_carga + m_H2_ini:.1f} kg")
    print("="*60)

    def dinamica(t, y):
        nonlocal ultimo_drop_anunciado
        v, h, m_H2, T_H2, Q_con, Q_gen = y
        m_H2 = max(m_H2, 0.1)

        # Lógica de liberación escalonada 
        drops_hechos = int(t // intervalo)
        drops_hechos = min(drops_hechos, n_entregas)
        
        m_carga_actual = p.m_carga - (drops_hechos * (toneladas_por_drop * 1000.0))
        M_total_nave = p.m_estructura + m_carga_actual + m_H2

        if drops_hechos > ultimo_drop_anunciado:
            print(f"[t={t:6.1f}s] DROP #{drops_hechos} ejecutado. Masa restante: {m_carga_actual/1000:4.1f} t")
            ultimo_drop_anunciado = drops_hechos

        # Sintonía Adaptativa de Kd 
        zeta = 1.35 
        Kd_auto = 2 * zeta * np.sqrt(p.Kp * (p.m_estructura + m_carga_actual))
        
        rho_aire = p.rho_aire_std * np.exp(-h / p.H_escala)
        V_H2 = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        V_eq = ((p.m_estructura + m_carga_actual) * p.R_H2 * p.T1) / (p.P_atmos * (rho_h0 / rho_H2_ini - 1))
        
        # Controlador PD [cite: 13]
        V_cmd = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        dm = p.K_mass_control * (V_cmd - V_H2)
        
        # Pulso proactivo (ajustado para drop en 1s)
        t_anticipo = min(2.0, intervalo)
        if (intervalo - t_anticipo) <= (t % intervalo) < intervalo and drops_hechos < n_entregas:
            dm = p.dm_H2_proactivo

        dm = np.clip(dm, -p.max_flujo_H2, p.max_flujo_H2)
        
        # Física de Vuelo [cite: 13]
        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * p.Area * rho_aire * v**2 * np.sign(v)
        dv = (F_emp - M_total_nave * p.g - F_arr) / M_total_nave

        dQgen = abs(dm) * p.Q_comb_H2 * p.rend_gen if dm < -0.05 else 0
        
        return [dv, v, dm, (p.T1 - T_H2)/5.0, 0, dQgen]

    # Integración RK45 
    y0 = [0.0, p.h0, m_H2_ini, p.T1, 0.0, 0.0]
    sol = solve_ivp(dinamica, [0, t_final], y0, method='RK45', max_step=0.1)
    
    t, v, h, m_H2, Q_gen = sol.t, sol.y[0], sol.y[1], sol.y[2], sol.y[5]
    a = np.gradient(v, t)
    dm_real = np.gradient(m_H2, t)
    
    if graficar:
        if n_entregas == 1:
            _graficar_mision_base_estetica(t, h, v, a, dm_real, m_H2, Q_gen/3.6e6, p, t_lib=intervalo)
        else:
            _graficar_super_mision(t, h, v, a, dm_real, m_H2, Q_gen/3.6e6, p, n_entregas, intervalo, toneladas_por_drop)

    return {"Energia": Q_gen[-1]/3.6e6}

def _graficar_mision_base_estetica(t, h, v, a, dm, m_h2, E, p, t_lib=1.0):
    """ Gráfica 3 paneles: Cinemática log-x, Flujo/Masa superpuestos y Energía. """
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True, 
                             gridspec_kw={'height_ratios': [2.5, 1.5, 1]})
    plt.subplots_adjust(hspace=0.15)
    
    # PANEL 1: Cinemática (Escala Log-X)
    ax_alt = axes[0]
    ax_vel, ax_acc = ax_alt.twinx(), ax_alt.twinx()
    ax_acc.spines["right"].set_position(("axes", 1.12))
    ax_acc.spines["right"].set_visible(True)
    ax_acc.yaxis.set_visible(True) 

    ax_alt.plot(t, h, color='steelblue', lw=2.5, label='Altitude')
    ax_vel.plot(t, v, color='forestgreen', lw=1.8, ls='--', label='Velocity')
    ax_acc.plot(t, a, color='indianred', lw=1.2, ls='-.', alpha=0.7, label='Acceleration')

    ax_alt.set_xscale('log')
    ax_alt.set_xlim(0.1, 250) 
    ax_alt.axhline(p.h0, color='steelblue', ls=':', alpha=0.6)
    ax_alt.axvline(t_lib, color='black', ls='--', lw=1.2, label=f'Release ({t_lib}s)')
    
    ax_alt.set_ylabel('Altitude [m]', color='steelblue', weight='bold')
    ax_vel.set_ylabel('Velocity [m/s]', color='forestgreen', weight='bold')
    ax_acc.set_ylabel('Acceleration [m/s²]', color='indianred', weight='bold')
    ax_alt.set_title(f'Dynamic Analysis (Log-X) — {p.nombre_dirigible}', fontsize=12, pad=10)
    ax_alt.grid(True, which='both', alpha=0.2)
    lines_alt, labels_alt = ax_alt.get_legend_handles_labels()
    lines_vel, labels_vel = ax_vel.get_legend_handles_labels()
    lines_acc, labels_acc = ax_acc.get_legend_handles_labels()
    ax_alt.legend(lines_alt + lines_vel + lines_acc,
              labels_alt + labels_vel + labels_acc,
              loc='upper right', fontsize=8)

    # PANEL 2: Flujo + Masa (SUPERPUESTOS)
    ax_flujo = axes[1]
    ax_masa = ax_flujo.twinx()
    
    lf, = ax_flujo.plot(t, dm, color='firebrick', lw=1.5, label='H₂ Flow')
    lm, = ax_masa.plot(t, m_h2, color='darkorchid', lw=0.8, ls=':', marker='.', 
                        markersize=2, markevery=15, alpha=0.5, label='H₂ Mass')
    
    ax_flujo.axvline(t_lib, color='black', ls='--', lw=1)
    ax_flujo.set_ylabel('H₂ Flow [kg/s]', color='firebrick', weight='bold')
    ax_masa.set_ylabel('H₂ Mass [kg]', color='darkorchid', weight='bold')
    ax_flujo.grid(True, which='both', alpha=0.2)
    ax_flujo.legend([lf, lm], [l.get_label() for l in [lf, lm]], loc='lower right', fontsize=8)

    # PANEL 3: Energía
    ax_energia = axes[2]
    ax_energia.plot(t, E, color='gold', lw=2.5)
    ax_energia.axvline(t_lib, color='black', ls='--', lw=1)
    ax_energia.set_ylabel('Energy [kWh]', weight='bold')
    ax_energia.set_xlabel('Time [s] (Logarithmic Scale)', weight='bold')
    ax_energia.grid(True, which='both', alpha=0.2)

    plt.savefig(f"mision_base_final_logx.png", dpi=200, bbox_inches='tight')
    plt.show()

def _graficar_super_mision(t, h, v, a, dm, m_h2, E, p, n, intervalo, tons_p):
    """ Gráfica estándar para misiones de múltiples drops. """
    fig, axes = plt.subplots(5, 1, figsize=(12, 16), sharex=True)
    # plt.subplots_adjust(hspace=0.25) # Lo reemplazaremos con tight_layout más abajo
    
    # 1. Altitud
    axes[0].plot(t, h, 'steelblue', lw=2.5)
    axes[0].set_ylabel("Altitude [m]", weight='bold')
    axes[0].set_title(f"Multi-Drop Mission: {n} releases of {tons_p}t each")
    
    # 2. Velocidad y Aceleración
    ax_v = axes[1]
    ax_v.plot(t, v, 'forestgreen')
    ax_v.set_ylabel("Velocity [m/s]", color='forestgreen', weight='bold') # <-- FALTABA ESTO
    
    ax_a = ax_v.twinx()
    ax_a.plot(t, a, 'firebrick', alpha=0.3, ls='-.')
    ax_a.set_ylabel("Acceleration [m/s²]", color='firebrick', weight='bold') # <-- FALTABA ESTO
    
    # 3. Flujo H2
    axes[2].plot(t, dm, 'firebrick', lw=1.2)
    axes[2].set_ylabel("H₂ Flow [kg/s]", weight='bold')
    
    # 4. Masa H2
    axes[3].plot(t, m_h2, 'darkorchid', lw=2)
    axes[3].set_ylabel("H₂ Mass [kg]", weight='bold')
    
    # 5. Energía
    axes[4].plot(t, E, 'gold', lw=3)
    axes[4].set_ylabel("Energy [kWh]", weight='bold')
    
    for ax in axes: 
        ax.grid(True, alpha=0.2)
    
    # Ajuste automático de márgenes para que los labels de los bordes no se corten
    fig.tight_layout()
    
    plt.savefig(f"mision_multi_{n}drops.png", dpi=150)
    plt.show()
class _Params:
    def __init__(self, overrides=None):
        import parametros as base
        for k, v in vars(base).items():
            if not k.startswith('_'): setattr(self, k, v)
        if overrides:
            for k, v in overrides.items(): setattr(self, k, v)

if __name__ == '__main__':
    # Configuración: 1 entrega de 1 tonelada
    correr(n_entregas=4, toneladas_por_drop=1, graficar=True)