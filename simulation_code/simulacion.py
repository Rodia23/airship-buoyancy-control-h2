# =============================================================================
# simulacion.py  (sim_v2)
# Simulación dinámica del dirigible H2 — maniobra de liberación de carga
#
# Uso directo:
#   python simulacion.py
#   python simulacion.py --sin-graficas
#
# Uso como módulo:
#   from simulacion import correr
#   resultado = correr({'max_flujo_H2': 5.0}, graficar=False)
#
# Retorna dict con métricas clave de la simulación.
# =============================================================================

import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import parametros as _p_default


def correr(overrides=None, graficar=True):
    """
    Ejecuta la simulación dinámica del dirigible.

    Parámetros
    ----------
    overrides : dict, opcional
        Variables de parametros.py a sobreescribir para esta ejecución.
    graficar : bool
        Si True, genera y guarda las gráficas.

    Retorna
    -------
    dict con: nombre, max_flujo_H2, a_max_g, v_max_ms, h_desp_m,
              dm_max_kgs, E_gen_kWh, t_lib_s, t_estab_s, exito
    """
    p = _Params(overrides)

    # =========================================================================
    # Cálculos iniciales
    # =========================================================================

    m_total_ini = p.m_estructura + p.m_carga
    m_total_fin = p.m_estructura

    rho_h0      = p.rho_aire_std * np.exp(-p.h0 / p.H_escala)
    rho_H2_ini  = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini    = m_total_ini / (rho_h0 / rho_H2_ini - 1)
    V_flot      = (m_H2_ini * p.R_H2 * p.T1) / p.P_atmos
    m_H2_T2     = (p.P_atmos * V_flot) / (p.R_H2 * p.T2)

    tau = p.tau_term * (p.escala_tiempo if p.acortar_tiempos else 1.0)

    t_calentar   = tau * 4.6
    t_ini_calor  = 1e-7
    t_fin_calor  = t_ini_calor + t_calentar
    t_fin_mante  = t_fin_calor + p.t_mantener
    t_ini_enfria = t_fin_mante

    rho_H2_post = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_post   = m_total_fin / (rho_h0 / rho_H2_post - 1)

    delta_m = m_H2_T2 - m_H2_post
    t_lib   = t_ini_enfria + (delta_m / abs(p.dm_H2_proactivo)) * 0.2
    t_lib   = max(t_lib, t_ini_enfria)

    V_eq_pre  = (m_H2_T2   * p.R_H2 * p.T2) / p.P_atmos
    V_eq_post = (m_H2_post  * p.R_H2 * p.T1) / p.P_atmos

    # =========================================================================
    # ODE
    # =========================================================================

    def dinamica(t, y):
        v, h, m_H2, T_H2, Q_con, Q_gen = y
        m_H2 = max(m_H2, 0.0)

        m_carga_act = m_total_ini if t < t_lib else m_total_fin
        M_total     = m_carga_act + m_H2_ini

        rho_aire = p.rho_aire_std * np.exp(-h / p.H_escala)
        V_H2     = (m_H2 * p.R_H2 * T_H2) / p.P_atmos

        # Temperatura
        if t < t_ini_calor:
            dT = 0.0;  dQcon = 0.0
        elif t < t_fin_calor:
            dT    = (p.T2 - T_H2) / tau
            dQcon = p.cp_H2 * m_H2 * abs(dT) / (1 - p.perdidas)
        elif t < t_fin_mante:
            dT = (p.T2 - T_H2) * 0.1;  dQcon = 0.0
        elif t < t_lib:
            dT = (p.T1 - T_H2) * p.cooling_factor / tau;  dQcon = 0.0
        else:
            dT = (p.T1 - T_H2) * 10.0 / tau;  dQcon = 0.0

        # Control PD de altitud
        V_eq  = V_eq_pre if t < t_lib else V_eq_post
        e_h   = p.h0 - h
        V_cmd = np.clip(V_eq + p.Kp * e_h - p.Kd * v,
                        V_eq - p.V_delta_max, V_eq + p.V_delta_max)
        dm_control = p.K_mass_control * (V_cmd - V_H2)

        if t_ini_enfria <= t < t_lib:
            dm = p.dm_H2_proactivo
        else:
            dm = dm_control

        dm = max(dm, -p.max_flujo_H2)
        if m_H2 <= 0.01 and dm < 0:
            dm = 0.0

        dQgen = abs(dm) * p.Q_comb_H2 * p.rend_gen if dm < 0 else 0.0

        # Dinámica de vuelo (Newton + Arquímedes)
        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * p.Area * rho_aire * v**2 * np.sign(v)
        dv    = (F_emp - M_total * p.g - F_arr) / M_total

        return [dv, v, dm, dT, dQcon, dQgen]

    # =========================================================================
    # Integración
    # =========================================================================

    y0   = [0.0, p.h0, m_H2_ini, p.T1, 0.0, 0.0]
    tmax = t_lib + p.t_post_liberacion

    try:
        sol = solve_ivp(
            dinamica, [0, tmax], y0,
            method='RK45', max_step=p.dt_max,
            rtol=1e-4, atol=1e-6
        )
        exito = sol.success
    except Exception as e:
        print(f"  [ERROR en ODE] {e}")
        return _resultado_vacio(p, t_lib)

    t    = sol.t
    v    = sol.y[0]
    h    = sol.y[1]
    m_H2 = sol.y[2]
    T_H2 = sol.y[3]
    Q_gen = sol.y[5]

    V_H2 = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
    a    = np.gradient(v, t)
    dm   = np.gradient(m_H2, t)
    P_el = np.gradient(Q_gen, t)

    # =========================================================================
    # Métricas
    # =========================================================================

    idx_post = t >= t_lib
    a_max    = np.max(np.abs(a[idx_post]))
    v_max    = np.max(np.abs(v[idx_post]))
    h_desp   = np.max(h[idx_post]) - p.h0
    dm_max   = np.max(np.abs(dm))
    E_kWh    = Q_gen[-1] / 3.6e6

    t_estab = float('nan')
    for i in np.where(idx_post)[0]:
        if abs(h[i] - p.h0) < 1.0 and abs(v[i]) < 0.05:
            t_estab = t[i] - t_lib
            break

    resultado = {
        'nombre'       : p.nombre_dirigible,
        'max_flujo_H2' : p.max_flujo_H2,
        'a_max_g'      : round(a_max / p.g, 5),
        'v_max_ms'     : round(v_max, 5),
        'h_desp_m'     : round(h_desp, 4),
        'dm_max_kgs'   : round(dm_max, 3),
        'E_gen_kWh'    : round(E_kWh, 3),
        't_lib_s'      : round(t_lib, 3),
        't_estab_s'    : round(t_estab, 2) if not np.isnan(t_estab) else 'N/A',
        'exito'        : exito,
    }

    if graficar:
        _graficar(t, v, h, m_H2, T_H2, V_H2, a, dm, P_el, Q_gen,
                  p, t_lib, V_eq_pre, V_eq_post, m_H2_ini, m_total_ini, m_total_fin)

    return resultado


# =============================================================================
# Auxiliares
# =============================================================================

def _resultado_vacio(p, t_lib):
    return {
        'nombre': p.nombre_dirigible, 'max_flujo_H2': p.max_flujo_H2,
        'a_max_g': None, 'v_max_ms': None, 'h_desp_m': None,
        'dm_max_kgs': None, 'E_gen_kWh': None,
        't_lib_s': round(t_lib, 3), 't_estab_s': None, 'exito': False,
    }


class _Params:
    """Copia los valores de parametros.py y aplica overrides."""
    def __init__(self, overrides=None):
        import parametros as base
        for k, v in vars(base).items():
            if not k.startswith('_'):
                setattr(self, k, v)
        if overrides:
            for k, v in overrides.items():
                setattr(self, k, v)


def _graficar(t, v, h, m_H2, T_H2, V_H2, a, dm, P_el, Q_gen,
              p, t_lib, V_eq_pre, V_eq_post, m_H2_ini, m_total_ini, m_total_fin):

    def marca(ax):
        ax.axvline(t_lib, color='red', linestyle='--', linewidth=1, label='Liberación')

    fig, axes = plt.subplots(5, 2, figsize=(14, 18))
    fig.suptitle(f'Simulación: {p.nombre_dirigible}  |  ṁ_H2_max = {p.max_flujo_H2} kg/s',
                 fontsize=13)

    pares = [
        (axes[0,0], a,           'Aceleración [m/s²]', 'Aceleración vertical', 'b'),
        (axes[0,1], T_H2,        'Temperatura [K]',    'Temperatura del H₂',   'b'),
        (axes[1,0], v,           'Velocidad [m/s]',    'Velocidad vertical',    'b'),
        (axes[1,1], m_H2,        'Masa H₂ [kg]',       'H₂ en ballonets',       'b'),
        (axes[2,0], h,           'Altura [m]',         'Altitud',               'b'),
        (axes[2,1], dm,          'Flujo H₂ [kg/s]',    'Flujo de H₂',           'r'),
        (axes[3,0], V_H2,        'Volumen H₂ [m³]',    'Volumen H₂',            'b'),
        (axes[3,1], P_el/1e3,    'Potencia [kW]',      'Potencia generada',     'g'),
        (axes[4,0], Q_gen/3.6e6, 'Energía [kWh]',      'Energía acumulada',     'g'),
    ]

    for ax, datos, ylabel, titulo, color in pares:
        ax.plot(t, datos, color, linewidth=1.5)
        marca(ax)
        ax.set_xscale('log')
        ax.set_ylabel(ylabel)
        ax.set_title(titulo)
        ax.legend(fontsize=8)
        ax.grid(True)

    axes[0,0].axhline( 0.2*p.g, color='r', linestyle=':', linewidth=1)
    axes[0,0].axhline(-0.2*p.g, color='r', linestyle=':', linewidth=1)
    axes[2,0].axhline(p.h0, color='green', linestyle='--', linewidth=1,
                      label=f'h_ref={p.h0} m')
    axes[2,0].legend(fontsize=8)
    axes[3,0].axhline(V_eq_pre,  color='green',  linestyle='--', linewidth=1, label='V_eq pre')
    axes[3,0].axhline(V_eq_post, color='orange', linestyle='--', linewidth=1, label='V_eq post')
    axes[3,0].legend(fontsize=8)

    ax = axes[4,1]
    rho_t   = p.rho_aire_std * np.exp(-h / p.H_escala)
    F_emp   = rho_t * p.g * V_H2
    m_car_t = np.where(t < t_lib, m_total_ini, m_total_fin)
    F_peso  = (m_car_t + m_H2_ini) * p.g
    ax.plot(t, F_emp/1e3, 'b', linewidth=1.5, label='Empuje [kN]')
    ax.plot(t, F_peso/1e3, 'r', linewidth=1.5, label='Peso [kN]')
    marca(ax)
    ax.set_xscale('log')
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('Fuerza [kN]')
    ax.set_title('Empuje vs Peso')
    ax.legend(fontsize=8)
    ax.grid(True)

    axes[4,0].set_xlabel('Tiempo [s]')

    nombre_arch = f"sim_{p.nombre_dirigible.replace(' ','_')}_flujo{p.max_flujo_H2:.0f}.png"
    plt.tight_layout()
    plt.savefig(nombre_arch, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfica guardada: {nombre_arch}")


if __name__ == '__main__':
    sin_graficas = '--sin-graficas' in sys.argv
    resultado = correr(graficar=not sin_graficas)
    print("\n--- Resultado ---")
    for k, v in resultado.items():
        print(f"  {k:<20}: {v}")
