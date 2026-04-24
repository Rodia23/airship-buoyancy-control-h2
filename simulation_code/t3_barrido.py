# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# t3_barrido.py  (v2 — corregido y mejorado)
# Análisis paramétrico: tasa de combustión vs desplazamiento/estabilización


import csv
import itertools
from datetime import datetime
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import parametros as p

# =============================================================================
# Configuración del barrido
# =============================================================================

DIRIGIBLES = [
    {'nombre': 'Pathfinder 3',         'Area': 4390,  'm_estructura': 99000,  'm_carga': 1000},
    {'nombre': 'Flying Whales LCA60T', 'Area': 7854,  'm_estructura': 129000, 'm_carga': 1000},
    {'nombre': 'Zeppelin NT',          'Area': 837,   'm_estructura': 5622,   'm_carga': 1000},
    {'nombre': 'LZ 129 Hindenburg',    'Area': 7950,  'm_estructura': 263600, 'm_carga': 1000},
]

# Tasas globales (todos los dirigibles)
TASAS_H2_BASE = [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 800]

# Tasas adicionales SOLO para Zeppelin NT (para acotar el mínimo real)
TASAS_H2_ZEPPELIN_EXTRA = [25, 30, 35, 40, 45]

FRANJAS_TECH = [
    (0.01,  0.1,  'royalblue',  0.15, 'Small PEM Cell'),
    (0.5,   5.0,  'seagreen',   0.15, 'H₂ Combustion Engine'),
    (5.0,   50.0, 'darkorange', 0.15, 'Industrial System'),
    (50.0,  800.0,'firebrick',  0.10, 'Theoretical Control Limit'),
]

# Criterios de estabilización (ajustados por dirigible en evaluar_caso)
CRITERIO_H_M    = 1.0   # [m]   umbral de altitud para t_estab
CRITERIO_V_PESADO = 0.05  # [m/s] para dirigibles pesados (Pathfinder, FW, Hindenburg)
CRITERIO_V_LIGERO = 0.10  # [m/s] para Zeppelin NT (oscilaciones más lentas)

# =============================================================================
# Motor de Física
# =============================================================================

def evaluar_caso(diri, tasa_max):
    m_est       = diri['m_estructura']
    Area        = diri['Area']
    m_carga_ini = diri['m_carga']
    es_ligero   = m_est < 10000  # Zeppelin NT: 5622 kg

    # Cálculos iniciales de flotabilidad
    rho_h0     = p.rho_aire_std * np.exp(-p.h0 / p.H_escala)
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini   = (m_est + m_carga_ini) / (rho_h0 / rho_H2_ini - 1)

    t_drop  = 30.0   # [s] instante de liberación
    t_final = 250.0  # [s] tiempo total de simulación

    def dinamica(t, y):
        v, h, m_H2, T_H2 = y
        m_H2 = max(m_H2, 0.1)

        ha_soltado     = t >= t_drop
        m_carga_actual = 0.0 if ha_soltado else m_carga_ini
        M_total_nave   = m_est + m_carga_actual + m_H2

        # Sintonía adaptativa de Kd
        zeta     = 1.35
        Kd_auto  = 2 * zeta * np.sqrt(p.Kp * (m_est + m_carga_actual))

        rho_aire = p.rho_aire_std * np.exp(-h / p.H_escala)
        V_H2     = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        V_eq     = ((m_est + m_carga_actual) * p.R_H2 * p.T1) / (
                    p.P_atmos * (rho_h0 / rho_H2_ini - 1))

        # Controlador PD
        V_cmd = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        dm    = p.K_mass_control * (V_cmd - V_H2)

        

        if (t_drop - 2.0) <= t < t_drop:
            dm = p.dm_H2_proactivo  # = -10 kg/s definido en parametros.py
        else:
            dm = np.clip(dm, -tasa_max, tasa_max)

        # Dinámica de vuelo
        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * Area * rho_aire * v**2 * np.sign(v)
        dv    = (F_emp - M_total_nave * p.g - F_arr) / M_total_nave
        dT    = (p.T1 - T_H2) / 5.0

        return [dv, v, dm, dT]

    y0  = [0.0, p.h0, m_H2_ini, p.T1]
    sol = solve_ivp(dinamica, [0, t_final], y0, method='RK45', max_step=0.2)

    t, v, h = sol.t, sol.y[0], sol.y[1]

    # Análisis post-drop
    idx_post = t >= t_drop
    if not np.any(idx_post):
        return None, None, None

    t_post = t[idx_post]
    h_post = h[idx_post]
    v_post = v[idx_post]

    h_desp = np.max(h_post) - p.h0
    a_post = np.gradient(v_post, t_post)
    a_max  = np.max(np.abs(a_post)) / p.g

    # FIX: criterio de velocidad relajado para dirigibles ligeros
    crit_v = CRITERIO_V_LIGERO if es_ligero else CRITERIO_V_PESADO

    
    t_estab = float('nan')
    for i in range(len(t_post)):
        if abs(h_post[i] - p.h0) < CRITERIO_H_M and abs(v_post[i]) < crit_v:
            t_estab = t_post[i] - t_drop
            break

    h_desp_r  = round(h_desp, 3)
    a_max_r   = round(a_max, 4)
    t_estab_r = round(t_estab, 1) if not np.isnan(t_estab) else 'N/A'
    return h_desp_r, a_max_r, t_estab_r

# =============================================================================
# Ejecutar Barrido
# =============================================================================

def ejecutar_barrido():
    timestamp   = datetime.now().strftime('%Y%m%d_%H%M%S')
    archivo_csv = f'resultados_t3_barrido_{timestamp}.csv'

    print(f"{'='*65}")
    print(f"  BARRIDO PARAMÉTRICO t3 (v2) - 1 Drop")
    print(f"  Zeppelin NT incluye tasas adicionales: {TASAS_H2_ZEPPELIN_EXTRA}")
    print(f"{'='*65}\n")

    # Construir lista de combinaciones: base para todos, extras solo para ZepNT
    combinaciones = []
    for diri in DIRIGIBLES:
        tasas = list(TASAS_H2_BASE)
        if diri['nombre'] == 'Zeppelin NT':
            tasas = sorted(set(tasas + TASAS_H2_ZEPPELIN_EXTRA))
        for tasa in tasas:
            combinaciones.append((diri, tasa))

    total      = len(combinaciones)
    resultados = []

    for i, (diri, tasa) in enumerate(combinaciones, 1):
        nombre = diri['nombre']
        print(f"  [{i:>3}/{total}] {nombre:<25} ṁ_max = {tasa:>6} kg/s ... ",
              end='', flush=True)

        h_d, a_m, t_e = evaluar_caso(diri, tasa)

        fila = {
            'nombre':       nombre,
            'max_flujo_H2': tasa,
            'h_desp_m':     h_d,
            'a_max_g':      a_m,
            't_estab_s':    t_e,
        }
        resultados.append(fila)
        print(f"|Δh|={h_d:>8} m  |  a_max={a_m:>7} g  |  t_estab={t_e}")

    # Guardar CSV
    with open(archivo_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(resultados[0].keys()))
        writer.writeheader()
        writer.writerows(resultados)

    print(f"\n  CSV guardado: {archivo_csv}\n")

    # Resumen de mínimos por dirigible
    _imprimir_resumen_minimos(resultados)

    # Gráficas
    _graficar_resumen(resultados, timestamp)


def _imprimir_resumen_minimos(resultados):
    """Imprime la tabla de tasas mínimas y óptimas al terminar el barrido."""
    print(f"\n{'='*65}")
    print("  RESUMEN: TASAS MÍNIMAS Y ÓPTIMAS")
    print(f"  Umbral seguro : |Δh| ≤ 5.0 m")
    print(f"  Umbral óptimo : |Δh| ≤ 1.0 m")
    print(f"{'='*65}")

    dirigibles = list(dict.fromkeys(r['nombre'] for r in resultados))
    for nombre in dirigibles:
        sub = sorted([r for r in resultados if r['nombre'] == nombre],
                     key=lambda r: float(r['max_flujo_H2']))
        min_seg = next((r for r in sub
                        if r['h_desp_m'] is not None and r['h_desp_m'] <= 5.0), None)
        min_opt = next((r for r in sub
                        if r['h_desp_m'] is not None and r['h_desp_m'] <= 1.0), None)
        print(f"\n  {nombre}")
        if min_seg:
            print(f"    ṁ_min (seguro) = {min_seg['max_flujo_H2']:>5} kg/s  "
                  f"→ |Δh|={min_seg['h_desp_m']} m, t_estab={min_seg['t_estab_s']} s")
        else:
            print(f"    ṁ_min (seguro) = > 800 kg/s  (no alcanzado)")
        if min_opt:
            print(f"    ṁ_min (óptimo) = {min_opt['max_flujo_H2']:>5} kg/s  "
                  f"→ |Δh|={min_opt['h_desp_m']} m, t_estab={min_opt['t_estab_s']} s")
        else:
            print(f"    ṁ_min (óptimo) = > 800 kg/s  (no alcanzado)")
    print()

# =============================================================================
# Gráficas
# =============================================================================

def _graficar_resumen(resultados, timestamp):
    dirigibles_unicos = list(dict.fromkeys(r['nombre'] for r in resultados))
    colores = plt.cm.tab10(np.linspace(0, 0.9, len(dirigibles_unicos)))

    # -------------------------------------------------------------------------
    # GRÁFICA 1: Desplazamiento Vertical
    # -------------------------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(7, 5))

    for x0, x1, col, alpha, _ in FRANJAS_TECH:
        ax1.axvspan(x0, x1, color=col, alpha=alpha, zorder=0)

    for diri, color in zip(dirigibles_unicos, colores):
        sub   = sorted([r for r in resultados if r['nombre'] == diri],
                       key=lambda r: float(r['max_flujo_H2']))
        tasas = [r['max_flujo_H2'] for r in sub]
        h_d   = [r['h_desp_m'] if r['h_desp_m'] is not None else float('nan')
                 for r in sub]
        ax1.plot(tasas, h_d, 'o-', color=color, linewidth=1.8, label=diri, zorder=3)

    ax1.axhline(5.0, color='black',   linestyle='--', linewidth=1.3, zorder=4,
                label='Safety threshold = 5.0 m')
    ax1.axhline(1.0, color='dimgray', linestyle=':',  linewidth=1.1, zorder=4,
                label='Precision threshold = 1.0 m')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(left=0.4,  right=1000)
    ax1.set_ylim(bottom=1e-3, top=1000)  # FIX: límite inferior 1e-3 (antes 0.1)
    ax1.grid(True, which='major', alpha=0.35)
    ax1.grid(True, which='minor', alpha=0.15)
    ax1.set_xlabel(r'$\dot{m}_{H_2,max}$ — Max. combustion rate [kg/s]', fontsize=10)
    ax1.set_ylabel('Max. vertical displacement after release [m]', fontsize=10)
    ax1.set_title('Vertical displacement vs combustion rate', fontsize=11)
    ax1.legend(fontsize=8, loc='upper right')

    plt.tight_layout()
    nombre_png1 = f't3_barrido_desplazamiento_{timestamp}.png'
    fig1.savefig(nombre_png1, dpi=150, bbox_inches='tight')
    plt.close(fig1)
    print(f"  Gráfica 1 guardada: {nombre_png1}")

    # -------------------------------------------------------------------------
    # GRÁFICA 2: Tiempo de Estabilización
    # -------------------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(7, 5))

    for x0, x1, col, alpha, _ in FRANJAS_TECH:
        ax2.axvspan(x0, x1, color=col, alpha=alpha, zorder=0)

    for diri, color in zip(dirigibles_unicos, colores):
        sub   = sorted([r for r in resultados if r['nombre'] == diri],
                       key=lambda r: float(r['max_flujo_H2']))
        tasas = [r['max_flujo_H2'] for r in sub]
        t_e   = [r['t_estab_s'] if isinstance(r['t_estab_s'], (int, float))
                 else float('nan') for r in sub]
        ax2.plot(tasas, t_e, '^-', color=color, linewidth=1.8, label=diri, zorder=3)

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(left=0.4,  right=1000)
    ax2.set_ylim(bottom=0.1, top=300)
    ax2.grid(True, which='major', alpha=0.35)
    ax2.grid(True, which='minor', alpha=0.15)
    ax2.set_xlabel(r'$\dot{m}_{H_2,max}$ — Max. combustion rate [kg/s]', fontsize=10)
    ax2.set_ylabel('Stabilization time [s]', fontsize=10)
    ax2.set_title('Stabilization time vs combustion rate', fontsize=11)

    lineas, labels = ax2.get_legend_handles_labels()
    leg_lineas = ax2.legend(handles=lineas, labels=labels,
                            loc='upper right', fontsize=8)
    ax2.add_artist(leg_lineas)

    parches = [mpatches.Patch(color=col, alpha=alpha + 0.2, label=lbl)
               for _, _, col, alpha, lbl in FRANJAS_TECH]
    ax2.legend(handles=parches, fontsize=8, loc='lower right', title='Technology')

    plt.tight_layout()
    nombre_png2 = f't3_barrido_tiempo_{timestamp}.png'
    fig2.savefig(nombre_png2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"  Gráfica 2 guardada: {nombre_png2}")


if __name__ == '__main__':
    ejecutar_barrido()