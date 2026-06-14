# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# t3_barrido_o2.py  (basado en t3_barrido.py v2)
# Igual que t3_barrido pero incluye absorción de masa por O₂ atmosférico.
#
# Física añadida:
#   Al quemar H₂ a tasa |ṁ| kg/s, el vehículo ingiere 8·|ṁ| kg/s de O₂
#   y produce 9·|ṁ| kg/s de H₂O que queda como lastre.
#   → M_total_nave += m_agua (estado adicional)
#   → Solo acumula cuando dm < 0 (combustión); no cuando dm > 0 (recarga).
# =============================================================================

import os
import csv
from datetime import datetime
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import parametros as p

OUTDIR = os.path.dirname(os.path.abspath(__file__))  # siempre junto al .py

# =============================================================================
# Configuración del barrido  (idéntica a t3_barrido.py)
# =============================================================================

DIRIGIBLES = [
    {'nombre': 'Pathfinder 3',         'Area': 4390,  'm_estructura': 99000,  'm_carga': 1000},
    {'nombre': 'Flying Whales LCA60T', 'Area': 7854,  'm_estructura': 129000, 'm_carga': 1000},
    {'nombre': 'Zeppelin NT',          'Area': 837,   'm_estructura': 5622,   'm_carga': 1000},
    {'nombre': 'LZ 129 Hindenburg',    'Area': 7950,  'm_estructura': 263600, 'm_carga': 1000},
]

# 6 puntos log-espaciados [1, 100] — mínimo 1 kg/s (por debajo es irrelevante)
TASAS_H2_BASE           = [1, 2, 5, 15, 50, 100]
# Extras solo para Zeppelin NT: su umbral está en ~20 kg/s, necesita resolución ahí
TASAS_H2_ZEPPELIN_EXTRA = [20, 30, 35]

FRANJAS_TECH = [
    (0.01,  0.1,   'royalblue',  0.15, 'Small PEM Cell'),
    (0.5,   5.0,   'seagreen',   0.15, 'H₂ Combustion Engine'),
    (5.0,   50.0,  'darkorange', 0.15, 'Industrial System'),
    (50.0,  120.0, 'firebrick',  0.10, 'Theoretical Control Limit'),
]

CRITERIO_H_M      = 1.0
CRITERIO_V_PESADO = 0.05
CRITERIO_V_LIGERO = 0.10

# =============================================================================
# Opciones de visualización
# =============================================================================
# True  → eje X = ṁ_max / m_est  [kg/(s·t)]  (compara plataformas a igual escala)
# False → eje X = ṁ_max           [kg/s]      (valores absolutos, igual que t3_barrido)
NORMALIZAR_X = False

# =============================================================================
# Motor de Física — con absorción de O₂
# =============================================================================

def evaluar_caso(diri, tasa_max):
    m_est       = diri['m_estructura']
    Area        = diri['Area']
    m_carga_ini = diri['m_carga']
    es_ligero   = m_est < 10000

    rho_h0     = p.rho_aire_std * np.exp(-p.h0 / p.H_escala)
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini   = (m_est + m_carga_ini) / (rho_h0 / rho_H2_ini - 1)

    t_drop  = 30.0
    t_final = 250.0

    # ── Controlador de dos fases (gain scheduling) ────────────────────────────
    # Fase 1 (emergencia, K_mass=60): respuesta rápida post-drop.
    # Fase 2 (trim, K_mass=0.5):      asentamiento suave una vez que m_H2
    #   llega al nuevo equilibrio — evita el overshooting prolongado.
    # Se resetea a Fase 1 en cada nuevo drop (aquí siempre es 1 drop).
    fase2 = [False]

    def dinamica(t, y):
        v, h, m_H2, T_H2, m_agua = y
        m_H2   = max(m_H2,  0.1)
        m_agua = max(m_agua, 0.0)

        ha_soltado     = t >= t_drop
        m_carga_actual = 0.0 if ha_soltado else m_carga_ini

        # Masa total incluye el agua producida por combustión
        M_total_nave = m_est + m_carga_actual + m_H2 + m_agua

        zeta     = 1.35
        rho_aire = p.rho_aire_std * np.exp(-h / p.H_escala)

        # Kd con fórmula dimensionalmente correcta: 2ζ√(Kp·M/(ρ_aire·g))
        # (idéntica a t2_gain_scheduled; m_agua incluida en M_eff)
        M_eff   = m_est + m_carga_actual + m_agua
        Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_eff / (rho_aire * p.g))

        V_H2 = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        V_eq = (M_eff * p.R_H2 * p.T1) / (p.P_atmos * (rho_h0 / rho_H2_ini - 1))

        V_cmd   = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        V_error = V_cmd - V_H2

        # Detección de fase: switch a Fase 2 cuando m_H2 cae al nuevo equilibrio
        if ha_soltado and not fase2[0]:
            m_H2_eq_nuevo = M_eff / (rho_h0 / rho_H2_ini - 1)
            if m_H2 <= m_H2_eq_nuevo:
                fase2[0] = True

        K_mass_eff = 0.5 if fase2[0] else p.K_mass_control   # trim vs emergencia
        dm = K_mass_eff * V_error

        # Pulso proactivo (overrides ambas fases)
        if (t_drop - 2.0) <= t < t_drop:
            dm = p.dm_H2_proactivo
        dm = np.clip(dm, -tasa_max, tasa_max)

        # Acumulación de agua: solo cuando hay combustión (dm < 0)
        # 1 kg H₂ + 8 kg O₂ → 9 kg H₂O
        d_agua = -9.0 * min(dm, 0.0)

        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * Area * rho_aire * v**2 * np.sign(v)
        dv    = (F_emp - M_total_nave * p.g - F_arr) / M_total_nave
        dT    = (p.T1 - T_H2) / 5.0

        return [dv, v, dm, dT, d_agua]

    y0  = [0.0, p.h0, m_H2_ini, p.T1, 0.0]
    sol = solve_ivp(dinamica, [0, t_final], y0, method='RK45', max_step=0.2)

    t, v, h, m_agua = sol.t, sol.y[0], sol.y[1], sol.y[4]

    idx_post = t >= t_drop
    if not np.any(idx_post):
        return None, None, None, None

    t_post  = t[idx_post]
    h_post  = h[idx_post]
    v_post  = v[idx_post]

    h_desp = np.max(h_post) - p.h0
    a_post = np.gradient(v_post, t_post)
    a_max  = np.max(np.abs(a_post)) / p.g

    crit_v = CRITERIO_V_LIGERO if es_ligero else CRITERIO_V_PESADO
    se_desestabilizo = (np.max(np.abs(h_post - p.h0)) >= CRITERIO_H_M) or \
                       (np.max(np.abs(v_post)) >= crit_v)

    if not se_desestabilizo:
        t_estab = 0.0
    else:
        # Criterio por velocidad: robusto para sobre-compensación (Δh negativo).
        # Busca el último instante en que |v| >= umbral → estabilizado justo después.
        mask_mov = np.abs(v_post) >= crit_v
        if not np.any(mask_mov):
            t_estab = 0.0
        else:
            last_idx = np.where(mask_mov)[0][-1]
            if last_idx + 1 < len(t_post):
                t_estab = t_post[last_idx + 1] - t_drop
            else:
                t_estab = float('250')   # sigue en movimiento al final de la sim

    m_agua_final = round(float(m_agua[-1]), 2)                   # <-- nuevo output
    return round(h_desp, 3), round(a_max, 4), \
           round(t_estab, 1) if t_estab < 250.0 else 'N/A', \
           m_agua_final

# =============================================================================
# Ejecutar Barrido
# =============================================================================

def ejecutar_barrido():
    timestamp   = datetime.now().strftime('%Y%m%d_%H%M%S')
    archivo_csv = os.path.join(OUTDIR, f'resultados_t3_barrido_o2_{timestamp}.csv')

    print(f"{'='*65}")
    print(f"  BARRIDO PARAMÉTRICO t3-O2 — con absorción de masa O₂")
    print(f"  Stoichiometry: 1 kg H₂ + 8 kg O₂ → 9 kg H₂O (lastre)")
    print(f"{'='*65}\n")

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

        h_d, a_m, t_e, m_ag = evaluar_caso(diri, tasa)

        fila = {
            'nombre':        nombre,
            'max_flujo_H2':  tasa,
            'h_desp_m':      h_d,
            'a_max_g':       a_m,
            't_estab_s':     t_e,
            'm_agua_kg':     m_ag,       # <-- nuevo campo
        }
        resultados.append(fila)
        print(f"|Δh|={h_d:>8} m  |  a_max={a_m:>7} g  |  t_estab={t_e}  |  H₂O={m_ag} kg")

    with open(archivo_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(resultados[0].keys()))
        writer.writeheader()
        writer.writerows(resultados)

    print(f"\n  CSV guardado: {archivo_csv}\n")
    _imprimir_resumen_minimos(resultados)
    _graficar_resumen(resultados, timestamp)


def _imprimir_resumen_minimos(resultados):
    print(f"\n{'='*65}")
    print("  RESUMEN: TASAS MÍNIMAS Y ÓPTIMAS  (con absorción O₂)")
    print(f"  Umbral seguro : |Δh| ≤ 5.0 m")
    print(f"  Umbral óptimo : |Δh| ≤ 1.0 m")
    print(f"{'='*65}")

    dirigibles = list(dict.fromkeys(r['nombre'] for r in resultados))
    for nombre in dirigibles:
        sub     = sorted([r for r in resultados if r['nombre'] == nombre],
                         key=lambda r: float(r['max_flujo_H2']))
        min_seg = next((r for r in sub if r['h_desp_m'] is not None
                        and r['h_desp_m'] <= 5.0), None)
        min_opt = next((r for r in sub if r['h_desp_m'] is not None
                        and r['h_desp_m'] <= 1.0), None)
        print(f"\n  {nombre}")
        if min_seg:
            print(f"    ṁ_min (seguro) = {min_seg['max_flujo_H2']:>5} kg/s  "
                  f"→ |Δh|={min_seg['h_desp_m']} m, t_estab={min_seg['t_estab_s']} s, "
                  f"H₂O={min_seg['m_agua_kg']} kg")
        else:
            print(f"    ṁ_min (seguro) = > 800 kg/s  (no alcanzado)")
        if min_opt:
            print(f"    ṁ_min (óptimo) = {min_opt['max_flujo_H2']:>5} kg/s  "
                  f"→ |Δh|={min_opt['h_desp_m']} m, t_estab={min_opt['t_estab_s']} s, "
                  f"H₂O={min_opt['m_agua_kg']} kg")
        else:
            print(f"    ṁ_min (óptimo) = > 800 kg/s  (no alcanzado)")
    print()

# =============================================================================
# Gráficas
# =============================================================================

def _graficar_resumen(resultados, timestamp):
    FIGSIZE   = (7, 5)
    FS_LABEL  = 10
    FS_TITLE  = 11
    FS_LEGEND = 8
    LW_MAIN   = 2.2
    LW_SEGURO = 1.3
    LW_PREC   = 1.1
    LOC_LEG1  = 'upper right'
    CMAP      = plt.cm.tab10

    dirigibles_unicos = list(dict.fromkeys(r['nombre'] for r in resultados))
    colores = CMAP(np.linspace(0, 0.9, len(dirigibles_unicos)))

    # Etiqueta y límites del eje X según flag NORMALIZAR_X
    if NORMALIZAR_X:
        xlabel_x = r'$\dot{m}_{H_2} / m_{est}$ — Normalized H$_2$ combustion rate [kg/(s·t)]'
        xlim     = (0.001, 200.0)
        titulo_x = 'normalized H$_2$ combustion rate'
    else:
        xlabel_x = r'$\dot{m}_{H_2}$ — H$_2$ combustion rate [kg/s]'
        xlim     = (0.8, 120.0)
        titulo_x = 'H$_2$ combustion rate'

    def _x_vals(sub, diri):
        tasas = [float(r['max_flujo_H2']) for r in sub]
        if NORMALIZAR_X:
            info  = next(d for d in DIRIGIBLES if d['nombre'] == diri)
            m_t   = info['m_estructura'] / 1000.0
            return [t / m_t for t in tasas]
        return tasas

    # ── Gráfica 1: Desplazamiento ──────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=FIGSIZE)

    # Franjas tecnológicas (fondo)
    for x0, x1, col, alpha, lbl in FRANJAS_TECH:
        x0_plot = x0 if NORMALIZAR_X else x0
        ax1.axvspan(x0_plot, x1, color=col, alpha=alpha, zorder=0)

    for diri, color in zip(dirigibles_unicos, colores):
        sub = sorted([r for r in resultados if r['nombre'] == diri],
                     key=lambda r: float(r['max_flujo_H2']))
        h_d = [r['h_desp_m'] if r['h_desp_m'] is not None else float('nan')
               for r in sub]
        ax1.plot(_x_vals(sub, diri), h_d, 'o-', color=color,
                 linewidth=LW_MAIN, markersize=5, label=diri, zorder=3)

    ax1.axhline(0.0, color='gray',    ls='-',  lw=0.8,      zorder=1)
    ax1.axhline(5.0, color='black',   ls='--', lw=LW_SEGURO, zorder=4,
                label='Safety threshold = 5.0 m')
    ax1.axhline(1.0, color='dimgray', ls=':',  lw=LW_PREC,   zorder=4,
                label='Precision threshold = 1.0 m')
    ax1.set_xscale('log')
    ax1.set_yscale('symlog', linthresh=0.1)   # linthresh=0.1: región lineal ±0.1 m
    ax1.set_xlim(*xlim)
    ax1.set_ylim(-0.15, 500)                  # recorta espacio negativo vacío
    ax1.grid(True, which='major', alpha=0.35)
    ax1.grid(True, which='minor', alpha=0.15)
    ax1.set_xlabel(xlabel_x, fontsize=FS_LABEL)
    ax1.set_ylabel(r'Max. vertical displacement $|\Delta h|$ [m]', fontsize=FS_LABEL)
    ax1.set_title(f'Vertical displacement vs {titulo_x}\n(with O₂ mass absorption)',
                  fontsize=FS_TITLE)

    # Leyenda: curvas + umbrales en primera leyenda, franjas en segunda
    leg_curvas = ax1.legend(fontsize=FS_LEGEND, loc='upper right')
    ax1.add_artist(leg_curvas)
    parches = [mpatches.Patch(color=col, alpha=alpha + 0.2, label=lbl)
               for _, _, col, alpha, lbl in FRANJAS_TECH]
    ax1.legend(handles=parches, fontsize=FS_LEGEND, loc='lower left',
               title='Technology')
    plt.tight_layout()
    f1 = os.path.join(OUTDIR, f't3_o2_desplazamiento_{timestamp}.png')
    fig1.savefig(f1, dpi=150, bbox_inches='tight')
    plt.close(fig1)
    print(f"  Gráfica 1 guardada: {f1}")

    # ── Gráfica 2: Tiempo de estabilización ───────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=FIGSIZE)

    # Franjas tecnológicas (fondo)
    for x0, x1, col, alpha, lbl in FRANJAS_TECH:
        ax2.axvspan(x0, x1, color=col, alpha=alpha, zorder=0)

    for diri, color in zip(dirigibles_unicos, colores):
        sub = sorted([r for r in resultados if r['nombre'] == diri],
                     key=lambda r: float(r['max_flujo_H2']))
        # Timeout → NaN → hueco en la línea (no contamina escala Y)
        t_e = []
        for r in sub:
            val = r['t_estab_s']
            if isinstance(val, (int, float)) and not np.isnan(val) and val < 250.0:
                t_e.append(val)
            else:
                t_e.append(float('nan'))
        t_masked = np.ma.masked_invalid(t_e)
        orden_z = 2 if diri == 'LZ 129 Hindenburg' else 3
        ax2.plot(_x_vals(sub, diri), t_masked, 'o-', color=color,
                 linewidth=LW_MAIN, markersize=5, label=diri, zorder=orden_z)

    ax2.set_xscale('log')
    ax2.set_xlim(*xlim)
    ax2.set_ylim(bottom=0)
    ax2.grid(True, which='major', alpha=0.35)
    ax2.grid(True, which='minor', alpha=0.15)
    ax2.set_xlabel(xlabel_x, fontsize=FS_LABEL)
    ax2.set_ylabel('Stabilization time [s]', fontsize=FS_LABEL)
    ax2.set_title(f'Stabilization time vs {titulo_x}\n(with O₂ mass absorption)',
                  fontsize=FS_TITLE)

    # Leyenda: curvas arriba-dcha, franjas abajo-izq
    leg_curvas2 = ax2.legend(fontsize=FS_LEGEND, loc=LOC_LEG1)
    ax2.add_artist(leg_curvas2)
    parches2 = [mpatches.Patch(color=col, alpha=alpha + 0.2, label=lbl)
                for _, _, col, alpha, lbl in FRANJAS_TECH]
    ax2.legend(handles=parches2, fontsize=FS_LEGEND, loc='upper right',
               title='Technology', bbox_to_anchor=(1.0, 0.62))
    plt.tight_layout()
    f2 = os.path.join(OUTDIR, f't3_o2_tiempo_{timestamp}.png')
    fig2.savefig(f2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"  Gráfica 2 guardada: {f2}")

    # ── Gráfica 3: Agua acumulada al final de la simulación ───────────────────
    fig3, ax3 = plt.subplots(figsize=FIGSIZE)

    for diri, color in zip(dirigibles_unicos, colores):
        sub  = sorted([r for r in resultados if r['nombre'] == diri],
                      key=lambda r: float(r['max_flujo_H2']))
        m_ag = [r['m_agua_kg'] if r['m_agua_kg'] is not None else float('nan')
                for r in sub]
        ax3.plot(_x_vals(sub, diri), m_ag, 'o-', color=color,
                 linewidth=LW_MAIN, markersize=5, label=diri, zorder=3)

    ax3.set_xscale('log')
    ax3.set_xlim(*xlim)
    ax3.grid(True, which='major', alpha=0.35)
    ax3.grid(True, which='minor', alpha=0.15)
    ax3.set_xlabel(xlabel_x, fontsize=FS_LABEL)
    ax3.set_ylabel('Water ballast accumulated [kg]', fontsize=FS_LABEL)
    ax3.set_title(f'H₂O ballast vs {titulo_x}\n(O₂ absorbed from atmosphere)',
                  fontsize=FS_TITLE)
    ax3.legend(fontsize=FS_LEGEND, loc='upper left')
    plt.tight_layout()
    f3 = os.path.join(OUTDIR, f't3_o2_agua_{timestamp}.png')
    fig3.savefig(f3, dpi=150, bbox_inches='tight')
    plt.close(fig3)
    print(f"  Gráfica 3 guardada: {f3}")

# =============================================================================
if __name__ == '__main__':
    ejecutar_barrido()
