# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# analisis_optimos.py  (sim_v2)
# Extrae tasas mínimas y óptimas del barrido paramétrico.
#
# Uso:
#   python analisis_optimos.py
#
# Entrada (busca automáticamente el CSV más reciente en la carpeta):
#   resultados_barrido_*.csv   — del barrido paramétrico
#   ciclo_condensacion_*.csv   — del balance del ciclo
#
# Salidas:
#   tabla_minimos_<ts>.csv         — tabla central del paper
#   grafica_minimos_<ts>.png       — curvas con umbrales anotados
#   grafica_condensacion_<ts>.png  — Q_cond vs t_cond con umbral
# =============================================================================

import csv
import glob
import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# CRITERIOS
# =============================================================================

H_DESP_SEGURO   = 5.0   # [m]  |Δh| máximo aceptable (umbral de seguridad)
H_DESP_OPTIMO   = 1.0   # [m]  |Δh| donde los rendimientos son decrecientes
T_COND_MAX      = 10.0  # [min] tiempo máximo de condensación operacional

# Rangos tecnológicos de ṁ_H2_max disponibles hoy [kg/s]
TECNOLOGIAS = [
    ('Celda PEM pequeña',      0.01,  0.1,  'lightblue'),
    ('Motor combustión H₂',    0.5,   5.0,  'lightyellow'),
    ('Sistema industrial',     5.0,   50.0, 'lightgreen'),
    ('Límite teórico control', 50.0,  800.0,'lightsalmon'),
]


# =============================================================================
# UTILIDADES
# =============================================================================

def csv_mas_reciente(patron):
    archivos = glob.glob(os.path.join(os.path.dirname(__file__), patron))
    if not archivos:
        raise FileNotFoundError(f"No se encontró: {patron}")
    return max(archivos, key=os.path.getmtime)


def leer_csv(path):
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def primera_tasa_bajo(subset, columna, umbral):
    """Retorna la primera tasa donde columna <= umbral. None si ninguna lo cumple."""
    for fila in subset:
        val = fila.get(columna)
        if val is None or val == 'N/A':
            continue
        try:
            if float(val) <= umbral:
                return float(fila['max_flujo_H2'])
        except ValueError:
            continue
    return None


def tecnologia_para(tasa):
    """Clasifica una tasa en la categoría tecnológica disponible."""
    if tasa is None:
        return 'No alcanzable con tecnología actual'
    for nombre, lo, hi, _ in TECNOLOGIAS:
        if lo <= tasa <= hi:
            return nombre
    if tasa < TECNOLOGIAS[0][1]:
        return 'Celda PEM pequeña'
    return 'Límite teórico control'


# =============================================================================
# ANÁLISIS PRINCIPAL
# =============================================================================

def analizar():
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')

    # --- Cargar datos ---
    path_barrido = csv_mas_reciente('resultados_barrido_*.csv')
    path_cond    = csv_mas_reciente('ciclo_condensacion_*.csv')
    print(f"  Barrido : {os.path.basename(path_barrido)}")
    print(f"  Cond.   : {os.path.basename(path_cond)}\n")

    barrido = leer_csv(path_barrido)
    cond    = leer_csv(path_cond)

    dirigibles = list(dict.fromkeys(r['nombre'] for r in barrido))

    # ==========================================================================
    # PARTE 1 — Tasas mínimas y óptimas por dirigible
    # ==========================================================================

    print(f"{'='*70}")
    print("  TASAS MÍNIMAS Y ÓPTIMAS POR DIRIGIBLE")
    print(f"  Criterio seguro : |Δh| ≤ {H_DESP_SEGURO} m")
    print(f"  Criterio óptimo : |Δh| ≤ {H_DESP_OPTIMO} m")
    print(f"{'='*70}\n")

    tabla_minimos = []

    for nombre in dirigibles:
        subset = [r for r in barrido if r['nombre'] == nombre]
        # Ordenar por tasa ascendente
        subset.sort(key=lambda r: float(r['max_flujo_H2']))

        tasa_min  = primera_tasa_bajo(subset, 'h_desp_m', H_DESP_SEGURO)
        tasa_opt  = primera_tasa_bajo(subset, 'h_desp_m', H_DESP_OPTIMO)
        tech_min  = tecnologia_para(tasa_min)
        tech_opt  = tecnologia_para(tasa_opt)

        # h_desp y t_estab en tasa_min y tasa_opt
        def get_val(tasa, col):
            for r in subset:
                if float(r['max_flujo_H2']) == tasa:
                    v = r.get(col, 'N/A')
                    return v if v != 'N/A' else 'N/A'
            return 'N/A'

        h_en_min = get_val(tasa_min, 'h_desp_m')  if tasa_min else 'N/A'
        h_en_opt = get_val(tasa_opt, 'h_desp_m')  if tasa_opt else 'N/A'
        t_en_min = get_val(tasa_min, 't_estab_s')  if tasa_min else 'N/A'
        t_en_opt = get_val(tasa_opt, 't_estab_s')  if tasa_opt else 'N/A'

        fila = {
            'dirigible'       : nombre,
            'tasa_min_kgs'    : tasa_min  if tasa_min  else '>800',
            'h_desp_en_min_m' : h_en_min,
            't_estab_en_min_s': t_en_min,
            'tasa_opt_kgs'    : tasa_opt  if tasa_opt  else '>800',
            'h_desp_en_opt_m' : h_en_opt,
            't_estab_en_opt_s': t_en_opt,
            'tecnologia_min'  : tech_min,
            'tecnologia_opt'  : tech_opt,
        }
        tabla_minimos.append(fila)

        print(f"  ── {nombre} ──")
        print(f"     ṁ_H2_min  = {fila['tasa_min_kgs']:>6} kg/s  "
              f"→ |Δh|={h_en_min} m, t_estab={t_en_min} s")
        print(f"     ṁ_H2_opt  = {fila['tasa_opt_kgs']:>6} kg/s  "
              f"→ |Δh|={h_en_opt} m, t_estab={t_en_opt} s")
        print(f"     Tecnología (mínimo) : {tech_min}")
        print(f"     Tecnología (óptimo) : {tech_opt}")
        print()

    # ==========================================================================
    # PARTE 2 — Q_cond mínimo
    # ==========================================================================

    print(f"{'='*70}")
    print(f"  Q_COND MÍNIMO (t_cond ≤ {T_COND_MAX} min)")
    print(f"{'='*70}\n")

    # Calcular exacto: t_cond = m_H2O * L_v / Q_cond
    
    m_H2O_ref = 667.9
    L_v       = 2257e3
    Q_cond_min_W  = (m_H2O_ref * L_v) / (T_COND_MAX * 60)
    Q_cond_min_kW = Q_cond_min_W / 1e3

    print(f"  m_H2O referencia  = {m_H2O_ref} kg")
    print(f"  Q_cond_min exacto = {Q_cond_min_kW:.0f} kW")
    print(f"  → Con {Q_cond_min_kW:.0f} kW, la condensación tarda exactamente {T_COND_MAX} min\n")

    # Tabla de condensación
    print(f"  {'Q_cond [kW]':>12}  {'t_cond [min]':>14}  {'Operacional?':>14}")
    for fila in cond:
        qc  = float(fila['Q_cond_kW'])
        tc  = float(fila['t_cond_min'])
        ok  = 'Sí ✓' if tc <= T_COND_MAX else 'No  ✗'
        mark = ' ← UMBRAL' if abs(qc - round(Q_cond_min_kW, -2)) < 200 else ''
        print(f"  {qc:>12.0f}  {tc:>14.1f}  {ok:>14}{mark}")
    print()

    # ==========================================================================
    # PARTE 3 — Guardar CSV tabla de mínimos
    # ==========================================================================

    path_csv = f'tabla_minimos_{ts}.csv'
    with open(path_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(tabla_minimos[0].keys()))
        writer.writeheader()
        writer.writerows(tabla_minimos)
    print(f"  CSV guardado: {path_csv}\n")

    # ==========================================================================
    # PARTE 4 — Gráficas
    # ==========================================================================

    _grafica_minimos(barrido, dirigibles, tabla_minimos, ts)
    _grafica_condensacion(cond, Q_cond_min_kW, ts)


# =============================================================================
# GRÁFICA 1 — Curvas de barrido con umbrales anotados
# =============================================================================

def _grafica_minimos(barrido, dirigibles, tabla_minimos, ts):
    colores = plt.cm.tab10(np.linspace(0, 0.6, len(dirigibles)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        'Análisis paramétrico: tasas mínimas y óptimas de combustión de H₂\n'
        'Maniobra de liberación de 1 tonelada de carga — 4 dirigibles representativos',
        fontsize=12
    )
    ax1, ax2 = axes

    for nombre, color in zip(dirigibles, colores):
        subset = sorted([r for r in barrido if r['nombre'] == nombre],
                        key=lambda r: float(r['max_flujo_H2']))
        tasas  = [float(r['max_flujo_H2']) for r in subset]
        h_desp = [float(r['h_desp_m'])     for r in subset]
        t_est  = [float(r['t_estab_s']) if r['t_estab_s'] not in ('N/A', '')
                  else float('nan') for r in subset]

        ax1.plot(tasas, h_desp, 'o-', color=color, linewidth=2, label=nombre)
        ax2.plot(tasas, t_est,  's-', color=color, linewidth=2, label=nombre)

    # Umbrales
    ax1.axhline(H_DESP_SEGURO, color='red',    linestyle='--', linewidth=1.5,
                label=f'Umbral seguro = {H_DESP_SEGURO} m')
    ax1.axhline(H_DESP_OPTIMO, color='orange', linestyle='--', linewidth=1.5,
                label=f'Umbral óptimo = {H_DESP_OPTIMO} m')

    # Anotar ṁ_H2_min por dirigible en ax1
    for fila, color in zip(tabla_minimos, colores):
        tasa_min = fila['tasa_min_kgs']
        if tasa_min == '>800':
            continue
        tasa_min = float(tasa_min)
        subset = sorted([r for r in barrido if r['nombre'] == fila['dirigible']],
                        key=lambda r: float(r['max_flujo_H2']))
        for r in subset:
            if float(r['max_flujo_H2']) == tasa_min:
                h = float(r['h_desp_m'])
                ax1.annotate(
                    f"  {tasa_min:.0f} kg/s",
                    xy=(tasa_min, h), fontsize=8, color=color,
                    xytext=(tasa_min * 1.3, h * 1.5),
                    arrowprops=dict(arrowstyle='->', color=color, lw=0.8)
                )
                break

    # Franjas tecnológicas en ax1
    for nombre_t, lo, hi, col in TECNOLOGIAS:
        ax1.axvspan(lo, hi, alpha=0.08, color=col)

    # Leyenda franjas
    parches = [mpatches.Patch(color=col, alpha=0.3, label=nombre_t)
               for nombre_t, lo, hi, col in TECNOLOGIAS]

    ax1.set_xscale('log')
    ax1.set_xlabel('ṁ_H2_max — Tasa máx. de combustión [kg/s]', fontsize=10)
    ax1.set_ylabel('Desplazamiento vertical máx. post-liberación |Δh| [m]', fontsize=10)
    ax1.set_title('Desplazamiento vertical vs tasa de combustión')
    ax1.legend(handles=ax1.get_legend_handles_labels()[0] + parches,
               labels=ax1.get_legend_handles_labels()[1] +
                      [t[0] for t in TECNOLOGIAS],
               fontsize=7, loc='upper right')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_ylim(bottom=0)

    ax2.set_xscale('log')
    ax2.set_xlabel('ṁ_H2_max — Tasa máx. de combustión [kg/s]', fontsize=10)
    ax2.set_ylabel('Tiempo de estabilización [s]\n(|Δh|<1 m, |v|<0.05 m/s)', fontsize=10)
    ax2.set_title('Tiempo de estabilización vs tasa de combustión')
    ax2.legend(fontsize=8)
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    nombre_png = f'grafica_minimos_{ts}.png'
    plt.savefig(nombre_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfica guardada: {nombre_png}")


# =============================================================================
# GRÁFICA 2 — Q_cond vs tiempo de condensación con umbral
# =============================================================================

def _grafica_condensacion(cond, Q_cond_min_kW, ts):
    qc_kw  = [float(f['Q_cond_kW'])  for f in cond]
    tc_min = [float(f['t_cond_min']) for f in cond]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(qc_kw, tc_min, 'o-', color='steelblue', linewidth=2, markersize=7)

    ax.axhline(T_COND_MAX, color='red', linestyle='--', linewidth=1.5,
               label=f'Umbral operacional = {T_COND_MAX} min')
    ax.axvline(Q_cond_min_kW, color='orange', linestyle='--', linewidth=1.5,
               label=f'Q_cond_min = {Q_cond_min_kW:.0f} kW')

    ax.fill_betweenx([0, T_COND_MAX], Q_cond_min_kW, max(qc_kw),
                     alpha=0.1, color='green', label='Zona operacional (< 10 min)')
    ax.fill_betweenx([T_COND_MAX, max(tc_min)*1.1], 0, Q_cond_min_kW,
                     alpha=0.08, color='red', label='Zona no operacional')

    ax.annotate(
        f'Q_cond_min\n= {Q_cond_min_kW:.0f} kW',
        xy=(Q_cond_min_kW, T_COND_MAX),
        xytext=(Q_cond_min_kW * 0.4, T_COND_MAX * 2.5),
        fontsize=9,
        arrowprops=dict(arrowstyle='->', color='orange', lw=1.2),
        color='darkorange'
    )

    ax.set_xscale('log')
    ax.set_xlabel('Potencia de condensación Q_cond [kW]', fontsize=11)
    ax.set_ylabel('Tiempo de condensación [min]', fontsize=11)
    ax.set_title('Tiempo de condensación del vapor de agua\n'
                 f'(m_H₂O = 667.9 kg, L_v = 2257 kJ/kg, m_carga = 1000 kg)',
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    nombre_png = f'grafica_condensacion_{ts}.png'
    plt.savefig(nombre_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfica guardada: {nombre_png}")


# =============================================================================
# Ejecución directa
# =============================================================================

if __name__ == '__main__':
    analizar()
