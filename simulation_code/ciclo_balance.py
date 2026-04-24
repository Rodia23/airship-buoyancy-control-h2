# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# ciclo_balance.py  (sim_v2)
# Balance del ciclo regenerativo de H2
#
# Calcula y grafica:
#   1. Masa de H2 consumida por maniobra y agua producida
#   2. Tiempo de condensación en función de Q_cond
#   3. Masa de H2 regenerada por día (solar + electrólisis)
#   4. Autonomía de ciclo: entregas por día antes de agotar el H2
#
# Uso:
#   python ciclo_balance.py
#
# Salidas:
#   ciclo_balance_<timestamp>.png
#   ciclo_balance_<timestamp>.csv
# =============================================================================

import csv
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import parametros as p


# =============================================================================
# FUNCIONES DE BALANCE
# =============================================================================

def balance_maniobra(m_carga, Q_comb_H2=p.Q_comb_H2, rend_gen=p.rend_gen):
    """
    Calcula la masa de H2 necesaria para compensar la liberación de m_carga [kg].

    Principio: la energía de combustión del H2 (a rend_gen de eficiencia)
    genera la electricidad para mantener el sistema. La masa de H2 quemada
    se estima por balance de flotabilidad: el H2 consumido ≈ m_carga
    (orden de magnitud), ajustado por la densidad relativa del H2 vs aire.

    Para el artículo: se usa la aproximación m_H2_quemado ≈ m_carga × f_masa,
    donde f_masa es la fracción de masa de H2 necesaria para compensar
    la pérdida de empuje debida a la carga liberada.

    Retorna
    -------
    dict con:
        m_H2_quemado  [kg]   masa de H2 consumida en la maniobra
        m_H2O_prod    [kg]   agua producida (9 × m_H2_quemado, reacción H2+O2→H2O)
        E_generada    [kWh]  energía eléctrica generada
    """
    # Aproximación: m_H2 quemado ~ m_carga × (rho_H2 / rho_aire)
    # Factor típico para p atmos ~85 kPa, T~281 K: rho_H2 ≈ 0.073 kg/m³
    # rho_aire ≈ 1.058 → ratio ≈ 0.069
    # Por balance de empuje: m_H2 ≈ m_carga × (rho_H2 / (rho_aire - rho_H2))
    rho_H2_bal = p.P_atmos / (p.R_H2 * p.T1)
    rho_aire   = p.rho_aire_std
    f_masa     = rho_H2_bal / (rho_aire - rho_H2_bal)
    m_H2_quemado = m_carga * f_masa

    # Agua producida: H2 + ½O2 → H2O; relación másica H2O/H2 = 18/2 = 9
    m_H2O_prod = 9.0 * m_H2_quemado

    # Energía eléctrica generada
    E_generada = (m_H2_quemado * Q_comb_H2 * rend_gen) / 3.6e6  # [kWh]

    return {
        'm_H2_quemado' : round(m_H2_quemado, 3),
        'm_H2O_prod'   : round(m_H2O_prod, 3),
        'E_generada'   : round(E_generada, 3),
    }


def tiempo_condensacion(m_H2_quemado, Q_cond_W, L_v=p.L_v_agua):
    """
    Calcula el tiempo para condensar el vapor de agua producido.

    t_cond = (m_H2O × L_v) / Q_cond

    Parámetros
    ----------
    m_H2_quemado : float  [kg]  H2 quemado en la maniobra
    Q_cond_W     : float  [W]   potencia de condensación del intercambiador
    L_v          : float  [J/kg] calor latente de vaporización del agua

    Retorna
    -------
    t_cond_s  [s]   tiempo de condensación
    t_cond_min [min]
    """
    m_H2O    = 9.0 * m_H2_quemado
    t_cond_s = (m_H2O * L_v) / Q_cond_W
    return t_cond_s, t_cond_s / 60.0


def regeneracion_solar(Area_total_m2,
                       rendimiento_efectivo=p.rendimiento_efectivo,
                       A_paneles_fraccion=p.A_paneles_fraccion,
                       eta_electrolisis=p.eta_electrolisis,
                       E_H2_electrolisis=p.E_H2_electrolisis):
    """
    Calcula la masa de H2 regenerada por día mediante paneles solares y electrólisis.

    H2_regen [kg/día] = (E_solar [kWh/día] × η_electrolisis) / E_H2 [kWh/kg]

    donde:
      E_solar = rendimiento_efectivo × A_paneles
      A_paneles = A_paneles_fraccion × Area_total_m2

    Parámetros
    ----------
    Area_total_m2 : float  [m²]  superficie total del dirigible (elipsoide)

    Retorna
    -------
    dict con: A_paneles, E_solar_kWh, H2_regen_kg_dia
    """
    A_paneles    = A_paneles_fraccion * Area_total_m2
    E_solar_kWh  = rendimiento_efectivo * A_paneles           # [kWh/día]
    H2_regen     = (E_solar_kWh * eta_electrolisis) / E_H2_electrolisis  # [kg/día]

    return {
        'A_paneles_m2'    : round(A_paneles, 1),
        'E_solar_kWh_dia' : round(E_solar_kWh, 1),
        'H2_regen_kg_dia' : round(H2_regen, 3),
    }


def autonomia_ciclo(m_H2_total, m_H2_quemado_por_entrega, H2_regen_kg_dia):
    """
    Calcula cuántas entregas de 1 tonelada puede hacer el dirigible por día
    antes de agotar el H2, considerando regeneración solar continua.

    Retorna
    -------
    dict con:
        entregas_sin_regen  → entregas posibles con stock inicial, sin regenerar
        entregas_con_regen  → entregas posibles con regeneración durante el día (8 h operativas)
        dias_ciclo_completo → días para recuperar el H2 consumido en 1 entrega
    """
    entregas_sin_regen = m_H2_total / m_H2_quemado_por_entrega if m_H2_quemado_por_entrega > 0 else 0

    # Asumimos 8 horas operativas por día → regen durante el vuelo
    H2_regen_8h  = H2_regen_kg_dia * (8.0 / 24.0)
    H2_neto_por_entrega = m_H2_quemado_por_entrega - H2_regen_8h
    if H2_neto_por_entrega <= 0:
        entregas_con_regen = float('inf')  # totalmente sostenible
    else:
        entregas_con_regen = m_H2_total / H2_neto_por_entrega

    dias_ciclo = m_H2_quemado_por_entrega / H2_regen_kg_dia if H2_regen_kg_dia > 0 else float('inf')

    return {
        'entregas_sin_regen'  : round(entregas_sin_regen, 1),
        'entregas_con_regen'  : round(entregas_con_regen, 1) if entregas_con_regen != float('inf') else 'inf',
        'dias_ciclo_completo' : round(dias_ciclo, 2),
    }


# =============================================================================
# ANÁLISIS COMPLETO
# =============================================================================

DIRIGIBLES = [
    {
        'nombre'        : 'Pathfinder 3',
        'Area_frontal'  : 4390,    # [m²] (para barrido dinámico)
        'Area_total'    : 13000,   # [m²] superficie elipsoide estimada
        'm_H2_total'    : 76.0,    # [kg] H2 disponible en ballonets para carga de 1t
        'm_carga'       : 1000,
    },
    {
        'nombre'        : 'Flying Whales LCA60T',
        'Area_frontal'  : 7854,
        'Area_total'    : 23000,
        'm_H2_total'    : 76.0,
        'm_carga'       : 1000,
    },
    {
        'nombre'        : 'Zeppelin NT',
        'Area_frontal'  : 837,
        'Area_total'    : 2800,
        'm_H2_total'    : 76.0,
        'm_carga'       : 1000,
    },
    {
        'nombre'        : 'LZ 129 Hindenburg',
        'Area_frontal'  : 7950,
        'Area_total'    : 24000,
        'm_H2_total'    : 76.0,
        'm_carga'       : 1000,
    },
]

# Rango de Q_cond para el análisis paramétrico [W]
Q_COND_RANGO = [10e3, 50e3, 100e3, 200e3, 500e3, 1000e3, 2000e3, 5000e3]


def analizar():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filas_csv = []

    print(f"\n{'='*70}")
    print("  BALANCE DEL CICLO REGENERATIVO — Dirigible H2")
    print(f"  rend_gen={p.rend_gen*100:.0f}%  |  eta_elec={p.eta_electrolisis*100:.0f}%  |  "
          f"E_H2={p.E_H2_electrolisis} kWh/kg  |  r_ef={p.rendimiento_efectivo} kWh/m2/dia")
    print(f"{'='*70}\n")

    for diri in DIRIGIBLES:
        nombre   = diri['nombre']
        m_carga  = diri['m_carga']
        m_H2_tot = diri['m_H2_total']
        A_total  = diri['Area_total']

        bal   = balance_maniobra(m_carga)
        regen = regeneracion_solar(A_total)
        aut   = autonomia_ciclo(m_H2_tot, bal['m_H2_quemado'], regen['H2_regen_kg_dia'])

        print(f"  ── {nombre} ──")
        print(f"     Maniobra (m_carga={m_carga} kg):")
        print(f"       H2 quemado     : {bal['m_H2_quemado']:>8.3f} kg")
        print(f"       Agua producida : {bal['m_H2O_prod']:>8.3f} kg")
        print(f"       Energía gen.   : {bal['E_generada']:>8.3f} kWh")
        print(f"     Regeneración solar (A_total={A_total} m²):")
        print(f"       Área paneles   : {regen['A_paneles_m2']:>8.1f} m²")
        print(f"       Energía solar  : {regen['E_solar_kWh_dia']:>8.1f} kWh/día")
        print(f"       H2 regenerado  : {regen['H2_regen_kg_dia']:>8.3f} kg/día")
        print(f"     Autonomía de ciclo (m_H2_stock={m_H2_tot} kg):")
        print(f"       Sin regen.     : {aut['entregas_sin_regen']:>8.1f} entregas")
        print(f"       Con regen.     : {aut['entregas_con_regen']} entregas/día (8 h operativas)")
        print(f"       Días ciclo     : {aut['dias_ciclo_completo']:>8.2f} días para recuperar 1 entrega")
        print()

        fila = {'dirigible': nombre, **bal, **regen, **aut}
        filas_csv.append(fila)

    # =========================================================================
    # Tabla de condensación (Q_cond paramétrico)
    # =========================================================================
    print(f"  ── Tiempo de condensación vs Q_cond (m_carga=1000 kg) ──")
    print(f"  {'Q_cond [kW]':>12}  {'t_cond [min]':>14}  {'t_cond [h]':>12}")
    bal_ref = balance_maniobra(1000)
    filas_cond = []
    for qc in Q_COND_RANGO:
        t_s, t_min = tiempo_condensacion(bal_ref['m_H2_quemado'], qc)
        print(f"  {qc/1e3:>12.0f}  {t_min:>14.1f}  {t_min/60:>12.2f}")
        filas_cond.append({'Q_cond_kW': qc/1e3, 't_cond_min': round(t_min, 1),
                           't_cond_h': round(t_min/60, 3)})
    print()

    # =========================================================================
    # CSV
    # =========================================================================
    csv_bal  = f'ciclo_balance_{timestamp}.csv'
    csv_cond = f'ciclo_condensacion_{timestamp}.csv'

    with open(csv_bal, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(filas_csv[0].keys()))
        writer.writeheader()
        writer.writerows(filas_csv)

    with open(csv_cond, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Q_cond_kW', 't_cond_min', 't_cond_h'])
        writer.writeheader()
        writer.writerows(filas_cond)

    print(f"  CSVs guardados: {csv_bal}, {csv_cond}")

    # =========================================================================
    # Gráficas
    # =========================================================================
    _graficar(filas_csv, filas_cond, bal_ref, timestamp)


def _graficar(filas_csv, filas_cond, bal_ref, timestamp):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle('Balance del ciclo regenerativo de H₂\n'
                 'Liberación de 1 tonelada de carga útil', fontsize=12)

    # --- Panel 1: Tiempo de condensación vs Q_cond ---
    ax = axes[0]
    qc_kw  = [f['Q_cond_kW']    for f in filas_cond]
    tc_min = [f['t_cond_min']   for f in filas_cond]
    ax.plot(qc_kw, tc_min, 'o-', color='steelblue', linewidth=2)
    ax.axhline(60, color='gray', linestyle='--', linewidth=1, label='1 hora')
    ax.axhline(10, color='green', linestyle='--', linewidth=1, label='10 min')
    ax.set_xscale('log')
    ax.set_xlabel('Potencia de condensación Q_cond [kW]')
    ax.set_ylabel('Tiempo de condensación [min]')
    ax.set_title('Condensación del vapor de agua')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.4)

    # Anotar punto de referencia (500 kW)
    for f in filas_cond:
        if f['Q_cond_kW'] == 500:
            ax.annotate(f"500 kW\n{f['t_cond_min']} min",
                        xy=(500, f['t_cond_min']),
                        xytext=(300, f['t_cond_min']*1.4),
                        fontsize=8, arrowprops=dict(arrowstyle='->', color='gray'))

    # --- Panel 2: H2 regenerado vs área de paneles ---
    ax = axes[1]
    areas = np.linspace(500, 15000, 100)
    escenarios = [
        ('Desfavorable (1.2 kWh/m²/día)', 1.2, 'red'),
        ('Medio (2.0 kWh/m²/día)',         2.0, 'steelblue'),
        ('Favorable (2.8 kWh/m²/día)',     2.8, 'green'),
    ]
    for label, r_ef, color in escenarios:
        E = r_ef * areas * p.A_paneles_fraccion
        H2 = (E * p.eta_electrolisis) / p.E_H2_electrolisis
        ax.plot(areas, H2, color=color, linewidth=1.8, label=label)

    ax.axhline(bal_ref['m_H2_quemado'], color='black', linestyle='--', linewidth=1.5,
               label=f"H₂ consumido/entrega ({bal_ref['m_H2_quemado']:.2f} kg)")
    ax.set_xlabel('Área superficial total del dirigible [m²]')
    ax.set_ylabel('H₂ regenerado por día [kg]')
    ax.set_title('Regeneración solar de H₂\n(fracción paneles = '
                 f'{p.A_paneles_fraccion*100:.0f}% del área total)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

    # --- Panel 3: Autonomía de ciclo por dirigible ---
    ax = axes[2]
    nombres = [f['dirigible'] for f in filas_csv]
    dias    = [f['dias_ciclo_completo'] for f in filas_csv]
    colores = ['steelblue', 'darkorange', 'green', 'red']
    bars    = ax.barh(nombres, dias, color=colores[:len(nombres)], height=0.5)
    ax.axvline(1.0, color='black', linestyle='--', linewidth=1, label='1 día')
    for bar, d in zip(bars, dias):
        ax.text(d + 0.05, bar.get_y() + bar.get_height()/2,
                f'{d:.1f} d', va='center', fontsize=9)
    ax.set_xlabel('Días para recuperar el H₂ de 1 entrega (regeneración solar)')
    ax.set_title('Autonomía de ciclo por dirigible')
    ax.legend(fontsize=9)
    ax.grid(True, axis='x', alpha=0.4)

    plt.tight_layout()
    nombre_png = f'ciclo_balance_{timestamp}.png'
    plt.savefig(nombre_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Gráfica guardada: {nombre_png}")


# =============================================================================
# Ejecución directa
# =============================================================================

if __name__ == '__main__':
    analizar()
