# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# comparacion_O2.py  (sim_v2)
#
# Comparación completa de dos escenarios de fuente de O₂ en el ciclo
# regenerativo de H₂ para dirigibles:
#
#   Escenario A — O₂ de tanques a bordo (ciclo de masa cerrado)
#   Escenario B — O₂ extraído del aire   (ciclo de masa abierto en O₂)
#
# Análisis por fases:
#   0. Equilibrio inicial (a altitud h0)
#   1. Combustión de m_H2_man kg de H₂
#   2. Liberación de carga (m_carga)
#   3. Condensación del vapor (H₂O vapor → líquido)
#   4. Electrólisis + retorno a estado inicial
#
# Salidas:
#   comparacion_O2_<timestamp>.png
#   comparacion_O2_<timestamp>.csv
#
# Uso:
#   python comparacion_O2.py
# =============================================================================

import csv
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import parametros as p


# =============================================================================
# CONSTANTES ADICIONALES
# =============================================================================

M_O2        = 0.032          # [kg/mol]  Masa molar del O₂
M_H2O       = 0.018          # [kg/mol]  Masa molar del H₂O
R_O2        = 8.314 / M_O2  # [J/kg·K]  ≈ 259.8
R_H2O       = 8.314 / M_H2O # [J/kg·K]  ≈ 461.9

# Densidades de referencia  —  a presión de operación (P_atmos) y T1
rho_H2_op   = p.P_atmos / (p.R_H2 * p.T1)          # [kg/m³] H₂ gaseoso en ballonets

# O₂ comprimido en tanques (350 bar, 25°C) — gas ideal, factor Z≈0.92
T_tanque    = 298.0          # [K]  temperatura tanque
Z_O2        = 0.92           # factor compresibilidad O₂ a 350 bar
rho_O2_tank = p.P_tanque / (Z_O2 * R_O2 * T_tanque)  # [kg/m³]

# H₂O vapor — 100°C, 1 atm (101325 Pa)
T_vapor     = 373.0          # [K]
P_vapor     = 101325.0       # [Pa]
rho_vapor   = P_vapor / (R_H2O * T_vapor)            # [kg/m³]  ≈ 0.587

# H₂O líquida
rho_liq     = 1000.0         # [kg/m³]

# Relación másica O₂/H₂ en la combustión: 2H₂ + O₂ → 2H₂O
#   m_O2 = 8 × m_H2    (32/2×2 = 8)
#   m_H2O = 9 × m_H2
RATIO_O2_H2 = 8.0
RATIO_H2O   = 9.0


# =============================================================================
# FUNCIONES CORE
# =============================================================================

def f_masa_equilibrio():
    """
    Fracción másica para equilibrio de flotabilidad:
        m_H2 = m_otros × ρ_H2 / (ρ_aire − ρ_H2)
    """
    return rho_H2_op / (p.rho_aire_std - rho_H2_op)


def m_H2_min_escenario_A(m_carga):
    """
    H₂ mínimo para déficit de flotabilidad = m_carga (solo efecto volumen H₂).

    Déficit_A = (ρ_aire − ρ_H2) × V_H2 × g = [(ρ_aire − ρ_H2)/ρ_H2] × m_H2 × g
    Requerimiento: Déficit_A ≥ m_carga × g
    """
    factor_A = (p.rho_aire_std - rho_H2_op) / rho_H2_op
    return m_carga / factor_A


def m_H2_min_escenario_B(m_carga):
    """
    H₂ mínimo para déficit de flotabilidad = m_carga (volumen H₂ + masa O₂ absorbida).

    Déficit_B = [(ρ_aire − ρ_H2)/ρ_H2 + 8] × m_H2 × g
    El término +8 viene del O₂ absorbido del aire (8 × m_H2).
    """
    factor_B = (p.rho_aire_std - rho_H2_op) / rho_H2_op + RATIO_O2_H2
    return m_carga / factor_B


def m_H2_penalizacion_ascenso_B(m_H2_man_B):
    """
    H₂ extra que debe existir en gas cells durante el ascenso post-entrega
    para soportar el peso del O₂ aún a bordo (antes de ventilar).

    m_H2_extra = m_O2 × f_masa_equilibrio()
    """
    m_O2 = RATIO_O2_H2 * m_H2_man_B
    return m_O2 * f_masa_equilibrio()


# =============================================================================
# ANÁLISIS FASE A FASE
# =============================================================================

def ciclo_completo(m_carga, m_H2_man_A, m_H2_man_B, m_estructura):
    """
    Devuelve tabla de fases para ambos escenarios.

    Cada fase: (label, m_total_A, lift_A, m_total_B, lift_B)
    lift = F_Archimedes/g - m_total  [kg equivalentes]
    """
    m_O2_man  = RATIO_O2_H2 * m_H2_man_B  # O₂ para el Escenario B

    # m_H2_ini: suficiente para flotar m_total_ini a h0
    fma = f_masa_equilibrio()
    m_total_ini = m_estructura + m_carga

    # H₂ inicial (equilibrio exacto antes de maniobra)
    m_H2_ini = m_total_ini * fma

    # Función de flotabilidad neta [kg equivalentes]
    # F_net/g = ρ_aire×V_H2 − m_total = m_H2×(ρ_aire/ρ_H2 − 1) − m_otros
    def lift(m_H2, m_total):
        return m_H2 * (p.rho_aire_std / rho_H2_op - 1) - (m_total - m_H2)

    fases = []

    # Fase 0: Equilibrio inicial
    # H₂ ini compartido por ambos (se calcula aparte; aquí usamos el real)
    m_H2_0 = m_H2_ini
    # A: O₂ ya a bordo → m_total incluye O₂ de tanques (pero lo compensó el diseño)
    # Simplificamos: ambos parten con lift≈0
    fases.append({
        'fase'         : '0. Equilibrio inicial',
        'm_total_A'    : round(m_total_ini, 1),
        'lift_A_kg'    : 0.0,
        'm_total_B'    : round(m_total_ini, 1),
        'lift_B_kg'    : 0.0,
        'nota'         : 'F_neta = 0 (diseñado para equilibrio)',
    })

    # Fase 1: Quema de H₂
    # A: m_total sin cambio (O₂ ya estaba contado), solo cae V_H2
    m_total_A1 = m_total_ini
    m_H2_A1    = m_H2_ini - m_H2_man_A
    lift_A1    = lift(m_H2_A1, m_total_A1)

    # B: m_total sube por O₂ absorbido del aire
    m_total_B1 = m_total_ini + m_O2_man
    m_H2_B1    = m_H2_ini - m_H2_man_B
    lift_B1    = lift(m_H2_B1, m_total_B1)

    fases.append({
        'fase'         : '1. Combustión (quema H₂)',
        'm_total_A'    : round(m_total_A1, 1),
        'lift_A_kg'    : round(lift_A1, 2),
        'm_total_B'    : round(m_total_B1, 1),
        'lift_B_kg'    : round(lift_B1, 2),
        'nota'         : f'A: -m_H2={m_H2_man_A:.2f} kg  |  B: -m_H2={m_H2_man_B:.2f} kg +{m_O2_man:.1f} kg O2 del aire',
    })

    # Fase 2: Liberación de carga
    m_total_A2 = m_total_A1 - m_carga
    lift_A2    = lift(m_H2_A1, m_total_A2)

    m_total_B2 = m_total_B1 - m_carga
    lift_B2    = lift(m_H2_B1, m_total_B2)

    fases.append({
        'fase'         : '2. Liberación de carga',
        'm_total_A'    : round(m_total_A2, 1),
        'lift_A_kg'    : round(lift_A2, 2),
        'm_total_B'    : round(m_total_B2, 1),
        'lift_B_kg'    : round(lift_B2, 2),
        'nota'         : f'-m_carga = -{m_carga} kg en ambos',
    })

    # Fase 3: Condensación (vapor → líquido)
    # Casco rígido: masa constante, ballonet de aire compensa volumen
    # → sin cambio adicional de flotabilidad
    fases.append({
        'fase'         : '3. Condensación (vapor→liquido)',
        'm_total_A'    : round(m_total_A2, 1),
        'lift_A_kg'    : round(lift_A2, 2),
        'm_total_B'    : round(m_total_B2, 1),
        'lift_B_kg'    : round(lift_B2, 2),
        'nota'         : 'Casco rigido: ballonet compensa; sin cambio neto de flotabilidad',
    })

    # Fase 4a: Electrólisis — H₂ recuperado (antes de ventilar O₂ en B)
    # H₂ recuperado = m_H2_man quemado (por estequiometría electrólisis)
    m_H2_A3 = m_H2_A1 + m_H2_man_A   # vuelve a m_H2_ini (aprox)
    m_H2_B3 = m_H2_B1 + m_H2_man_B   # ídem, pero O₂ aún a bordo en B

    lift_A3 = lift(m_H2_A3, m_total_A2)   # A: O₂ ya vuelve a tanques
    lift_B3 = lift(m_H2_B3, m_total_B2)   # B: O₂ aún a bordo

    fases.append({
        'fase'         : '4a. Electrolisis (H₂ restaurado)',
        'm_total_A'    : round(m_total_A2, 1),
        'lift_A_kg'    : round(lift_A3, 2),
        'm_total_B'    : round(m_total_B2, 1),
        'lift_B_kg'    : round(lift_B3, 2),
        'nota'         : 'B: O₂ aun a bordo antes de ventilar',
    })

    # Fase 4b: A — O₂ retorna a tanques (sin cambio); B — O₂ se ventila
    m_total_B4 = m_total_B2 - m_O2_man   # ventila O₂
    lift_B4    = lift(m_H2_B3, m_total_B4)

    fases.append({
        'fase'         : '4b. O₂ retorna/ventila',
        'm_total_A'    : round(m_total_A2, 1),
        'lift_A_kg'    : round(lift_A3, 2),
        'm_total_B'    : round(m_total_B4, 1),
        'lift_B_kg'    : round(lift_B4, 2),
        'nota'         : f'A: O₂ a tanques (sin cambio)  |  B: ventila {m_O2_man:.1f} kg O₂',
    })

    return fases


# =============================================================================
# TRANSICIONES DE VOLUMEN
# =============================================================================

def tabla_volumenes(m_H2_man):
    """
    Volúmenes en cada transición del ciclo para m_H2_man kg quemados.
    Común a ambos escenarios (misma cantidad de H₂O producida).
    """
    m_O2       = RATIO_O2_H2 * m_H2_man
    m_H2O      = RATIO_H2O  * m_H2_man

    V_H2_gas   = m_H2_man / rho_H2_op          # [m³] H₂ en gas cells
    V_O2_tank  = m_O2 / rho_O2_tank            # [m³] O₂ en tanques (350 bar)
    V_vapor    = m_H2O / rho_vapor             # [m³] H₂O vapor a 100°C
    V_liquido  = m_H2O / rho_liq              # [m³] H₂O líquida

    return {
        'V_H2_m3'      : round(V_H2_gas, 1),
        'V_O2_tank_m3' : round(V_O2_tank, 3),
        'V_vapor_m3'   : round(V_vapor, 1),
        'V_liq_m3'     : round(V_liquido, 3),
        'ratio_H2_vapor': round(V_vapor / V_H2_gas, 2),   # H₂→vapor
        'ratio_vap_liq' : round(V_vapor / V_liquido, 0),  # vapor→líquido
    }


# =============================================================================
# BARRIDO: déficit de flotabilidad vs m_H2 quemado
# =============================================================================

def barrido_deficit(m_carga, m_estructura, n_puntos=200):
    """
    Calcula déficit de flotabilidad para ambos escenarios barriendo m_H2_man.
    """
    fma    = f_masa_equilibrio()
    m_total_ini = m_estructura + m_carga
    m_H2_ini    = m_total_ini * fma

    m_H2_A_min = m_H2_min_escenario_A(m_carga)
    m_H2_B_min = m_H2_min_escenario_B(m_carga)
    rango_max  = m_H2_A_min * 2.2

    m_H2_arr = np.linspace(0, rango_max, n_puntos)

    deficit_A = []
    deficit_B = []
    for mq in m_H2_arr:
        m_O2 = RATIO_O2_H2 * mq
        # A: sin cambio de masa, cae V_H2
        dA = (p.rho_aire_std - rho_H2_op) / rho_H2_op * mq
        # B: cae V_H2 + añade masa O₂
        dB = dA + m_O2
        deficit_A.append(dA)
        deficit_B.append(dB)

    return m_H2_arr, np.array(deficit_A), np.array(deficit_B), m_H2_A_min, m_H2_B_min


# =============================================================================
# COMPARACIÓN POR DIRIGIBLES
# =============================================================================

DIRIGIBLES = [
    {'nombre': 'Zeppelin NT',         'm_carga': 1000,  'm_estructura': 5600},
    {'nombre': 'Pathfinder 3',        'm_carga': 1000,  'm_estructura': 99000},
    {'nombre': 'Flying Whales LCA60T','m_carga': 60000, 'm_estructura': 180000},
    {'nombre': 'LZ 129 Hindenburg',   'm_carga': 1000,  'm_estructura': 264000},
]


# =============================================================================
# ANÁLISIS PRINCIPAL
# =============================================================================

def analizar():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    print(f"\n{'='*72}")
    print("  COMPARACIÓN: O₂ de tanques (A) vs O₂ del aire (B)")
    print(f"  ρ_H2_op={rho_H2_op:.4f} kg/m³  ρ_aire={p.rho_aire_std:.3f} kg/m³")
    print(f"  ρ_O2_tank={rho_O2_tank:.1f} kg/m³ (350 bar, Z={Z_O2})")
    print(f"  ρ_vapor_H2O={rho_vapor:.3f} kg/m³ (100°C, 1 atm)")
    print(f"{'='*72}\n")

    # -------------------------------------------------------------------------
    # 1. m_H2_min por escenario y dirigible
    # -------------------------------------------------------------------------
    print("  [1] m_H2_min requerido para maniobra completa:\n")
    print(f"  {'Dirigible':<25} {'m_carga':>8}  {'A [kg]':>8}  {'B [kg]':>8}  {'B/A %':>7}  {'Ahorro B':>9}")
    print(f"  {'-'*68}")
    filas_csv = []
    for d in DIRIGIBLES:
        mA = m_H2_min_escenario_A(d['m_carga'])
        mB = m_H2_min_escenario_B(d['m_carga'])
        pen_B = m_H2_penalizacion_ascenso_B(mB)
        ciclo_B = mB + pen_B   # costo ciclo completo B (H₂ total involucrado)
        ratio = mB / mA * 100
        ahorro_pct = (1 - mB / mA) * 100
        print(f"  {d['nombre']:<25} {d['m_carga']:>8,}  {mA:>8.3f}  {mB:>8.3f}  {ratio:>7.1f}%  {ahorro_pct:>8.1f}%")
        filas_csv.append({
            'dirigible'      : d['nombre'],
            'm_carga'        : d['m_carga'],
            'm_H2_min_A_kg'  : round(mA, 4),
            'm_H2_min_B_kg'  : round(mB, 4),
            'ratio_B_A_pct'  : round(ratio, 2),
            'ahorro_B_pct'   : round(ahorro_pct, 2),
            'penalizacion_ascenso_B_kg': round(pen_B, 4),
            'costo_ciclo_B_kg': round(ciclo_B, 4),
        })
    print()

    # -------------------------------------------------------------------------
    # 2. Ciclo fase a fase (Pathfinder 3, m_carga=1000 kg)
    # -------------------------------------------------------------------------
    m_carga_ref = 1000
    m_est_ref   = 99000
    mA_ref = m_H2_min_escenario_A(m_carga_ref)
    mB_ref = m_H2_min_escenario_B(m_carga_ref)
    fases  = ciclo_completo(m_carga_ref, mA_ref, mB_ref, m_est_ref)

    print(f"  [2] Ciclo fase a fase — Pathfinder 3 (m_carga={m_carga_ref} kg):\n")
    print(f"  {'Fase':<35} {'m_total_A':>10} {'lift_A':>8} || {'m_total_B':>10} {'lift_B':>8}")
    print(f"  {'-'*80}")
    for f in fases:
        print(f"  {f['fase']:<35} {f['m_total_A']:>10.1f} {f['lift_A_kg']:>+8.2f} || "
              f"{f['m_total_B']:>10.1f} {f['lift_B_kg']:>+8.2f}")
    print()

    # -------------------------------------------------------------------------
    # 3. Transiciones de volumen
    # -------------------------------------------------------------------------
    vol = tabla_volumenes(mA_ref)
    vol_B = tabla_volumenes(mB_ref)

    print(f"  [3] Transiciones de volumen (m_carga=1000 kg):\n")
    print(f"  Escenario A (quema {mA_ref:.3f} kg H₂):")
    print(f"    H₂ gas cells      : {vol['V_H2_m3']:>8.1f} m³")
    print(f"    O₂ tanques 350 bar: {vol['V_O2_tank_m3']:>8.3f} m³   (negligible)")
    print(f"    H₂O vapor (100°C) : {vol['V_vapor_m3']:>8.1f} m³   (ratio H₂→vap = {vol['ratio_H2_vapor']:.2f})")
    print(f"    H₂O líquida       : {vol['V_liq_m3']:>8.3f} m³   (ratio vap→liq = {vol['ratio_vap_liq']:.0f})")
    print()
    print(f"  Escenario B (quema {mB_ref:.3f} kg H₂):")
    print(f"    H₂ gas cells      : {vol_B['V_H2_m3']:>8.1f} m³")
    print(f"    O₂ del aire       :        — m³   (no ocupa volumen interno)")
    print(f"    H₂O vapor (100°C) : {vol_B['V_vapor_m3']:>8.1f} m³   (ratio H₂→vap = {vol_B['ratio_H2_vapor']:.2f})")
    print(f"    H₂O líquida       : {vol_B['V_liq_m3']:>8.3f} m³   (ratio vap→liq = {vol_B['ratio_vap_liq']:.0f})")
    print()

    # -------------------------------------------------------------------------
    # 4. Balance energético del ciclo completo
    # -------------------------------------------------------------------------
    print(f"  [4] Balance del ciclo completo (Pathfinder 3, m_carga=1000 kg):\n")
    pen_A = 0.0   # Escenario A no tiene penalización de ascenso por O₂ extra
    pen_B_val = m_H2_penalizacion_ascenso_B(mB_ref)
    ciclo_A = mA_ref
    ciclo_B = mB_ref + pen_B_val
    E_A = mA_ref * p.Q_comb_H2 * p.rend_gen / 3.6e6
    E_B = mB_ref * p.Q_comb_H2 * p.rend_gen / 3.6e6

    print(f"  {'Concepto':<40} {'Escen. A':>10}  {'Escen. B':>10}")
    print(f"  {'-'*63}")
    print(f"  {'H₂ quemado en descenso [kg]':<40} {mA_ref:>10.3f}  {mB_ref:>10.3f}")
    print(f"  {'Energia generada en descenso [kWh]':<40} {E_A:>10.2f}  {E_B:>10.2f}")
    print(f"  {'H₂ extra para ascenso (O₂ a bordo) [kg]':<40} {pen_A:>10.3f}  {pen_B_val:>10.3f}")
    print(f"  {'H₂ total involucrado ciclo completo [kg]':<40} {ciclo_A:>10.3f}  {ciclo_B:>10.3f}")
    print(f"  {'Diferencia A-B (descenso) [kg]':<40} {mA_ref - mB_ref:>10.3f}  (B usa menos)")
    print(f"  {'Costo ciclo completo A ≈ B?':<40} {'~igual' if abs(ciclo_A - ciclo_B) < 0.5 else 'diferente':>10}")
    print()

    # CSV
    csv_path = f'comparacion_O2_{timestamp}.csv'
    with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(filas_csv[0].keys()))
        writer.writeheader()
        writer.writerows(filas_csv)
    print(f"  CSV guardado: {csv_path}")

    # Gráficas
    _graficar(fases, filas_csv, mA_ref, mB_ref, m_carga_ref, m_est_ref, timestamp)


# =============================================================================
# GRÁFICAS
# =============================================================================

def _graficar(fases, filas_csv, mA_ref, mB_ref, m_carga, m_estructura, timestamp):

    COLOR_A = '#1f77b4'   # azul
    COLOR_B = '#d62728'   # rojo

    fig = plt.figure(figsize=(17, 12))
    fig.suptitle('O₂ source comparison: tanks (A) vs air (B)\n'
                 'Regenerative H₂ cycle in airship', fontsize=13, fontweight='bold')

    gs  = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.38)
    ax1 = fig.add_subplot(gs[0, 0])   # déficit de flotabilidad vs m_H2
    ax2 = fig.add_subplot(gs[0, 1])   # flotabilidad por fase
    ax3 = fig.add_subplot(gs[0, 2])   # m_H2_min por dirigible (barras)
    ax4 = fig.add_subplot(gs[1, 0])   # volúmenes transición
    ax5 = fig.add_subplot(gs[1, 1])   # masa total por fase
    ax6 = fig.add_subplot(gs[1, 2])   # balance ciclo completo

    # -------------------------------------------------------------------------
    # Panel 1: Déficit de flotabilidad vs m_H2 quemado
    # -------------------------------------------------------------------------
    m_arr, dA, dB, mAmin, mBmin = barrido_deficit(m_carga, m_estructura)

    ax1.plot(m_arr, dA, color=COLOR_A, linewidth=2, label='A — O₂ from tanks')
    ax1.plot(m_arr, dB, color=COLOR_B, linewidth=2, label='B — O₂ from air')
    ax1.axhline(m_carga, color='black', linestyle='--', linewidth=1.2, label=f'm_carga = {m_carga} kg')
    ax1.axvline(mAmin, color=COLOR_A, linestyle=':', linewidth=1.5,
                label=f'ṁ_min A = {mAmin:.3f} kg')
    ax1.axvline(mBmin, color=COLOR_B, linestyle=':', linewidth=1.5,
                label=f'ṁ_min B = {mBmin:.3f} kg')
    ax1.fill_betweenx([0, m_carga * 1.5], mBmin, mAmin, alpha=0.12, color='green',
                      label=f'Savings B: {(mAmin - mBmin):.3f} kg')
    ax1.set_xlabel('H₂ burned [kg]')
    ax1.set_ylabel('Buoyancy deficit [kg equiv.]')
    ax1.set_title('Buoyancy deficit\nvs H₂ burned')
    ax1.legend(fontsize=7.5, loc='upper left')
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.4)

    # -------------------------------------------------------------------------
    # Panel 2: Flotabilidad neta por fase
    # -------------------------------------------------------------------------
    etiquetas = [f['fase'].split('.')[1].strip()[:22] for f in fases]
    lift_A_vals = [f['lift_A_kg'] for f in fases]
    lift_B_vals = [f['lift_B_kg'] for f in fases]
    x_idx = np.arange(len(fases))

    ax2.plot(x_idx, lift_A_vals, 'o-', color=COLOR_A, linewidth=2, markersize=6, label='A')
    ax2.plot(x_idx, lift_B_vals, 's--', color=COLOR_B, linewidth=2, markersize=6, label='B')
    ax2.axhline(0, color='gray', linewidth=0.8, linestyle='-')
    ax2.set_xticks(x_idx)
    ax2.set_xticklabels(etiquetas, rotation=35, ha='right', fontsize=7.5)
    ax2.set_ylabel('Net buoyancy [kg equiv.]')
    ax2.set_title('Net buoyancy per phase')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.4)

    # -------------------------------------------------------------------------
    # Panel 3: m_H2_min por dirigible
    # -------------------------------------------------------------------------
    nombres  = [d['dirigible'].replace(' ', '\n') for d in filas_csv]
    vals_A   = [d['m_H2_min_A_kg'] for d in filas_csv]
    vals_B   = [d['m_H2_min_B_kg'] for d in filas_csv]
    x3       = np.arange(len(filas_csv))
    w        = 0.35

    bars_A = ax3.bar(x3 - w/2, vals_A, w, color=COLOR_A, alpha=0.85, label='A — tanks')
    bars_B = ax3.bar(x3 + w/2, vals_B, w, color=COLOR_B, alpha=0.85, label='B — air')
    for bar, v in zip(bars_A, vals_A):
        ax3.text(bar.get_x() + bar.get_width()/2, v * 1.02, f'{v:.2f}',
                 ha='center', va='bottom', fontsize=7)
    for bar, v in zip(bars_B, vals_B):
        ax3.text(bar.get_x() + bar.get_width()/2, v * 1.02, f'{v:.2f}',
                 ha='center', va='bottom', fontsize=7)
    ax3.set_xticks(x3)
    ax3.set_xticklabels(nombres, fontsize=7.5)
    ax3.set_ylabel('m_H₂_min [kg]')
    ax3.set_title('Minimum H₂ for descent\nper airship and scenario')
    ax3.legend(fontsize=9)
    ax3.grid(True, axis='y', alpha=0.4)

    # -------------------------------------------------------------------------
    # Panel 4: Transiciones de volumen (A vs B lado a lado)
    # -------------------------------------------------------------------------
    vol_A = tabla_volumenes(mA_ref)
    vol_B = tabla_volumenes(mB_ref)

    categorias = ['H₂\n(gas cells)', 'O₂\n(tank/air)', 'H₂O\nvapor', 'H₂O\nliquid']
    vols_A = [vol_A['V_H2_m3'], vol_A['V_O2_tank_m3'], vol_A['V_vapor_m3'], vol_A['V_liq_m3']]
    vols_B = [vol_B['V_H2_m3'], 0.0,                  vol_B['V_vapor_m3'], vol_B['V_liq_m3']]

    x4 = np.arange(len(categorias))
    ba = ax4.bar(x4 - w/2, vols_A, w, color=COLOR_A, alpha=0.85, label='A')
    bb = ax4.bar(x4 + w/2, vols_B, w, color=COLOR_B, alpha=0.85, label='B')
    for bar, v in zip(ba, vols_A):
        if v > 0.01:
            ax4.text(bar.get_x() + bar.get_width()/2, v * 1.03,
                     f'{v:.1f}' if v >= 1 else f'{v:.3f}',
                     ha='center', va='bottom', fontsize=7)
    for bar, v in zip(bb, vols_B):
        if v > 0.01:
            ax4.text(bar.get_x() + bar.get_width()/2, v * 1.03,
                     f'{v:.1f}' if v >= 1 else f'{v:.3f}',
                     ha='center', va='bottom', fontsize=7)
    ax4.set_xticks(x4)
    ax4.set_xticklabels(categorias, fontsize=9)
    ax4.set_ylabel('Volume [m³]')
    ax4.set_title('Volume transitions\nby H₂O state')
    ax4.legend(fontsize=9)
    ax4.grid(True, axis='y', alpha=0.4)
    # Nota: escala log para ver líquido y vapor juntos
    ax4.set_yscale('log')
    ax4.set_ylim(bottom=0.0001)

    # -------------------------------------------------------------------------
    # Panel 5: Masa total por fase
    # -------------------------------------------------------------------------
    m_total_A = [f['m_total_A'] for f in fases]
    m_total_B = [f['m_total_B'] for f in fases]
    ax5.plot(x_idx, m_total_A, 'o-', color=COLOR_A, linewidth=2, markersize=6, label='A')
    ax5.plot(x_idx, m_total_B, 's--', color=COLOR_B, linewidth=2, markersize=6, label='B')
    ax5.fill_between(x_idx, m_total_A, m_total_B, alpha=0.15, color='orange',
                     label='Difference (O₂ on board)')
    ax5.set_xticks(x_idx)
    ax5.set_xticklabels(etiquetas, rotation=35, ha='right', fontsize=7.5)
    ax5.set_ylabel('Total airship mass [kg]')
    ax5.set_title('Total airship mass\nper phase')
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.4)

    # -------------------------------------------------------------------------
    # Panel 6: Resumen cuantitativo del ciclo completo
    # -------------------------------------------------------------------------
    ax6.axis('off')
    pen_B = m_H2_penalizacion_ascenso_B(mB_ref)
    ciclo_A = mA_ref
    ciclo_B = mB_ref + pen_B
    ahorro_desc = (mA_ref - mB_ref) / mA_ref * 100
    m_O2_B = RATIO_O2_H2 * mB_ref

    texto = (
        f"Pathfinder 3  |  m_carga = {m_carga} kg\n"
        f"{'─'*44}\n"
        f"{'Quantity':<32} {'A':>6}  {'B':>6}\n"
        f"{'─'*44}\n"
        f"{'H₂ burned descent [kg]':<32} {mA_ref:>6.3f}  {mB_ref:>6.3f}\n"
        f"{'O₂ involved [kg]':<32} {RATIO_O2_H2*mA_ref:>6.1f}  {m_O2_B:>6.1f}\n"
        f"{'H₂O produced [kg]':<32} {RATIO_H2O*mA_ref:>6.2f}  {RATIO_H2O*mB_ref:>6.2f}\n"
        f"{'V_H₂O vapor [m³]':<32} {RATIO_H2O*mA_ref/rho_vapor:>6.1f}  {RATIO_H2O*mB_ref/rho_vapor:>6.1f}\n"
        f"{'V_H₂O liq [m³]':<32} {RATIO_H2O*mA_ref/rho_liq:>6.4f}  {RATIO_H2O*mB_ref/rho_liq:>6.4f}\n"
        f"{'Ascent penalty B [kg H₂]':<32} {'0':>6}  {pen_B:>6.3f}\n"
        f"{'─'*44}\n"
        f"{'Full cycle cost [kg H₂]':<32} {ciclo_A:>6.3f}  {ciclo_B:>6.3f}\n"
        f"{'─'*44}\n"
        f"Descent savings: {ahorro_desc:.1f}%  (B uses less H₂)\n"
        f"Cycle cost A ≈ B: {'YES' if abs(ciclo_A - ciclo_B) < 0.5 else 'NO'} "
        f"(diff = {abs(ciclo_A-ciclo_B):.3f} kg)\n"
        f"ratio vap→liq: 1:{vol_A['ratio_vap_liq']:.0f}  (ambos escenarios igual)"
    )

    ax6.text(0.05, 0.97, texto, transform=ax6.transAxes,
             fontsize=8.5, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax6.set_title('Quantitative summary', fontsize=10)

    plt.savefig(f'comparacion_O2_{timestamp}.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Grafica guardada: comparacion_O2_{timestamp}.png")


# =============================================================================
# EJECUCIÓN
# =============================================================================

if __name__ == '__main__':
    analizar()
