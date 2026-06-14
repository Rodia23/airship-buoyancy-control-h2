# -*- coding: utf-8 -*-
"""
comparacion_O2_corrected.py
============================
Comparison of O2 source (on-board tanks vs. atmospheric air) for descent-phase
H2 consumption after a 1-tonne payload release.

PHYSICAL MODEL — key distinction between the two scenarios:

  Scenario A — O2 from on-board tanks (closed system):
    When H2 burns with on-board O2:
      - H2 leaves ballonets       →  volume decreases → buoyancy ↓
      - O2 leaves storage tank    →  becomes H2O water ballast (stays on board)
      - H2O retained on board     →  mass conserved, ΔM_vehicle = 0
    Net downward force per kg H2 burned = (ρ_air / ρ_H2) × g
    → factor_A = ρ_air / ρ_H2

  Scenario B — O2 from atmosphere (open system, this work):
    When H2 burns with atmospheric O2:
      - H2 leaves ballonets       →  volume decreases → buoyancy ↓
      - 8 kg O2 absorbed from air →  becomes H2O water ballast → ΔM_vehicle = +8·ṁ
    Net downward force per kg H2 burned = (ρ_air / ρ_H2 + 8) × g  =  L_eff × g
    → factor_B = ρ_air / ρ_H2 + 8  =  L_eff

NOTE — what the original comparacion_O2.py computed (incorrectly for Scenario A):
    factor_A_wrong = (ρ_air - ρ_H2) / ρ_H2 = ρ_air/ρ_H2 - 1
    This corresponds to VENTING H2 to atmosphere (both buoyancy and H2 weight are
    lost), not to combustion where the H2 mass is conserved as on-board water.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(__file__))
import parametros as p

# =============================================================================
# CONDITIONS
# Use ISA sea-level (STP) so that Scenario B gives 44.7 kg, consistent with
# Eq. mh2comb in the paper (operational conditions T=281 K give 44.5 kg).
# =============================================================================
P_STP        = 101325.0   # Pa    — ISA sea-level pressure
T_STP        = 288.0      # K     — ISA sea-level temperature (~15 °C)
rho_air_stp  = 1.225      # kg/m³ — ISA sea-level air density
rho_H2_stp   = P_STP / (p.R_H2 * T_STP)   # ≈ 0.0853 kg/m³

RATIO_O2_H2  = 8.0        # stoichiometric O2/H2 mass ratio

# =============================================================================
# EFFECTIVE FORCE FACTORS  (net downward force per kg H2 burned / g)
# =============================================================================
factor_A = rho_air_stp / rho_H2_stp           # O2 from tanks: buoyancy only
factor_B = rho_air_stp / rho_H2_stp + RATIO_O2_H2   # air: buoyancy + absorbed O2 mass

# =============================================================================
# MINIMUM H2 REQUIRED  (to cancel m_load × g)
# =============================================================================
def m_H2_required(m_load, factor):
    return m_load / factor

m_load = 1000.0   # kg payload

mA = m_H2_required(m_load, factor_A)
mB = m_H2_required(m_load, factor_B)

reduction_pct = (mA - mB) / mA * 100.0
savings_kg    = mA - mB

sep = "=" * 50
print(sep)
print("O2 SOURCE COMPARISON - corrected model")
print(sep)
print(f"  rho_air / rho_H2    = {rho_air_stp/rho_H2_stp:.4f}")
print(f"  L_eff  (air)        = {factor_B:.4f}")
print()
print(f"  Scenario A (O2 tanks):  factor = {factor_A:.4f}")
print(f"    m_H2 required = {mA:.2f} kg")
print()
print(f"  Scenario B (air):       factor = {factor_B:.4f}")
print(f"    m_H2 required = {mB:.2f} kg")
print()
print(f"  Reduction (A->B): {reduction_pct:.1f}%  (saves {savings_kg:.1f} kg H2 per delivery)")
print()
print(f"  Water produced (A): {9*mA:.1f} kg/delivery")
print(f"  Water produced (B): {9*mB:.1f} kg/delivery")
print()
print(f"  E_gen (A): {mA * 120e6 * p.rend_gen / 3.6e6:.0f} kWh/delivery")
print(f"  E_gen (B): {mB * 120e6 * p.rend_gen / 3.6e6:.0f} kWh/delivery")
print()

platforms = [
    ("LZ 129 Hindenburg", 24000),
    ("Pathfinder 3",       13000),
    ("Flying Whales",      23000),
    ("Zeppelin NT",         2800),
]
print("  Solar regeneration (Lima, 2.0 kWh/m2/day, target = Scenario B):")
print(f"  {'Platform':<20} {'reg(kg/d)':>10} {'days':>8} {'hours':>8}")
for name, A_total in platforms:
    reg  = 2.0 * 0.45 * A_total * p.eta_electrolisis / p.E_H2_electrolisis
    days = mB / reg
    print(f"  {name:<20} {reg:>10.0f} {days:>8.2f} {days*24:>8.1f}")

print(sep)

# =============================================================================
# FIGURE
# =============================================================================
FIGSIZE   = (7.5, 5.0)
DPI       = 300
FS_LABEL  = 11
FS_LEGEND = 10
FS_TITLE  = 13
LW_MAIN   = 2.5
LW_DASH   = 1.5
COLOR_A   = '#1f77b4'   # blue  — O2 tanks
COLOR_B   = '#d62728'   # red   — air

n_pts    = 300
x_max    = mA * 1.5
m_H2_arr = np.linspace(0, x_max, n_pts)

deficit_A = factor_A * m_H2_arr
deficit_B = factor_B * m_H2_arr

fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
ax.set_facecolor('#fcfcfc')
ax.grid(True, linestyle='-', alpha=0.3, color='gray')

ax.plot(m_H2_arr, deficit_A, color=COLOR_A, lw=LW_MAIN,
        label=f'Scenario A — O$_2$ from tanks   ($\\dot{{m}}_{{\\mathrm{{min}}}}$ = {mA:.1f} kg)')
ax.plot(m_H2_arr, deficit_B, color=COLOR_B, lw=LW_MAIN,
        label=f'Scenario B — O$_2$ from air      ($\\dot{{m}}_{{\\mathrm{{min}}}}$ = {mB:.1f} kg)')

ax.axhline(m_load, color='black', ls='--', lw=LW_DASH,
           label=f'Payload $m_{{\\mathrm{{load}}}}$ = {m_load:.0f} kg')

ax.axvline(mA, color=COLOR_A, ls=':', lw=1.8)
ax.axvline(mB, color=COLOR_B, ls=':', lw=1.8)

# double-headed arrow showing savings
arrow_y = m_load * 1.12
ax.annotate('', xy=(mA, arrow_y), xytext=(mB, arrow_y),
            arrowprops=dict(arrowstyle='<->', color='green', lw=2.0))
ax.text((mA + mB) / 2, arrow_y + m_load * 0.04,
        f'{savings_kg:.1f} kg saved ({reduction_pct:.1f}%)',
        color='green', fontsize=FS_LEGEND, ha='center', va='bottom', fontweight='bold')

ax.set_xlabel('H$_2$ burned [kg]', fontsize=FS_LABEL, labelpad=10)
ax.set_ylabel('Net buoyancy deficit [kg equiv.]', fontsize=FS_LABEL, labelpad=10)

ax.set_xlim(0, x_max)
ax.set_ylim(0, m_load * 1.8)
ax.legend(fontsize=FS_LEGEND, loc='upper left', framealpha=0.95, edgecolor='gray')

plt.tight_layout()

output_path = os.path.join(os.path.dirname(__file__),
                           '..', 'paper_manuscript', 'buoyancy deficit.png')
plt.savefig(output_path, bbox_inches='tight')
print(f"\nFigure saved: {os.path.abspath(output_path)}")
plt.show()
