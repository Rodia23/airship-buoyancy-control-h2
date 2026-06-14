# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# apendice_dlim_flujo.py
# D_lim vs flujo óptimo mínimo para los 4 dirigibles
# =============================================================================

import os
import numpy as np
import matplotlib.pyplot as plt

OUTDIR   = os.path.dirname(os.path.abspath(__file__))
BASENAME = 'apendice_dlim_flujo'

def siguiente_nombre(base, ext='png'):
    """Devuelve el próximo nombre incremental que no exista."""
    i = 1
    while os.path.exists(os.path.join(OUTDIR, f'{base}_{i:03d}.{ext}')):
        i += 1
    return os.path.join(OUTDIR, f'{base}_{i:03d}.{ext}')

# ── Parámetros físicos (parametros.py, Lima ~1500 m s.n.m.) ─────────────────
g            = 9.81
R_H2         = 8.314 / 0.002016
T1           = 281.0
P_atmos      = 84700.0
rho_aire_std = 1.058
H_escala     = 8500.0
h0           = 50.0
h_max        = 5.0
m_cargo      = 1000.0

rho_h0 = rho_aire_std * np.exp(-h0 / H_escala)
rho_H2 = P_atmos / (R_H2 * T1)
denom  = rho_h0 / rho_H2 - 1
ratio  = rho_H2 / rho_h0
L_eff  = rho_h0 / rho_H2 + 8.0

# ── Dirigibles ───────────────────────────────────────────────────────────────
DIRIGIBLES = [
    {'nombre': 'LZ 129 Hindenburg',    'm_est': 263600, 'color': '#d62728'},
    {'nombre': 'Flying Whales LCA60T', 'm_est': 129000, 'color': '#1f77b4'},
    {'nombre': 'Pathfinder 3',         'm_est':  99000, 'color': '#2ca02c'},
    {'nombre': 'Zeppelin NT',          'm_est':   5622, 'color': '#ff7f0e'},
]
for d in DIRIGIBLES:
    m_H2_ini = (d['m_est'] + m_cargo) / denom
    d['M']   = d['m_est'] + m_cargo + m_H2_ini
    d['dm0'] = ratio * np.sqrt(m_cargo**3 * g / (2.0 * h_max * d['M']))

# ── Funciones analíticas ─────────────────────────────────────────────────────
def alpha_star(C):
    if C <= 0: return 0.0
    roots = np.roots([1.0, C, 0.0, -C])
    cands = [r.real for r in roots
             if abs(r.imag) < 1e-9 and 1e-9 < r.real < 1.0 - 1e-9]
    return min(cands)**2 if cands else 0.0

def C_de_D(D):
    return ratio * L_eff * np.sqrt(3.0 * D / h_max) if D > 0 else 0.0

# Grilla D_lim: escala log de 0.005 m hasta h_max = 5 m
D_arr  = np.logspace(np.log10(5e-3), np.log10(h_max), 800)
a_arr  = np.array([alpha_star(C_de_D(D)) for D in D_arr])

# ── Figura ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8.0, 5.2))

# ── Curvas ṁ_opt(D_lim) por dirigible ───────────────────────────────────────
for d in DIRIGIBLES:
    dm_arr = d['dm0'] * (1.0 - a_arr)
    ax.plot(D_arr, dm_arr, color=d['color'], lw=2.3, label=d['nombre'], zorder=3)

# ── Línea vertical en D_lim = h_max ─────────────────────────────────────────
ax.axvline(h_max, color='black', lw=1.5, ls='--', alpha=0.5, zorder=1)
ax.text(h_max * 0.92, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 20,
        f'$D_\\mathrm{{lim}} = h_\\mathrm{{max}} = {h_max:.0f}$ m',
        fontsize=8, ha='right', va='top', color='black', alpha=0.6)

# ── Decoración ───────────────────────────────────────────────────────────────
ax.set_xscale('log')
ax.set_xlim(5e-3, h_max)
ax.set_xticks([0.01, 0.05, 0.1, 0.5, 1.0, 5.0])
ax.set_xticklabels(['0.01', '0.05', '0.1', '0.5', '1.0', '5.0'])
ax.set_ylim(0, None)
ax.set_xlabel(r'Max. preconditioning depth  $D_\mathrm{lim}$  [m]', fontsize=10)
ax.set_ylabel(r'Minimum optimal flow  $\dot{m}_\mathrm{opt}$  [kg/s]', fontsize=10)
ax.set_title('Minimum required flow vs. preconditioning depth', fontsize=11, pad=7)
ax.grid(True, alpha=0.25)
ax.tick_params(labelsize=9)
ax.legend(fontsize=8.5, loc='upper right', framealpha=0.93)

plt.tight_layout()
fname = siguiente_nombre(BASENAME)
plt.savefig(fname, dpi=180, bbox_inches='tight')
plt.show()
print(f"Guardado: {fname}")

# ── Tabla consola ────────────────────────────────────────────────────────────
print(f"\n  {'Dirigible':<28} {'ṁ₀ (reactivo)':>15} {'ṁ_opt @ 5m':>12}")
for d in DIRIGIBLES:
    dm_5 = d['dm0'] * (1 - alpha_star(C_de_D(h_max)))
    print(f"  {d['nombre']:<28} {d['dm0']:>14.2f}  {dm_5:>11.2f} kg/s")
