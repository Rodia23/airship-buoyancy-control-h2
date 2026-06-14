"""
Comparación exacta vs. aproximación cúbica para el hundimiento del dirigible.

Fórmula exacta (sin aproximar M(t) ≈ M_total):
  |Δh| = (L_eff·g·M²) / (512·ṁ²) · φ(S)
  φ(S) = S²/2 + S - (1+S)·ln(1+S)
  S    = 8·ṁ·t / M_total   (aumento fraccional de masa por O₂)

Fórmula aproximada (M constante):
  |Δh| = (L_eff·g·ṁ) / (6·M_total) · t³
       = (L_eff·g·M²) / (512·ṁ²) · S³/6

El error nace cuando S = 8ṁt/M_total ya no es ≪ 1.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec

# ── Funciones ──────────────────────────────────────────────────────────────────

def phi(S):
    """Función de corrección exacta. phi(S) ≈ S³/6 para S pequeño."""
    return 0.5*S**2 + S - (1 + S)*np.log1p(S)   # log1p = ln(1+S), más estable

def h_exacto(S, L_eff, g, M, mdot):
    """Desplazamiento exacto en función del parámetro adimensional S."""
    return (L_eff * g * M**2) / (512 * mdot**2) * phi(S)

def h_aprox(S, L_eff, g, M, mdot):
    """Desplazamiento aproximado (fórmula cúbica)."""
    return (L_eff * g * M**2) / (512 * mdot**2) * (S**3 / 6)

def error_relativo(S):
    """Error relativo = (h_aprox - h_exacto) / h_exacto = S³/6/φ(S) - 1."""
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = (S**3 / 6) / phi(S)
        ratio = np.where(S < 1e-6, 1.0, ratio)
    return (ratio - 1) * 100   # en %

# ── Parámetros físicos de referencia ──────────────────────────────────────────
L_eff_0  = 14.0    # L_eff = ρ_air/ρ_H2 + 8  (~14 a baja altitud)
g        = 9.81    # m/s²
M_0      = 5000.0  # kg  — masa total del vehículo
mdot_0   = 0.05    # kg/s — caudal másico de H₂

# ── Grilla S (parámetro adimensional) ─────────────────────────────────────────
S_vals = np.linspace(1e-4, 3.0, 1000)

# ── Figura ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 7))
gs  = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)
plt.subplots_adjust(bottom=0.28)

ax1 = fig.add_subplot(gs[0, 0])   # φ(S) exacta vs S³/6
ax2 = fig.add_subplot(gs[0, 1])   # Error relativo (%)
ax3 = fig.add_subplot(gs[1, :])   # |Δh| vs tiempo real

# ── Ax1: φ(S) ─────────────────────────────────────────────────────────────────
ax1.plot(S_vals, phi(S_vals),    color='darkorchid', lw=2.5, label=r'$\phi(S)$ exacta')
ax1.plot(S_vals, S_vals**3 / 6, color='tomato',     lw=1.8, ls='--', label=r'$S^3/6$ (aprox. cúbica)')
ax1.axvline(1.0, color='gray', ls=':', lw=1)
ax1.set_xlabel('$S = 8\\dot{m}t / M_{total}$', weight='bold')
ax1.set_ylabel(r'$\phi(S)$', weight='bold')
ax1.set_title('Función de corrección', fontsize=10)
ax1.set_xlim(0, 3); ax1.set_ylim(0, None)
ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)

# ── Ax2: Error relativo ────────────────────────────────────────────────────────
err_vals = error_relativo(S_vals)
ax2.plot(S_vals, err_vals, color='steelblue', lw=2)
ax2.axhline(5,  color='orange', ls='--', lw=1, label='5% error')
ax2.axhline(20, color='red',    ls='--', lw=1, label='20% error')
# S crítico al 5%
S_5pct = S_vals[np.argmin(np.abs(err_vals - 5))]
ax2.axvline(S_5pct, color='orange', ls=':', lw=1)
ax2.text(S_5pct + 0.05, 25, f'$S \\approx {S_5pct:.2f}$', color='orange', fontsize=8)
ax2.set_xlabel('$S = 8\\dot{m}t / M_{total}$', weight='bold')
ax2.set_ylabel('Sobreestimación (%)', weight='bold')
ax2.set_title('Error de la aprox. cúbica', fontsize=10)
ax2.set_xlim(0, 3); ax2.set_ylim(0, 60)
ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

# ── Ax3: |Δh| vs tiempo real ──────────────────────────────────────────────────
t_max_s = M_0 * 3.0 / (8 * mdot_0)          # tiempo para S=3
t_vals  = np.linspace(0, t_max_s, 500)
S_t     = 8 * mdot_0 * t_vals / M_0

line_ex,  = ax3.plot(t_vals, h_exacto(S_t, L_eff_0, g, M_0, mdot_0),
                     color='darkorchid', lw=2.5, label='Exacto')
line_ap,  = ax3.plot(t_vals, h_aprox(S_t,  L_eff_0, g, M_0, mdot_0),
                     color='tomato', lw=1.8, ls='--', label='Aprox. cúbica')
dot_ex,  = ax3.plot([], [], 'o', color='darkorchid', ms=9, zorder=5)
dot_ap,  = ax3.plot([], [], 'o', color='tomato',     ms=9, zorder=5)

ax3.set_xlabel('Tiempo $t$ (s)', weight='bold')
ax3.set_ylabel('Hundimiento $|\\Delta h|$ (m)', weight='bold')
ax3.set_title('Hundimiento real vs. tiempo', fontsize=10)
ax3.legend(fontsize=9); ax3.grid(True, alpha=0.3)

caja = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
info = ax3.text(0.65, 0.15, '', transform=ax3.transAxes, fontsize=10, bbox=caja)

# ── Sliders ───────────────────────────────────────────────────────────────────
sl_ax_M    = plt.axes([0.10, 0.17, 0.35, 0.03])
sl_ax_mdot = plt.axes([0.10, 0.12, 0.35, 0.03])
sl_ax_t    = plt.axes([0.55, 0.12, 0.35, 0.03])

sl_M    = Slider(sl_ax_M,    '$M_{total}$ (kg)',  500,  20000, valinit=M_0,    color='mediumpurple')
sl_mdot = Slider(sl_ax_mdot, '$\\dot{m}$ (kg/s)', 0.005, 0.5, valinit=mdot_0, color='steelblue')
sl_t    = Slider(sl_ax_t,    'Tiempo $t$ (s)',      0,  2000,  valinit=300,    color='darkorange')

def update(_=None):
    M    = sl_M.val
    mdot = sl_mdot.val
    t_s  = sl_t.val

    # Actualizar curvas
    t_max_new = M * 3.0 / (8 * mdot)
    t_new  = np.linspace(0, t_max_new, 500)
    S_new  = 8 * mdot * t_new / M
    line_ex.set_data(t_new, h_exacto(S_new, L_eff_0, g, M, mdot))
    line_ap.set_data(t_new, h_aprox( S_new, L_eff_0, g, M, mdot))
    ax3.relim(); ax3.autoscale_view()

    # Punto en t_slider
    S_ts = 8 * mdot * t_s / M
    he   = h_exacto(np.array([S_ts]), L_eff_0, g, M, mdot)[0]
    ha   = h_aprox( np.array([S_ts]), L_eff_0, g, M, mdot)[0]
    dot_ex.set_data([t_s], [he])
    dot_ap.set_data([t_s], [ha])

    err = abs(ha - he) / he * 100 if he > 0 else 0
    zona = "⚠️ aprox. dudosa" if S_ts > S_5pct else "✓ aprox. válida"
    info.set_text(
        f"t = {t_s:.0f} s  |  S = {S_ts:.3f}\n"
        f"Exacto:  {he:.3f} m\n"
        f"Aprox.:  {ha:.3f} m\n"
        f"Error:   {err:.1f}%  {zona}"
    )
    fig.canvas.draw_idle()

sl_M.on_changed(update)
sl_mdot.on_changed(update)
sl_t.on_changed(update)
update()
plt.show()
