# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# t_apendice.py
# Simulación con parámetros de preacondicionamiento calculados analíticamente
# según las derivaciones del Apéndice A (apendice.tex).
#
# Reemplaza los valores hardcodeados de t2_gain_scheduled.py:
#   ANTES: dm_H2_proactivo = -10 kg/s  (fijo)
#          t_anticipo      =  2.0 s    (fijo)
#   AHORA: dm_pre = η · ṁ_max          [A.1: L_eff + flujo proporcional]
#          t_anticipo = τ · t_max       [A.2: fórmula cúbica]
#
# Parámetros operacionales que el usuario fija:
#   DM_MAX    — capacidad del actuador [kg/s]
#   ETA       — fracción de DM_MAX usada en el precondicionamiento (0 < η ≤ 1)
#   D_MAX     — hundimiento máximo permitido antes del drop [m]
#   TAU       — fracción de t_max que se usa efectivamente (0 < τ ≤ 1)
#   DIRIGIBLE — plataforma a simular
# =============================================================================

import os
from datetime import datetime
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import parametros as p_base

# Directorio del propio .py — las imágenes se guardan aquí siempre,
# independientemente del directorio de trabajo desde donde se ejecute.
OUTDIR = os.path.dirname(os.path.abspath(__file__))


def _siguiente_path(stem):
    """
    Devuelve la próxima ruta disponible con sufijo incremental:
    stem_001.png, stem_002.png, ...
    Nunca sobreescribe un archivo existente.
    """
    for i in range(1, 10000):
        path = os.path.join(OUTDIR, f"{stem}_{i:03d}.png")
        if not os.path.exists(path):
            return path
    raise RuntimeError(f"No se encontró un nombre libre para {stem}")

# =============================================================================
# Selección de dirigible y parámetros operacionales
# =============================================================================

DIRIGIBLE = {
    'nombre':        'Pathfinder 3',
    'Area':          4390,
    'm_estructura':  99000,
    'm_carga':       1000,
}

DM_MAX  = 30.0    # [kg/s]  capacidad del actuador
ETA     = 1/3     # [—]     fracción de DM_MAX para el pulso proactivo
D_MAX   = 0.10    # [m]     hundimiento máximo operacional antes del drop
TAU     = 1.0     # [—]     fracción de t_max realmente usada (1.0 = usar todo t_max)

# =============================================================================
# A.1 — Factor de Sustentación Efectiva
# =============================================================================

def calcular_L_eff(rho_air, rho_H2):
    """
    L_eff = ρ_air/ρ_H2 + 8
    Cada kg de H2 quemado genera L_eff·g N de fuerza neta hacia abajo:
      - ρ_air/ρ_H2: pérdida de empuje aeroestático
      - 8:          ganancia de masa por O2 absorbido (1 kg H2 + 8 kg O2 → 9 kg H2O)
    """
    return rho_air / rho_H2 + 8.0


# =============================================================================
# A.2 — Tiempo máximo de anticipación (fórmula cúbica)
# =============================================================================

def calcular_t_max(M_total, D_max, L_eff, g, eta, dm_max):
    """
    t_max = (6·M·D_max / (L_eff·g·η·ṁ_max))^(1/3)
    Tiempo máximo que el dirigible puede precondicionarse sin superar D_max.
    """
    return (6.0 * M_total * D_max / (L_eff * g * eta * dm_max)) ** (1.0/3.0)


def sinking_cubico(t, M_total, L_eff, g, dm_pre):
    """
    |Δh(t)| = (L_eff·g·|ṁ|) / (6·M) · t³   [Ec. A.1]
    Desplazamiento vertical durante el precondicionamiento (sin drag).
    """
    return (L_eff * g * abs(dm_pre)) / (6.0 * M_total) * t**3


# =============================================================================
# A.3 — Hundimiento óptimo D* para compensación perfecta
# =============================================================================

def calcular_D_star(m_cargo, M_total, L_eff, g, eta, dm_max, tau):
    """
    D* = m_cargo³·g / (6·M·L_eff²·η²·ṁ_max²·τ³)
    Hundimiento al que la masa extraída exactamente compensa la carga liberada.
    """
    return (m_cargo**3 * g) / (6 * M_total * L_eff**2 * eta**2 * dm_max**2 * tau**3)


# =============================================================================
# A.4 — Fracción de compensación
# =============================================================================

def calcular_alpha(D, D_star):
    """
    α(D) = (D / D*)^(1/3)
    Fracción de la perturbación pre-compensada para un hundimiento D.
    """
    return (D / D_star) ** (1.0/3.0)


# =============================================================================
# A.5 — Flujo mínimo para compensación objetivo
# =============================================================================

def calcular_dm_min(m_cargo, M_total, L_eff, g, eta, tau, D_limit, alpha_target=1.0):
    """
    ṁ_min = m_cargo^(3/2)·g^(1/2)·α^(3/2) / (η·τ^(3/2)·√(6·M·L_eff²·D_limit))
    Flujo mínimo para lograr compensación alpha_target con hundimiento ≤ D_limit.
    """
    return (m_cargo**1.5 * g**0.5 * alpha_target**1.5 /
            (eta * tau**1.5 * np.sqrt(6 * M_total * L_eff**2 * D_limit)))


# =============================================================================
# A.6 — Solución Minimax: fracción óptima α*
# =============================================================================

def _cbrt(x):
    """Raíz cúbica con signo (evita NaN para negativos)."""
    return np.sign(x) * abs(x) ** (1.0/3.0)


def calcular_alpha_star(C):
    """
    Resuelve x³ + C·x² - C = 0, devuelve α* = (x*)²
    donde C = constante del sistema (Ec. A.6).

    Régimen 1 (Δ > 0, alta inercia): solución de Cardano.
    Régimen 2 (Δ ≤ 0, baja inercia): sustitución trigonométrica.
    """
    Delta = C**2 / 4.0 - C**4 / 27.0

    if Delta > 0:
        # Régimen 1: plataformas de alta inercia / gran escala
        q = 2*C**3/27.0 - C
        x_star = _cbrt(-q/2 + Delta**0.5) + _cbrt(-q/2 - Delta**0.5) - C/3.0
    else:
        # Régimen 2: Casus Irreducibilis — sustitución trigonométrica de Vieta
        arg = np.clip((27.0 - 2*C**2) / (2*C**2), -1.0, 1.0)
        x_star = 2*C/3.0 * np.cos(np.arccos(arg) / 3.0) - C/3.0

    alpha_star = x_star**2
    return float(np.clip(alpha_star, 0.0, 1.0))


def calcular_C(dm_post0, eta, tau, M_total, L_eff, D_limit, m_cargo, g):
    """
    C = (ṁ_post0·η·τ^(3/2)·√(6·M·L_eff²·D_limit)) / (m_cargo^(3/2)·√g)
    Constante adimensional del sistema que determina el régimen óptimo.
    """
    return (dm_post0 * eta * tau**1.5 * np.sqrt(6*M_total*L_eff**2*D_limit) /
            (m_cargo**1.5 * g**0.5))


def calcular_dm_post0(m_cargo, M_total, rho_H2, rho_air, g, h_max=5.0):
    """
    Demanda reactiva base (sin precondicionamiento) según la ley de escala del paper:
    ṁ_post0 = (ρ_H2/ρ_air) · √(m_cargo³·g / (2·h_max·M_total))
    """
    return (rho_H2 / rho_air) * np.sqrt(m_cargo**3 * g / (2*h_max*M_total))


# =============================================================================
# Motor de física (adaptado de t2_gain_scheduled.py)
# =============================================================================

class _Params:
    def __init__(self, overrides=None):
        for k, v in vars(p_base).items():
            if not k.startswith('_'):
                setattr(self, k, v)
        if overrides:
            for k, v in overrides.items():
                setattr(self, k, v)


def correr(n_entregas=1, toneladas_por_drop=1, graficar=True, verbose=True):
    p = _Params({'m_carga': n_entregas * toneladas_por_drop * 1000.0,
                 'Area':         DIRIGIBLE['Area'],
                 'm_estructura': DIRIGIBLE['m_estructura'],
                 'nombre_dirigible': DIRIGIBLE['nombre']})

    # Siempre 30 s — intervalo de 1 s mutilaba t_anticipo (≈3 s) en la gráfica base.
    # El xlim del gráfico log-X se ajusta por separado para la vista visual.
    intervalo = 30.0
    t_final   = 160.0 if n_entregas == 1 else (n_entregas + 3) * intervalo

    # ── Condiciones iniciales ─────────────────────────────────────────────────
    rho_h0     = p.rho_aire_std * np.exp(-p.h0 / p.H_escala)
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini   = (p.m_estructura + p.m_carga) / (rho_h0 / rho_H2_ini - 1)
    M_ini      = p.m_estructura + p.m_carga + m_H2_ini

    # ── Parámetros analíticos A.1-A.6 (orden correcto) ───────────────────────
    L_eff     = calcular_L_eff(rho_h0, rho_H2_ini)
    m_drop_kg = toneladas_por_drop * 1000.0

    # A.3-A.6: calculados primero, basados en M_ini y drop 1 como referencia
    D_star     = calcular_D_star(m_drop_kg, M_ini, L_eff, p.g, ETA, DM_MAX, TAU)
    alpha_cond = calcular_alpha(D_MAX, D_star)
    dm_min_pre = calcular_dm_min(m_drop_kg, M_ini, L_eff, p.g, ETA, TAU, D_MAX)
    dm_post0   = calcular_dm_post0(m_drop_kg, M_ini, rho_H2_ini, rho_h0, p.g)
    C_sys      = calcular_C(dm_post0, ETA, TAU, M_ini, L_eff, D_MAX, m_drop_kg, p.g)
    alpha_star = calcular_alpha_star(C_sys)
    K_pre      = calcular_dm_min(m_drop_kg, M_ini, L_eff, p.g, ETA, TAU, D_MAX)
    dm_opt     = K_pre * alpha_star**1.5

    # A.1: flujo óptimo Minimax como flujo de precondicionamiento [revisor punto 2]
    # dm_opt balancea perfectamente la demanda proactiva y reactiva (α*)
    # En lugar de η·ṁ_max (fracción fija), se usa el flujo teóricamente óptimo.
    dm_pre = dm_opt

    # ── Por-drop: t_anticipo recalculado por masa usando dm_pre = dm_opt ──────
    # t_max usa dm_pre en denominador (≡ η·ṁ_max), ahora coherente con Minimax.
    # D_MAX es relativo a la posición del dirigible AL INICIO del pulso proactivo
    # (no necesariamente h0 si aún no recuperó de un drop anterior).
    t_drops   = [(k + 1) * intervalo for k in range(n_entregas)]

    t_anticipos = []
    for k in range(n_entregas):
        m_carga_k = p.m_carga - k * m_drop_kg
        M_k       = p.m_estructura + m_carga_k + m_H2_ini
        # t_max con dm_pre (= dm_opt): (6·M·D_max/(L_eff·g·dm_pre))^(1/3)
        t_ant_k   = TAU * (6 * M_k * D_MAX / (L_eff * p.g * dm_pre)) ** (1.0/3.0)
        t_ant_k   = min(t_ant_k, intervalo * 0.95)   # no solapar con drop anterior
        t_anticipos.append(t_ant_k)

    t_pre_starts = [t_drops[k] - t_anticipos[k] for k in range(n_entregas)]

    if verbose:
        print("\n" + "="*66)
        print(f"  {DIRIGIBLE['nombre']} — Parámetros analíticos (Apéndice A)")
        print("="*66)
        print(f"  L_eff              = {L_eff:.3f}  (ρ_air/ρ_H2 + 8)")
        print(f"  dm_pre  [A.1]      = {dm_pre:.2f} kg/s  (η={ETA:.3f}·ṁ_max={DM_MAX})")
        print(f"  D_max operacional  = {D_MAX*100:.1f} cm  (hundimiento máx. antes del drop)")
        print(f"  D_star  [A.3]      = {D_star*100:.4f} cm")
        print(f"  α_cond  [A.4]      = {alpha_cond:.3f}  (para D_max={D_MAX} m)")
        print(f"  ṁ_min   [A.5]      = {dm_min_pre:.3f} kg/s")
        print(f"  ṁ_post0 (reactivo) = {dm_post0:.3f} kg/s")
        regime = "Alta inercia (Cardano)" if (C_sys**2/4 - C_sys**4/27) > 0 \
                 else "Baja inercia (Casus Irreducibilis)"
        print(f"  C_sys / Régimen    = {C_sys:.4f}  [{regime}]")
        print(f"  α*      [A.6]      = {alpha_star:.4f}  (fracción óptima Minimax)")
        print(f"  ṁ_opt   [A.6]      = {dm_opt:.3f} kg/s")
        print(f"\n  t_anticipo por drop (recalculado por masa):")
        for k in range(n_entregas):
            M_k = p.m_estructura + (p.m_carga - k*m_drop_kg) + m_H2_ini
            print(f"    Drop {k+1}: M≈{M_k:.0f} kg  "
                  f"→ t_anticipo={t_anticipos[k]:.3f} s  "
                  f"(precond. inicia t={t_pre_starts[k]:.2f} s)")
        print(f"\n  → Comparación vs hardcoded:")
        print(f"     dm_H2_proactivo : -10.0 kg/s  →  -{dm_pre:.2f} kg/s")
        print(f"     t_anticipo      :   2.0 s     →    {t_anticipos[0]:.2f} s (drop 1)")
        print("="*66 + "\n")

    # ── Simulación ODE ────────────────────────────────────────────────────────
    fase2               = [False]
    fase2_at_drop       = [0]
    ultimo_drop_anun    = -1
    precond_notif       = [False] * n_entregas   # evitar prints repetidos por paso ODE

    def dinamica(t, y):
        nonlocal ultimo_drop_anun
        v, h, m_H2, T_H2, Q_con, Q_gen = y
        m_H2 = max(m_H2, 0.1)

        drops_hechos  = min(int(t // intervalo), n_entregas)
        m_carga_actual = p.m_carga - drops_hechos * (toneladas_por_drop * 1000.0)
        M_total_nave   = p.m_estructura + m_carga_actual + m_H2

        if drops_hechos > ultimo_drop_anun:
            if verbose and drops_hechos > 0:   # ignorar el "drop 0" en t=0
                print(f"  [t={t:6.1f}s] DROP #{drops_hechos}  "
                      f"carga restante: {m_carga_actual/1000:.1f} t")
            ultimo_drop_anun = drops_hechos

        rho_aire = p.rho_aire_std * np.exp(-h / p.H_escala)
        V_H2     = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        V_eq     = ((p.m_estructura + m_carga_actual) * p.R_H2 * p.T1) / \
                   (p.P_atmos * (rho_h0 / rho_H2_ini - 1))

        # Controlador PD con ganancia adaptativa (idéntico a t2_gain_scheduled)
        zeta    = 1.35
        Kd_auto = 2 * zeta * np.sqrt(p.Kp * (p.m_estructura + m_carga_actual)
                                      / (rho_aire * p.g))
        V_cmd   = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        V_error = V_cmd - V_H2

        # Gain scheduling: Fase 1 (emergencia) → Fase 2 (trim)
        if fase2[0] and drops_hechos > fase2_at_drop[0]:
            fase2[0] = False
        if not fase2[0] and drops_hechos > 0:
            carga_rest  = (n_entregas - drops_hechos) * toneladas_por_drop * 1000.0
            m_H2_eq_now = (p.m_estructura + carga_rest) / (rho_h0 / rho_H2_ini - 1)
            if m_H2 <= m_H2_eq_now:
                fase2[0]        = True
                fase2_at_drop[0] = drops_hechos

        K_mass_eff = 0.5 if fase2[0] else p.K_mass_control
        dm = K_mass_eff * V_error
        if drops_hechos == 0:
            dm = 0.0

        # Pulso proactivo — ventanas de tiempo puro (sin depender de drops_hechos)
        # Las ventanas [t_pre_starts[k], t_drops[k]) son no solapantes por diseño,
        # por lo que el loop activa a lo sumo un pulso por instante de tiempo.
        for k in range(n_entregas):
            if t_pre_starts[k] <= t < t_drops[k]:
                dm = -dm_pre
                if verbose and not precond_notif[k]:
                    print(f"  [t={t:6.2f}s] PRECOND #{k+1} activado  "
                          f"dm={-dm_pre:.1f} kg/s  t_anticipo={t_anticipos[k]:.3f}s  "
                          f"→ drop en t={t_drops[k]:.1f}s")
                    precond_notif[k] = True
                break

        # Saturar en dm_opt, NO en DM_MAX — coherencia con el Minimax [A.6]:
        # dm_opt es la capacidad mínima instalada que balancea proactivo y reactivo.
        # Usar DM_MAX=30 kg/s haría que el controlador reactivo ignorase la optimización
        # y saturase ahí, desconectando la simulación de la teoría del apéndice.
        dm = np.clip(dm, -dm_opt, dm_opt)

        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * DIRIGIBLE['Area'] * rho_aire * v**2 * np.sign(v)
        dv    = (F_emp - M_total_nave * p.g - F_arr) / M_total_nave
        dQgen = abs(dm) * p.Q_comb_H2 * p.rend_gen if dm < -0.05 else 0.0

        return [dv, v, dm, (p.T1 - T_H2)/5.0, 0.0, dQgen]

    y0  = [0.0, p.h0, m_H2_ini, p.T1, 0.0, 0.0]
    sol = solve_ivp(dinamica, [0, t_final], y0, method='RK45', max_step=0.05)

    t_arr = sol.t
    v_arr, h_arr, m_H2_arr = sol.y[0], sol.y[1], sol.y[2]
    Q_gen_arr = sol.y[5]
    a_arr     = np.gradient(v_arr, t_arr)
    dm_arr    = np.gradient(m_H2_arr, t_arr)

    if graficar:
        if n_entregas == 1:
            _graficar_base(t_arr, h_arr, v_arr, a_arr, dm_arr, m_H2_arr,
                           Q_gen_arr/3.6e6, p, t_lib=intervalo,
                           dm_pre=dm_pre, t_anticipo=t_anticipos[0],
                           dm_opt=dm_opt, alpha_star=alpha_star)
        else:
            _graficar_multi(t_arr, h_arr, v_arr, a_arr, dm_arr, m_H2_arr,
                            Q_gen_arr/3.6e6, p, n_entregas=n_entregas,
                            intervalo=intervalo, tons_p=toneladas_por_drop,
                            dm_pre=dm_pre, t_anticipo=t_anticipos[0],
                            t_pre_starts=t_pre_starts)

    return {'Energia_kWh': Q_gen_arr[-1]/3.6e6,
            'dm_pre': dm_pre, 't_anticipo': t_anticipos[0],
            'alpha_star': alpha_star, 'dm_opt': dm_opt}


# =============================================================================
# Gráfica
# =============================================================================

def _graficar_base(t, h, v, a, dm, m_h2, E, p, t_lib,
                   dm_pre, t_anticipo, dm_opt, alpha_star):

    FIGSIZE  = (10, 12)
    FS_TITLE = 14
    FS_LEG   = 10
    LW       = 2.2

    COL_ALT   = 'steelblue'
    COL_VEL   = 'forestgreen'
    COL_ACC   = 'indianred'
    COL_FLUJO = 'firebrick'
    COL_MASA  = 'darkorchid'
    COL_ENERG = 'darkorange'
    COL_LIB   = 'black'

    plt.rc('axes', labelsize=11)

    fig, axes = plt.subplots(3, 1, figsize=FIGSIZE, sharex=True,
                             gridspec_kw={'height_ratios': [2.5, 1.5, 1]})
    plt.subplots_adjust(hspace=0.15)

    # Panel 1: Cinemática
    ax_alt = axes[0]
    ax_vel = ax_alt.twinx()
    ax_acc = ax_alt.twinx()
    ax_acc.spines["right"].set_position(("axes", 1.12))
    ax_acc.spines["right"].set_visible(True)

    ax_alt.plot(t, h, color=COL_ALT, lw=LW,       label='Altitude [m]')
    ax_vel.plot(t, v, color=COL_VEL, lw=LW, ls='--', label='Velocity [m/s]')
    ax_acc.plot(t, a, color=COL_ACC, lw=1.8, ls='-.', alpha=0.7, label='Accel [m/s²]')

    ax_alt.set_xscale('linear')
    ax_alt.set_xlim(0, t[-1])   # simulación completa desde t=0
    ax_alt.axhline(p.h0, color=COL_ALT, ls=':', alpha=0.5)
    ax_alt.axvline(t_lib, color=COL_LIB, ls='--', lw=1.5, label=f'Release ({t_lib} s)')

    subtitle = (
        f'$\\dot{{m}}_{{\\mathrm{{pre}}}}={dm_pre:.1f}$ kg/s'
        f'  |  $t_{{\\mathrm{{anticip}}}}={t_anticipo:.2f}$ s'
        f'  |  $\\alpha^*={alpha_star:.3f}$'
        f'  |  $\\dot{{m}}_{{\\mathrm{{opt}}}}={dm_opt:.2f}$ kg/s'
        f'  [Appendix A]'
    )
    ax_alt.set_title(f'Dynamic Analysis — {p.nombre_dirigible}\n{subtitle}',
                     fontsize=FS_TITLE, pad=8)
    ax_alt.set_ylabel('Altitude [m]',      color=COL_ALT, weight='bold')
    ax_vel.set_ylabel('Velocity [m/s]',    color=COL_VEL, weight='bold')
    ax_acc.set_ylabel('Accel [m/s²]',      color=COL_ACC, weight='bold')
    ax_alt.grid(True, which='both', alpha=0.2)

    lines  = (ax_alt.get_legend_handles_labels()[0] +
              ax_vel.get_legend_handles_labels()[0] +
              ax_acc.get_legend_handles_labels()[0])
    labels = (ax_alt.get_legend_handles_labels()[1] +
              ax_vel.get_legend_handles_labels()[1] +
              ax_acc.get_legend_handles_labels()[1])
    ax_alt.legend(lines, labels, loc='upper right', fontsize=FS_LEG)

    # Panel 2: Flujo + Masa
    ax_f = axes[1]
    ax_m = ax_f.twinx()
    ax_f.plot(t, dm,   color=COL_FLUJO, lw=LW,  label='H₂ Flow [kg/s]')
    ax_m.plot(t, m_h2, color=COL_MASA,  lw=1.8, ls=':', marker='.', markersize=2,
              markevery=15, alpha=0.5, label='H₂ Mass [kg]')
    ax_f.axvline(t_lib, color=COL_LIB, ls='--', lw=1.2)
    ax_f.set_ylabel('H₂ Flow [kg/s]',  color=COL_FLUJO, weight='bold')
    ax_m.set_ylabel('H₂ Mass [kg]',    color=COL_MASA,  weight='bold')
    ax_f.grid(True, which='both', alpha=0.2)
    ax_f.legend([ax_f.lines[0], ax_m.lines[0]],
                ['H₂ Flow [kg/s]', 'H₂ Mass [kg]'],
                loc='lower right', fontsize=FS_LEG)

    # Panel 3: Energía
    axes[2].plot(t, E, color=COL_ENERG, lw=LW)
    axes[2].axvline(t_lib, color=COL_LIB, ls='--', lw=1.2)
    axes[2].set_ylabel('Energy [kWh]',  weight='bold')
    axes[2].set_xlabel('Time [s]', weight='bold')
    axes[2].grid(True, which='both', alpha=0.2)

    plt.tight_layout()
    nombre = DIRIGIBLE['nombre'].replace(' ', '_')
    path   = _siguiente_path(f"apendice_{nombre}_base")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Figura guardada: {path}")
    plt.show()


def _graficar_multi(t, h, v, a, dm, m_h2, E, p, n_entregas, intervalo,
                    tons_p, dm_pre, t_anticipo, t_pre_starts=None):
    """5-panel linear-scale figure for multi-drop missions."""

    FIGSIZE   = (12, 16)
    LW        = 2.2
    COL_ALT   = 'steelblue'
    COL_VEL   = 'forestgreen'
    COL_ACC   = 'indianred'
    COL_FLUJO = 'firebrick'
    COL_MASA  = 'darkorchid'
    COL_ENERG = 'darkorange'
    COL_LIB   = 'black'

    plt.rc('axes', labelsize=12)
    plt.rc('axes', titlesize=14)

    fig, axes = plt.subplots(5, 1, figsize=FIGSIZE, sharex=True)

    # Líneas verticales en cada drop
    t_drops = [intervalo * (k + 1) for k in range(n_entregas)]

    # Panel 1: Altitud
    axes[0].plot(t, h, COL_ALT, lw=LW)
    axes[0].axhline(p.h0, color=COL_ALT, ls=':', alpha=0.5)
    axes[0].set_ylabel('Altitude [m]', weight='bold')
    axes[0].set_title(
        f"Multi-Drop Mission: {n_entregas} releases of {tons_p} t  —  "
        f"{DIRIGIBLE['nombre']}\n"
        f"$\\dot{{m}}_{{\\mathrm{{pre}}}}={dm_pre:.1f}$ kg/s"
        f"  |  $t_{{\\mathrm{{anticip}}}}={t_anticipo:.2f}$ s"
        f"  [Appendix A]",
        fontsize=13)

    # Panel 2: Velocidad + Aceleración
    ax_v = axes[1]
    ax_a = ax_v.twinx()
    l_v, = ax_v.plot(t, v, COL_VEL, lw=LW)
    l_a, = ax_a.plot(t, a, COL_ACC, lw=2.0, ls='--')
    ax_v.set_ylabel('Velocity [m/s]',    color=COL_VEL, weight='bold')
    ax_a.set_ylabel('Accel [m/s²]',      color=COL_ACC, weight='bold')
    vl = max(abs(ax_v.get_ylim()[0]), abs(ax_v.get_ylim()[1]))
    al = max(abs(ax_a.get_ylim()[0]), abs(ax_a.get_ylim()[1]))
    ax_v.set_ylim(-vl*1.1, vl*1.1)
    ax_a.set_ylim(-al*2.2, al*2.2)
    ax_v.legend([l_v, l_a], ['Velocity [m/s]', 'Accel [m/s²]'],
                loc='upper right', fontsize=9)

    # Panel 3: Flujo H2 — con marcadores de inicio de precondicionamiento
    axes[2].plot(t, dm, COL_FLUJO, lw=LW)
    axes[2].set_ylabel('H₂ Flow [kg/s]', weight='bold')
    if t_pre_starts:
        for k, tps in enumerate(t_pre_starts):
            axes[2].axvline(tps, color='royalblue', ls=':', lw=1.4, alpha=0.8)
            axes[2].annotate(f'P{k+1}\n{tps:.1f}s', xy=(tps, 0),
                             xytext=(tps + 0.3, -dm_pre * 0.35),
                             fontsize=7, color='royalblue', alpha=0.9)

    # Panel 4: Masa H2
    axes[3].plot(t, m_h2, COL_MASA, lw=LW)
    axes[3].set_ylabel('H₂ Mass [kg]', weight='bold')

    # Panel 5: Energía
    axes[4].plot(t, E, COL_ENERG, lw=LW)
    axes[4].set_ylabel('Energy [kWh]',   weight='bold')
    axes[4].set_xlabel('Time [s]',       weight='bold')

    # Líneas de drop en todos los paneles
    for ax in axes:
        for td in t_drops:
            ax.axvline(td, color=COL_LIB, ls='--', lw=1.2, alpha=0.7)
        ax.grid(True, alpha=0.2)

    fig.tight_layout()
    nombre = DIRIGIBLE['nombre'].replace(' ', '_')
    path   = _siguiente_path(f"apendice_{nombre}_multi{n_entregas}")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Figura guardada: {path}")
    plt.show()


# =============================================================================
# Entrada principal
# =============================================================================

if __name__ == '__main__':
    # Misión base (1 drop)
    resultado = correr(n_entregas=1, toneladas_por_drop=1, graficar=True)
    print(f"\n  Energía generada: {resultado['Energia_kWh']:.1f} kWh")

    # Descomenta para probar multi-drop:
    # correr(n_entregas=4, toneladas_por_drop=1, graficar=True)
