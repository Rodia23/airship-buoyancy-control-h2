# =============================================================================

# Para simular otro dirigible, edita la sección "DIRIGIBLE" y "CONTROL".
# =============================================================================


# =============================================================================
# CONSTANTES FÍSICAS
# =============================================================================

g             = 9.81       # [m/s²]      Aceleración gravitacional
R_univ        = 8.314      # [J/mol·K]   Constante universal de gases ideales
M_H2          = 0.002016   # [kg/mol]    Masa molar del hidrógeno
R_H2          = R_univ / M_H2   # [J/kg·K]  ≈ 4124
cp_H2         = 14300      # [J/kg·K]    Calor específico del H2
Q_comb_H2     = 120e6      # [J/kg]      PCI del H2 (~33.33 kWh/kg)
L_v_agua      = 2257e3     # [J/kg]      Calor latente de vaporización del agua
                            


# =============================================================================
# ATMÓSFERA
# =============================================================================

rho_aire_std  = 1.058      # [kg/m³]   Densidad del aire a h0 (~500 m s.n.m.)
H_escala      = 8500       # [m]       Altura de escala: rho(h) = rho_std·exp(-h/H_escala)
P_atmos       = 84700      # [Pa]      Presión atmosférica inicial (~1500 m)


# =============================================================================
# DIRIGIBLE — Pathfinder 3
# =============================================================================

nombre_dirigible = 'Pathfinder 3'

Area          = 4390       # [m²]   Área frontal efectiva (para cálculo de arrastre)
m_estructura  = 99000      # [kg]   Masa de estructura y equipos (sin carga)
m_carga       = 1000       # [kg]   Masa de la carga útil a liberar
h0            = 50         # [m]    Altura inicial y referencia del controlador


# =============================================================================
# AERODINÁMICA
# =============================================================================

Cd            = 0.5        # [—]    Coeficiente de arrastre (cuerpo elipsoidal)


# =============================================================================
# HIDRÓGENO — almacenamiento y ballonets
# =============================================================================

P_tanque      = 350e5      # [Pa]     Presión de almacenamiento (350 bar)
rho_H2_tanque = 27.8       # [kg/m³]  Densidad H2 comprimido a 350 bar, 25°C

m_descomp     = 2.0        # [kg/s]   Tasa de descompresión: flujo tanque → ballonets
                          # Ref: 0.1–0.5 kg/s (PEM), 0.5–5 kg/s (mediano), 5–20 (industrial)


# =============================================================================
# TÉRMICA — preacondicionamiento de ballonets
# =============================================================================

T1            = 281.0      # [K]   Temperatura inicial del H2 (~8°C)
T2            = 281.2      # [K]   Temperatura objetivo del calentamiento (+0.2°C)
                            # NOTA: ΔT = 0.2 K → efecto de flotabilidad casi nulo.
                            #   El efecto de control principal viene del flujo proactivo.

tau_term      = 35.0       # [s]   Constante de tiempo térmica del calentamiento de ballonets
                            # Reemplaza h_conv, n_ballonets, T_pared, A_superficie
                            # Ref: 10–120 s según tamaño del sistema


# =============================================================================
# ENERGÍA — generador
# =============================================================================

rend_gen      = 0.30       # [—]   Rendimiento del generador eléctrico de H2 (30%)
                            # Contexto: E_generada = flujo_H2 · Q_comb_H2 · rend_gen
                            # NOTA: unificado a 30% (valor del código).
                            #   El artículo declarará este valor explícitamente.

perdidas      = 0.30       # [—]   Fracción de pérdidas térmicas en el calentamiento


# =============================================================================
# CONDENSACIÓN — parámetro de capacidad del intercambiador
# =============================================================================

Q_cond        = 500e3      # [W]   Potencia de condensación asumida del intercambiador
                            # Reemplaza U, A, ΔTm (diseño del equipo queda fuera del artículo)
                            # Contexto: t_cond = m_H2O · L_v_agua / Q_cond
                            #           donde m_H2O = 9 · m_H2_quemado
                            # Ref: compacto pequeño 10–50 kW, auto 50–200 kW,
                            #      industrial 200–2000 kW
                            # Con Q_cond=500 kW y m_carga=1t → t_cond ≈ 51 min


# =============================================================================
# CONTROL — controlador PD de altitud
# =============================================================================

Kp              = 15      # [m³/m]      Ganancia proporcional
Kd              = 2400     # [m³·s/m]    Ganancia derivativa
K_mass_control  = 60       # [kg/s/m³]   Conversión volumen→flujo másico

max_flujo_H2    = 30      # [kg/s]   Límite teórico del controlador (no físicamente realizable)
                            # ADVERTENCIA: solo define el techo del controlador PD.
                            #   Para un generador real usar 0.1–5 kg/s.
                            #   El barrido paramétrico evalúa el rango completo.

V_delta_max     = 1000     # [m³]    Saturación del volumen comandado


# =============================================================================
# MANIOBRA PROACTIVA
# =============================================================================

dm_H2_proactivo = -10      # [kg/s]   Consumo de H2 antes de la liberación (proactivo)
cooling_factor  = 12       # [—]      Factor de aceleración del enfriamiento proactivo
t_mantener      = 2.0      # [s]      Duración de la fase de mantenimiento térmico


# =============================================================================
# SIMULACIÓN
# =============================================================================

acortar_tiempos = False    # Si True, escala tiempos térmicos por escala_tiempo
escala_tiempo   = 0.2
dt_max          = 0.05     # [s]   Paso máximo del integrador ODE
t_post_liberacion = 400    # [s]   Duración de la sim. después de liberar la carga


# =============================================================================
# ENERGÍA SOLAR — modelo paramétrico (dos parámetros)
# =============================================================================

rendimiento_efectivo = 2.0    # [kWh/m²/día]  HSP·η_panel·(1-pérdidas_sistema)
                               # Desfavorable (invierno, Lima): ~1.2
                               # Medio (Lima):                  ~2.0
                               # Favorable (verano, Lima):      ~2.8

A_paneles_fraccion   = 0.45   # [—]   Fracción del área superficial cubierta con paneles
                               # Rango típico: 0.30 – 0.60


# =============================================================================
# REGENERACIÓN — electrólisis
# =============================================================================

eta_electrolisis     = 0.65   # [—]       Eficiencia del banco de electrólisis (PEM: 0.55–0.75)
E_H2_electrolisis    = 45.0   # [kWh/kg]  Energía para producir 1 kg de H2 (39–50 kWh/kg)
