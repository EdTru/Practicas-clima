"""
Simulación de la Histéresis de las Capas de Hielo de la Antártida
Basado en los mecanismos descritos en el informe:
  1. Retroalimentación Elevación-Balance de Masa (Height-SMB Feedback)
  2. Inestabilidad de las capas de hielo marinas (MISI)
  3. Retroalimentación hielo-albedo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# PALETA Y ESTILO
# ─────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.facecolor': '#0a0f1e',
    'figure.facecolor': '#060b18',
    'text.color': '#e8eaf6',
    'axes.labelcolor': '#b0bec5',
    'xtick.color': '#78909c',
    'ytick.color': '#78909c',
    'axes.edgecolor': '#1e3a5f',
    'grid.color': '#1e3a5f',
    'grid.alpha': 0.4,
})

ICE_CMAP = LinearSegmentedColormap.from_list(
    'ice_loss',
    ['#cce7ff', '#7ec8e3', '#2196f3', '#0d47a1', '#880000', '#3e0000'],
    N=256
)

# ─────────────────────────────────────────────────────────────────────────────
# MODELO FÍSICO DE HISTÉRESIS
# ─────────────────────────────────────────────────────────────────────────────

class AntarcticIceModel:
    """
    Modelo simplificado de la Antártida Occidental (WAIS) que integra:
      - Height-SMB Feedback  (ecuación BM = SMB − D)
      - MISI: flujo q ∝ h^4.7
      - Retroalimentación hielo-albedo
    """
    def __init__(self, nx=80, ny=80):
        self.nx, self.ny = nx, ny
        self.dx = 50e3          # 50 km por celda
        self.dt = 0.5           # años por paso

        # Cuenca batimétrica de la Antártida Occidental (pendiente retrógrada)
        x = np.linspace(-1, 1, nx)
        y = np.linspace(-1, 1, ny)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)

        # Topografía: cuenca con forma de "cuenco" en el oeste
        self.bed = -800 + 1200 * R**2 - 400 * np.exp(-((X+0.3)**2 + Y**2) / 0.3)
        # Zona oriental más alta
        self.bed += 600 * np.clip(X, 0, 1)
        # Máscara "marina": celdas con lecho por debajo del nivel del mar
        self.marine_mask = (self.bed < 0).astype(float)

        # Condición inicial: capa de hielo completa
        self.h = np.where(R < 0.95, 2500 * (1 - R**2), 0).clip(0)
        self.h_init = self.h.copy()

        # Estado físico
        self.time = 0.0
        self.T_forcing = 0.0         # forzamiento de temperatura (ΔT sobre preindustrial)
        self.grounding_line = R.copy()
        self.X, self.Y, self.R = X, Y, R

    # ── Mecanismo 1: Height-SMB Feedback ──────────────────────────────────
    def smb(self, h):
        """Balance de masa superficial: disminuye al bajar la elevación."""
        elevation = self.bed + h
        # Gradiente térmico: –6.5 °C / km
        T_local = self.T_forcing + 6.5e-3 * (3000 - elevation)  # T relativa al nivel medio
        # SMB positivo si frío, negativo si cálido
        smb_base = 0.5 - 0.35 * T_local   # m/año equivalente
        return smb_base

    # ── Mecanismo 2: MISI  q ∝ h^4.7 ──────────────────────────────────────
    def dynamic_discharge(self, h):
        """Descarga dinámica en la línea de apoyo (proporcional a h^4.7)."""
        h_grounding = np.where(self.bed + h < 0, h, 0)  # sólo zona marina
        D = 2e-11 * h_grounding**4.7 * self.marine_mask
        return D

    # ── Mecanismo 3: Retroalimentación hielo-albedo ────────────────────────
    def albedo_feedback(self, h):
        """Reducción del albedo al perder hielo → más absorción → más derretimiento."""
        ice_fraction = np.clip(h / 500, 0, 1)
        alpha = 0.25 + 0.55 * ice_fraction          # albedo: 0.25 (océano) → 0.80 (nieve)
        Q_absorbed = (1 - alpha) * 340              # W/m² absorbidos
        melt_albedo = np.clip((Q_absorbed - 200) * 0.002, 0, 2)  # m/año extra de fusión
        return melt_albedo

    # ── Paso temporal ──────────────────────────────────────────────────────
    def step(self):
        h = self.h
        SMB  = self.smb(h)
        D    = self.dynamic_discharge(h)
        M_alb = self.albedo_feedback(h)

        dh = (SMB - D - M_alb) * self.dt
        self.h = np.clip(h + dh, 0, None)
        self.time += self.dt

    # ── Métricas globales ──────────────────────────────────────────────────
    @property
    def ice_volume(self):
        return float(np.sum(self.h) * self.dx**2 * 1e-9)  # km³ × 10⁻⁹ → normalizado

    @property
    def ice_fraction(self):
        return float(np.sum(self.h > 10)) / (self.nx * self.ny)

    @property
    def grounding_line_position(self):
        """Radio medio de la línea de apoyo (donde h > 0 y zona marina)."""
        marine_ice = (self.h > 50) & (self.marine_mask > 0.5)
        if marine_ice.any():
            return float(np.mean(self.R[marine_ice]))
        return 0.0

# ─────────────────────────────────────────────────────────────────────────────
# CURVA DE HISTÉRESIS (calentamiento → enfriamiento)
# ─────────────────────────────────────────────────────────────────────────────

def compute_hysteresis_curve(n_steps=400):
    """Calcula la curva de histéresis completa avanzando y retrocediendo T."""
    T_warm = np.linspace(0, 4.5, n_steps // 2)
    T_cool = np.linspace(4.5, 0, n_steps // 2)
    T_path = np.concatenate([T_warm, T_cool])

    model = AntarcticIceModel(nx=50, ny=50)
    volumes_w, volumes_c = [], []

    for i, T in enumerate(T_path):
        model.T_forcing = T
        for _ in range(8):
            model.step()
        if i < n_steps // 2:
            volumes_w.append(model.ice_volume)
        else:
            volumes_c.append(model.ice_volume)

    return T_warm, T_cool, volumes_w, volumes_c

# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZACIÓN PRINCIPAL
# ─────────────────────────────────────────────────────────────────────────────

def build_figure():
    fig = plt.figure(figsize=(20, 14), facecolor='#060b18')
    fig.suptitle(
        'Histéresis de las Capas de Hielo Antárticas\n'
        'Simulación de retroalimentaciones físicas: Height-SMB · MISI · Albedo',
        fontsize=17, fontweight='bold', color='#e8eaf6', y=0.97,
        fontfamily='DejaVu Sans'
    )

    gs = gridspec.GridSpec(
        3, 4,
        figure=fig,
        hspace=0.45, wspace=0.35,
        left=0.05, right=0.96, top=0.92, bottom=0.06
    )

    # ── Mapa de la Antártida (proyección polar manual) ────────────────────
    ax_map = fig.add_subplot(gs[0:2, 0:2], projection='polar')
    ax_map.set_facecolor('#0d2137')
    ax_map.set_theta_zero_location('N')
    ax_map.set_theta_direction(-1)

    # Círculo polar antártico base (hielo)
    theta_full = np.linspace(0, 2*np.pi, 500)
    r_ice = 0.85   # radio del casquete polar
    ax_map.fill(theta_full, np.full_like(theta_full, r_ice),
                color='#b3e5fc', alpha=0.30, zorder=2)

    # Contorno de la Antártida (forma simplificada)
    # Coordenadas aproximadas del contorno antártico en proyección polar
    theta_coast = np.linspace(0, 2*np.pi, 200)
    r_coast = 0.70 + 0.12 * np.sin(3 * theta_coast) + 0.06 * np.cos(5 * theta_coast)
    ax_map.fill(theta_coast, r_coast, color='#e0e8f0', alpha=0.85, zorder=3)

    # Interior del continente (tierra/roca bajo el hielo)
    r_interior = 0.55 + 0.08 * np.sin(3 * theta_coast) + 0.04 * np.cos(5 * theta_coast)
    ax_map.fill(theta_coast, r_interior, color='#8d9e8d', alpha=0.6, zorder=4)

    # WAIS: Antártida Occidental (sector ~90°W a 180°, aprox theta pi/2 a pi)
    theta_wais = np.linspace(np.pi * 0.5, np.pi * 1.25, 80)
    r_wais_out = 0.72 + 0.10 * np.sin(3 * theta_wais)
    r_wais_in  = np.zeros_like(theta_wais)
    ax_map.fill_between(theta_wais, r_wais_in, r_wais_out,
                        color='#ff5722', alpha=0.35, zorder=5, label='WAIS')

    # Línea de apoyo WAIS (zona costera marina)
    r_gl = 0.58 + 0.07 * np.sin(3 * theta_wais)
    ax_map.plot(theta_wais, r_gl, color='#ffee58', linewidth=1.8,
                linestyle='--', zorder=6, label='Línea de apoyo')

    # Meridiano de referencia y etiquetas
    for ang_deg, label in [(0, 'N'), (90, 'O'), (180, 'S'), (270, 'E')]:
        ang = np.radians(ang_deg)
        ax_map.plot([ang, ang], [0, 1.0], color='#1e3a5f',
                    linewidth=0.5, linestyle=':', alpha=0.6, zorder=1)

    # Círculos de latitud
    for r_lat, lat_label in [(0.5, '80°S'), (0.75, '70°S'), (1.0, '60°S')]:
        ax_map.plot(theta_full, np.full_like(theta_full, r_lat),
                    color='#1e3a5f', linewidth=0.5, linestyle=':', alpha=0.6)
        ax_map.text(0.1, r_lat, lat_label, color='#546e7a', fontsize=6.5, va='bottom')

    # Punto del Polo Sur
    ax_map.plot(0, 0, 'o', color='#4dd0e1', markersize=5, zorder=7)
    ax_map.text(0.3, 0.05, 'Polo Sur', color='#4dd0e1', fontsize=7)

    # Etiqueta WAIS
    ax_map.text(np.radians(230), 0.62, 'WAIS', color='#ff8a65',
                fontsize=9, fontweight='bold', ha='center')
    ax_map.text(np.radians(60),  0.62, 'EAIS', color='#b0bec5',
                fontsize=9, ha='center')

    ax_map.set_rticks([])
    ax_map.set_xticks([])
    ax_map.set_ylim(0, 1.05)
    ax_map.set_title('Antártida — Vista Polar Sur\n(proyección estereográfica)', 
                     color='#b0bec5', fontsize=10, pad=12)

    # Leyenda del mapa
    wais_patch = mpatches.Patch(color='#ff5722', alpha=0.5, label='WAIS (pendiente retrógrada)')
    ice_patch  = mpatches.Patch(color='#e0e8f0', alpha=0.85, label='Casquete de hielo')
    gl_line    = plt.Line2D([0], [0], color='#ffee58', linestyle='--', linewidth=1.8,
                            label='Línea de apoyo')
    ax_map.legend(handles=[wais_patch, ice_patch, gl_line], loc='lower left',
                  bbox_to_anchor=(0.0, -0.05),
                  fontsize=7.5, facecolor='#0a0f1e', edgecolor='#1e3a5f',
                  labelcolor='#e8eaf6')

    return fig, gs


def run_simulation_and_plot():
    print("⏳ Calculando curva de histéresis...")
    T_warm, T_cool, vol_w, vol_c = compute_hysteresis_curve(n_steps=300)

    print("⏳ Ejecutando simulación espacio-temporal...")
    model = AntarcticIceModel(nx=70, ny=70)

    n_snapshots = 6
    T_values    = [0.0, 0.8, 1.5, 2.5, 3.4, 4.0]
    snapshots   = []
    times_snap  = []
    gl_snap     = []

    for T_target in T_values:
        model.T_forcing = T_target
        for _ in range(200):
            model.step()
        snapshots.append(model.h.copy())
        times_snap.append(model.time)
        gl_snap.append(model.grounding_line_position)

    # ── Layout principal ──────────────────────────────────────────────────
    fig, gs = build_figure()

    # ── Curva de histéresis ───────────────────────────────────────────────
    ax_hyst = fig.add_subplot(gs[0:2, 2:4])
    ax_hyst.set_facecolor('#0a0f1e')

    v_max = max(max(vol_w), max(vol_c)) * 1.05

    # Relleno entre ramas
    T_common = np.linspace(0, 4.5, 150)
    vol_w_interp = np.interp(T_common, T_warm, vol_w)
    vol_c_interp = np.interp(T_common[::-1], T_cool, vol_c)[::-1]
    ax_hyst.fill_between(T_common, vol_w_interp, vol_c_interp,
                         alpha=0.12, color='#ff7043', label='Área de histéresis')

    ax_hyst.plot(T_warm, vol_w, color='#ef5350', linewidth=2.5,
                 label='Calentamiento (pérdida de hielo)')
    ax_hyst.plot(T_cool, vol_c, color='#42a5f5', linewidth=2.5,
                 linestyle='--', label='Enfriamiento (recuperación imposible)')

    # Flechas de dirección
    mid = len(T_warm) // 2
    ax_hyst.annotate('', xy=(T_warm[mid+5], vol_w[mid+5]),
                     xytext=(T_warm[mid-5], vol_w[mid-5]),
                     arrowprops=dict(arrowstyle='->', color='#ef5350', lw=2))
    ax_hyst.annotate('', xy=(T_cool[mid+5], vol_c[mid+5]),
                     xytext=(T_cool[mid-5], vol_c[mid-5]),
                     arrowprops=dict(arrowstyle='->', color='#42a5f5', lw=2))

    # Umbrales críticos (Robinson et al. 2012)
    for T_crit, label, col in [(0.8, '0.8°C umbral inferior', '#ffb74d'),
                                (3.4, '3.4°C umbral superior', '#ff7043')]:
        ax_hyst.axvline(T_crit, color=col, linestyle=':', linewidth=1.5, alpha=0.8)
        ax_hyst.text(T_crit + 0.05, v_max * 0.95, label,
                     color=col, fontsize=8, rotation=90, va='top')

    ax_hyst.set_xlabel('ΔT (°C sobre niveles preindustriales)', fontsize=11)
    ax_hyst.set_ylabel('Volumen de hielo (u.a.)', fontsize=11)
    ax_hyst.set_title('Curva de Histéresis: WAIS', color='#e8eaf6', fontsize=12, pad=8)
    ax_hyst.legend(fontsize=8.5, facecolor='#0a0f1e', edgecolor='#1e3a5f',
                   labelcolor='#e8eaf6', loc='upper right')
    ax_hyst.set_xlim(-0.1, 4.7)
    ax_hyst.set_ylim(0, v_max)
    ax_hyst.grid(True, alpha=0.25)

    # ── Snapshots espaciales ──────────────────────────────────────────────
    h_global_max = max(s.max() for s in snapshots)
    snap_titles  = [f'ΔT = {t}°C' for t in T_values]
    colors_snap  = ['#1565c0', '#1976d2', '#ffa726', '#f57c00', '#d32f2f', '#b71c1c']

    for idx in range(n_snapshots):
        row = 2
        col = idx
        if idx >= 4:
            row = 3 - 1   # solo 4 en la fila 2, los últimos van abajo — usamos los 4 primeros
        ax_s = fig.add_subplot(gs[2, idx if idx < 4 else idx - 4])
        ax_s.set_facecolor('#0a0f1e')

        snap_norm = snapshots[idx] / (h_global_max + 1e-6)
        im = ax_s.imshow(snap_norm, cmap=ICE_CMAP, vmin=0, vmax=1,
                         origin='lower', interpolation='bilinear')

        # Línea de apoyo aproximada
        gl_radius = gl_snap[idx]
        nx, ny = model.nx, model.ny
        cx, cy = nx / 2, ny / 2
        radius_px = gl_radius * min(nx, ny) / 2
        circle = plt.Circle((cx, cy), radius_px, fill=False,
                             color='#ffee58', linewidth=1.2, linestyle='--', alpha=0.8)
        ax_s.add_patch(circle)

        ax_s.set_title(snap_titles[idx], color=colors_snap[idx],
                       fontsize=9.5, fontweight='bold', pad=4)
        ax_s.set_xticks([]); ax_s.set_yticks([])
        ax_s.spines[:].set_edgecolor('#1e3a5f')

        # Porcentaje de hielo restante
        frac = np.sum(snapshots[idx] > 10) / (nx * ny) * 100
        ax_s.text(0.05, 0.05, f'{frac:.0f}% hielo', transform=ax_s.transAxes,
                  color='white', fontsize=7.5, fontweight='bold',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='#00000088', edgecolor='none'))

    # Colorbar para snapshots
    cbar_ax = fig.add_axes([0.97, 0.06, 0.012, 0.22])
    sm = plt.cm.ScalarMappable(cmap=ICE_CMAP, norm=mcolors.Normalize(vmin=0, vmax=1))
    cb = plt.colorbar(sm, cax=cbar_ax)
    cb.set_label('Espesor relativo', color='#b0bec5', fontsize=8)
    cb.ax.yaxis.set_tick_params(color='#78909c')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='#78909c', fontsize=7)

    # Línea guía: línea de apoyo
    ax_hyst.axhline(0, color='#546e7a', linewidth=0.7, alpha=0.5)
    leyenda_gl = mpatches.Patch(facecolor='none', edgecolor='#ffee58',
                                linestyle='--', linewidth=1.5,
                                label='Línea de apoyo (MISI)')

    # Anotación de mecanismos en el gráfico de histéresis
    ax_hyst.annotate('MISI activa\n(q ∝ h⁴·⁷)',
                     xy=(2.8, vol_w[int(len(vol_w)*0.62)]),
                     xytext=(3.3, vol_w[int(len(vol_w)*0.62)] * 1.3),
                     color='#ffcc80', fontsize=8,
                     arrowprops=dict(arrowstyle='->', color='#ffcc80', lw=1.2),
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a1a2e', alpha=0.8))

    ax_hyst.annotate('Retroalim.\nalbedo',
                     xy=(1.5, vol_c[int(len(vol_c)*0.35)]),
                     xytext=(0.4, vol_c[int(len(vol_c)*0.35)] * 0.7),
                     color='#80deea', fontsize=8,
                     arrowprops=dict(arrowstyle='->', color='#80deea', lw=1.2),
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a1a2e', alpha=0.8))

    # ── Nota sobre mecanismos ─────────────────────────────────────────────
    note = (
        "Mecanismos modelados  ①  Height-SMB: BM = SMB − D, con gradiente térmico −6.5°C/km  "
        "②  MISI: q ∝ h⁴·⁷ en pendiente retrógrada (Schoof 2007)  "
        "③  Albedo: α∈[0.25, 0.80] → forzamiento adicional de fusión\n"
        "Umbrales críticos (Robinson et al. 2012): 0.8°C – 3.4°C de calentamiento global "
        "│  Línea amarilla discontinua = línea de apoyo (grounding line)"
    )
    fig.text(0.50, 0.005, note, ha='center', va='bottom',
             color='#78909c', fontsize=7.5, style='italic',
             wrap=True)

    plt.savefig('/mnt/user-data/outputs/antarctic_hysteresis.png',
                dpi=160, bbox_inches='tight', facecolor='#060b18')
    print("✅ Figura guardada en /mnt/user-data/outputs/antarctic_hysteresis.png")

    # ── Guardar también el script limpio ──────────────────────────────────
    plt.close()


if __name__ == '__main__':
    run_simulation_and_plot()
