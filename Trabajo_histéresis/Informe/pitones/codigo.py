import numpy as np
import matplotlib.pyplot as plt

# 1. Crear una rejilla de coordenadas (Mapa de la Antártida simplificado)
n = 200
x = np.linspace(-1, 1, n)
y = np.linspace(-1, 1, n)
X, Y = np.meshgrid(x, y)

# 2. Definir la forma del continente (Radio aproximado)
radio_antartida = np.sqrt(X**2 + Y**2)
mascara_hielo = radio_antartida < 0.9

# 3. Simular la Tasa de Derretimiento basada en la Altitud y el Océano
# Los bordes (radio cercano a 0.9) tienen menor altitud y más contacto oceánico
tasa_derretimiento = np.exp(radio_antartida * 3) * mascara_hielo

# 4. Simular puntos críticos de Inestabilidad (MISI/MICI) en la Antártida Occidental
# Añadimos un "hotspot" en el cuadrante inferior izquierdo (Mar de Amundsen)
hotspot = np.exp(-((X + 0.5)**2 + (Y + 0.5)**2) / 0.1)
tasa_derretimiento += hotspot * 10 * mascara_hielo

# 5. Visualización del Mapa
plt.figure(figsize=(10, 8))
mapa = plt.imshow(tasa_derretimiento, cmap='YlOrRd', extent=[-1, 1, -1, 1])
plt.colorbar(mapa, label='Tasa de pérdida de masa (Simulada)')
plt.title('Simulación de Puntos Críticos de Deshielo en la Antártida', fontsize=14)
plt.axis('off') # Quitar ejes para que parezca un mapa real

# Añadir etiquetas de zonas críticas
plt.text(-0.8, -0.7, 'Antártida Occidental\n(Inestabilidad Marina)', color='black', fontweight='bold')
plt.text(0.3, 0, 'Meseta Central\n(Estable)', color='blue')

plt.show()