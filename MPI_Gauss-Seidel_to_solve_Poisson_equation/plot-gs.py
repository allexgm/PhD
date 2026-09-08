import pandas as pd
import matplotlib.pyplot as plt

# 1. Cargar el archivo usando el comportamiento por defecto (separado por comas).
# Al no poner header=None, pandas usa la primera fila como encabezado 
# y convierte el resto a números automáticamente.
df = pd.read_csv('pot_gs.dat', sep=r'\s+')

# Extraer las columnas 1, 2 y 3
x = df.iloc[:, 0]
y = df.iloc[:, 1]
z = df.iloc[:, 2]

# 2. Configurar la figura
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 3. Crear el gráfico de superficie (malla de triángulos)
grafico = ax.plot_trisurf(x, y, z, cmap='viridis', edgecolor='none')

# Configurar etiquetas
ax.set_xlabel('Eje X')
ax.set_ylabel('Eje Y')
ax.set_zlabel('Eje Z')

# 4. Barra de color
barra_color = fig.colorbar(grafico, ax=ax, pad=0.1)
barra_color.set_label('Intensidad Z')

ax.view_init(elev=0, azim=-90)
plt.title('Superficie 3D')

plt.show()
plt.savefig('pot_gs.png', dpi=300)
