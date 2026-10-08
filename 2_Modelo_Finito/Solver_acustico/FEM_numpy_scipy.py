import meshio
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# 1. LECTURA DE LA MALLA Y GRUPOS FÍSICOS
malla = meshio.read("placa_64_cilindros.msh")
nodos = malla.points[:, :2]
elementos = malla.cells_dict["triangle"]
num_nodos = len(nodos)
num_elementos = len(elementos) # <--- Variable definida desde el principio

lineas = malla.cells_dict["line"]
inlet_lines = lineas[malla.cell_sets_dict["Inlet"]["line"]]
nodos_entrada = np.unique(inlet_lines)

# 2. PARÁMETROS FÍSICOS
frecuencia = 2490.0 
c0 = 343.0 
k0 = (2 * np.pi * frecuencia) / c0 

# 3. ENSAMBLAJE GLOBAL VECTORIZADO
x = nodos[elementos, 0]
y = nodos[elementos, 1]

x21 = x[:, 1] - x[:, 0]
x31 = x[:, 2] - x[:, 0]
y21 = y[:, 1] - y[:, 0]
y31 = y[:, 2] - y[:, 0]

area = 0.5 * np.abs(x21 * y31 - x31 * y21)

b = np.column_stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]]) / (2 * area[:, None])
c = np.column_stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]]) / (2 * area[:, None])

# --- INICIO ZONA PML (Absorbing Layer) ---
xc = np.mean(x, axis=1)
yc = np.mean(y, axis=1)

xmin, xmax = np.min(nodos[:,0]), np.max(nodos[:,0])
ymin, ymax = np.min(nodos[:,1]), np.max(nodos[:,1])

# El PML en X solo se aplica fuertemente al final del dominio (Outlet)
# para no ahogar la fuente que inyecta en X=0
dist_x = np.maximum(0, xmax - xc) 
dist_y = np.minimum(yc - ymin, ymax - yc)

grosor_pml = 0.20  
sigma_max = 50.0   

sigma_x = np.zeros(num_elementos)
sigma_y = np.zeros(num_elementos)

mask_x = dist_x < grosor_pml
sigma_x[mask_x] = sigma_max * ((grosor_pml - dist_x[mask_x]) / grosor_pml)**2

mask_y = dist_y < grosor_pml
sigma_y[mask_y] = sigma_max * ((grosor_pml - dist_y[mask_y]) / grosor_pml)**2

k_local = k0 - 1j * (sigma_x + sigma_y)
# --- FIN ZONA PML ---

num_elementos = len(elementos)
filas = np.zeros((num_elementos, 3, 3), dtype=int)
columnas = np.zeros((num_elementos, 3, 3), dtype=int)
# IMPORTANTE: valores_M debe ser complejo para no perder el término iterativo del PML
valores_K = np.zeros((num_elementos, 3, 3), dtype=float)
valores_M = np.zeros((num_elementos, 3, 3), dtype=complex) 

M_ref = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]]) / 12.0

for i in range(3):
    for j in range(3):
        filas[:, i, j] = elementos[:, i]
        columnas[:, i, j] = elementos[:, j]
        valores_K[:, i, j] = (b[:, i] * b[:, j] + c[:, i] * c[:, j]) * area
        valores_M[:, i, j] = M_ref[i, j] * area * (k_local**2)

filas = filas.flatten()
columnas = columnas.flatten()
valores_K = valores_K.flatten()
valores_M = valores_M.flatten()

K_global = sp.coo_matrix((valores_K, (filas, columnas)), shape=(num_nodos, num_nodos)).tocsr()
M_global_mod = sp.coo_matrix((valores_M, (filas, columnas)), shape=(num_nodos, num_nodos)).tocsr()

A = K_global - M_global_mod
b_vector = np.zeros(num_nodos, dtype=complex)

# 4. CONDICIONES DE FRONTERA (Método de Penalización)
penalizacion = 1e15

for nodo in nodos_entrada:
    A[nodo, nodo] += penalizacion 
    b_vector[nodo] = 1.0 * penalizacion

# 5. SOLUCIÓN DIRECTA DEL SISTEMA
presion_compleja = spla.spsolve(A, b_vector)

# 6. POST-PROCESAMIENTO Y VISUALIZACIÓN TERMOGRÁFICA
presion_real = np.real(presion_compleja)

plt.figure(figsize=(10, 5))
mapa = plt.tripcolor(nodos[:, 0], nodos[:, 1], elementos, presion_real, 
                     shading='gouraud', cmap='RdBu_r')

plt.colorbar(mapa, label='Presión Acústica (Pa)')
plt.title(f'Campo de Presión Acústica a {frecuencia} Hz con PML')
plt.xlabel('Eje X (m)')
plt.ylabel('Eje Y (m)')
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()