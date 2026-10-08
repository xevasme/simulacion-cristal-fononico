import meshio
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# 1. LECTURA DE LA MALLA
malla = meshio.read("placa_64_cilindros.msh")
nodos = malla.points[:, :2]           # Extraemos coordenadas X e Y
elementos = malla.cells_dict["triangle"] # Nodos que componen cada triángulo
num_nodos = len(nodos)

# 2. PARÁMETROS FÍSICOS (Ecuación de Helmholtz)
frecuencia = 2490.0   # Frecuencia de prueba en Hz
c0 = 343.0            # Velocidad del sonido en el aire (m/s)
# Se añade una pequeñísima parte imaginaria para simular leve atenuación del aire 
# y evitar resonancias numéricas infinitas al no tener una frontera PML perfecta
k = (2 * np.pi * frecuencia) / c0 + 0.02j 

# 3. INICIALIZACIÓN DE LA MATRIZ DISPERSA
filas = []
columnas = []
valores_K = [] # Rigidez (gradientes de presión)
valores_M = [] # Masa (término de inercia acústica)

# 4. ENSAMBLAJE GLOBAL
for elemento in elementos:
    # Coordenadas de los 3 vértices del triángulo
    x = nodos[elemento, 0]
    y = nodos[elemento, 1]
    
    # Cálculo del área (Determinante del Jacobiano)
    J = np.array([[x[1]-x[0], x[2]-x[0]], 
                  [y[1]-y[0], y[2]-y[0]]])
    area = 0.5 * abs(np.linalg.det(J))
    
    # Derivadas espaciales (Funciones de forma)
    b = np.array([y[1]-y[2], y[2]-y[0], y[0]-y[1]]) / (2*area)
    c = np.array([x[2]-x[1], x[0]-x[2], x[1]-x[0]]) / (2*area)
    
    # Matrices elementales de FEM
    K_local = (np.outer(b, b) + np.outer(c, c)) * area
    M_local = np.array([[2, 1, 1], 
                        [1, 2, 1], 
                        [1, 1, 2]]) * (area / 12.0)
    
    # Mapeo a las listas de la matriz global dispersa
    for i in range(3):
        for j in range(3):
            filas.append(elemento[i])
            columnas.append(elemento[j])
            valores_K.append(K_local[i, j])
            valores_M.append(M_local[i, j])

# Construcción eficiente de matrices dispersas (evita el colapso de memoria)
K_global = sp.coo_matrix((valores_K, (filas, columnas)), shape=(num_nodos, num_nodos)).tocsr()
M_global = sp.coo_matrix((valores_M, (filas, columnas)), shape=(num_nodos, num_nodos)).tocsr()

# Sistema de Helmholtz: [K - k^2 M] * P = F
A = K_global - (k**2) * M_global
b_vector = np.zeros(num_nodos, dtype=complex)

# 5. CONDICIONES DE FRONTERA
# Identificar nodos en el margen izquierdo (x cercano a 0) para inyectar la onda
tolerancia = 1e-3
nodos_entrada = np.where(nodos[:, 0] <= tolerancia)[0]

# Aplicar condición de Dirichlet (Presión = 1 Pa) en el emisor
for nodo in nodos_entrada:
    # Desconectar el nodo del resto del dominio vaciando su fila
    A.data[A.indptr[nodo]:A.indptr[nodo+1]] = 0 
    A[nodo, nodo] = 1.0                         
    b_vector[nodo] = 1.0 + 0j

# Nota: Los cilindros ya tienen condición de Neumann (velocidad normal nula) 
# de forma implícita gracias a la formulación variacional de FEM.

# 6. SOLUCIÓN DIRECTA DEL SISTEMA
# Inversión matemática mediante solver subyacente en C/Fortran
presion_compleja = spla.spsolve(A, b_vector)

# 7. POST-PROCESAMIENTO Y VISUALIZACIÓN TERMOGRÁFICA
presion_real = np.real(presion_compleja)

plt.figure(figsize=(10, 5))
# tripcolor interpola los valores nodales sobre la malla triangular
mapa = plt.tripcolor(nodos[:, 0], nodos[:, 1], elementos, presion_real, 
                     shading='gouraud', cmap='RdBu_r')

plt.colorbar(mapa, label='Presión Acústica (Pa)')
plt.title(f'Campo de Presión Acústica a {frecuencia} Hz')
plt.xlabel('Eje X (m)')
plt.ylabel('Eje Y (m)')
plt.gca().set_aspect('equal') # Mantiene la proporción geométrica real
plt.tight_layout()
plt.show()