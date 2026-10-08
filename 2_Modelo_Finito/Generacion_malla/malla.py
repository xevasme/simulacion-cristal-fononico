import gmsh

# 1. Inicializar el motor geométrico
gmsh.initialize()
gmsh.model.add("placa_64_cilindros")

# 2. Definición de parámetros físicos y geométricos
radio = 0.012          
a = 0.035              
margen_izq = 0.5     
margen_der = 0.5     
margen_y = 0.05       

# 8 columnas en X = 7 separaciones 'a' + radios de los extremos
longitud_x = margen_izq + (7 * a) + (2 * radio) + margen_der
# 8 filas en Y = 7 separaciones 'a' + radios de los extremos
longitud_y = margen_y + (7 * a) + (2 * radio) + margen_y

# 3. Crear la superficie base (El dominio continuo de aire)
# occ = OpenCASCADE, el núcleo geométrico de Gmsh
placa_aire = gmsh.model.occ.addRectangle(0, 0, 0, longitud_x, longitud_y)

## 4. Construir la matriz rectangular de 16 cilindros
cilindros = []
x_inicial = margen_izq + radio
y_inicial = margen_y + radio

for fila in range(8):         # AUMENTO A 4 FILAS
    for columna in range(8):  # MANTIENE 4 COLUMNAS
        centro_x = x_inicial + (columna * a)
        centro_y = y_inicial + (fila * a)
        
        disco = gmsh.model.occ.addDisk(centro_x, centro_y, 0, radio, radio)
        cilindros.append((2, disco))

gmsh.model.occ.synchronize()

# 5. Operación Booleana: Restar los cilindros a la placa de aire
# cut(Objeto_principal, Herramientas_de_corte)
gmsh.model.occ.cut([(2, placa_aire)], cilindros)
gmsh.model.occ.synchronize()

# --- INICIO DEL NUEVO BLOQUE(Sebastian): DEFINICIÓN DE GRUPOS FÍSICOS ---

# Obtener la superficie resultante después de restar los cilindros
superficies = gmsh.model.getEntities(2)
tag_superficie = superficies[0][1] # Extraemos el ID de la superficie principal

# Obtener todas las líneas (fronteras de dimensión 1) que limitan esta superficie
fronteras = gmsh.model.getBoundary(superficies, oriented=False)

inlet_tags = []
outlet_tags = []
paredes_tags = []
cilindros_tags = []

eps = 1e-4 # Tolerancia para errores de precisión de punto flotante

# Clasificar cada línea analizando sus coordenadas espaciales
for dim, tag in fronteras:
    # getBoundingBox devuelve (xmin, ymin, zmin, xmax, ymax, zmax)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    
    if abs(xmin) < eps and abs(xmax) < eps:
        # La curva está pegada al eje X=0
        inlet_tags.append(tag)
    elif abs(xmin - longitud_x) < eps and abs(xmax - longitud_x) < eps:
        # La curva está en el límite derecho (X = longitud_x)
        outlet_tags.append(tag)
    elif (abs(ymin) < eps and abs(ymax) < eps) or (abs(ymin - longitud_y) < eps and abs(ymax - longitud_y) < eps):
        # La curva está en el piso (Y=0) o en el techo (Y=longitud_y)
        paredes_tags.append(tag)
    else:
        # Si no es ninguna de las fronteras externas, es un recorte interno (cilindro)
        cilindros_tags.append(tag)

# Crear los Grupos Físicos (Dimensión 1)
inlet_group = gmsh.model.addPhysicalGroup(1, inlet_tags)
gmsh.model.setPhysicalName(1, inlet_group, "Inlet")

outlet_group = gmsh.model.addPhysicalGroup(1, outlet_tags)
gmsh.model.setPhysicalName(1, outlet_group, "Outlet")

paredes_group = gmsh.model.addPhysicalGroup(1, paredes_tags)
gmsh.model.setPhysicalName(1, paredes_group, "Paredes_Exteriores")

cilindros_group = gmsh.model.addPhysicalGroup(1, cilindros_tags)
gmsh.model.setPhysicalName(1, cilindros_group, "Cilindros")

# Crear el Grupo Físico para el fluido (Dimensión 2)
fluido_group = gmsh.model.addPhysicalGroup(2, [tag_superficie])
gmsh.model.setPhysicalName(2, fluido_group, "Dominio_Aire")

# --- FIN DEL NUEVO BLOQUE(sebastian) ---

# 6. Configurar la resolución de la malla
# Una regla en acústica es tener al menos 6-10 elementos por longitud de onda
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.002)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.008)

# 7. Generar y exportar la malla
gmsh.model.mesh.generate(2)
gmsh.write("placa_64_cilindros.msh")
gmsh.finalize()

print("Malla exportada exitosamente como 'placa_64_cilindros.msh'")