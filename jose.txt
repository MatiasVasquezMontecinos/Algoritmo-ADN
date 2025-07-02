import pandas as pd
#import cupy as cp
import random
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

random.seed(42)
np.random.seed(42)

def calcular_matriz_distancias(posiciones):
    return np.linalg.norm(posiciones[:, None] - posiciones[None, :], axis=2)

def generar_poblacion_inicial(pob_max, n_diputados, quorum_q):
    poblacion = np.zeros((pob_max, n_diputados), dtype=np.int32)
    for i in range(pob_max):
        indices = random.sample(range(n_diputados), quorum_q)
        poblacion[i, indices] = 1
    return poblacion

def EvalCromo_batch(poblacion, D):
    # Calcula la suma de distancias dentro de cada coalición
    N = poblacion.shape[0]
    Z = np.zeros(N, dtype=np.float32)
    for i in range(N):
        idx = np.where(poblacion[i] == 1)[0]
        if len(idx) < 2:
            continue
        submat = D[idx[:, None], idx]
        Z[i] = np.sum(np.triu(submat, k=1))
    return Z

def cruzamiento_consenso(padre1, padre2):
    mask_iguales = padre1 == padre2
    hijo = np.copy(padre1)
    dif_pos = np.where(mask_iguales == 0)[0]
    for pos in dif_pos:
        hijo[pos] = padre1[pos] if random.random() < 0.5 else padre2[pos]
    return hijo

def mutacion_adaptativa(individuo, posiciones, quorum_q, prob_mut=0.1):
    if random.random() < prob_mut:
        idx_presentes = np.where(individuo == 1)[0]
        idx_ausentes = np.where(individuo == 0)[0]
        centro = np.mean(posiciones[idx_presentes], axis=0)
        dist_presentes = np.linalg.norm(posiciones[idx_presentes] - centro, axis=1)
        idx_peor = idx_presentes[np.argmax(dist_presentes)]
        dist_ausentes = np.linalg.norm(posiciones[idx_ausentes] - centro, axis=1)
        idx_mejor = idx_ausentes[np.argmin(dist_ausentes)]
        individuo[idx_peor] = 0
        individuo[idx_mejor] = 1
    return individuo

def seleccionar_padres_torneo_fuerte(poblacion, fitness, k=3):
    topk = int(len(poblacion) * 0.2)
    mejores_idx = (np.argsort(fitness))[:topk]
    seleccionados = []
    for _ in range(2):
        candidatos = random.sample(mejores_idx.tolist(), k)
        mejor = candidatos[int(np.argmin(fitness[candidatos]))]
        seleccionados.append(poblacion[mejor])
    return seleccionados

def nueva_generacion(poblacion, D, posiciones, quorum_q, prob_mut=0.1, elitismo=True):
    fitness = EvalCromo_batch(poblacion, D)
    elite = poblacion[np.argmin(fitness)]
    
    nueva_pob = []
    # élite reforzada
    nueva_pob.append(np.copy(elite))
    nueva_pob.append(np.copy(elite))
    
    while len(nueva_pob) < len(poblacion):
        padre1, padre2 = seleccionar_padres_torneo_fuerte(poblacion, fitness)
        hijo = cruzamiento_consenso(padre1, padre2)
        hijo = mutacion_adaptativa(hijo, posiciones, quorum_q, prob_mut=prob_mut)
        if np.sum(hijo) == quorum_q:
            nueva_pob.append(hijo)
    nueva_pob = np.stack(nueva_pob)
    return nueva_pob

def encontrar_Gc(posiciones, quorum_q):
    n = posiciones.shape[0]
    mejor_z = np.inf
    mejor_G = None
    for i in range(n):
        distancias = np.linalg.norm(posiciones - posiciones[i], axis=1)
        vecinos = np.argsort(distancias)[1:quorum_q]
        Gi = np.append(i, vecinos)
        submat = calcular_matriz_distancias(posiciones[Gi])
        z = np.sum(np.triu(submat, 1))
        if z < mejor_z:
            mejor_z = z
            mejor_G = Gi
    return mejor_G, mejor_z

def evolucionar_generaciones(poblacion, D, posiciones, quorum_q, max_gen=10000):
    mejor_fitness = np.inf
    mejor_individuo = None
    
    for gen in range(max_gen):
        poblacion = nueva_generacion(poblacion, D, posiciones, quorum_q, prob_mut=0.15, elitismo=True)
        fitness = EvalCromo_batch(poblacion, D)
        f_min = np.min(fitness)
        if f_min < mejor_fitness:
            mejor_fitness = f_min
            idx = np.argmin(fitness)
            mejor_individuo = poblacion[idx]
        
        if gen % 1 == 0:
            print(f"Gen {gen}: Z = {float(mejor_fitness):.2f}")
        
    return mejor_individuo, mejor_fitness

if __name__ == "__main__":
    start_time = time.time()
    quorum_q = 216
    pob_max = 30
    
    df = pd.read_csv("H094_members.csv")
    df = df[(df["congress"] == 94) & (df["chamber"] == "House")]
    columnas_utiles = ['nominate_dim1', 'nominate_dim2']
    df_numerico = df[columnas_utiles].dropna()
    datos_np = np.asarray(df_numerico.values.astype(np.float32))
    n_diputados = datos_np.shape[0]
    
    print("\n🔧 Calculando matriz de distancias...")
    D = calcular_matriz_distancias(datos_np)
    
    print("🔎 Buscando heurística Gc...")
    mejor_G, mejor_z = encontrar_Gc(datos_np, quorum_q)
    print("Mejor subconjunto Gc:", mejor_G)
    print("Costo total z:", float(mejor_z))
    
    print("🧬 Generando población inicial...")
    poblacion = generar_poblacion_inicial(pob_max, n_diputados, quorum_q)
    
    # inyectar mejor_G
    individuo_gc = np.zeros(n_diputados, dtype=np.int32)
    individuo_gc[mejor_G] = 1
    poblacion = np.vstack([poblacion, individuo_gc])
    
    # variantes mutadas
    for _ in range(3):
        mutado = mutacion_adaptativa(individuo_gc.copy(), datos_np, quorum_q, prob_mut=0.25)
        if np.sum(mutado) == quorum_q:
            poblacion = np.vstack([poblacion, mutado])
    if poblacion.shape[0] > pob_max:
        poblacion = poblacion[:pob_max]
    
    print("📊 Fitness inicial:", (EvalCromo_batch(poblacion, D)))
    
    mejor_individuo, mejor_fitness = evolucionar_generaciones(
        poblacion, D, datos_np, quorum_q, max_gen=10000
    )
    indices_mejor = np.where(mejor_individuo == 1)[0]
    
    print("🎯 Mejor individuo final:", (indices_mejor))
    print("🎯 Mejor fitness:", float(mejor_fitness))
    print(f"⏱ Tiempo total: {time.time() - start_time:.2f} segundos")

        # Añadir columnas necesarias al DataFrame
    df_filtrado = df[columnas_utiles + ["party_code"]].dropna().copy()
    df_filtrado["Pertenece"] = mejor_individuo.astype(bool)

    # Mapeo colores y formas
    colores_partido = {100: "blue", 200: "red"}
    formas_cgm = {True: "o", False: "x"}

    # Crear figura
    plt.figure(figsize=(10, 8))
    plt.style.use("ggplot")

    # Graficar cada punto según partido y pertenencia
    for partido in df_filtrado["party_code"].unique():
        for pertenece in [True, False]:
            subset = df_filtrado[(df_filtrado["party_code"] == partido) & (df_filtrado["Pertenece"] == pertenece)]
            plt.scatter(
                subset["nominate_dim1"],
                subset["nominate_dim2"],
                c=colores_partido[partido],
                label=f"{partido} - {'Pertenece' if pertenece else 'No pertenece'}",
                alpha=0.7,
                marker=formas_cgm[pertenece],
                edgecolor="black" if pertenece else "none"
            )

    # Envolvente convexa (CGM)
    puntos_cgm = datos_np[mejor_individuo.astype(bool)]
    if len(puntos_cgm) >= 3:
        hull = ConvexHull(puntos_cgm)
        for simplex in hull.simplices:
            plt.plot(puntos_cgm[simplex, 0], puntos_cgm[simplex, 1], 'purple')
        plt.fill(puntos_cgm[hull.vertices, 0], puntos_cgm[hull.vertices, 1], 'purple', alpha=0.2)

    # Leyendas y etiquetas
    plt.xlabel("Dimensión 1")
    plt.ylabel("Dimensión 2")
    plt.title("Visualización de CGM y Partidos")
    plt.legend()
    plt.grid(True)
    plt.show()