import pandas as pd
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
    # genera una matriz binaria (pob_max x n_diputados) para representar la población inicial de cromosomas
    poblacion = np.zeros((pob_max, n_diputados), dtype=np.int32)
    for i in range(pob_max):
        # selecciona aleatoriamente quorum_q diputados para activarlos (ponerlos en 1)
        indices = random.sample(range(n_diputados), quorum_q)
        # activa (pone en 1) las posiciones correspondientes en el cromosoma
        poblacion[i, indices] = 1
    return poblacion

def EvalCromo_batch(poblacion, D):
    # tamaño de la población
    N = poblacion.shape[0]
    # vector que almacena el fitness de cada cromosoma
    Z = np.zeros(N, dtype=np.float32)
    for i in range(N):
        # toma las posiciones de los diputados "activos"
        idx = np.where(poblacion[i] == 1)[0]
        if len(idx) < 2:
            continue
        # crea la submatriz de distancias entre los diputados activos
        submat = D[idx[:, None], idx]
        # toma solamente la parte superior triangular sin contar la diagonal ni duplicar las distancias
        Z[i] = np.sum(np.triu(submat, k=1))
    # Guada el fitness
    return Z

def cruzamiento_consenso(padre1, padre2):
    # Compara elemento a elemento entre los padres seleccionados, devolviendo un bector booleano
    mask_iguales = padre1 == padre2
    # Se crea una copia del primer padre para usuarla como punto de partida del nuevo hijo.
    hijo = np.copy(padre1)
    # verfica el punto en donde los padres difieren, para empezar el procesos de seleccion aleatorio
    dif_pos = np.where(mask_iguales == 0)[0]
    # recorre las pocsisciones en donde difieren para decidir en cada una de esas posiciones que gen hereda el hijo
    for pos in dif_pos:
        """
            Genera un nuemro aleatorio que puede ser 0 o 1, tiene un 50% de que hijo herede el gene de padre1
            si no, tiene el gen del padre2 
        """
        hijo[pos] = padre1[pos] if random.random() < 0.5 else padre2[pos]
    return hijo

def mutacion_adaptativa(individuo, posiciones, quorum_q, prob_mut=0.1):
    # con probabildad se activa la mutacion, sino el individou se va tal cual
    if random.random() < prob_mut:
        # busca lo indices de la individuos "activo" o "no activos" 
        idx_presentes = np.where(individuo == 1)[0]
        idx_ausentes = np.where(individuo == 0)[0]
        # Cañcula el centroide de los individuo activos.
        centro = np.mean(posiciones[idx_presentes], axis=0)
        # Calcula la distancia de cada individuo activo
        dist_presentes = np.linalg.norm(posiciones[idx_presentes] - centro, axis=1)
        # Buscamos el individuo mas alejado de coalicion, para poder reemplazarlo 
        idx_peor = idx_presentes[np.argmax(dist_presentes)]
        # Calcula la distancia entre cada individuo inactivo respecto al centroide, para ver quien puede mejorar la cohecion
        dist_ausentes = np.linalg.norm(posiciones[idx_ausentes] - centro, axis=1)
        # Busca loe mejores de o los mas cercanos al centroide, o el integracion para mejorar la cohesion
        idx_mejor = idx_ausentes[np.argmin(dist_ausentes)]
        # Se remueve o activa con respecto a los resultados del mejor o peor
        individuo[idx_peor] = 0
        individuo[idx_mejor] = 1
        # Se devuelve el mejor
    return individuo

def seleccionar_padres_torneo_fuerte(poblacion, fitness, k=3):
    # Se toman solo el 20% mejor del grupo de elite
    topk = int(len(poblacion) * 0.2)
    # Se ordena los indices de la poblacion por fitness de menor a mayor y se queda con los indices de los mejores individuos
    mejores_idx = (np.argsort(fitness))[:topk]
    # lista en donde se guardaran los dos padres seleccionados
    seleccionados = []
    for _ in range(2):
        # Se elegien aleatoreamente individuos dentro de los mejores indiceys para formar el torneo, siendo posibles futuros padres
        candidatos = random.sample(mejores_idx.tolist(), k)
        # Como se minimiza se coje segun el bajo desempeño y se seleciona como padres
        mejor = candidatos[int(np.argmin(fitness[candidatos]))]
        # SE agrega ese nuevo padre a la lista de candidatos
        seleccionados.append(poblacion[mejor])
    return seleccionados

def nueva_generacion(poblacion, D, posiciones, quorum_q, prob_mut=0.1, elitismo=True):
    #  Calcula el fitness de toda la población actual, usando la matriz de distancias.
    fitness = EvalCromo_batch(poblacion, D)
    # Identifica el mejor individuo (el de menor fitness) de la población actual, para preservarlo como élite.
    elite = poblacion[np.argmin(fitness)]
    
    #  Inicializa la lista para almacenar la nueva generación.
    nueva_pob = []
    # Introduce dos copias del mejor individuo al comienzo de la nueva población, reforzando la élite.
    # Esto aumenta la probabilidad de mantener la calidad del mejor cromosoma.
    nueva_pob.append(np.copy(elite))
    nueva_pob.append(np.copy(elite))
    

    # Repite el proceso hasta que la nueva generación alcance el mismo tamaño que la población original.
    while len(nueva_pob) < len(poblacion):
        # Selecciona dos padres mediante torneo fuerte, priorizando los mejores
        padre1, padre2 = seleccionar_padres_torneo_fuerte(poblacion, fitness)
        #  Genera un hijo haciendo cruzamiento (crossover) entre esos dos padres.
        hijo = cruzamiento_consenso(padre1, padre2)
        # Aplica una mutación adaptativa al hijo, con probabilidad prob_mut, para mejorar su cohesión.
        hijo = mutacion_adaptativa(hijo, posiciones, quorum_q, prob_mut=prob_mut)
        #  Comprueba que el hijo resultante sigue respetando la cantidad de diputados activos exactamente igual al quorum_q.(si es o no factible)
        if np.sum(hijo) == quorum_q:
            # Si pasa el control, agrega el hijo a la nueva población.
            nueva_pob.append(hijo)
    # Convierte la lista nueva_pob a un array de numpy, para mantener la misma estructura que la población original.
    nueva_pob = np.stack(nueva_pob)
    # Devuelve la nueva generación lista para la próxima iteración.
    return nueva_pob

def encontrar_Gc(posiciones, quorum_q):
    # guarda el número total de diputados, es decir el número de filas de la matriz posiciones.
    n = posiciones.shape[0]
    # Inicializa mejor_z con infinito, de modo que cualquier solución encontrada será mejor que esto en la primera iteración.
    mejor_z = np.inf
    # Inicializa mejor_G como None. Aquí se guardará el conjunto óptimo provisional que vaya encontrando.
    mejor_G = None
    # Empieza un bucle que recorre cada diputado i como posible centro del grupo inicial.
    for i in range(n):
        #  Calcula un vector de distancias euclidianas entre el diputado i y todos los demás diputados.
        distancias = np.linalg.norm(posiciones - posiciones[i], axis=1)
        # Ordena las distancias de menor a mayor con argsort, exceptuando al primero ([1:quorum_q]) porque el primero es el propio i.
        vecinos = np.argsort(distancias)[1:quorum_q]
        # Construye el grupo Gi formado por el diputado central i más sus vecinos más cercanos.
        # En total quorum_q diputados.
        Gi = np.append(i, vecinos)
        #  Calcula la submatriz de distancias interna de este grupo Gi
        # Es decir, la matriz de distancias completa solo para esos diputados seleccionados.
        submat = calcular_matriz_distancias(posiciones[Gi])
        #  Suma la parte triangular superior de la submatriz de distancias, quitando la diagonal y sin duplicar las distancias (por simetría).
        # Esto equivale a calcular la suma de todas las distancias entre pares en el grupo.
        z = np.sum(np.triu(submat, 1))
        # Si la suma de distancias del grupo actual z es menor que la mejor conocida hasta ahora mejor_z, la actualiza.
        if z < mejor_z:
            #  Guarda este Gi como el mejor grupo encontrado hasta ahora, y su coste total z como el mejor coste.
            mejor_z = z
            mejor_G = Gi
    return mejor_G, mejor_z

def evolucionar_generaciones(poblacion, D, posiciones, quorum_q, max_gen=10000):
    # Se maraca el mejor fitness cono infinito ya que buscamos minimizar
    mejor_fitness = np.inf
    # Y none come le mejor individuo
    mejor_individuo = None
    
    # Comienza a iterar durante max_gen generaciones.
    for gen in range(max_gen):
        # enera la siguiente población usando la función nueva_generacion con: probabilidad de mutación 0.15 y elitismo activado
        poblacion = nueva_generacion(poblacion, D, posiciones, quorum_q, prob_mut=0.1, elitismo=True)
        #  Evalúa el fitness de toda la población actual (menor valor significa mejor cohesión).
        fitness = EvalCromo_batch(poblacion, D)
        # Extrae el mejor fitness de esta generación.
        f_min = np.min(fitness)
        # Si este nuevo mejor fitness es mejor que el mejor encontrado hasta ahora, actualiza el registro:
        if f_min < mejor_fitness:
            # Guarda el nuevo mejor fitness y su cromosoma correspondiente.
            mejor_fitness = f_min
            idx = np.argmin(fitness)
            mejor_individuo = poblacion[idx]
        # Imprime el progreso de cada generación (en realidad % 1 imprime en todas las generaciones, es redundante, podrías poner solo if True:).
        # Muestra la generación actual y el mejor fitness hasta ahora.
        if gen % 1 == 0:
            print(f"Gen {gen}: Z = {float(mejor_fitness):.2f}")
        
    return mejor_individuo, mejor_fitness

if __name__ == "__main__":
    start_time = time.time()
    quorum_q = 216
    pob_max = 30
    max_gen = 10000
    repeticiones = 4

    # === Datos ===
    df = pd.read_csv("H094_members.csv")
    df = df[(df["congress"] == 94) & (df["chamber"] == "House")]
    columnas_utiles = ["nominate_dim1", "nominate_dim2"]
    df_numerico = df[columnas_utiles].dropna()
    datos_np = df_numerico.to_numpy(dtype=np.float32)
    n_diputados = datos_np.shape[0]

    print("\n🔧 Calculando matriz de distancias...")
    D = calcular_matriz_distancias(datos_np)

    print("🔎 Buscando heurística Gc...")
    mejor_G, mejor_z = encontrar_Gc(datos_np, quorum_q)
    print("Mejor heurístico Gc:", mejor_G)
    print("Costo total z:", float(mejor_z))

    # === Población inicial ===
    print("🧬 Generando población inicial...")
    poblacion = generar_poblacion_inicial(pob_max, n_diputados, quorum_q)

    # añade heurística Gc
    individuo_gc = np.zeros(n_diputados, dtype=np.int32)
    individuo_gc[mejor_G] = 1
    poblacion = np.vstack([poblacion, individuo_gc])

    # perturbaciones
    for _ in range(3):
        mutado = mutacion_adaptativa(individuo_gc.copy(), datos_np, quorum_q, prob_mut=0.25)
        if np.sum(mutado) == quorum_q:
            poblacion = np.vstack([poblacion, mutado])
    if poblacion.shape[0] > pob_max:
        poblacion = poblacion[:pob_max]

    print("📊 Fitness inicial:", EvalCromo_batch(poblacion, D))

    # === Evolución principal ===
    mejor_individuo, mejor_fitness = evolucionar_generaciones(
        poblacion, D, datos_np, quorum_q, max_gen=max_gen
    )
    indices_mejor = np.where(mejor_individuo == 1)[0]

    print("🎯 Mejor individuo final:", indices_mejor)
    print("🎯 Mejor fitness:", float(mejor_fitness))
    print(f"⏱ Tiempo total: {time.time() - start_time:.2f} segundos")

    # === Repeticiones estadísticas ===
    print("\n▶ Ejecutando repeticiones estadísticas...")

    fitness_finales = []
    tiempos = []
    iteraciones = []

    for rep in range(repeticiones):
        print(f"\n🔁 Repetición {rep + 1}/{repeticiones}")
        start_rep = time.time()

        poblacion = generar_poblacion_inicial(pob_max, n_diputados, quorum_q)
        individuo_gc = np.zeros(n_diputados, dtype=np.int32)
        individuo_gc[mejor_G] = 1
        poblacion = np.vstack([poblacion, individuo_gc])
        for _ in range(3):
            mutado = mutacion_adaptativa(individuo_gc.copy(), datos_np, quorum_q, prob_mut=0.25)
            if np.sum(mutado) == quorum_q:
                poblacion = np.vstack([poblacion, mutado])
        if poblacion.shape[0] > pob_max:
            poblacion = poblacion[:pob_max]

        mejor_fitness_rep = np.inf
        iteracion_encontrado = 0
        for gen in range(max_gen // 10):  # menos generaciones en repetición
            poblacion = nueva_generacion(poblacion, D, datos_np, quorum_q, prob_mut=0.15)
            fitness = EvalCromo_batch(poblacion, D)
            f_min = np.min(fitness)
            if f_min < mejor_fitness_rep:
                mejor_fitness_rep = f_min
                iteracion_encontrado = gen

        tiempo_total = time.time() - start_rep
        fitness_finales.append(mejor_fitness_rep)
        tiempos.append(tiempo_total)
        iteraciones.append(iteracion_encontrado)

    # === Estadísticas ===
    fitness_finales = np.array(fitness_finales)
    tiempos = np.array(tiempos)
    iteraciones = np.array(iteraciones)

    precision_array = 100 * mejor_z / fitness_finales
    precision_mean = np.mean(precision_array)
    precision_std = np.std(precision_array)

    fitness_mean = np.mean(fitness_finales)
    fitness_std = np.std(fitness_finales)

    iter_mean = np.mean(iteraciones)
    iter_std = np.std(iteraciones)

    tiempo_mean = np.mean(tiempos)
    tiempo_std = np.std(tiempos)

    print("\n=== 📈 Resultados Finales ===")
    print(f"Resultado esperado (Z):")
    print(f"  ➤ Media: {fitness_mean:.5f}")
    print(f"  ➤ Desviación estándar: {fitness_std:.5f}")
    print(f"\nPrecisión (%):")
    print(f"  ➤ Media: {precision_mean:.2f} %")
    print(f"  ➤ Desviación estándar: {precision_std:.2f} %")
    print(f"\nIteraciones hasta el mejor fitness:")
    print(f"  ➤ Media: {iter_mean:.2f}")
    print(f"  ➤ Desviación estándar: {iter_std:.2f}")
    print(f"\nTiempo de ejecución por repetición (segundos):")
    print(f"  ➤ Media: {tiempo_mean:.5f}")
    print(f"  ➤ Desviación estándar: {tiempo_std:.5f}")

    # === Graficar ===
    df_filtrado = df[columnas_utiles + ["party_code"]].dropna().copy()
    df_filtrado["Pertenece"] = mejor_individuo.astype(bool)

    colores_partido = {100: "blue", 200: "red"}
    formas_cgm = {True: "o", False: "x"}

    plt.figure(figsize=(10, 8))
    plt.style.use("ggplot")

    for partido in df_filtrado["party_code"].unique():
        for pertenece in [True, False]:
            subset = df_filtrado[
                (df_filtrado["party_code"] == partido)
                & (df_filtrado["Pertenece"] == pertenece)
            ]
            plt.scatter(
                subset["nominate_dim1"],
                subset["nominate_dim2"],
                c=colores_partido.get(partido, "gray"),
                label=f"{partido} - {'Pertenece' if pertenece else 'No pertenece'}",
                alpha=0.7,
                marker=formas_cgm[pertenece],
                edgecolor="black" if pertenece else "none",
            )

    puntos_cgm = datos_np[mejor_individuo.astype(bool)]
    if len(puntos_cgm) >= 3:
        hull = ConvexHull(puntos_cgm)
        for simplex in hull.simplices:
            plt.plot(puntos_cgm[simplex, 0], puntos_cgm[simplex, 1], "purple")
        plt.fill(
            puntos_cgm[hull.vertices, 0],
            puntos_cgm[hull.vertices, 1],
            "purple",
            alpha=0.2,
        )

    plt.xlabel("Dimensión 1")
    plt.ylabel("Dimensión 2")
    plt.title("Visualización de CGM y Partidos")
    plt.legend()
    plt.grid(True)
    plt.show()