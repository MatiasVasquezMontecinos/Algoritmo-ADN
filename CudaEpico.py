import pandas as pd
import cupy as cp
import random


def calcular_matriz_distancias(posiciones):
    return cp.linalg.norm(posiciones[:, None] - posiciones[None, :], axis=2)

def generar_poblacion_inicial(pob_max, n_diputados, quorum_q):
    poblacion = cp.zeros((pob_max, n_diputados), dtype=cp.int32)
    for i in range(pob_max):
        indices = random.sample(range(n_diputados), quorum_q)
        poblacion[i, indices] = 1
    return poblacion

def EvalCromo_batch(poblacion, D):
    N = poblacion.shape[0]
    Z = cp.zeros(N, dtype=cp.float32)

    for i in range(N):
        idx = cp.where(poblacion[i] == 1)[0]
        if len(idx) < 2:
            continue
        submat = D[idx[:, None], idx]
        Z[i] = cp.sum(cp.triu(submat, k=1))
    return Z

def SortCromo(poblacion, D):
    fitness = EvalCromo_batch(poblacion, D)
    idx_sorted = cp.argsort(fitness)
    return poblacion[idx_sorted]

def selectDad(poblacion):
    idx = random.sample(range(len(poblacion)), 2)
    return poblacion[idx[0]], poblacion[idx[1]]

def cruzamiento(Dad1, Dad2):
    punto = random.randint(1, len(Dad1) - 2)
    Son1 = cp.concatenate((Dad1[:punto], Dad2[punto:]))
    Son2 = cp.concatenate((Dad2[:punto], Dad1[punto:]))
    return Son1, Son2

def Mutacion(individuo, prob_mut=0.1):
    if random.random() < prob_mut:
        unos = cp.where(individuo == 1)[0]
        ceros = cp.where(individuo == 0)[0]
        if len(unos) > 0 and len(ceros) > 0:
            i = random.choice(unos.tolist())
            j = random.choice(ceros.tolist())
            individuo[i], individuo[j] = individuo[j], individuo[i]
    return individuo

def verificar(individuo, quorum_q):
    return cp.sum(individuo) == quorum_q

def nueva_generacion(poblacion, D, quorum_q, prob_mutacion=0.1):
    nueva_pob = []
    while len(nueva_pob) < len(poblacion):
        padre1, padre2 = selectDad(poblacion)
        hijo1, hijo2 = cruzamiento(padre1, padre2)
        hijo1 = Mutacion(hijo1, prob_mutacion)
        hijo2 = Mutacion(hijo2, prob_mutacion)

        if verificar(hijo1, quorum_q):
            nueva_pob.append(hijo1)
        if len(nueva_pob) < len(poblacion) and verificar(hijo2, quorum_q):
            nueva_pob.append(hijo2)

    nueva_pob = cp.stack(nueva_pob)
    return SortCromo(nueva_pob, D)


def encontrar_Gc(posiciones, quorum_q):
    n = posiciones.shape[0]
    mejor_z = cp.inf
    mejor_G = None

    for i in range(n):
        distancias = cp.linalg.norm(posiciones - posiciones[i], axis=1)
        vecinos = cp.argsort(distancias)[1:quorum_q]  
        Gi = cp.append(i, vecinos)
        submat = calcular_matriz_distancias(posiciones[Gi])
        z = cp.sum(cp.triu(submat, 1))
        if z < mejor_z:
            mejor_z = z
            mejor_G = Gi

    return mejor_G, mejor_z




df = pd.read_csv("Tabla.csv")
datos_cp = cp.asarray(df.values)