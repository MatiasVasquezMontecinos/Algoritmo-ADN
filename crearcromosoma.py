import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull



def calcular_matriz_distancias(posiciones):
    return np.linalg.norm(posiciones[:, None] - posiciones[None, :], axis=2)


def generar_poblacion_inicial(pob_max, n_diputados, quorum_q):
    poblacion = []
    while len(poblacion) < pob_max:
        indices = random.sample(range(n_diputados), quorum_q)
        individuo = np.zeros(n_diputados, dtype=int)
        individuo[indices] = 1
        poblacion.append(individuo)
    return np.array(poblacion)


def EvalCromo(individuo, D):
    index = np.where(individuo == 1)[0]
    z = 0
    for i in range(len(index) - 1):
        for j in range(i + 1, len(index)):
            z += D[index[i], index[j]]
    return z


def SortCromo(poblacion, D):
    fitness = [(ind, EvalCromo(ind, D)) for ind in poblacion]
    fitness.sort(key=lambda x: x[1])
    return np.array([ind for ind, _ in fitness])


def selectDad(poblacion):
    index = random.sample(range(len(poblacion)), 2)
    return poblacion[index[0]], poblacion[index[1]]


def cruzamiento(Dad1, Dad2):
    punto = random.randint(1, len(Dad1) - 2)
    Son1 = np.concatenate([Dad1[:punto], Dad2[punto:]])
    Son2 = np.concatenate([Dad2[:punto], Dad1[punto:]])
    return Son1, Son2


def Mutacion(individuo, prob_mut=0.1):
    if random.random() < prob_mut:
        unos = np.where(individuo == 1)[0]
        ceros = np.where(individuo == 0)[0]
        if len(unos) > 0 and len(ceros) > 0:
            i = random.choice(unos)
            j = random.choice(ceros)
            individuo[i], individuo[j] = individuo[j], individuo[i]
    return individuo


def verificar(individuo, quorum_q):
    return np.sum(individuo) == quorum_q


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

    return SortCromo(np.array(nueva_pob), D)



# -------- Implementacion del algoritmo determinista Gc del paper -------- #

def encontrar_Gc(posiciones, quorum_q):
    n = posiciones.shape[0]
    mejor_z = float('inf')
    mejor_G = None

    for i in range(n):
        distancias = np.linalg.norm(posiciones - posiciones[i], axis=1)
        vecinos = np.argsort(distancias)[1:quorum_q]  # q-1 más cercanos
        Gi = np.append(i, vecinos)
        submat = calcular_matriz_distancias(posiciones[Gi])
        z = np.sum(np.triu(submat, 1))
        if z < mejor_z:
            mejor_z = z
            mejor_G = Gi
    return mejor_G, mejor_z


