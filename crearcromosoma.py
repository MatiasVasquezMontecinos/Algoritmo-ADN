import numpy as np
import random


## Generar la poblacion inicial de tamaño IMPAR(Tam. poblacion) respetando las restricciones
def generar_poblacion_inicial(pob_max, n_diputados, quorum_q):
    poblacion = []

    while len(poblacion) < pob_max:
        indices = random.sample(range(n_diputados), quorum_q)  # q diputados únicos
        individuo = np.zeros(n_diputados, dtype=int)
        individuo[indices] = 1
        poblacion.append(individuo)

    return np.array(poblacion)

##  Evaluar cada cromosoma de acuerdo a la Funcion objetivo
def EvalCromo(Individuos, fun_obj):
    return 0 # la funcion

## Ordenar los cromosomas de acuerdo a la Funcion objetivo (Fitness)
def SortCromo(Individuos, fun_obj):

    fitness = [(i, EvalCromo(i)) for i in Individuos]
    fitness.sort(key=lambda x: x[1], reverse=True)
    return np.array([i for i, fit in fitness])


