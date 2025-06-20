import numpy as np
import random

def generar_poblacion_inicial(pob_max, n_diputados, quorum_q):
    poblacion = []

    while len(poblacion) < pob_max:
        indices = random.sample(range(n_diputados), quorum_q)  # q diputados únicos
        individuo = np.zeros(n_diputados, dtype=int)
        individuo[indices] = 1
        poblacion.append(individuo)

    return np.array(poblacion)
