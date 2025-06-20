import numpy as np

def generar_poblacion_inicial(pob_max, N, pesos, quorum_min):
    poblacion = []

    while len(poblacion) < pob_max:
        individuo = np.random.randint(0, 2, N)  # genera vector binario aleatorio
        suma_votos = np.sum(individuo * pesos)

        if suma_votos >= quorum_min:
            poblacion.append(individuo)

    return np.array(poblacion)
