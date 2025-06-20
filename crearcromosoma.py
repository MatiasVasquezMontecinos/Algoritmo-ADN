import numpy as np

def generar_poblacion_inicial(pob_max, n_diputados, quorum_min):
    poblacion = []
    while len(poblacion) < pob_max:
        cromosoma = np.random.randint(0, 2, size=n_diputados)
        if cromosoma.sum() >= quorum_min:
            poblacion.append(cromosoma)

    return np.array(poblacion)
