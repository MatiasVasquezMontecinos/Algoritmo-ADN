# Algoritmo Genetico

## Generar_poblacion_inicial

-  genera una cantidad de cromosomas en este caso 30

-  Cada cromosoma es un vector binario de longitud especifica (en este caso los 435 diputados)

-  Da una cantidad de espeficas de posiciones con un numero 1 (desde 0 a 216)

-  A los que sobran se marcan como no elegidos


## EvalCromo_batch

-  Calcula el fitness de cada cromosoma: siendo la suma de las distancias entre los diputados elegidos

-  Usa la matriz de distancia

- Solo cuenta la parte superior del triangulo de la matriz (np.triu) para no duplicar distancias


## Cruzamento_consenso

-  Combina dos cromosomas (Los padres)

-  Donde los padres tiene el mismo bit, copia el mismo bit

-  Donde difieren, elige aletoreamente entrse el gen del padre1 o del padre2


## Mutacion_adaptiva

-  Con probabilidad, reemplaza un diputado muy alejado del centro del grupo por uno extreno mas cercano

-  Calcula
	
	-  el centro de massas del grupo elegido

	-  identifica el diputado dentro del grupo que esta mas lejos del centro

	-  identifica el diputado fuera del grupo mas cerca del centro
	
	-  y los intercambia
(refuerza cohesion)

## Seleccionar_padre_torneo_fuerte

-  selecciona padre de la elite (20% mejora cromosomas)

-  Entre esos mejores, hece un torneo de tamaño especficico, y toma el mejor

-  Esto aumenta la presión selectiva para mejorar


## Nueva_generacion

-  Calcula el fittnes del la población incial

-  salva al mejor (elitismo) dos veces

-  completa la nueva generación cruzando y mutando padres seleccionados

-  garantiza que todos los hijos tengan exactamente el quorum de diputados


## Evolucionar_generaciones

- Para cada diputado

	-  lo considera como "centro"
	
	-  toma sus vecinos mas cercanos 
	
	-  calcula la suma de distancias entre ellos

-  se queda con el grupo con mejor suma y lo devuelve como heristica incial  
