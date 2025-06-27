## Contexto

El objetivo de esta prueba práctica es diseñar e implementar un algoritmo genético
capaz de resolver el problema de la coalición ganadora mínima (MWC), que consiste en
identificar un subconjunto de diputados cuya suma de distancias políticas a pares sea
mínima, sujeto a alcanzar el quórum requerido para la aprobación de una propuesta.

Este modelo se basa en la formulación presentada por Lincolao-Venegas y cols. (2023),
donde se detallan las ecuaciones (1)-(5) y el esquema de evaluación mediante distancias
en el espacio político.

Adicionalmente, el algoritmo deberá identificar todos los miembros de la MWC, así
como los vértices del polígono convexo que la contiene en el espacio político
multidimensional.

## Ecuaciones

- 1. P = {p1, p2, p3, ...,pn} ----> Conjunto de n-congresistas
- 2. C(n,q) = (n q) = n! / (q! * (n - q)!) ----> Convinaciones posibles para formar coaliciones de size q
- 3. G = {g1, g2, g3, ..., gq} ----> Coalicion Candidata
- 4. Z(G) = Sum(i=1 -> q-1) Sum(j= i+1 -> q) d(gi, gj) ------> funcion obj
- 5. G* = arg_G min Z(G)  --------> Resultado optimo (el minimo o mejor coalicion)

## Algoritmo

- [X] Generar la poblacion inicial de tamaño IMPAR(Tam. poblacion) respetando las restricciones

- [] Evaluar cada cromosoma de acuerdo a la Funcion objetivo
- [] Ordenar los cromosomas de acuerdo a la Funcion objetivo (Fitness)
- [] Seleccionar dos "Padres"
- [] Seleccionar el punto de corte
- [] Cruzar los cromosomas
- [] Realizar Mutaciones
- [] Verificar hijos contra restricciones
- [] Mover hijos a una nueva poblacion

## Elementos visules

- [] Crear una tabla con los resultados
- [] Explicaciones
- [] Dibujo de la grafica

## Dudas

- [X] Hay que saber alguna funcion?

- Si,  esta referenciadas en paper del profesor ignacio lincolao

- [X] ¿Se debe usar un .CSV?

- Si, se debe usar el .CSV ya que apartir de ahi se toma los puntos 

- [X] Se debe elaborar un Parser?

- No es necesario ya que se puede trabajar con la funcion de Distancia entre puntos es decir  D(x,y) = sqrt((x2-x1)^2 + (y2-y1)^2) de manera matematica.

- [X] El grafico de punto como debe ser,¿Igual al Paper o diferente?

-  Si, debe ser como el grafico que esta al final del paper

- [X] ¿Como los resultados obtenidos deben basarse en los del paper?

- Franco le sabe

- [] ¿Que explicacion deberiamos llegar?

- 

## Objetivos

- [] Si de aqui al sabado el codigo no es mas rapido o mas optimo, se recomienda cambiar de lenguaje

- [] Ver si se puede usar la libreria de CuPy (Para RTX)



