
import numpy as np
#Regresión polinomial para una ecuación de segundo grado
def RP_2(x,y):
	# x : vector de la variable independiente
	# y : vector de la variable dependiente
	n=len(x)
	x=np.asarray(x)
	y=np.asarray(y)
	#Sumatorias
	s_xi=np.sum(x)
	s_xi2=np.sum(x**2)
	s_xi3=np.sum(x**3)
	s_xi4=np.sum(x**4)
	s_yi=np.sum(y)
	s_xiyi=np.sum(x*y)
	s_xi2yi=np.sum(x**2*y)
	#Se define sistema de ec a resolver para encontrar los coef ai
	X=[[n,s_xi,s_xi2],
	  [s_xi,s_xi2,s_xi3],
	  [s_xi2,s_xi3,s_xi4]] #Será la matriz de los terminos a la izq de la igualdad
	Y=[[s_yi],
	  [s_xiyi],
	  [s_xi2yi]] 	  #Será la matriz de los terminos de la derecha de la igualdad
	#y=x*c-> y*1/c=x por tanto se debe hacer la inversa y multiplicar por y
	#Se encuentran valores de constastes ai
	Inv=np.linalg.inv(X)
	#Vector [ao,a1,a2]
	v_Ca=np.dot(Inv,Y) #Como son de distintas dimensiones se debe multiplicar con la función dot
	#Se crea una función para evaluar
	#Yi=lambda x: v_Ca[0]+v_Ca[1]*x+v_Ca[2]*x**2
	Yi=v_Ca[0]+v_Ca[1]*x+v_Ca[2]*x**2 #evaluar todos los valores ingresados
	return Yi #Yi(xi)

x=[0,1,2,3,4,5]	
y=[2.1,7.7,13.6,27.2,40.9,61.1]

print(RP_2(x,y))



