#Método de Muller
#----------------------
import numpy as np
from math import *
import sympy as sp
import cmath #modulo que permite trabajar con complejos de una mejor manera que math

def Coef(Pol,X_0,X_1,X_2):
	'''Pol: Lista con los coeficientes del polinomio
	X_0,X_1,X_2: Valores iniciales del proceso iterativo '''
	h_0=X_1-X_0
	h_1=X_2-X_1
	#Se evalua el polinomio
	P=np.poly1d(Pol)   
	S_0=(P(X_1)-P(X_0))/(X_1-X_0)
	S_1=(P(X_2)-P(X_1))/(X_2-X_1)
	#Se calculan los coeficientes
	a=(S_1-S_0)/(h_1+h_0)
	b=a*h_1+S_1
	c=P(X_2)
	return a,b,c,P

def Iter_Muller(Pol,X_0,X_1,X_2):
	#Se calculan las constantes a,b y c
	Cons=Coef(Pol,X_0,X_1,X_2)
	while (abs((X_2-X_1)/X_2)>=0.0001):
		a=Cons[0]
		b=Cons[1]
		c=Cons[2]
		#Condición para tomar la raiz positiva o negativa
		if (b>0):
			X_3=X_2+((-2*c)/(b+cmath.sqrt(b**2-4*a*c)))
		else:
			X_3=X_2+((-2*c)/(b-cmath.sqrt(b**2-4*a*c)))
		#Se actualizan los valores para una nueva iteración
		X_0=X_1
		X_1=X_2
		X_2=X_3
		#Se calculan nuevas constantes a,b y c
		Cons=Coef(Pol,X_0,X_1,X_2)	
	return X_3

#x^3-13*x-12
print(Iter_Muller([1,0,-13,-12],4.5,5.5,5))



