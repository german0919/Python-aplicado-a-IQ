#German Hernández 
#Método de Newton-Raphson
print("-----------------------------------------------------------")
print("         Este programa resulve ecuaciones no lineales ")
print("         usando el Método de Newton-Raphson, usando x ")
print("         como variable, use ** para representar exponente ^")
print("-----------------------------------------------------------")
#Se importan bibliotecas
from math import *
from sympy import * #esta vez importamos todas porque no se usan funciones en especifico
#Se define variable simbolica
x=symbols("x")
y=str(input("Introduzca la ecuacion a resolver: "))
#Se deriva la funcion
dy=str(diff(y))
#Se linealiza funcion similar a inline de Matlab
fy=lambda x: float(eval(y))
fdy=lambda x: float(eval(dy))
#El valor de x se define como float ya que sino se registra como str y no se evalua en la funcion solo se sustituye
x=float(input("Indique el valor en el que se iniciara el calculo: "))
tol=0.00001 #Se indica la tolerancia
#Se evaluan las funciones
fy0=fy(x)
dfy0=fdy(x)
er=tol+1 #Se inicia el error
n=0 #Este es el contador de iteraciones
#Inicia el ciclo
while (er>=tol):
	x0=x-(fy0/dfy0) #Formula de NR
	fy0=fy(x0)
	dfy0=fdy(x0)
	er=abs(x0-x)
	x=x0
	n=n+1
	print(er)
#Fin del ciclo	
print("La ecuacion se resolvia a las %u itereaciones" %n)
print("La respuesta es %g" %x)

