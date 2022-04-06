#German Hernández
#Programa para simular condiciones de operacion de destilacion
#Flash isotermico utilizando el procedimiento de Rachford-Rice
import numpy as np 
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
from FUG import TBD
#---------------------------------------------------------------------------
#Se define clase que contendra las funciones para calcular condiciones de operacion de flash isotermico

class FLASH_RR(object):
	#Se procede a definir funcion para calcular parametro FI
	def FI(Componentes,Composiciones,P):
		"""Componentes: lista con nombres de componente de acuerdo a los que se tienen en el prog FUG
		[Componente 1, Componente 2....Componente n]
		Composiciones: lista de composiciones en la alimentacion del proceso, [Z1,Z2,Z3,Z4...Zn]
		P: Presion de operacion dada en atm
		T: Temperatura del sistema dada en K """
		#Se procede a calcular parametros Ki para la presion indicada
		#Se asumira condicion de equilibrio para el calculo de estos parametros
		Datos=TBD.TB(Componentes,Composiciones,[P])
		ki=np.asarray(Datos[1][0])
		T=Datos[0][0]
		#Se define funcion anonima para resolver funcion de fi
		#Transformamos listas a vectores de numpy
		Zi=np.asarray(Composiciones)
		f_FI=lambda PHI: np.sum((Zi*(1-ki))/(1+PHI*(ki-1)))
		#Se encuentra el valor de phi
		PHI_v=fsolve(f_FI,0) #El valor de phi es pequeño entonces se recomienda inciar a con 0
		return PHI_v,T,ki

	#Funcion para encontrar las composiciones de liquido y vapor
	def XI_YI(Componentes,Zi,P):
		#Se calcula parametros PHI
		DATA=FLASH_RR.FI(Componentes,Zi,P)
		#Se calculan composiciones de liq y vap
		PHI=DATA[0]
		Ki=DATA[2]
		Xi=np.asarray(Zi)/(1+PHI*(Ki-1))
		Yi=Ki*Xi
		T=DATA[1]
		return Xi,Yi,T
	#Se procede a definir la funcion para calcular flash a distintas P
	def FLASH_C(Componentes,Zi,Pmin,Pmax,n=10):
		"""Pmin: Presion minima de operacion
		Pmax: Presion maxima de operacion
		n: numero de puntos en los que se dividira el intervalo de presiones
		Este algoritmo se puede utilizar para encontrar una presion optima de operacion para el sistema flash isotermico"""
		P=np.linspace(Pmin,Pmax,n)
		#Se inicia ciclo for para guardar los resultados
		v_Xi=[]
		v_Yi=[]
		v_Ti=[]
		for k in list(range(0,len(P))):
			Datos=FLASH_RR.XI_YI(Componentes,Zi,P[k])
			v_Xi.append(list(Datos[0]))
			v_Yi.append(list(Datos[1]))
			v_Ti.append(Datos[2])
		#Finaliza el ciclo for
		#Se muestran los resultados	 
		Composiciones_liquido=np.transpose(np.asarray(v_Xi)) #Se ordenan asi para poder graficarlos 
		Composiciones_vapor=np.transpose(np.asarray(v_Yi))   #Se tiene inicialmente xi=[x1,x2,x3,x4],[x1,x2,x3,x4],[x1,x2,x3,x4],[x1,x2,x3,x4] y eso se convierte en xi=[x1,x1,x1],[x2,x2,x2],[x3,x3,x3],[x4,x4,x4] lo cual facilita seleccionar por elementos
		return Composiciones_liquido,Composiciones_vapor,v_Ti
	#Finalmente funcion para mostrar  graficos
	def FLASH_G(Componentes,Zi,Pmin,Pmax,n=10):
		Datos=FLASH_RR.FLASH_C(Componentes,Zi,Pmin,Pmax,n)
		#Grafico de composiciones vs temperatura
		#Debemos acceder a cada componente de la matriz general 
		for k in list(range(0,len(Datos[0]))): 
			plt.figure(0)
			plt.plot(Datos[2],Datos[0][k],"--x",label="Xi-"+Componentes[k])
			plt.plot(Datos[2],Datos[1][k],"--1",label="Yi-"+Componentes[k])
		#lo que no tiene el iterador queda fuera del ciclo
		plt.legend(loc="best")	
		plt.ylabel("Xi-Yi")
		plt.xlabel("Temperatura (K)")
		plt.title("Composiciones Xi-Yi\n de flash isotermico")
		plt.grid(True)
		#Grafico de Composiciones vs Presion
		for k in list(range(0,len(Datos[0]))): 
			plt.figure(1)
			plt.plot(np.linspace(Pmin,Pmax,n)
				       ,Datos[0][k],"--x",label="Xi-"+Componentes[k])
			plt.plot(np.linspace(Pmin,Pmax,n)
				       ,Datos[1][k],"--1",label="Yi-"+Componentes[k])
		plt.legend(loc="best")	
		plt.ylabel("Xi-Yi")
		plt.xlabel("Presion (atm)")
		plt.title("Composiciones Xi-Yi\n de flash isotermico")
		plt.grid(True)
		return plt.show()

#from Flash_iso import FLASH_RR
#FLASH_RR.FLASH_C(["Propano","n-Butano","Pentano","Hexano"],[0.1,0.2,0.3,0.4],0.1,10,6)

FLASH_RR.FLASH_G(["Propano","n-Butano","Pentano","Hexano"],[0.1,0.2,0.3,0.4],0.1,10,100)







