#Programa para el calculo de LMTD
#En posteriores sesiones se utilizara
#Algoritmo para calcular numero optimo de corazas
from math import *
class FIC(object):
	#Funcion para calculo de parametros P y R
	def P_R(T1,T2,t1,t2):
		#Se calcula diferencia mayor de temp
		dt1=T1-t2
		#Se calcula diferencia menor de temp
		dt2=T2-t1
		#Se calcula P
		P=(t2-t1)/(T1-t1)
		#Se calcula R
		R=(T1-T2)/(t2-t1)
		return dt1,dt2,P,R
	def LMTD(T1,T2,t1,t2):
		#Se reciben valores de P y R
		CTES=FIC.P_R(T1,T2,t1,t2)
		dt1=CTES[0]
		dt2=CTES[1]
		P=CTES[2]
		R=CTES[3]
		#Ciclo while
		A=-1
		N=3 #inicializarlo
		B=A
		F=-1
		#Inicia el ciclo while
		while ((A<0) or (B<0) & (F<0.75)):
			#Inicia estructura condicional
			if (R!=1):
				P_p=(1-((P*R-1)/(P-1))**(1/N))/(R-((P*R-1)/(P-1))**(1/N))
				A=((2/P_p)-1-R+sqrt(R**2+1))/((2/P_p)-1-R-sqrt(R**2+1))
				#Factor de correccion
				if (A<0): #Estamos garantizando las condiciones de arriba 
					F=0
				else:
					F=((sqrt(R**2+1)/(R-1))*log10((1-P_p)/(1-R*P_p)))/log10(A)
			else:
				P_bp=P/(N-P*(N-1))
				B=((1/P_bp)-2+sqrt(2))/((2/P_bp)-2-sqrt(2))
				#Factor de correccion
				if (B<0):
					F=0  #Estamos garantizando las condiciones de arriba 
				else:
					F=((sqrt(R**2+1)*P_bp)/(log(10*(1-P_bp))))/log10(B)
		N=N+1 #Que vaya aumentando el numero de corazas para que siga optimizando			
		#Se va a detener cuando las condiciones dejen de ser ciertas
		#Fin de la estructura condicional
		dLMTD=F*((dt1-dt2)/log(dt1/dt2))
		return dLMTD,N

#print(FIC.LMTD(250,100,80,120))