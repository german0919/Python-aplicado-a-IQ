#German Hernández
"""Programa con funciones para determinar condiciones de operacion
de una columna de destilación multicomponente, utilizando ecuaciones
Fenske-Underwood-Gilliland"""
#-----------------------------------------------------------------------
import numpy as np 
from scipy.optimize import fsolve #fsolve resuelve las ecuaciones objetivo
from math import ceil,floor  
#--------------------------------------------------------------------------
#Se define clase con funciones para calculo de temp de burbuja y rocío
#Condiciones ideales
class TBD(object): 
	"""Se define un diccionario con las constantes de Antoine
	y Temperaturas de ebullición de cada componente"""
	Dicc_ANT={ 
	           "Etano":([15.6637,1511.42,-17.16], 184.45),
	           "Propano":([15.726,1872.46,-25.16], 231.05),
	           "Isobutano":([15.5381,2032.73,-33.15], 261.25),
	           "n-Butano":([15.6782,2154.9,-34.42],272.65),
	           "Pentano":([15.8333,2477.07,-39.94],309.15),
	           "Hexano":([15.8366,2697.55,-48.78],341.85)}
	def TB(Componentes,Composiciones,v_P):           
		"""Funcion para calcular la temp de burbuja, junto con otros parametros de interes (ki,Yi..)
		Componentes: [C1,C2,C3...Cn],se indican los nombres del diccionario.
		Composiciones: [X1,X2,X3...Xn], fracciones molares de componentes.
		v_P:[P1,P2,P3,P4...Pi], lista de presiones del sistema dadas en atm"""
		#Se toman datos del diccionario
		#Inicia ciclo for
		v_Ai=[]
		v_Bi=[]
		v_Ci=[]
		v_TBi=[]
		for k in list(range(0,len(Componentes))):
			#Se ordenan datos a utilizar para ello se crearon las listas vacias anteriormente
			Datos_c=TBD.Dicc_ANT[Componentes[k]]
			v_Ai.append(Datos_c[0][0])
			v_Bi.append(Datos_c[0][1])
			v_Ci.append(Datos_c[0][2])
			v_TBi.append(Datos_c[1])
		#Fin de ciclo for
		#Transformamos a array de numpy
		v_Ai=np.asarray(v_Ai)
		v_Bi=np.asarray(v_Bi)
		v_Ci=np.asarray(v_Ci)
		Composiciones=np.asarray(np.abs(Composiciones))
		#Se define Temperatura para iniciar el proceso iterativo
		TB_mix0=np.sum(np.asarray(v_TBi))/len(v_TBi)
		#Se normalizan valores de xi si su suma es diferente de 1
		if(np.sum(Composiciones)!=1):
			Composiciones=Composiciones/np.sum(Composiciones)
		#Se inicia ciclo for para resolver el sistema con cada presión
		v_TBubble=[]	
		v_ki=[]
		v_Yi=[]
		for k in list(range(0,len(v_P))):
			#Se define funcion anonima a resolver
			f_TB=lambda T: np.sum(
				           np.exp(v_Ai-(v_Bi/(T+v_Ci))
				           	)*Composiciones)-v_P[k]*760
			#Se resuelve función 
			T_Bubble=fsolve(f_TB,TB_mix0)
			#Se calcula ki
			ki=np.exp(v_Ai-(v_Bi/(T_Bubble+v_Ci)))/(v_P[k]*760)
			#Se calcula Yi
			Yi=ki*Composiciones
			#Se guardan resultados/ se pone 0 en T_Bubble porque la funcion fsolve me lo arroja [213]
			v_TBubble.append(T_Bubble[0])
			v_ki.append(list(ki))
			v_Yi.append(list(Yi))
		#Fin del ciclo for
		#Se procede a mostrar resultados
		return np.asarray(v_TBubble),np.asarray(v_ki),np.asarray(v_Yi)	
	def TD(Componentes,Composiciones,v_P):	
			"""Funcion para calcular la temp de rocío, junto con otros parametros de interes (ki,Yi..)
			Componentes: [C1,C2,C3...Cn],se indican los nombres del diccionario.
			Composiciones: [X1,X2,X3...Xn], fracciones molares de componentes.
			v_P:[P1,P2,P3,P4...Pi], lista de presiones del sistema dadas en atm"""
			#Se toman datos del diccionario
			#Inicia ciclo for
			v_Ai=[]
			v_Bi=[]
			v_Ci=[]
			v_TBi=[]
			for k in list(range(0,len(Componentes))): 
				#Se ordenan datos a utilizar para ello se crearon las listas vacias anteriormente
				Datos_c=TBD.Dicc_ANT[Componentes[k]]
				v_Ai.append(Datos_c[0][0])
				v_Bi.append(Datos_c[0][1])
				v_Ci.append(Datos_c[0][2])
				v_TBi.append(Datos_c[1])
			#Fin de ciclo for
			#Transformamos a array de numpy
			v_Ai=np.asarray(v_Ai)
			v_Bi=np.asarray(v_Bi)
			v_Ci=np.asarray(v_Ci)
			Composiciones=np.asarray(np.abs(Composiciones))
			#Se define Temperatura para iniciar el proceso iterativo
			TD_mix0=np.min(v_TBi)
			#Se normalizan valores de xi si su suma es diferente de 1
			if(np.sum(Composiciones)!=1):
				Composiciones=Composiciones/np.sum(Composiciones)
			#Se inicia ciclo for para resolver el sistema con cada presión
			v_TDew=[]	
			v_ki=[]
			v_Xi=[]
			for k in list(range(0,len(v_P))):  
				#Se define funcion anonima a resolver
				f_TD=lambda T: (
				 np.sum(Composiciones/np.exp(v_Ai-(v_Bi/(T+v_Ci))))
				                                         )*v_P[k]*760-1
				#Se resuelve función 
				T_Dew=fsolve(f_TD,TD_mix0)
				#Se calcula ki
				ki=np.exp(v_Ai-(v_Bi/(T_Dew+v_Ci)))/(v_P[k]*760)
				#Se calcula Xi
				Xi=Composiciones/ki
				#Se guardan resultados/ se pone 0 en T_Bubble porque la funcion fsolve me lo arroja [213]
				v_TDew.append(T_Dew[0])
				v_ki.append(list(ki))
				v_Xi.append(list(Xi))
			#Fin del ciclo for
			#Se procede a mostrar resultados
			return np.asarray(v_TDew),np.asarray(v_ki),np.asarray(v_Xi)	

#Clase con funciones de Feske-Underwood-Gilliland
class F_FUG(object):

	def Fenske(Componentes,XiD,XiB,LK,HK,P):
		"""Funcion para calculo de etapas minimas usando ecuaciones de Fenske
		Componentes: [C1,C2,C3...Cn], se indiacan conforme el diccionario
		XiD: Composiciones en el domo
		XiB: Composiciones en el fondo
		LK: Entero que se indica componente clave ligero, segun su posicion en la lista
		de Componentes, recordar que el conteo empieza en 0.
		HK: Entero que se indica componente clave pesado, segun su posicion en la lista
		de Componentes, recordar que el conteo empieza en 0.
		P: Presión de operacion de la columna dada en atm """

		#Se procede a calcular la volatilidad relativa con valores de ki
		#de componentes clave ligero y pesado en domo y fondo
		kiD=TBD.TD(Componentes,XiD,[P])[1][0] #P es solo un valor pero se pone como lista dado a las funciones TB y TD
		kiB=TBD.TB(Componentes,XiB,[P])[1][0]
		alfa_TD=kiD[LK]/kiD[HK]
		alfa_TB=kiB[LK]/kiB[HK]
		#Se calcula numero minimo de etapas con ecuacion de Fenske
		#ceil es util para que redondee el numero de etapas ya que por lo general no da enteros(redondeo hacia arriba) 0.89 es 0.9
		Nmin=ceil(np.log(
			  (XiD[LK]/XiD[HK])*(
			  	XiB[HK]/XiB[LK]))/np.log(np.sqrt(alfa_TD*alfa_TB)))

		return Nmin

	def Underwood(Componentes,XiD,XiF,LK,HK,P,q=1): 
		"""Funcion para calculo de reflujo minimo usando la ec de Underwood
		Componentes: [C1,C2,C3...Cn], se indiacan conforme el diccionario
		XiD: Composiciones en el domo
		XiB: Composiciones en el fondo
		LK: Entero que se indica componente clave ligero, segun su posicion en la lista
		de Componentes, recordar que el conteo empieza en 0.
		HK: Entero que se indica componente clave pesado, segun su posicion en la lista
		de Componentes, recordar que el conteo empieza en 0.
		P: Presión de operacion de la columna dada en atm 
		q: Condicion de alimentacion"""
		#Se procede a calcular la volatilidad relativa respecto a componente HK con valores de ki
		#de componentes en etapa de alimentacion y domo
		kiD=TBD.TB(Componentes,XiD,[P])[1][0]
		kiF=TBD.TB(Componentes,XiF,[P])[1][0]
		alfa_TD=kiD/kiD[HK]
		alfa_TF=kiF/kiF[HK]
		#Se define funcion para calcular tetha
		f_Tetha=lambda Tetha: 1-q-np.sum((alfa_TF*XiF)/(alfa_TF-Tetha))
		#Se resuelve funcion 
		Tetha0=1.5 #valor inicial
		Tetha=fsolve(f_Tetha,Tetha0)[0] #esta funcion retorna una lista por tanto le pido que me de el valor tomandolo[0]
		#Se calcula reflujo minimo
		Rmin=np.sum((alfa_TD*XiD)/(alfa_TD-Tetha))-1
		return Rmin

	def Gilliland(Componentes,XiD,XiF,XiB,LK,HK,P,B,D,q=1,F_Rmin=1.5): 	     
		"""Funcion para calcular etapas teoricas, se utiliza ecuacion de Gilliland.
		Mismos argumentos que la funcion de Underwood
		F_Rmin: constante que multiplica a Rmin(se puede modificar aqui se deja como 1.5)
		B: flujo molar en fondo
		D: flujo molar en domo """
		#Se calcula Rmin
		Rmin=F_FUG.Underwood(Componentes,XiD,XiF,LK,HK,P,q) #no es necesario poner q=1 en esto porque ya se llamo al inicio en Underwood, cuando vamos a llamar funciones en otra solo se escriben sus parametros y no los predefinidos
		#Se calcula parametro X
		R=F_Rmin*Rmin 
		X=(R-Rmin)/(R+1)
		"""Se utilizara la correlacion de McCromick para encontrar valor de Y.
		McCromick, J. E., A correlation for destillation stages and reflux, Chemical Engineering,
		95, 75, September 26, 1988."""
		B=0.105*np.log(X)+0.44
		Y=1-X**B
		#Se calcula Nmin
		Nmin=F_FUG.Fenske(Componentes,XiD,XiB,LK,HK,P)
		#Se calculan etapas teoricas
		N=ceil((-Nmin-Y)/(Y-1))
		#Relacion de etapas por encima y por debajo de alimentacion
		Ratio=np.exp(0.206*np.log((B/D)*(XiF[HK]/XiF[LK])*(XiB[LK]/XiD[HK])**2))
		#Se calculan etapas bajo la alimentacion
		Bajo_F=N/(1+Ratio)
		#Se calculan etapas sobre la alimentacion
		Sobre_F=Bajo_F*Ratio
		return Nmin,Rmin,N,floor(Bajo_F),ceil(Sobre_F) #floor redondea hacia abajo 7.9=7

	def M_FUG(Componentes,XiD,XiF,XiB,LK,HK,P,B,D,q=1,F_Rmin=1.5):
		"""Funcion que retorna numero de etapas de columna de destilacion y distintos parametros de operacion"""
		DATA=F_FUG.Gilliland(Componentes,XiD,XiF,XiB,LK,HK,P,B,D,q,F_Rmin)
		#Temperaturas en domo y fondo respectivamente
		T_Domo=TBD.TD(Componentes,XiD,[P])[0][0]
		T_Fondo=TBD.TB(Componentes,XiB,[P])[0][0]
		return T_Domo,T_Fondo,DATA

TBD.TB(["Etano","Propano","Isobutano","n-Butano"],[0.15,0.18,0.18,0.49],[7])

#from FUG import F_FUG
#print(F_FUG.M_FUG(["Etano","Propano","Isobutano","n-Butano"],[0.02,0.96,0.01,0.01],[0.15,0.18,0.18,0.49],[0.02,0.02,0.95,0.01],1,2,5,300,200))









 






