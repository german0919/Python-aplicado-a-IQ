#German Hernández
import numpy as np 
from matplotlib import pyplot as plt
class F_CSTR(object):
	#Se define funcion para obtener constantes K y alfa para nuestra ecuacion de velocidad
	def C_rA(t,CA,n):
		#Se genera un polinmio CA(t)
		F_CA=np.polyfit(t,CA,n)*(-1) #Esta regresion se puede optimizar con una funcion para obtener el mejor polinomio
		#Se deriva el polinomio
		d_FCA_t=np.polyder(F_CA)
		#Se evalua la derivada
		F_dfCA=np.polyval(d_FCA_t,t)
		#Se obtiene logaritmo de la derivada y de la CA
		L_CA_dCA=np.polyfit(np.log(CA),np.log(F_dfCA),1)
		#Se encuentra alfa y k
		ALFA=L_CA_dCA[0]
		K=np.exp(L_CA_dCA[1])
		return ALFA,K
	#Se define una funcion para determinar volumenes de reactor para datos de conversión
	def V_r(X,r_A,FA_0):
		#Se calcula FA_0/r_A
		inv_rA=FA_0/np.asarray(r_A) #Para hacer la grafica de levenspiel donde el area bajo la curva es el vol
		#Se procede a calcular el area bajo la curva(volumen del reactor)
		V1=[np.trapz(X,x=[inv_rA[0],inv_rA[k]]) for k in list(range(0,len(inv_rA)))] 
		#Se debe calcular el volumen sobrante  para restarlo
		V2=inv_rA*X-np.asarray(V1)
		#Se suman los volumenes
		V=np.asarray(V1)+np.asarray(V2)
		return V,X,r_A
	#Se define funcion para calcular volumenes de reactor en funcion de alfa y K
	def V_Ctes(CA_0,FA_0,t,CA,n):	
		#Se calculan alfa y K
		Ctes=F_CSTR.C_rA(t,CA,n)
		#Se calcula conversion y velocidad de reaccion 
		X=1-(np.asarray(CA)**Ctes[0]/CA_0)
		r_A=Ctes[1]*np.asarray(CA)**Ctes[0]
		Datos=F_CSTR.V_r(X,r_A,FA_0)
		#Se calcula volumen
		return Datos[0],Datos[1],Datos[2],t,CA
	#Se define una funcion para generar un grafico dinamico v=1 es la velocidad de la animación
	def G_D(X,Y,ejes,titulo,v=1):
		#Se genera una ventana para redibujar el grafico
		ax=plt.axes()
		#Se define ventana de grafico
		plt.title(titulo)
		plt.xlabel(ejes[0])
		plt.ylabel(ejes[1])
		plt.grid(True)
		#inicia for para redibujar el grafico
		for k in list(range(0,len(X))):
			#Se actualizan datos
			ax.plot(X[k], Y[k], "r>")
			#Se dibuja datos respecto a el tiempo
			plt.draw()
			plt.pause(v)
		#fin del ciclo for
		plt.plot(X,Y,"--r",linewidth=.5)
		return plt.show()	
	#Grafico V vs t
	def GD_Vt(CA_0,FA_0,t,CA,n,velocidad=1):
		Datos=F_CSTR.V_Ctes(CA_0,FA_0,t,CA,n)
		return F_CSTR.G_D(Datos[0],Datos[3],("Volumen","Tiempo"),"Perfil de Volumen de CSTR \n respecto al tiempo",velocidad)
	#Grafico de V vs rA
	def GD_VrA(CA_0,FA_0,t,CA,n,velocidad=1):	
		Datos=F_CSTR.V_Ctes(CA_0,FA_0,t,CA,n)
		return F_CSTR.G_D(Datos[0],Datos[2],("Volumen","-rA"),"Perfil de Volumen de CSTR \n respecto a -rA",velocidad)		
	#Grafico de V vs X
	def GD_VX(CA_0,FA_0,t,CA,n,velocidad=1):			
		Datos=F_CSTR.V_Ctes(CA_0,FA_0,t,CA,n)
		return F_CSTR.G_D(Datos[0],Datos[1],("Volumen","X"),"Perfil de Volumen de CSTR \n respecto a X",velocidad)				
	#Grafico de V vs CA
	def GD_VCA(CA_0,FA_0,t,CA,n,velocidad=1):
		Datos=F_CSTR.V_Ctes(CA_0,FA_0,t,CA,n)
		return F_CSTR.G_D(Datos[0],Datos[4],("Volumen","CA"),"Perfil de Volumen de CSTR \n respecto a CA",velocidad)						


#print(F_CSTR.GD_Vt(.05,.4,[0,1,2,3,4,5,6],[.050,.038,.0306,.0256,.0222,.0195,.0174],4,0.5))

#print(F_CSTR.GD_VrA(.05,.4,[0,1,2,3,4,5,6],[.050,.038,.0306,.0256,.0222,.0195,.0174],4,0.5))

#print(F_CSTR.GD_VX(.05,.4,[0,1,2,3,4,5,6],[.050,.038,.0306,.0256,.0222,.0195,.0174],4,0.5))

#print(F_CSTR.GD_VCA(.05,.4,[0,1,2,3,4,5,6],[.050,.038,.0306,.0256,.0222,.0195,.0174],4,0.5))








#print(F_CSTR.C_rA([0,50,100,150,200,250,300],[.050,.038,.0306,.0256,.0222,.0195,.0174],5))		
#print(F_CSTR.V_r([0,.1,.2,.4,.6,.7,.8],[.45,.37,.30,.195,.113,.079,.05],.5))
#print(F_CSTR.V_Ctes(.05,.4,[0,1,2,3,4,5,6],[.05,.038,.0306,.0256,.0222,.0195,.0174],4))

''' 	#Se define una funcion para determinar volumenes de reacotr para datos de conversión
	def V_r(X,r_A,FA_0):
		#Se calcula FA_0/r_A
		inv_rA=FA_0/np.asarray(r_A) #Para hacer la grafica de levenspiel donde el area bajo la curva es el vol
		#Se procede a calcular el area bajo la curva(volumen del reactor)
		V1=[np.trapz(X,x=[inv_rA[0],inv_rA[k]]) for k in list(range(0,len(inv_rA)))] 
		#Se debe calcular el volumen sobrante  para restarlo
		V2=inv_rA*X-np.asarray(V1)
		#Se suman los volumenes
		V=np.asarray(V1)+np.asarray(V2)
		plt.figure(1)
		plt.plot(V,X,"o-")
		plt.grid(True)
		plt.figure(2)
		plt.plot(V,r_A,"o-")
		plt.grid(True)
		return plt.show()
	#Se define funcion para calcular volumenes de reactor en funcion de alfa y K
	def V_Ctes(CA_0,FA_0,t,CA,n):	
		#Se calculan alfa y K
		Ctes=F_CSTR.C_rA(t,CA,n)
		#Se calcula conversion y velocidad de reaccion 
		X=1-(np.asarray(CA)**Ctes[0]/CA_0)
		r_A=Ctes[1]*np.asarray(CA)**Ctes[0]
		#Se calcula volumen
		return F_CSTR.V_r(X,r_A,FA_0)'''