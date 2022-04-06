#German Hernández

#-----------------------------------------------------------------------------
import numpy as np 
from matplotlib import pyplot as plt
from scipy import integrate
from  sklearn.metrics import r2_score
import seaborn as sns
from tabulate import tabulate
#-------------------------------------------------------------------------------
class CPFR(object):
	'''Funcion para encontrar el mejor ajuste a los datos ingresados, basandose en su R^2.
	   Devuelve coeficientes del polinomio, orden y  R^2. '''
	def Best_reg(x,y):
		x=np.asarray(x)
		y=np.asarray(y)
		pol=[]
		r2=[]
		for n in list(range(0,len(x))):
			poli=np.polyfit(x,y,n)
			pol.append(poli)
			model=np.poly1d(pol[n]) #Hacemos que todos los polinomios se puedan evaluar
			R2=r2_score(y,model(x))	#Encontramos el R^2 segun los datos 
			r2.append(R2)
			if 0.90<r2[n]<1:
				best_pol=pol[n]
				r_s=r2[n]
		return best_pol,len(best_pol)-1,r_s

	#Se define funcion para obtener constantes K y alfa para nuestra ecuacion de velocidad
	def CrA(t,CA):
		'''t=tiempo(list) unidades correpondientes
		   CA=Concentracion(list) unidades correspondientes
        Se genera el mejor ajuste segun su R^2  para obtener un polinmio CA(t)'''
		P_CA=(-1)*CPFR.Best_reg(t,CA)[0] #La rxn es de la forma A->B #Se multiplica por (-1)
		#Se deriva el polinomio 
		dCA_t=np.polyder(P_CA)
		#Se evalua la derivada
		dfCA_et=np.polyval(dCA_t,t)
		#Se obtiene logaritmo de la derivada y de la CA
		L_CA_dCA=np.polyfit(np.log(CA),np.log(dfCA_et),1)
		#Se encuentra alfa y k
		alpha=L_CA_dCA[0]
		K=np.exp(L_CA_dCA[1])
		return alpha,K

	def VR(X,r_A,FA_0):
		# V= dm^3
		r_A=np.asarray(r_A)
		#Se calcula FA_0/r_A
		FA_rA=FA_0/r_A #Para hacer la grafica de levenspiel donde el area bajo la curva es el vol
		#Se procede a calcular el area bajo la curva
		V=[]
		for k in range(0,len(X)):
			v=integrate.simps(FA_rA[0:k+1],X[0:k+1])
			V.append(v)
		return np.asarray(V),X,r_A		
	
	def VC(t,CA,CA_0,FA_0):	
		CA=np.asarray(CA)
		#Se calculan alfa y K
		A_K=CPFR.CrA(t,CA)
		#Se calcula la conversion y velocidad de reaccion  
		X=1-((CA**A_K[0])/CA_0)
		r_A=A_K[1]*(CA**A_K[0])
		rta=CPFR.VR(X,r_A,FA_0)
		#Se calcula volumen
		return rta[0],rta[1],rta[2],t,CA

	#Comportamiento de la conversion y la velocidad de reacción respecto al volumen
	def GVR(X,r_A,FA_0):	
		Vr=CPFR.VR(X,r_A,FA_0)[0]
		X=CPFR.VR(X,r_A,FA_0)[1]
		r_A=CPFR.VR(X,r_A,FA_0)[2]
		plt.style.use('seaborn-ticks')
		fig,axs=plt.subplots(1,2)
		axs[0].plot(Vr,X,"bo--",markersize=4)
		axs[0].set_xlabel("V(dm^3)")
		axs[0].set_ylabel("X")
		axs[0].grid(True)
		axs[1].plot(Vr,r_A,"ro--",markersize=4)
		axs[1].set_xlabel("V(dm^3)")
		axs[1].set_ylabel("-rA (mol*m^3/s)")
		axs[1].grid(True)
		#Tabla
		Datos=np.transpose([X,r_A,Vr])
		Encabezados=["X","r_A","V"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="fancy_grid",stralign="center"))				
		return plt.show()

	def GVC(t,CA,CA_0,FA_0):	
		VC=CPFR.VC(t,CA,CA_0,FA_0,)[0]
		X=CPFR.VC(t,CA,CA_0,FA_0,)[1]
		r_A=CPFR.VC(t,CA,CA_0,FA_0,)[2]
		#plt.style.use('seaborn-ticks')
		#sns.set_context("talk")
		sns.set_context("notebook", font_scale=0.9, rc={"lines.linewidth": 2})
		fig,axs=plt.subplots(1,2)
		axs[0].plot(VC,X,"bo-",markersize=4)
		axs[0].set_xlabel("V(dm^3)")
		axs[0].set_ylabel("X")
		axs[0].grid(True)
		axs[1].plot(VC,r_A,"ro-",markersize=4)
		axs[1].set_xlabel("V(dm^3)")
		axs[1].set_ylabel("-rA (mol*m^3/s)")
		axs[1].grid(True)
		Datos=np.transpose([t,CA,X,r_A,VC])
		Encabezados=["t","CA","X","r_A","V"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="fancy_grid",stralign="center"))						
		return plt.show()
	'''Funcion para calculo de las dimensiones del reactor, que grafique la longitud y el diametro.
	   Se conoce X y r_A '''
	def GLDR(X,r_A,FA_0):
		V=CPFR.VR(X,r_A,FA_0)[0]
		D=((4/3)*np.pi*V)**0.5
		L=3*D
		sns.set_context("notebook", font_scale=0.9, rc={"lines.linewidth": 2})
		fig,axs=plt.subplots(2,2)
		axs[0][0].plot(D,X,"bo-",markersize=4)
		axs[0][0].set_xlabel("D(dm)")
		axs[0][0].set_ylabel("X")
		axs[0][0].grid(True)
		axs[0][1].plot(D,r_A,"ro-",markersize=4)
		axs[0][1].set_xlabel("D(dm)")
		axs[0][1].set_ylabel("-rA (mol*m^3/s)")
		axs[0][1].grid(True)		
		axs[1][1].plot(L,r_A,"ro-",markersize=4)
		axs[1][1].set_xlabel("L(dm)")
		axs[1][1].set_ylabel("-rA (mol*m^3/s)")
		axs[1][1].grid(True)				
		axs[1][0].plot(L,X,"bo-",markersize=4)
		axs[1][0].set_xlabel("L(dm)")
		axs[1][0].set_ylabel("X")
		axs[1][0].grid(True)		
		Datos=np.transpose([X,r_A,V,D,L])
		Encabezados=["X","r_A","V","D","L"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="fancy_grid",stralign="center"))										
		return plt.show()
	'''Funcion para calculo de las dimensiones del reactor, que grafique la longitud y el diametro.
	   Se conocen datos de concentracion y tiempo'''
	def GLDC(t,CA,CA_0,FA_0):
		V=CPFR.VC(t,CA,CA_0,FA_0)[0]
		X=CPFR.VC(t,CA,CA_0,FA_0,)[1]
		r_A=CPFR.VC(t,CA,CA_0,FA_0,)[2]
		D=((4/3)*np.pi*V)**0.5
		L=3*D
		sns.set_context("notebook", font_scale=0.9, rc={"lines.linewidth": 2})
		fig,axs=plt.subplots(2,2)
		axs[0][0].plot(D,X,"bo-",markersize=4)
		axs[0][0].set_xlabel("D(dm)")
		axs[0][0].set_ylabel("X")
		axs[0][0].grid(True)
		axs[0][1].plot(D,r_A,"ro-",markersize=4)
		axs[0][1].set_xlabel("D(dm)")
		axs[0][1].set_ylabel("-rA (mol*m^3/s)")
		axs[0][1].grid(True)		
		axs[1][1].plot(L,r_A,"ro-",markersize=4)
		axs[1][1].set_xlabel("L(dm)")
		axs[1][1].set_ylabel("-rA (mol*m^3/s)")
		axs[1][1].grid(True)				
		axs[1][0].plot(L,X,"bo-",markersize=4)
		axs[1][0].set_xlabel("L(dm)")
		axs[1][0].set_ylabel("X")
		axs[1][0].grid(True)	
		Datos=np.transpose([t,CA,X,r_A,V,D,L])
		Encabezados=["t","CA","X","r_A","V","D","L"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="fancy_grid",stralign="center"))						
		return plt.show()


#Estructurar un programa donde se indiquen opciones para seleccionar el tipo de funcion a utilizar
#Escribir un encabezado

#print(CPFR.Best_reg([0,50,100,150,200,250,300],[.050,.038,.0306,.0256,.0222,.0195,.0174]))
#print(CPFR.CrA([0,50,100,150,200,250,300],[.050,.038,.0306,.0256,.0222,.0195,.0174]))
#print(CPFR.VR([0,0.1,0.2,0.4,0.6,0.7,0.8],[0.4494382 , 0.37037037, 0.30075188, 0.19512195, 0.11299435,0.07905138, 0.05],0.4))
#print(CPFR.GVR([0,0.1,0.2,0.4,0.6,0.7,0.8],[0.45,0.37,0.30,0.195, 0.113,0.079, 0.05],0.4))
#print(CPFR.VC([0,3,6,9,12,15,18,21],[0.0635,0.0597,.0564,.0531,.0502,.0474,.0446,.0422],0.03,0.1))
#print(CPFR.GVC([0,3,6,9,12,15,18,21],[0.00635,0.00597,.00564,.00531,.00502,.00474,.00446,.00422],0.03,0.1))
#print(CPFR.GLDR([0,0.1,0.2,0.4,0.6,0.7,0.8],[0.4494382 , 0.37037037, 0.30075188, 0.19512195, 0.11299435,0.07905138, 0.05],0.4))
#print(CPFR.GLDC([0,3,6,9,12,15,18,21],[0.00635,0.00597,.00564,.00531,.00502,.00474,.00446,.00422],0.03,0.1))

#Otra forma de hacer la integral es generando un polinomio e integrandolo
'''
	def VRI(X,r_A,FA_0):
		#Se calcula FA_0/r_A
		FA_rA=FA_0/np.asarray(r_A) #Para hacer la grafica de levenspiel donde el area bajo la curva es el vol
		#Se calcula un polinomio para los datos
		polnom=CPFR.Best_reg(X,FA_rA)[0]
		#Se procede a calcular el area bajo la curva
		Vi=np.polyint(polnom)
		V=np.polyval(Vi,X)
		return np.asarray(V)*1000,X,r_A				'''