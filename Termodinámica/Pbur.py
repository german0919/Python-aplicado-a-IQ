from thermo import Chemical
from thermo import VaporPressure
from thermo import PRMIX, SRKMIX, VDWMIX
from chemicals import search_chemical
import thermo.unifac
from thermo.unifac import UFIP, UFSG, UNIFAC
#UFIP parametros de interacción para el modelo UNIFAC original 
#UFSG o UFMG son los datos para el modelo original de UNIFAC
import numpy as np 
from tabulate import tabulate
from matplotlib import pyplot as plt

class P_b(object): 
	#T[=]K ; P[=]Pa
	def Psat(Componente1,Componente2,Ti):
		T_b=[Chemical(Componente1).Tb, Chemical(Componente2).Tb]
		P_c=[Chemical(Componente1).Pc, Chemical(Componente2).Pc]
		T_c=[Chemical(Componente1).Tc, Chemical(Componente2).Tc]
		Omega=[Chemical(Componente1).omega, Chemical(Componente2).omega]
		Comp1 = VaporPressure(Tb=T_b[0], Tc=T_c[0], Pc=P_c[0], omega=Omega[0])
		Comp2= VaporPressure(Tb=T_b[1], Tc=T_c[1], Pc=P_c[1], omega=Omega[1])
		Psati=[Comp1(Ti),Comp2(Ti)]
		return Psati
	#Cálculo de gammas mediante el uso del método UNIFAC
	def Ufac(Componente1,Componente2,Ti,xs1):
		#Se llama la función los datos se almacenan en el formato InChI key bool bool bool recuento de subgrupos ... recuento de subgrupos recuento de subgrupos ... donde los bools se refieren a si las asignaciones originales de UNIFAC, UNIFAC modificado y PSRK se completaron correctamente o no. Los subgrupos y su recuento tienen una longitud indefinida.
		thermo.unifac.load_group_assignments_DDBST()
		A=search_chemical(Componente1).InChI_key
		B=search_chemical(Componente2).InChI_key
		ka=thermo.unifac.DDBST_UNIFAC_assignments[A]
		kb=thermo.unifac.DDBST_UNIFAC_assignments[B]
		GE = UNIFAC.from_subgroups(chemgroups=[ka, kb], T=Ti, xs=[xs1, 1-xs1], version=0, interaction_data=UFIP, subgroups=UFSG)
		gamma_i=GE.gammas()
		return gamma_i
	#Cálculo de coeficientes de fugacidad para la fase vapor mediante EOS	 	
	def PHI_eos(Componente1,Componente2,Ti,Pi,xs1):
		P_c=[Chemical(Componente1).Pc, Chemical(Componente2).Pc]
		T_c=[Chemical(Componente1).Tc, Chemical(Componente2).Tc]
		Omega=[Chemical(Componente1).omega, Chemical(Componente2).omega]
		eos_mix=PRMIX(T=Ti, P=Pi, Tcs=T_c, Pcs=P_c, omegas=Omega, zs=[xs1, 1-xs1], kijs=[[0.0, 0.0], [0.0, 0.0]])
		philv=eos_mix.phis_g
		return philv	
	#Función para el cálculo de Presión de burbuja en Kpa
	def f_pb(Componente1,Componente2,Ti,xs1):
		Pi_sat=P_b.Psat(Componente1,Componente2,Ti)
		gammai=P_b.Ufac(Componente1,Componente2,Ti,xs1)
		phi_i=np.asarray([1,1])
		Pi=((xs1*gammai[0]*Pi_sat[0])/phi_i[0])+(((1-xs1)*gammai[1]*Pi_sat[1])/phi_i[1])
		tol=1e-5
		dif_p=tol+1
		n=0
		while (dif_p>=tol): 
			phi=P_b.PHI_eos(Componente1,Componente2,Ti,Pi,xs1)
			phi_i=phi
			yi=[((xs1*gammai[0]*Pi_sat[0])/(phi[0]*Pi)),(((1-xs1)*gammai[1]*Pi_sat[1])/(phi[1]*Pi))]
			P=((xs1*gammai[0]*Pi_sat[0])/phi[0])+(((1-xs1)*gammai[1]*Pi_sat[1])/phi[1])
			Pi=P					
			dif_p=np.sqrt(((P-Pi)**2))	
			n=n+1
			#print("Se llega a la convergencia en la iteración número:",n)
		return np.array(P)/1000,yi	
	#Gráfico de EVL	
	def G_pb(Componente1,Componente2,Ti):
		x=np.linspace(0,1,50)
		P=[P_b.f_pb(Componente1,Componente2,Ti,x[k])[0] for k in list(range(0,len(x)))]		
		y=[P_b.f_pb(Componente1,Componente2,Ti,x[k])[1][0] for k in list(range(0,len(x)))]		
		plt.plot(x,P,'b--',label="Líquido")
		plt.plot(y,P,'r--',label="Vapor")
		plt.title("Diagrama EVL \""+Componente1+ "\"-\""+ Componente2)
		plt.xlabel("x-y")
		plt.ylabel("Presión (kPa)")
		plt.legend(loc="best")
		plt.grid(True)
		#Tabla
		Datos=np.transpose([x,y,P])
		Encabezados=["x1","y1","P_b(kPa)"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="grid"))		
		return plt.show()
print(P_b.G_pb("chloroform","ethanol",55+273.15))
#print(P_b.Psat("123-35-3","5392-40-5",300))->tocaria ingresar los valores Tb,Pc,Tc,Omega


	