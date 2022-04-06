#German Hernández
import numpy as np 
from thermo import Chemical
from matplotlib import pyplot as plt
from tabulate import tabulate
class PB(object):
	Dic_Propc={"BENCENO":{"W":0.210,"Tc":562.2,"Pc":48.98,"Vc":259,"Zc":0.271},
	          "ACETONA":{"W":0.307,"Tc":508.2,"Pc":47.01,"Vc":209,"Zc":0.233},
	          "HEXANO":{"W":0.301,"Tc":507.6,"Pc":30.25,"Vc":371,"Zc":0.266},
	          "METANOL":{"W":0.564,"Tc":512.6,"Pc":80.97,"Vc":118,"Zc":0.224},
	          "ETANOL":{"W":0.645,"Tc":513.9,"Pc":61.48,"Vc":167,"Zc":0.240},
	          "1-BUTANOL":{"W":0.594,"Tc":563.1,"Pc":44.23,"Vc":275,"Zc":0.260},
	          "1-PROPANOL":{"W":0.622,"Tc":536.8,"Pc":51.75,"Vc":219,"Zc":0.254},	          	          	       
	          "2-PROPANOL":{"W":0.668,"Tc":508.3,"Pc":47.62,"Vc":220,"Zc":0.248},	
	          "TETRAHIDROFURANO":{"W":0.226,"Tc":540.15,"Pc":51.88,"Vc":223.9,"Zc":0.259},	 
	          "TETRACLORURO DE CARBONO":{"W":0.193,"Tc":556.4,"Pc":45.60,"Vc":276,"Zc":0.272},	#R_Temp:23-100 
	          "ACETATO DE ETILO":{"W":0.366,"Tc":523.3,"Pc":38.80,"Vc":286,"Zc":0.255},	          
	          "ACIDO FORMICO":{"W":0.473,"Tc":580,"Pc":73.90,"Vc":125,"Zc":0.192},	          	          
	          "ACETATO DE METILO":{"W":0.331,"Tc":506.6,"Pc":47.50,"Vc":228,"Zc":0.257},	          	          
	          "CLOROFORMO":{"W":0.222,"Tc":536.4,"Pc":54.72,"Vc":239,"Zc":0.293},
	          "AGUA":{"W":0.345,"Tc":647.1,"Pc":220.55,"Vc":55.9,"Zc":0.229}} #Ref:Smith Van Ness Introduccion a la termodinamica en IQ/ Yaws
#T[=]°C ; P[=]torr
	Dic_Ant={"BENCENO":{"A":6.87987,"B":1196.760,"C":219.161}, #R_Temp:8-80
	          "ACETONA":{"A":7.11714,"B":1210.595,"C":229.664}, #R_Temp:(-13)-55
	          "HEXANO":{"A":6.91058,"B":1189.640,"C":226.280}, #R_Temp:(-30)-170
	          "METANOL":{"A":8.08097,"B":1582.271,"C":239.726}, #R_Temp:15-84
	          "ETANOL":{"A":8.11220,"B":1592.864,"C":226.184}, #R_Temp:20-93
	          "1-BUTANOL":{"A":7.36366,"B":1305.198,"C":173.427}, #R_Temp:89-126
	          "1-PROPANOL":{"A":8.37895,"B":1788.020,"C":227.438},	 #R_Temp:(-15)-98         	          	       
	          "2-PROPANOL":{"A":8.87829,"B":2010.320,"C":252.636},	#R_Temp:(-26)-83
	          "TETRAHIDROFURANO":{"A":6.99515,"B":1202.290,"C":226.254},	#R_Temp:23-100 
	          "TETRACLORURO DE CARBONO":{"A":6.84083,"B":1177.910,"C":220.576},	#R_Temp:(-20)-77
	          "ACETATO DE ETILO":{"A":7.10179,"B":1244.951,"C":217.881},	#R_Temp:16-76          
	          "ACIDO FORMICO":{"A":6.94459,"B":1295.260,"C":218.00},	  #R_Temp:36-108      	          
	          "ACETATO DE METILO":{"A":7.06524,"B":1157.630,"C":219.726},	#R_Temp:2-56          	          
	          "CLOROFORMO":{"A":6.95465,"B":1170.966,"C":226.232}, #R_Temp:(-10)-60
	          "AGUA":{"A":8.07131,"B":1730.630,"C":233.426}} #R_Temp:1-100
	          #Ref:Perry tomo IV-13.22
	def C2VIRIAL(Componentes,T):
		R=83.14
		Dic_c=PB.Dic_Propc
		w=[Dic_c[Componentes[0]]["W"],Dic_c[Componentes[1]]["W"]]
		Tc=[Dic_c[Componentes[0]]["Tc"],Dic_c[Componentes[1]]["Tc"]]
		Pc=[Dic_c[Componentes[0]]["Pc"],Dic_c[Componentes[1]]["Pc"]]
		Vc=[Dic_c[Componentes[0]]["Vc"],Dic_c[Componentes[1]]["Vc"]]
		Zc=[Dic_c[Componentes[0]]["Zc"],Dic_c[Componentes[1]]["Zc"]]
		wij=[]
		Tcij=[]
		Zcij=[]
		Vcij=[]
		Pcij=[]
		for i in list(range(0,len(Componentes))):
			for j in list(range(0,len(Componentes))):
				w_ij=(w[i]+w[j])/2
				T_cij=(Tc[i]*Tc[j])**0.5
				Z_cij=(Zc[i]+Zc[j])/2
				V_cij=((Vc[i]**(1/3)+Vc[j]**(1/3))/2)**3
				if w_ij not in wij:
					wij.append(w_ij)
				if T_cij not in Tcij:
					Tcij.append(T_cij)    
				if Z_cij not in Zcij:   
					Zcij.append(Z_cij)    
				if V_cij not in Vcij:   
					Vcij.append(V_cij)        
		Zcij=np.asarray(Zcij)
		Vcij=np.asarray(Vcij)
		Tcij=np.asarray(Tcij)
		wij=np.asarray(wij)
		Pcij=(Zcij*R*Tcij)/Vcij
		Trij=T/Tcij
		B0=0.083-(0.422/(Trij**1.6))
		B1=0.139-(0.172/(Trij**4.2))
		Bcij=B0+wij*B1
		Bij=(Bcij*R*Tcij)/(Pcij)
		#si quiero calcular esto se debe poner la presión y xs
		#B=xs[0]**2*Bij[0]+2*xs[0]*xs[1]*Bij[1]+xs[1]**2*Bij[2]
		#Z=1+((B*P)/(R*T))
		#v=(Z*R*T)/P	
		#El resultado se pasara de cm^3 a m^3	
		return Bij*1e-6 
	def Margules(Componentes,x1):
		x1=np.asarray(x1)
		if Componentes == ["ACETONA","CLOROFORMO"]:
			A12=-0.8404
			A21=-0.5610
		if Componentes == ["ACETONA","METANOL"]:
			A12=0.6184
			A21=0.5788
		if Componentes == ["ACETONA","AGUA"]:
			A12=2.0400
			A21=1.5461					
		if Componentes == ["TETRACLORURO DE CARBONO","BENCENO"]:
			A12=0.0948
			A21=0.0922
		if Componentes == ["CLOROFORMO","METANOL"]:
			A12=0.8320
			A21=1.7365
		if Componentes == ["ETANOL","BENCENO"]:
			A12=1.8362
			A21=1.4717
		if Componentes == ["ETANOL","AGUA"]:
			A12=1.6022
			A21=0.7947
		if Componentes == ["ACETATO DE ETILO","ETANOL"]:
			A12=0.8557
			A21=0.7476
		if Componentes == ["HEXANO","ETANOL"]:
			A12=1.9398
			A21=2.7054
		if Componentes == ["METANOL","BENCENO"]:
			A12=2.1411
			A21=1.7905
		if Componentes == ["METANOL","ACETATO DE ETILO"]:
			A12=1.0016
			A21=1.0517
		if Componentes == ["METANOL","AGUA"]:
			A12=0.7923
			A21=0.5434			
		if Componentes == ["ACETATO DE METILO","METANOL"]:
			A12=0.9605
			A21=1.0120
		if Componentes == ["1-PROPANOL","AGUA"]:
			A12=2.7070
			A21=0.7172
		if Componentes == ["2-PROPANOL","AGUA"]:
			A12=2.3319
			A21=0.8976
		if Componentes == ["TETRAHIDROFURANO","AGUA"]:
			A12=2.8258
			A21=1.9450	
		if Componentes == ["AGUA","ACIDO ACETICO"]:
			A12=0.4178
			A21=0.9533										
		if Componentes == ["AGUA","1-BUTANOL"]:
			A12=0.8608
			A21=3.2051				
		if Componentes == ["AGUA","ACIDO FORMICO"]:
			A12=-0.2966
			A21=-0.2715			
		if Componentes == ["CLOROFORMO","ETANOL"]:
			A12=0.59
			A21=1.42						
		gamma1=np.exp((A12 + 2*(A21 - A12)*x1)*(1-x1)**2)
		gamma2=np.exp((A21 + 2*(A12 - A21)*(1-x1))*x1**2)
		gammas=[gamma1,gamma2]
		return gammas
	def Antoine(Componentes,T):	
		T=T-273.15 #Esta funcion trabaja con °C
		Dic_ctes=PB.Dic_Ant
		A=np.asarray([Dic_ctes[Componentes[0]]["A"],Dic_ctes[Componentes[1]]["A"]])
		B=np.asarray([Dic_ctes[Componentes[0]]["B"],Dic_ctes[Componentes[1]]["B"]])
		C=np.asarray([Dic_ctes[Componentes[0]]["C"],Dic_ctes[Componentes[1]]["C"]])
		Pi_s=10**(A-(B/(T+C)))
		return Pi_s/7.501 #Pasarlos a kPa
	def PHI2V(Componentes,T,P,x1):
		R=83.14
		Bij=PB.C2VIRIAL(Componentes,T)
		Pi_s=PB.Antoine(Componentes,T)
		d12=2*Bij[1]-Bij[2]-Bij[0]
		phi_1=np.exp((Bij[0]*(P-Pi_s[0])+P*(1-x1)**2*d12)/(R*T))
		phi_2=np.exp((Bij[2]*(P-Pi_s[1])+P*x1**2*d12)/(R*T))
		phi_ij=np.asarray([phi_1,phi_2])
		return phi_ij
	def f_pb(Componentes,T,x1):
		Pi_sat=PB.Antoine(Componentes,T)
		gammai=PB.Margules(Componentes,x1)
		phi_i=np.asarray([1,1])
		Pi=((x1*gammai[0]*Pi_sat[0])/phi_i[0])+((((1-x1))*gammai[1]*Pi_sat[1])/phi_i[1])
		tol=1e-5
		dif_p=tol+1
		n=0
		while (dif_p>=tol): 
			phi=PB.PHI2V(Componentes,T,Pi,x1)
			phi_i=phi
			yi=[((x1*gammai[0]*Pi_sat[0])/(phi[0]*Pi)),((((1-x1))*gammai[1]*Pi_sat[1])/(phi[1]*Pi))]
			P=((x1*gammai[0]*Pi_sat[0])/phi[0])+((((1-x1))*gammai[1]*Pi_sat[1])/phi[1])
			Pi=P					
			dif_p=np.sqrt(((P-Pi)**2))	
			n=n+1
			#print("Se llega a la convergencia en la iteración número:",n)
		return np.array(P),yi
	#Gráfico de EVL	
	def G_pb(Componentes,T):
		x1=np.linspace(0,1,30)
		P=[PB.f_pb(Componentes,T,x1[k])[0] for k in list(range(0,len(x1)))]
		y1=[PB.f_pb(Componentes,T,x1[k])[1][0] for k in list(range(0,len(x1)))]
		plt.plot(x1,P,'b--',label="Líquido")
		plt.plot(y1,P,'r--',label="Vapor")
		plt.title("Diagrama EVL")
		plt.xlabel("x-y")
		plt.ylabel("Presión (kPa)")
		plt.legend(loc="best")
		plt.grid(True)
		#Tabla
		Datos=np.transpose([x1,y1,P])
		Encabezados=["x1","y1","P_b(kPa)"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="fancy_grid",stralign="center"))		
		return plt.show()		

#print(PB.Margules(["CLOROFORMO","ETANOL"],0.35))		
#print(PB.C2VIRIAL(["CLOROFORMO","ETANOL"],400))
#print(PB.Antoine(["CLOROFORMO","ETANOL"],55+273.15))
#print(PB.PHI2V(["CLOROFORMO","ETANOL"],55+273.15,72.83,0.35))
#print(PB.f_pb(["CLOROFORMO","ETANOL"],55+273.15,0.35))
#print(PB.G_pb(["ACETONA","AGUA"],65+273.15))
print(PB.G_pb(["CLOROFORMO","ETANOL"],55+273.15))

