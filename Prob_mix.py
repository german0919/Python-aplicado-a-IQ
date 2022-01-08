#German Hernández

#Factor de compresibilidad y volumen molar mediante Peng-Robinson para mezclas
import numpy as np 
from matplotlib import pyplot as plt
from scipy.optimize import fsolve 
from tabulate import tabulate
class ZVPR(object): 
	#T[=]K; P[=]bar ;Vc[=] cm^3/mol
	Dic_Prop={"METANO":{"W":0.012,"Tc":190.6,"Pc":45.99,"Vc":98.6},
	          "ETANO":{"W":0.100,"Tc":305.3,"Pc":48.72,"Vc":145.5},
	          "PROPANO":{"W":0.152,"Tc":369.8,"Pc":42.48,"Vc":200},
	          "PROPILENO":{"W":0.140,"Tc":365.6,"Pc":46.65,"Vc":188.4},
	          "ACETILENO":{"W":0.187,"Tc":308.3,"Pc":61.39,"Vc":113},
	          "BENCENO":{"W":0.210,"Tc":562.2,"Pc":48.98,"Vc":259},
	          "TOLUENO":{"W":0.262,"Tc":591.8,"Pc":41.06,"Vc":316},
	          "ETILBENCENO":{"W":0.303,"Tc":617.2,"Pc":36.06,"Vc":374},
	          "CUMENO":{"W":0.326,"Tc":631.1,"Pc":32.09,"Vc":427},
	          "ACETONA":{"W":0.307,"Tc":508.2,"Pc":47.01,"Vc":209},
	          "ESTIRENO":{"W":0.297,"Tc":636.0,"Pc":38.40,"Vc":352},
	          "NAFTALENO":{"W":0.302,"Tc":748.4,"Pc":40.51,"Vc":413},
	          "BUTANO":{"W":0.200,"Tc":425.1,"Pc":37.96,"Vc":255},
	          "PENTANO":{"W":0.252,"Tc":469.7,"Pc":33.70,"Vc":313},
	          "HEXANO":{"W":0.301,"Tc":507.6,"Pc":30.25,"Vc":371},
	          "METANOL":{"W":0.564,"Tc":512.6,"Pc":80.97,"Vc":118},
	          "ETANOL":{"W":0.645,"Tc":513.9,"Pc":61.48,"Vc":167},
	          "CLOROFORMO":{"W":0.222,"Tc":536.4,"Pc":54.72,"Vc":239},
	          "AGUA":{"W":0.345,"Tc":647.1,"Pc":220.55,"Vc":55.9},
	          "AMONIACO":{"W":0.253,"Tc":405.7,"Pc":112.80,"Vc":72.5},
	          "ACIDO NITRICO":{"W":0.714,"Tc":520.0,"Pc":68.90,"Vc":145}} #Ref:Smith Van Ness Introduccion a la termodinamica en IQ
	

	def ZV(Componentes,T,P,zs):
		R=83.14
		T=np.asarray(T)
		P=np.asarray(P)
		zs=np.asarray(zs)
		Dic_pc=ZVPR.Dic_Prop
		#Parametros sustancia pura
		w=[]
		Tc=[]
		Pc=[]
		Vc=[]
		ai=[]
		bi=[]
		Pr=[]
		Tr=[]
		alpha_i=[]
		Bi=[]
		Aii=[]
		for k in list(range(0,len(zs))):
			tc=Dic_pc[Componentes[k]]["Tc"]
			w_a=Dic_pc[Componentes[k]]["W"]
			pc=Dic_pc[Componentes[k]]["Pc"]
			vc=Dic_pc[Componentes[k]]["Vc"]
			w.append(w_a)
			Tc.append(tc)
			Pc.append(pc)
			Vc.append(vc)
			pr=(P/Pc[k])
			tr=(T/Tc[k])
			a=0.45724*(R**2)*((Tc[k]**2)/(Pc[k]))
			b=0.07780*R*(Tc[k]/Pc[k])
			Tr.append(tr)
			Pr.append(pr)
			ai.append(a)
			bi.append(b)
			alpha=(1+((0.37464+1.54226*w[k]-0.26992*w[k]**2)*(1-Tr[k]**0.5)))**2
			alpha_i.append(alpha)
			B_i=bi[k]*(P/(R*T))
			Bi.append(B_i)
			A_i=0.45724*(Pr[k]/Tr[k]**2)*alpha_i[k]
			Aii.append(A_i)
		#Parametros para mezclas
		kij=[]
		for i in list(range(0,len(zs))):
			for j in list(range(0,len(zs))):
				k_ij=1-((8*(Vc[i]*Vc[j])**0.5)/(Vc[i]**(1/3)+Vc[j]**(1/3))**3)
				kij.append(k_ij)
		kij=np.asarray(kij).reshape(-1,len(zs)) 
		Aij=[]
		for i in list(range(0,len(zs))):
			for j in list(range(0,len(zs))):
				A_ij=(1-kij[i][j])*((Aii[i]*Aii[j]))**0.5
				Aij.append(A_ij)
		Aij=np.asarray(Aij).reshape(-1,len(zs))         
		s=0
		A=[]
		for i in list(range(0,len(zs))):
			for j in list(range(0,len(zs))):
				s += zs[i]*zs[j]*Aij[i][j]       
		A.append(s)        
		B=0        
		for m in list(range(0,len(zs))):
			B += zs[m]*(Bi[m]) 

		f= lambda z: z**3-(1-B)*z**2+(A-3*B**2-2*B)*z-(A*B-B**2-B**3)
		Z=fsolve(f,1)		
		v=(np.asarray(Z)*R*T)/(P)			
		return Z,v,A,B	

	#Grafico a presion constante(Z y V vs T	)
	def GZP(Componentes,P,zs,Tmin,Tmax,n):
		T=np.linspace(Tmin,Tmax,n)
		Z=[]
		Vm=[]
		for k in list(range(0,len(T))):
			z=ZVPR.ZV(Componentes,T[k],P,zs)[0]
			Z.append(z)
			v=ZVPR.ZV(Componentes,T[k],P,zs)[1]
			Vm.append(v)		
		fig, axs=plt.subplots(1,2)
		axs[0].plot(T,Z,"r*--")
		axs[0].set_title("Z vs T")
		axs[0].set_ylabel("Z")
		axs[0].set_xlabel("T")
		axs[0].grid(True)
		axs[1].plot(T,Vm,"c*--")
		axs[1].set_title("Vm vs T")
		axs[1].set_ylabel("Vm ")
		axs[1].set_xlabel("T")
		axs[1].grid(True)
		#Tabla
		Datos=np.transpose([T,Z,Vm])
		Encabezados=["T(K)","Z","Vm(m^3)"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="grid"))					
		return plt.show()
	#Grafico a temperatura constante(Z y V vs P	)		
	def GZT(Componentes,T,zs,Pmin,Pmax,n):
		P=np.linspace(Pmin,Pmax,n)
		Z=[]
		Vm=[]
		for k in list(range(0,len(P))):
			z=ZVPR.ZV(Componentes,T,P[k],zs)[0]
			Z.append(z)
			v=ZVPR.ZV(Componentes,T,P[k],zs)[1]
			Vm.append(v)				
		fig, axs=plt.subplots(1,2)
		axs[0].plot(P,Z,"r*--")
		axs[0].set_title("Z vs P")
		axs[0].set_ylabel("Z")
		axs[0].set_xlabel("P")
		axs[0].grid(True)
		axs[1].plot(P,Vm,"c*--")
		axs[1].set_title("Vm vs P")
		axs[1].set_ylabel("Vm ")
		axs[1].set_xlabel("P")
		axs[1].grid(True)	
		#Tabla
		Datos=np.transpose([P,Z,Vm])
		Encabezados=["P(bar)","Z","Vm(m^3)"]
		print(tabulate(Datos,headers=Encabezados,tablefmt="grid"))							
		return plt.show()					

#print(ZVPR.GZP(["PROPANO","BUTANO","PENTANO"],[40],[0.50,0.20,0.30],393,500,50))
#print(ZVPR.ZV(["PROPANO","BUTANO","PENTANO"],[393],[40],[0.50,0.20,0.30]))
#print(ZVPR.ZV(["PROPANO","BUTANO","PENTANO"],[300],[1],[0.50,0.20,0.30]))
#print(ZVPR.ZV(["PROPANO","PENTANO"],[393],[40],[0.50,0.50]))
#print(ZVPR.ZV(["PROPANO","PENTANO"],[393],[30],[0.50,0.50]))
#print(ZVPR.ZV(["PROPANO","PENTANO","BUTANO"],[393],[40],[0.50,0.3,0.2]))
#print(ZVPR.GZT(["PROPANO","BUTANO"],[395],[0.3,0.7],40,50,10))




'''def ZV(Componentes,T,P,zs):
		R=83.14
		T=np.asarray(T)
		P=np.asarray(P)
		zs=np.asarray(zs)
		Dic_pc=ZVPR.Dic_Prop
		#Parametros sustancia pura
		w=[]
		Tc=[]
		Pc=[]
		Vc=[]
		ai=[]
		bi=[]
		Pr=[]
		Tr=[]
		alpha_i=[]
		Bi=[]
		Aii=[]
		for k in list(range(0,len(zs))):
			tc=Dic_pc[Componentes[k]]["Tc"]
			w_a=Dic_pc[Componentes[k]]["W"]
			pc=Dic_pc[Componentes[k]]["Pc"]
			vc=Dic_pc[Componentes[k]]["Vc"]
			w.append(w_a)
			Tc.append(tc)
			Pc.append(pc)
			Vc.append(vc)
			pr=(P/Pc[k])
			tr=(T/Tc[k])
			a=0.45724*(R**2)*((Tc[k]**2)/(Pc[k]))
			b=0.07780*R*(Tc[k]/Pc[k])
			Tr.append(tr)
			Pr.append(pr)
			ai.append(a)
			bi.append(b)
			alpha=(1+((0.37464+1.54226*w[k]-0.26992*w[k]**2)*(1-Tr[k]**0.5)))**2
			alpha_i.append(alpha)
			B_i=bi[k]*(P/(R*T))
			Bi.append(B_i)
			A_i=0.45724*(Pr[k]/Tr[k]**2)*alpha_i[k]
			Aii.append(A_i)
		#Parametros para mezclas
		#Mezcla binaria	
		if len(zs)==2:
			k12=1-((8*(Vc[0]*Vc[1])**0.5)/(Vc[0]**(1/3)+Vc[1]**(1/3))**3)
			k21=1-((8*(Vc[1]*Vc[0])**0.5)/(Vc[1]**(1/3)+Vc[0]**(1/3))**3)
			kij=np.asarray([k12,k21])
		#Mezcla ternaria	
		if len(zs)==3:
			k12=1-((8*(Vc[0]*Vc[1])**0.5)/(Vc[0]**(1/3)+Vc[1]**(1/3))**3)
			k21=k12
			k13=1-((8*(Vc[0]*Vc[2])**0.5)/(Vc[0]**(1/3)+Vc[2]**(1/3))**3)
			k31=k13
			k23=1-((8*(Vc[1]*Vc[2])**0.5)/(Vc[1]**(1/3)+Vc[2]**(1/3))**3)
			k32=k23
			kij=np.asarray([k12,k13,k23])    
		B=0       
		for m in list(range(0,len(zs))): 
			B += zs[m]*(Bi[m]) 	
		if len(zs)==2:
			A_12=(1-kij[0])*(Aii[0]*Aii[1])**0.5
			A_21=(1-kij[0])*(Aii[1]*Aii[0])**0.5
			Aij=np.asarray([A_12,A_21])
		if len(zs)==3:
			A_12=(1-kij[0])*(Aii[0]*Aii[1])**0.5
			A_21=(1-kij[0])*(Aii[1]*Aii[0])**0.5
			A_13=(1-kij[1])*(Aii[0]*Aii[2])**0.5
			A_31=(1-kij[1])*(Aii[2]*Aii[0])**0.5
			A_23=(1-kij[2])*(Aii[1]*Aii[2])**0.5
			A_32=(1-kij[2])*(Aii[2]*Aii[1])**0.5
			Aij=np.asarray([A_12,A_13,A_23])
		#A=zs[0]**2*Aij[0]+2*zs[0]*zs[1]*Aij[1]+zs[1]**2*Aij[2]->para binaria
		A=0
		for i in list(range(0,len(zs))):  #Es probable que haya un error en esta sumatoria
			for j in list(range(0,len(zs))):
				A += zs[i]*zs[j]*Aij[i]		
		if len(T)>=1: 
			A=list(A)
			B=list(B)
			Z=[]
			zi=np.ones(len(T))
			for k in list(range(0,len(T))):
				f= lambda z: z**3-(1-B[k])*z**2+(A[k]-3*B[k]**2-2*B[k])*z-(A[k]*B[k]-B[k]**2-B[k]**3)
				sol=fsolve(f,zi)
				Z.append(sol[0])
		else:
			f= lambda z: z**3-(1-B)*z**2+(A-3*B**2-2*B)*z-(A*B-B**2-B**3)
			Z=fsolve(f,1)		
		if len(P)>=1: 
			A=list(A)
			B=list(B)
			Z=[]
			zi=np.ones(len(P))
			for k in list(range(0,len(P))):
				f= lambda z: z**3-(1-B[k])*z**2+(A[k]-3*B[k]**2-2*B[k])*z-(A[k]*B[k]-B[k]**2-B[k]**3)
				sol=fsolve(f,zi)
				Z.append(sol[0])
		else:
			f= lambda z: z**3-(1-B)*z**2+(A-3*B**2-2*B)*z-(A*B-B**2-B**3)
			Z=fsolve(f,1)									
		v=(np.asarray(Z)*R*T)/(P)
		return Z,v '''






	          
    
