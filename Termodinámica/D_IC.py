#German Hernández
from scipy import interpolate
import numpy as np   
#---------Diseño de intercambiadores de calor_ Metodo NTU ------------#
'''Se realizara una clase con las funciones necesarias para estimar la temperatura a la salida
   para distintos tipos de IC suponiendo una temperatura promedio para hallar
   los cp e iterando hasta que el calculo llega a la convergencia'''


class MNTU(object):
	#T[=]°C ; Cp[=] J/kg*K ref:cengel
	Dic_TCP={"ACEITE PARA MOTOR":{"T":[0,20,40,60,80,100,120,140,150],"Cp":[1797,1881,1964,2048,2132,2220,2308,2395,2441]},
         "AGUA LIQUIDA":{"T":[0.01,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360],"Cp":[4217,4105,4194,4185,4182,4180,4178,4178,4179,4180,4181,4183,4185,4187,4190,4193,4197,4201,4206,4212,4217,4229,4244,4263,4286,4311,4340,4370,4410,4460,4500,4610,4760,4970,5280,5750,6540,8240,14690]},
         "AIRE":{"T":[0,5,10,15,20,25,30,80,90,100,120,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000],"Cp":[1006,1006,1006,1007,1007,1007,1007,1008,1008,1009,1011,1013,1016,1019,1023,1033,1044,1056,1069,1081,1093,1115,1135,1153,1169,1184]}}


	def Efectividad(Tipo,F_c,F_h,Tci,Thi,m_c,m_h,U,D,L,Tprom_c,Tprom_h,N):
		#N: numero de pasos por los tubos 
		Dic=MNTU.Dic_TCP
		T_fc=np.asarray(Dic[F_c]["T"])
		T_fh=np.asarray(Dic[F_h]["T"])
		cp_fc=np.asarray(Dic[F_c]["Cp"])
		cp_fh=np.asarray(Dic[F_h]["Cp"])
		in_cpc=interpolate.interp1d(T_fc,cp_fc)
		Cp_c=in_cpc(Tprom_c)
		in_cph=interpolate.interp1d(T_fh,cp_fh)
		Cp_h=in_cph(Tprom_h)
		#Calculo de c
		C_m=[m_c*Cp_c,m_h*Cp_h]
		Cmin=np.min(C_m)
		Cmax=np.max(C_m)
		c=Cmin/Cmax
		#Calculo de NTU
		if Tipo=="Doble Tubo-Flujo paralelo": 
			As=np.pi*D*L
			Ntu=(U*As)/(Cmin)
			eff=(1-np.exp(-Ntu*(1+c)))/(1+c)
		
		elif 	Tipo=="Doble Tubo-Contraflujo": 
			As=np.pi*D*L
			Ntu=(U*As)/(Cmin)			
			eff=(1-np.exp(-Ntu*(1-c)))/(1-(c*np.exp(-Ntu*(1-c))))
		
		elif 	Tipo=="Tubo y coraza":  #un paso por la coraza y 2n pasos por los tubos
			As=N*np.pi*D*L
			Ntu=(U*As)/(Cmin)	
			eff=2*((1+c+((np.sqrt(1+c**2)*(1+np.exp(-Ntu*(np.sqrt(1+c**2))))/(1-(np.exp(-Ntu*(np.sqrt(1+c**2))))))))**(-1))

		elif 	Tipo=="Flujo cruzado-No mix": 
			#As=1
			As=N*np.pi*D*L
			Ntu=(U*As)/(Cmin)
			eff=1-np.exp((Ntu**0.22/c)*(np.exp(-c*Ntu**0.78)-1))

		elif 	Tipo=="Flujo cruzado-Cmax mix":  #Cmax mix, Cmin no mix
			As=N*np.pi*D*L
			Ntu=(U*As)/(Cmin)
			eff=(1/c)*(1-np.exp((1-c)*(1-np.exp(-Ntu))))

		elif 	Tipo=="Flujo cruzado-Cmin mix":  #Cmax no mix, Cmin mix
			As=N*np.pi*D*L
			Ntu=(U*As)/(Cmin)
			eff=1-np.exp((-1/c)*(1-np.exp(-c*Ntu)))

		elif c==0:
			As=N*np.pi*D*L
			Ntu=(U*As)/(Cmin)	
			eff=1-np.exp(-Ntu)

		Q=eff*Cmin*(Thi-Tci)
		Tco=Tci+(Q/(m_c*Cp_c))
		Tho=Thi-(Q/(m_h*Cp_h))
		return Tco,Tho,Cp_h,Cp_c,Ntu,eff,Q

	def NTU(Tipo,F_c,F_h,Tci,Thi,m_c,m_h,U,D,L,Tprom_c,Tprom_h,N):
		Tco_0=MNTU.Efectividad(Tipo,F_c,F_h,Tci,Thi,m_c,m_h,U,D,L,Tprom_c,Tprom_h,N)[0]
		Tho_0=MNTU.Efectividad(Tipo,F_c,F_h,Tci,Thi,m_c,m_h,U,D,L,Tprom_c,Tprom_h,N)[1]
		tol=1e-5
		er1=tol+1
		er2=tol+1
		n=0
		while (er1>=tol and er2>=tol):
			Tpromc=(Tci+Tco_0)/2	
			Tprom_c=Tpromc
			Tpromh=(Thi+Tho_0)/2		
			Tpromh=Tprom_h
			er1=np.sqrt(((Tpromc-Tprom_c)**2))	
			er2=np.sqrt(((Tpromh-Tprom_h)**2))	
			n=n+1 
		efec=MNTU.Efectividad(Tipo,F_c,F_h,Tci,Thi,m_c,m_h,U,D,L,Tpromc,Tpromh,N)				
		return efec
	def InterCP(F_c,F_h,Tprom_c,Tprom_h):
		#N: numero de pasos por los tubos 
		Dic=MNTU.Dic_TCP
		T_fc=np.asarray(Dic[F_c]["T"])
		T_fh=np.asarray(Dic[F_h]["T"])
		cp_fc=np.asarray(Dic[F_c]["Cp"])
		cp_fh=np.asarray(Dic[F_h]["Cp"])
		in_cpc=interpolate.interp1d(T_fc,cp_fc)
		Cp_c=in_cpc(Tprom_c)
		in_cph=interpolate.interp1d(T_fh,cp_fh)
		Cp_h=in_cph(Tprom_h)
		#Calculo de c
		#C_m=[m_c*Cp_c,m_h*Cp_h]
		#Cmin=np.min(C_m)
		#Cmax=np.max(C_m)
		#c=Cmin/Cmax
		return Cp_c,Cp_h

#print(MNTU.NTU("Doble Tubo-Contraflujo","AGUA LIQUIDA","AIRE",22,90,0.1,0.3,80,0.012,12,40,90,None))

#print(MNTU.NTU("Flujo cruzado-No mix","AIRE","ACEITE PARA MOTOR",30,75,0.21,0.026,53,None,None,45,63,None)) #nos dan el  As=1 #11-114

#print(MNTU.NTU("Tubo y coraza","AGUA LIQUIDA","ACEITE PARA MOTOR",18,160,0.1,0.2,340,0.018,3,40,100,12))

print(MNTU.InterCP("AGUA LIQUIDA","AGUA LIQUIDA",25.25,42.8))
print(MNTU.InterCP("AGUA LIQUIDA","AGUA LIQUIDA",25.35,42.7))
print(MNTU.InterCP("AGUA LIQUIDA","AGUA LIQUIDA",25.4,42.7))
print(MNTU.InterCP("AGUA LIQUIDA","AGUA LIQUIDA",25.45,43.2))
print(MNTU.InterCP("AGUA LIQUIDA","AGUA LIQUIDA",25.5,43.7))



