#German Hernández
import numpy as np 
from P_Bur import PB
from P_roc import PRO
from scipy.optimize import fsolve 

def VI(Componentes,P,T,zi):
	y1=zi
	x1=zi
	z=np.asarray([zi,1-zi])
	Pr=PRO.f_pr(Componentes,T,y1)[0]
	Pb=PB.f_pb(Componentes,T,x1)[0]
	Psat=PB.Antoine(Componentes,T)
	if Pr<P<Pb:
		gammai_r=PRO.Margules(Componentes,x1)
		gammai_b=PB.Margules(Componentes,y1)
		gammas=[(((P-Pr)/(Pb-Pr))*(gammai_b[0]-gammai_r[0]))+gammai_r[0],(((P-Pr)/(Pb-Pr))*(gammai_b[1]-gammai_r[1]))+gammai_r[1]]		
		phi_r=PRO.PHI2V(Componentes,T,P,x1)
		phi_b=PB.PHI2V(Componentes,T,P,y1)
		#valores iniciales para comenzar a iterar
		phi=[(((P-Pr)/(Pb-Pr))*(phi_b[0]-phi_r[0]))+phi_r[0],(((P-Pr)/(Pb-Pr))*(phi_b[1]-phi_r[1]))+phi_r[1]]		
		Vi=(Pb-P)/(Pb-Pr)
		tol=1e-5
		er1=tol+1
		er2=tol+1
		er3=tol+1
		n=0
		while er1>=tol and er2>=tol and er3>=tol:  
			ki=np.asarray([(gammas[0]*Psat[0])/(phi[0]*P),(gammas[1]*Psat[1])/(phi[1]*P)])
			f=lambda V: np.sum((z*(ki-1)/(1+V*(ki-1))))
			#f=lambda V: ((zi*(ki[0]-1)/(1+V*(ki[0]-1)))+(1-zi)*(ki[1]-1)/(1+V*(ki[1]-1)))
			V=fsolve(f,Vi)
			Vi=V
			er1=np.sqrt(((V-Vi)**2))	
			L=1-V
			xi=[zi/(1+V*(ki[0]-1)),(1-zi)/(1+V*(ki[1]-1))]
			x1=xi[0]
			yi=[xi[0]*ki[0],xi[1]*ki[1]]
			y1=yi[0]
			gammai=[PRO.Margules(Componentes,x1),PB.Margules(Componentes,y1)]
			gammas=gammai
			er2=np.sqrt(((np.asarray(gammai)-np.asarray(gammas))**2))	
			phis=[PRO.PHI2V(Componentes,T,P,x1),PB.PHI2V(Componentes,T,P,y1)]
			phi=phis
			er3=np.sqrt(((np.asarray(phis)-np.asarray(phi))**2))	
			n+1

		print("El sistema se encuentra en dos fases")

	else:
		print("El sistema no se encuentra en dos fases")
	return V,L,xi,yi


print(VI(["ACETONA","METANOL"],479,105+273.15,0.65))


