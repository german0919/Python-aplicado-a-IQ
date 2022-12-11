import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
#German Hernández
'''                   OPTIMIZACIÓN NO LINEAL 
Estimación de parámetros de isotermas de sorción con el modelo BET  
utilizando un algoritmo de optimización minimizando la suma del error.
Modelo: q_eq=(qs*kla*Ceq)/((1-klb*Ceq)*(1-klb*Ceq+kla*Ceq))
'''
def BET_opti(v_i):
    #Matriz con lo datos experimentales
    M=np.asarray([[0.17089701,0.085324916],
                [0.066582226,0.069637725]])
    Ceq=M[:,0]
    qexp=M[:,1]
    v_i=np.asarray(v_i)
    #Parámetros que voy a optimizar
    kla=v_i[0]
    klb=v_i[1]
    qs=v_i[2]
    #Modelo de isoterma BET
    q_eq=(qs*kla*Ceq)/((1-klb*Ceq)*(1-klb*Ceq+kla*Ceq))
    n=len(v_i)
    #Vector para almacenar la diferencia **2
    Sres=np.zeros(n)
    for i in list(range(0,n-1)):
        Sres[i]=(qexp[i]-q_eq[i])**2  
    #Vamos a calcular la suma del error    
    f=np.sum(Sres)   
    return f

def G_opti(kla,klb,qs,Ceqexp,qexp):
    #Con el modelo BET
    Ceq=np.linspace(0,0.26)
    q_eq=(qs*kla*Ceq)/((1-klb*Ceq)*(1-klb*Ceq+kla*Ceq))
    plt.plot(Ceq,q_eq,"b--")
    #Experimental
    plt.plot(Ceqexp,qexp,"o",color="red")
    plt.xlabel("Concentración del aceite\n en la fase líquida en el equilibrio Ce")
    plt.ylabel("Concentración del soluto\n en la fase sólida en el equilibrio qe")
    plt.grid(True)
    return plt.show()
    
v_i=[30,1e-3,0] #valores supuestos
Opt=minimize(BET_opti,v_i,method="Nelder-Mead")        
print(Opt)    
G_opti(34.64,-8.51e-3,9.99e-2,[0.17089701,0.066582226],[0.085324916,0.069637725])    

