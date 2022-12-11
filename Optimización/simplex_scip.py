#German Hernández
#Optimización lineal por el método Simplex usando scipy
from scipy.optimize import linprog
'''-------------------------------------------------------------------------
    Minimizar los costos de producir gasolina con las especificaciones dadas
    OBJ Min: 43*x1 + 31x2 + 47*x3 + 37*x4  donde xi:composición  
    Restricciones: 80x1+30x2+70x3+40x4>=20  ;   x1<=0.3   
                   10x1+30x2+10x3+50x4>=30  ;   x2<=0.4 ; x1+x2+x3+x4=1
                   10x1+40x2+20x3+10x4>=20  ; 0<=x3,x4<=1
-----------------------------------------------------------------------'''
c=[43,31,47,37] #Coef fobj
#Coef restricciones que son inecuaciones TODAS POR DEFECTO <= toca * (-)
A=[[-80,-30,-70,-40],[-10,-30,-10,-50],[-10,-40,-20,-10]]
#Vector b (lado derecho de las restricciones)
b=[-20,-30,-20]
#Restricciones que son ec ==
Aeq=[[1,1,1,1]]
beq=[1]
#Rangos (intervalos) para cada solución
x1_bns=(0,0.3) ; x2_bns=(0,0.4) ; x3_bns=(0,1) ; x4_bns=(0,1)
res=linprog(c,A,b,Aeq,beq,bounds=(x1_bns,x2_bns,x3_bns,x4_bns),method="simplex")
print(res)















'''
c=[-1,4] #coeficientes fobj Min: -1*x1+4*x2
#Restricciones 
# -3*x1+x2<=6
A=[[-3,1],[1,2]]
#x1+2*x2<=4
b=[6,4]
x1_bns=(None,None) #rango de las soluciones esta entre para x1 (-inf, inf)
# x2>=-3 
x2_bns=(-3,None) #rango de las soluciones de x2 (-3, inf)
res=linprog(c,A,b,bounds=(x1_bns,x2_bns),method="simplex")
print(res) '''

'''import numpy as np
from scipy.optimize import minimize
#Optimización sin restricciones usando minimize
x0 = np.array([0.5])
def fun(x):
    return 3*x[0]**2+(12/x[0]**3)-5
#res = minimize(fun, x0, method='nelder-mead',options={'xatol': 1e-8, 'disp': True})
res = minimize(fun, x0)
print(res)
'''

'''
c=[-4,-5] #coeficientes fobj Min: -1*x1+4*x2
#Restricciones 
# -3*x1+x2<=6
A=[[2,1]]
#x1+2*x2<=4
b=[8]
x1_bns=(0,None) #rango de las soluciones esta entre para x1 (-inf, inf)
# x2>=-3 
x2_bns=(0,5) #rango de las soluciones de x2 (-3, inf)
res=linprog(c,A,b,bounds=(x1_bns,x2_bns),method="simplex")
print(res) 
'''