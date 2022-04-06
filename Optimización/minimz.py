#German Hernández
import numpy as np
from scipy.optimize import minimize
#Ejemplo de minimizar una funcion con varias condiciones
def objective(x):
	x1 = x[0]
	x2 = x[1]
	x3 = x[2]
	x4 = x[3]
	return x1*x4*(x1+x2+x3)+x3
def constraint1(x):
	return x[0]*x[1]*x[2]*x[3]-25.0
#Condicion sum(xi**2)=40
def constraint2(x):	
	sum_sq=40
	for i in range(4):
		sum_sq=sum_sq-x[i]**2
	return sum_sq	

x0=[1,5,5,1]
b=(1.0,5.0)
bnds=(b,b,b,b) #rango de las soluciones 1<=x1,x2,x3,x4<=5
#Identificar de que tipo son las condiciones
con1={'type':'ineq','fun':constraint1}
con2={'type':'eq','fun':constraint2}
cons=[con1,con2]
sol=minimize(objective,x0,method='SLSQP', bounds=bnds,constraints=cons)
print(sol)
print(sol.x)

'''
def objective(x):
	x1 = x[0]
	x2 = x[1]
	x3 = x[2]
	return 19.4*x1**-1.47+16.8*x2**-1.66+91.5*x3**-0.30
def constraint1(x):
	return x[2]
def constraint2(x):
	return x[1]
def constraint3(x):
	return x[0]-x[1]
def constraint4(x):
	return x[0]-x[2]	
def constraint5(x):
	return x[1]-x[2]	



x0=[1,0.5,0.3]
b=(0,0.05)
b2=(0,1)
bnds=(b2,b2,b) #rango de las soluciones 1<=x1,x2,x3,x4<=5
#Identificar de que tipo son las condiciones
con1={'type':'ineq','fun':constraint1}
con2={'type':'ineq','fun':constraint2}
con3={'type':'ineq','fun':constraint3}
con4={'type':'ineq','fun':constraint4}
con5={'type':'ineq','fun':constraint5}
cons=[con1,con2,con3,con4,con5]
sol=minimize(objective,x0,method='SLSQP', bounds=bnds,constraints=cons)
print(sol)
print(sol.x) '''

