#CAIQ
#Programa para resolver sistema 
#de ecuaciones multivaribale por metodo de
#N-R
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
def NR_multi(V_i,X,v_F,error):
	'''V_i: Vector con valores iniciales.
	   X: Tupla con simbolos a utilizar, cada simbolo 
		  se debe definir entre comillas.
	   v_F: Vector que contiene funciones que definen 
			el sistema de ecuaciones no lineales a resolver,
			cada funcion debe definirse entre comillas.'''
	#Se procede a convertir simbolos str a simbolos de SymPy
	v_Var=sp.symbols(X)
	#Se convierten funciones de str a funciones operables por SymPy
	#Se utiliza una lista comprimida
	Func=[parse_expr(v_F[k]) for k in list(range(0,len(v_F)))] 
	#Se define el Jacobiano
	m_J=[]# Matriz donde se guardaran derivadas
	#Inicia ciclo for para crear Jacobiano y convertir 
	#lista de funciones a columna
	c_F=[]
	for k in list(range(0,len(v_F))):
		v_Deriv=[]
		funcion=Func[k]
		#Ciclo for anidado para ingresar derivadas por filas
		for m in list(range(0,len(v_Var))):
			d_F=sp.diff(funcion,v_Var[m])
			v_Deriv.append(d_F)
		#Fin de for anidado
		c_F.append([Func[k]])
		m_J.append(v_Deriv)
	#Fin de ciclo for
	#Se convierten funciones a expresiones evaluables, que reciben
	#la lista de valores iniciales
	Func_eval=sp.lambdify([v_Var],sp.Matrix(c_F))
	mJ_eval=sp.lambdify([v_Var],sp.Matrix(m_J))
	#Se definen condiciones de iteracion
	V_ic=np.asarray(V_i).reshape(-1,1)#Se transforma vector a matriz
	#con una columna
	Xi=V_ic-np.dot(np.linalg.inv(mJ_eval(V_i)),Func_eval(V_i))
	sum_e=np.sum(np.abs(Xi-V_ic))
	#Inicia ciclo while para solucionar problema
	n=0
	while (sum_e>error):
		n=n+1
		V_i=Xi
		V_if=list(V_i.reshape(1,-1)[0])
		Xi=V_i-np.dot(np.linalg.inv(mJ_eval(V_if)),Func_eval(V_if))
		sum_e=np.sqrt(np.sum((Xi-V_i)**2))
	#Fin de while
	return Xi,sum_e,n



'''print(
	NR_multi(
		[1,1,1],("x1","x2","x3"),
		["3*x1-cos(x2*x3)-0.5","x1**2-625*x2**2"
		 ,"exp(-x1*x2)+20*x3+(10*pi-3)/3"],0.00001))'''
"""print(
	NR_multi(
		[1,1,1,1],("x1","x2","x3","x4"),
		["3*x1-cos(x2*x3)-0.5*x4","x1**2-625*x2**2+x4**2"
		 ,"exp(-x1*x2)+20*x3-15*x4+(10*pi-3)/3","x4**2+5*x3"],0.00001))"""

'''		
print(
	NR_multi(
		[0.011808,0.034412],("x","y"),
		["((0.333+x)*(0.167+3*x+4*y)**3/(0.167-x-y)*(0.167-x-2*y))-1.30","((0.167+y)*(0.167+3*x+4*y)**4/(0.167-x-y)*(0.167-x-2*y)**2)-2.99"
		 ],0.00001))		 '''
'''
print(
	NR_multi(
		[1,1],("Xm","Nm"),
		["0.825*Xm**3-1.255*Xm**2-0.078*Xm+2-Nm"
		 ,"5.348*Xm-Nm"],0.00001))'''


		 