#German Hernández

import numpy as np 
from matplotlib import pyplot as plt
from tabulate import tabulate
###### Euler #######
def Euler(y0,xi,xf,n):
	func=lambda x,y: y+2*x-x**2
	#crear tanto x como el tamaño de paso en ellos h
	x,h=np.linspace(xi,xf,n, retstep=True) #retstep->obtener paso se debe aumentar n si quiero un tamaño de paso mas peq
	#Se crea un vector donde almacenar los resultados
	y=np.zeros(x.shape) #shape forma de x
	y[0]=y0 
	for n in range(0,len(x)-1):
		y[n+1]=y[n]+func(x[n],y[n])*h
	#fin del ciclo for
	plt.plot(x,y,"ro--", label="númerico")	
	plt.plot(x,x**2+np.exp(x),"bx--",label="analítico")
	plt.title("Método de Euler")
	plt.grid(True)
	plt.legend()
	return plt.show()

############# RK4 ###########
def RK4(y0,xi,xf,n):
	func=lambda x,y: y+2*x-x**2
	#crear tanto x como el tamaño de paso en ellos h
	x,h=np.linspace(xi,xf,n, retstep=True) #retstep->obtener paso se debe aumentar n si quiero un tamaño de paso mas peq
	#Se crea un vector donde almacenar los resultados
	y=np.zeros(x.shape) #shape forma de x
	y[0]=y0 
	for n in range(0,len(x)-1):
		k1=func(x[n],y[n])
		k2=func(x[n]+h/2,y[n]+h*(k1/2))
		k3=func(x[n]+h/2,y[n]+h*(k2/2))
		k4=func(x[n]+h,y[n]+h*k3)
		y[n+1]=y[n]+(h*(k1+2*k2+2*k3+k4)/6)
	#fin del ciclo for
	#plt.style.use('dark_background')
	plt.plot(x,y,"ro--", label="numerico")	
	plt.plot(x,x**2+np.exp(x),"bx--",label="analitico")
	plt.title("Método de RK4")
	plt.grid(True)
	plt.legend()
	Datos=np.transpose([x,y])
	Encabezados=["x","y"]
	print(tabulate(Datos,headers=Encabezados,tablefmt="fancy_grid",stralign="center"))							
	return plt.show()

#print(Euler(1,0,1.5,30))
print(RK4(1,0,1.5,30))

'''
def RK4(fun,y0,xi,xf,n):
	#crear tanto x como el tamaño de paso en ellos h
	x,h=np.linspace(xi,xf,n, retstep=True) #retstep->obtener paso se debe aumentar n si quiero un tamaño de paso mas peq
	#Se crea un vector donde almacenar los resultados
	y=np.zeros(x.shape) #shape forma de x
	y[0]=y0 
	for n in range(0,len(x)-1):
		k1=func(x[n],y[n])
		k2=func(x[n]+h/2,y[n]+h*(k1/2))
		k3=func(x[n]+h/2,y[n]+h*(k2/2))
		k4=func(x[n]+h,y[n]+h*k3)
		y[n+1]=y[n]+(h*(k1+2*k2+2*k3+k4)/6)
	#fin del ciclo for
	return y
'''