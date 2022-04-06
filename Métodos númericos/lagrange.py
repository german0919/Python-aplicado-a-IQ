#German Hernández
from scipy import interpolate
#Interpolacion utilizando Lagrage
T=[26.67,93.33,148.89,315.56]
v=[1.35,0.085,0.012,0.00075]
#Datos a interpolar
x=[80,100,200,300]
f=interpolate.lagrange(T,v)
print(f(x))