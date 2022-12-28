#German Hernández
#Optimización lineal:SISTEMA BOILER/TURBO GENERADOR
import numpy as np
from pyomo.environ import *
#Creamos el modelo
model=ConcreteModel()
#Agregamos las variables al modelo
model.HPS=Var(domain=NonNegativeReals); model.PP=Var(domain=NonNegativeReals);model.EP=Var(domain=NonNegativeReals)
model.P1=Var(domain=NonNegativeReals); model.HE1=Var(domain=NonNegativeReals);model.C=Var(domain=NonNegativeReals)
model.I1=Var(domain=NonNegativeReals); model.I2=Var(domain=NonNegativeReals);model.LE2=Var(domain=NonNegativeReals)
model.P2=Var(domain=NonNegativeReals); model.BF1=Var(domain=NonNegativeReals);model.BF2=Var(domain=NonNegativeReals)
model.MPS=Var(domain=NonNegativeReals); model.LPS=Var(domain=NonNegativeReals);model.HE2=Var(domain=NonNegativeReals)
model.LE1=Var(domain=NonNegativeReals) 
#*Función objetivo
model.cost=Objective(expr = 0.00261*model. HPS + 0.0239*model.PP + 0.00983*model.EP, sense= minimize)
#*Restricciones
#Turbina 1
model.c1=Constraint(expr = model.P1 <= 6250); model.c2=Constraint(expr = model.P1 >= 2500)
model.c3=Constraint(expr = model.HE1 <= 192000); model.c4=Constraint(expr = model.C <= 62000)
model.c5=Constraint(expr = model.I1 - model.HE1 <= 132000)
#Turbina 2
model.c6=Constraint(expr = model.P1 <= 9000); model.c7=Constraint(expr = model.P2 >= 3000)
model.c8=Constraint(expr = model.I2 <= 244000) ;model.c9=Constraint(expr = model.LE2 <= 142000)
#Balances de materia
model.c10=Constraint(expr = 0 == model.HPS - model.I1 - model.I2 - model.BF1 ) 
model.c11=Constraint(expr = 0 == model.I1 + model.I2 + model.BF1 - model.C - model.MPS - model.LPS)
model.c12=Constraint(expr = 0 == model.I1 - model.HE1 - model.LE1 - model.C ) 
model.c13=Constraint(expr = 0 == model.I2 - model.HE2 - model.LE2)
model.c14=Constraint(expr = 0 == model.HE1 + model.HE2 + model.BF1 - model.BF2 - model.MPS )
model.c15=Constraint(expr = 0 == model.LE1 + model.LE2 + model.BF2 - model.LPS )
#Energía comprada
model.c16=Constraint(expr = model.EP + model.PP >= 12000)
#Demanda
model.c17=Constraint(expr = model.MPS >= 271536); model.c18=Constraint(expr = model.LPS >= 100623)
model.c19=Constraint(expr = model.P1 + model.P2 + model.PP >= 24550)
#Balances de energía
model.c20=Constraint(expr = 0 == 1359.8*model.I1 - 1267.8*model.HE1 - 1251.4*model.LE1 - 192*model.C - 3413*model.P1) 
model.c21=Constraint(expr = 0 == 1359.8*model.I2 - 1267.8*model.HE2 - 1251.4*model.LE2 - 3413*model.P2) 
#*SOLVER
sol=SolverFactory("cbc")
sol.solve(model)
#print(model.display())
print("HPS:", round(model.HPS(),3), "[]")
print("PP:", round(model.PP(),3), "[]")
print("EP:", round(model.EP(),3), "[]")
print("El costo mínimo es:", round(model.cost(),3), "[$/h]")
