#German Hernández
from PFR import CPFR
import time  

print("---------------------------------------------------------")
print("Programa para calcular las dimensiones de un reactor PFR")
print("---------------------------------------------------------")
#------------------------------------------------------------------------------
print("¿Que desea graficar?, Ingrese si desea graficar volumen(V) o longitud y diametro(LD)")
Datos=input() #upper convertir a mayuscula.
if Datos=="V":
    respuesta=input("¿Tiene datos de conversión y velocidad de reacción?, escriba si :")
    if respuesta=="si":
        X=[]
        while True:
            element1=float(input("Ingrese los valores de X :"))
            X.append(element1)
            choice1=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice1 == "ok":
                break
        r_A=[]
        while True:
            element2=float(input("Ingrese los valores de rA :"))
            r_A.append(element2)
            choice2=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice2 == "ok":
                break        
        print("Ingrese el flujo molar inicial FA0 :")        
        FA_0=float(input())
        print(CPFR.GVR(X,r_A,FA_0))
    else:    
        t=[]
        while True:
            element1=float(input("Ingrese los valores de t :"))
            t.append(element1)
            choice1=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice1 == "ok":
                break
        CA=[]
        while True:
            element2=float(input("Ingrese los valores de CA :"))
            CA.append(element2)
            choice2=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice2 == "ok":
                break        
        print("Ingrese el concentracion  inicial CA0 :")        
        CA_0=float(input())                
        print("Ingrese el flujo molar inicial FA0 :")        
        FA_0=float(input())
        print(CPFR.GVC(t,CA,CA_0,FA_0))        
else:
    print("¿Tiene datos de conversión y velocidad de reacción?, escriba si :")
    respuesta=input()
    if respuesta=="si":
        X=[]
        while True:
            element1=float(input("Ingrese los valores de X :"))
            X.append(element1)
            choice1=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice1 == "ok":
                break
        r_A=[]
        while True:
            element2=float(input("Ingrese los valores de rA :"))
            r_A.append(element2)
            choice2=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice2 == "ok":
                break        
        print("Ingrese el flujo molar inicial FA0 :")        
        FA_0=float(input())
        print(CPFR.GLDR(X,r_A,FA_0))
    else:    
        t=[]
        while True:
            element1=float(input("Ingrese los valores de t :"))
            t.append(element1)
            choice1=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice1 == "ok":
                break
        CA=[]
        while True:
            element2=float(input("Ingrese los valores de CA :"))
            CA.append(element2)
            choice2=input("¿Quiere detener el ingreso de datos?, escriba ok, presione enter si desea continuar :")
            if choice2 == "ok":
                break        
        print("Ingrese el concentracion  inicial CA0 :")        
        CA_0=float(input())                    
        print("Ingrese el flujo molar inicial FA0 :")        
        FA_0=float(input())
        print(CPFR.GLDC(t,CA,CA_0,FA_0))            

time.sleep(6000)
