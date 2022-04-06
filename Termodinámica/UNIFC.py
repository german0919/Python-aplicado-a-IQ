

'''#Para ejecutar el programa que calcule gammai con el metodo UNIFAC
from thermo.unifac import UFIP, UFSG, UNIFAC
GE = UNIFAC.from_subgroups(chemgroups=[{1:2, 2:4}, {1:1, 2:1, 18:1}], T=60+273.15, xs=[0.5, 0.5], version=0, interaction_data=UFIP, subgroups=UFSG)
print(GE.gammas())'''

'''#Para buscar la key de cada componente de la mezcla en la base de datos
from chemicals import search_chemical
import thermo.unifac
thermo.unifac.load_group_assignments_DDBST()
A=search_chemical('diethylamine').InChI_key
B=search_chemical('benzene').InChI_key
ka=thermo.unifac.DDBST_UNIFAC_assignments[A]
kb=thermo.unifac.DDBST_UNIFAC_assignments[B]
print("Componente1:",ka)
print("Componente2:",kb)'''


#-----------------------------------------------------------------------------------
#Ejemplo del libro
#heptano={1: 2, 2: 5}
#dietilamina={1: 2, 2: 1, 32: 1}

'''from thermo.unifac import UFIP, UFSG, UNIFAC
GE = UNIFAC.from_subgroups(chemgroups=[{1: 2, 2: 5}, {1: 2, 2: 1, 32: 1}], T=308.15, xs=[0.4, 0.6], version=0, interaction_data=UFIP, subgroups=UFSG)
print(GE.gammas())'''

'''#--------------------------------------------------------------------------------------
#Ejemplo del parcial 
#La desviacion es pequeña puede ser por ser de diferentes fuentes
from chemicals import search_chemical
import thermo.unifac
thermo.unifac.load_group_assignments_DDBST()
A=search_chemical('ethylamine').InChI_key
B=search_chemical('diethylamine').InChI_key
ka=thermo.unifac.DDBST_MODIFIED_UNIFAC_assignments[A]
kb=thermo.unifac.DDBST_MODIFIED_UNIFAC_assignments[B]
print("Componente1:",ka)
print("Componente2:",kb)

from thermo.unifac import UFIP, UFSG, UNIFAC
GE = UNIFAC.from_subgroups(chemgroups=[ka, kb], T=385, xs=[0.25, 1-0.25], version=0, interaction_data=UFIP, subgroups=UFSG)
print(GE.gammas())'''
#------------------------------------------------------------------------------

#Otra forma es aprovechar la función
'''Calcula los coeficientes de actividad utilizando el modelo UNIFAC (opcionalmente modificado),
 dada la temperatura de una mezcla, las fracciones molares líquidas y, opcionalmente,
 los datos del subgrupo y los datos de los parámetros de interacción de su elección. 
 El valor predeterminado es utilizar el modelo UNIFAC original, con los últimos parámetros publicados por DDBST.
  El modelo admite formas modificadas (Dortmund, NIST) cuando el parámetro modificado es Verdadero.   '''

#from thermo import UNIFAC_gammas  
#UNIFAC_gammas(T=308.15, xs=[0.4, 0.6], chemgroups=[{1: 2, 2: 5}, {1: 2, 2: 1, 32: 1}])
#------------------------------------------------------------------------------------------------


#Ejemplo del parcial pero ahora con una funcion para generalizar aún más.
from chemicals import search_chemical
import thermo.unifac
from thermo.unifac import UFIP, UFSG, UNIFAC
#UFIP parametros de interacción para el modelo UNIFAC original 
#UFSG o UFMG son los datos para el modelo original de UNIFAC

def Ufac(Componente1,Componente2,Ti,xs1):
	#Se llama la función los datos se almacenan en el formato InChI key bool bool bool recuento de subgrupos ... recuento de subgrupos recuento de subgrupos ... donde los bools se refieren a si las asignaciones originales de UNIFAC, UNIFAC modificado y PSRK se completaron correctamente o no. Los subgrupos y su recuento tienen una longitud indefinida.
	thermo.unifac.load_group_assignments_DDBST()
	A=search_chemical(Componente1).InChI_key
	B=search_chemical(Componente2).InChI_key
	ka=thermo.unifac.DDBST_UNIFAC_assignments[A]
	kb=thermo.unifac.DDBST_UNIFAC_assignments[B]
	GE = UNIFAC.from_subgroups(chemgroups=[ka, kb], T=Ti, xs=[xs1, 1-xs1], version=0, interaction_data=UFIP, subgroups=UFSG)
	gamma_i=GE.gammas()
	return gamma_i
print(Ufac("18172-67-3","5392-40-5",300,0.5))
#Es mejor usar el CAS que la formula para no confundir con isomeros
#α-Caryophyllene->C15H24
#Geranial->C10H16O 141-27-5
#β-Bisabolene->	C10H16O CAS:76-22-2
#Acetic acid nonyl ester->C11H22O2	CAS:143-13-5 
#Limoneno->138-86-3
#Vainillina->121-33-5
#(S)-(+)-Carvona CAS:2244-16-8
#β-mirceno CAS:123-35-3
#(1S)-(-)-α-Pineno CAS 7785-26-4
#ß-Pineno CAS:18172-67-3
#(S) -β-bisaboleno(CAS:495-61-4)
#(+)-Canfeno (CAS 79-92-5)
#Citral (CAS No. 5392-40-5)

#varias opciones según el metodo	
#thermo.unifac.DDBST_MODIFIED_UNIFAC_assignments
#thermo.unifac.DDBST_UNIFAC_assignments
#thermo.unifac.DDBST_PSRK_assignments