#Código en Python 3.7 para la determinación de la solubilidad mediante el método UNIFAC 
import pandas as pd
import numpy 
import csv
print("Solubilidad mediante UNIFAC por David Rivera Ramírez")
#	Introducción de las constantes desde una base de datos v, R, Q y a
print(""" Este programa predice la solubilidad en fracción mol de un soluto frente a un disolvente. Tomando en cuenta lo anterior y de la lista mostrada a continuación, se considera que el índice 1 está reservado para el disolvente y el resto a un conjunto de solutos que interaccionan 1 a 1 con el disolvente, por lo cual debes digitar un número mayor a 1.""")
compuestos=pd.read_csv('BBDDCompuestosColillasA.csv') 
print(compuestos.iloc[1:,0:1])
ielegido=int(input("Seleccione hasta que número será el rango del cálculo de la mezcla de compuestos""\n"))  
indice=ielegido+1
listC0=numpy.array(compuestos.iloc[1:indice,:1])
matdis=numpy.array(compuestos.iloc[1:2,:1])
Disolvente=str(matdis).lstrip('[').rstrip(']')
#	Parámetros del disolvente
df0d=compuestos.iloc[1,1:109]
df0xd=df0d.values 
paradis=numpy.matrix(df0xd,dtype=float)
#	Parámetros de los solutos
df0s=compuestos.iloc[2:indice,1:109]
df0xs=df0s.values 
parasol=numpy.matrix(df0xs,dtype=float)
listsolutos=numpy.delete(listC0,0,0)
#	Entrada de datos de área y volumen 
arevol=pd.read_csv('Area y Volumen.csv') 
#	Obtención de la matriz R 
df1=arevol.iloc[0:108,2:3] 
df1x=df1.values 
Ri_1=numpy.matrix(df1x,dtype=float)
#	Obtención de la matriz Q 
df2=arevol.iloc[0:108,3:4] 
df2x=df2.values 
Qi_1=numpy.matrix(df2x,dtype=float) 
#	Obtención de matrices ri
ri_d=numpy.mat(paradis) * numpy.mat(Ri_1)
ri_s=numpy.mat(parasol) * numpy.mat(Ri_1)
#	Obtención de matrices qi
qi_d=numpy.mat(paradis) * numpy.mat(Qi_1)
qi_s=numpy.mat(parasol) * numpy.mat(Qi_1)
#	Obtención de la matriz Gk=Vk*Qi (Área parcial por grupo estructural)
qi2=numpy.reshape(Qi_1,(1,108))
Gk_d=numpy.multiply(paradis,qi2)
Gk_s=numpy.multiply(parasol,qi2)
#	Definiendo la temperatura para obtener la matriz "τ" 
Temp=float(input("Introduce la temperatura de la mezcla en grados Kelvin ""\n"))
#	Cálculo de la energía funcional de la mezcla
interac=pd.read_csv('interacciones.csv')
dfi=interac.iloc[1:109,2:110]
dfix=dfi.values 
itr=numpy.matrix(dfix,dtype=float)
itrx=itr*(-1)
it=itrx/Temp
thau=numpy.exp(it)
#	Energía "parcial molar" por especie (S)
Sk_d=numpy.mat(Gk_d) * numpy.mat(thau)
Sk_s=numpy.mat(Gk_s) * numpy.mat(thau)
#	Entalpias de fusión para cálculo de solubilidad
dfe=compuestos.iloc[2:indice,109:110]
dfent=dfe.values 
entalpias=numpy.matrix(dfent,dtype=float)
#	Temperaturas de fusión para cálculo de solubilidad
dfe=compuestos.iloc[2:indice,110:111]
dfTx=dfe.values 
Tempfusion=numpy.matrix(dfTx,dtype=float)
#	Cálculo de las solubilidades ideales
ener=numpy.multiply((-entalpias/(0.008314472*Tempfusion)),((Tempfusion/Temp)-1))
solub_1=numpy.exp(ener)
print("Las solubilidades ideales de referencia para cada soluto son""\n",solub_1)
indicesolutos=ielegido-1
totsln=numpy.ones([(indicesolutos),1], dtype=float)
disolvente=totsln-solub_1
a=[]
er=[]
#Inicia el ciclo determinado para el número de compuestos i
for i in range(0,indicesolutos):
	no_iteraciones=0
			
	while solub_1[i,0]>=0:
# Se obtiene el arreglo ri e qi, añadiendo los valores de r y q tanto del compuesto i, #como del disolvente
		ri=numpy.vstack((ri_d,ri_s[i,0]))
		qi=numpy.vstack((qi_d,qi_s[i,0]))
		fracciniciales=numpy.vstack((disolvente[i,0],solub_1[i,0]))
		print("Las fracciones iniciales en el bucle son""\n",fracciniciales) 
		fraccresh=numpy.reshape(fracciniciales,(1,2))
		
		#Calculando Ji
		c0xri=numpy.mat(fraccresh) * numpy.mat(ri)
		Ji=numpy.mat(ri)/numpy.mat(c0xri)
				
		#Calculando Li
		c0_xqi=numpy.mat(fraccresh) * numpy.mat(qi)
		Li=numpy.mat(qi)/numpy.mat(c0_xqi)
		
        		#Calculando Gk		
		Gk=numpy.vstack((Gk_d,Gk_s[i,]))
		#Calculando teta
		teta=numpy.mat(fraccresh) * numpy.mat(Gk)
		#Calculando Sk
		Sk=numpy.vstack((Sk_d,Sk_s[i,]))
		#Calculando nuk
		nuk=numpy.mat(fraccresh) * numpy.mat(Sk)
		
		# Calculando LnRi
		a1=numpy.log(Sk/nuk)
		a2=numpy.multiply(Gk,a1)
		a22=numpy.array(a2)
		a3=numpy.sum(a22,axis=1, keepdims=True)
		
		b1=(Sk/nuk)
		b2=numpy.multiply(teta,b1)
		b22=numpy.array(b2)
		b3=numpy.sum(b22,axis=1, keepdims=True)
		
		c1=numpy.log(Li)
		c2=1-c1
		c21=numpy.multiply(qi,c2)
		c22=numpy.array(c21)
		c3=numpy.sum(c22,axis=1, keepdims=True)
		
		LnRi=c3-(b3-a3)
		#Calculando LnCi
		d1=numpy.log(Ji/Li)
		d2=Ji/Li
		d3=1-d2+d1
		d4=5*numpy.multiply(qi,d3)
		d5=numpy.log(Ji)
		LnCi=1-Ji+d5-d4
				
		e=LnRi+LnCi
		alfa=numpy.exp(e)
		alfasoluto=numpy.delete(alfa,0,0)
		print("El coeficiente de actividad del soluto es""\n",alfasoluto)
		xid=numpy.multiply((entalpias[i,0]/(0.008314472*Tempfusion[i,0])),((Tempfusion[i,0]/Temp)-1))		
		xideal=numpy.exp(-xid)
		fracnueva=xideal/alfasoluto
		print("La fraccion nueva calculada es ""\n",fracnueva)
		
		errmol=(solub_1[i,0]-fracnueva)/solub_1[i,0]				
		print("El error es  ""\n",errmol)
		
		if errmol<=0.001 and fracnueva<=1 :
			print("La fraccion mol definitiva es ""\n",fracnueva)
			print("Numero de iteraciones hechas fue " ,no_iteraciones)
			er.append(errmol)
			a.append(fracnueva)
			print(a)
			break;
		elif errmol<=0.001 and fracnueva>=1 :
			print("La fraccion mol definitiva es ""\n",solub_1[i,0])
			print("Numero de iteraciones hechas fue " ,no_iteraciones)
			er.append(errmol)
			a.append(solub_1[i,0])
			print(a)
			break;
		else:
			solub_1[i,0]=fracnueva
			disolvente[i,0]=1-fracnueva	
			no_iteraciones=no_iteraciones+1
		
res1=numpy.array(a)
res2=res1.flatten()
res3=numpy.reshape(res2,(indicesolutos,1))

resul = numpy.append([listsolutos], [res3],axis=2)
resultados=numpy.reshape(resul,(indicesolutos,2))

erro1=numpy.array(er)
erro2=erro1.flatten()
erro3=numpy.reshape(erro2,(indicesolutos,1))

extrares=numpy.append(resultados, erro3,axis=1)
print("Los resultados fueron""\n",resultados)
print("El numero de iteraciones fué :   ",no_iteraciones)

encabezado= ('  Disolvente:', Disolvente,'  Temperatura ',Temp)
subencabezados = (' Solutos', ' Fraccion','  % de Error')
dfresultados=pd.DataFrame(listsolutos,res3)
#Generando un archivo para resguardar los resultados
with open ('ResultadosUNIFAC.csv','a',newline='', encoding='utf-8') as archivocsv:
	escritor = csv.writer(archivocsv)
	escritor.writerow(encabezado)
	escritor.writerow(subencabezados)
	escritor.writerows(extrares)

