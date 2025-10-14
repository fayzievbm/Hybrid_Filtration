# -*- coding: utf-8 -*-
"""
Created on Thu May 23 18:26:55 2019

@author: Komiljon
"""
import numpy as np
import matplotlib.pyplot as plt
#import copy
import modul25
#boshlang'ich ma'lumotlar
#beta=0.005

m0=0.35
m1=0.3

v=0.001
c00=0.01
c10=0.1

D=1e-5
D1=1e-5

r0=0.1
ra0=0.1
ra01=0.08
ra1=0.01
ra11=0.007

betar=20
betaa=32
betad=betaa*v*c00/ra01*0.8 #0.002
betar1=12
betaa1=22#32#22
betad1=betad#a1*v*c00/ra01 #0.002

nn=41
tm=1201 

#2700
#N=range(nn)
#print(N[0])
tmax=range(1,tm) #vaqt
tau=1
h=0.01
h1=c10*v/(1-m1)*tau
print(h1)

#print(h2)
#Procedure Nach_Dan;{Boshlang'ich va chegaraviy shartlar}
roa=np.zeros((nn,tm), dtype=float)
c=np.zeros((nn,tm), dtype=float)
c1=np.zeros((tm+nn,tm), dtype=float)
roa1=np.zeros((tm+nn,tm), dtype=float)
alpha=np.zeros(nn, dtype=float)
bet=np.zeros(nn, dtype=float)
alpha1=np.zeros(tm+nn, dtype=float)
bet1=np.zeros(tm+nn, dtype=float)

#c[nn-1,:]=c00

#  /// bettaa ni fisoblash funksiyasi
def bettaa(ra,cij):
   if  ra<=ra1: 
     return betar*v*cij
   elif ra<ra0:
     return betaa*v*cij-betad*ra
   else: 
     return 0.0
 
def bettaa1(ra,cij):
   if  ra<=ra11: 
     return betar1*v*cij
   elif ra<ra01:
     return betaa1*v*cij-betad1*ra
   else: 
     return 0.0

def progonka():
    #A=tau*D/(h*h)+tau*v/h
    #B=2*tau*D/(h*h)+tau*v/h+m0
    #E=tau*D/(h*h)
    A=tau*D/(h*h)
    B=2*tau*D/(h*h)+tau*v/h+m0
    E=tau*D/(h*h)+tau*v/h
    alpha[1]=1#0
    bet[1]=0#c00
    for i in range(1,nn-1):
        F=m0*c[i,j-1]-(roa[i,j]-roa[i,j-1])
        bet[i+1]=(F+A*bet[i])/(B-A*alpha[i])
        alpha[i+1]=E/(B-A*alpha[i])
    c[nn-1,j]=c00#(F+A*bet[nn-1])/(B-A*alpha[nn-1])
    for i in range(nn-2,-1,-1):
        c[i,j]=alpha[i+1]*c[i+1,j]+bet[i+1]



#Boshlandi 
for j in tmax:
     for i in range(nn):
        
         roa[i,j]=roa[i,j-1]+tau*bettaa(roa[i,j-1],c[i,j-1])
         if  roa[i,j]>ra0: 
            roa[i,j]=ra0
     #govak muhit uchun progonka
     A=tau*D/(h*h)
     B=2*tau*D/(h*h)+tau*v/h+m0
     E=tau*D/(h*h)+tau*v/h
     alpha[1]=1#0
     bet[1]=0#c00
     for i in range(1,nn-1):
         F=m0*c[i,j-1]-(roa[i,j]-roa[i,j-1])
         bet[i+1]=(F+A*bet[i])/(B-A*alpha[i])
         alpha[i+1]=E/(B-A*alpha[i])
    
     #keyk qatlam uchun progonka
     for i in range(nn-1,nn+j-1):
        
         roa1[i,j]=roa1[i,j-1]+tau*bettaa1(roa1[i,j-1],c1[i,j-1])
         if  roa1[i,j]>ra01: 
            roa1[i,j]=ra01
     A1=tau*D1/(h1*h1)
     B1=2*tau*D1/(h1*h1)+tau*v/h1+m1
     E1=tau*D1/(h1*h1)+tau*v/h1
     
     alpha1[nn-1]=(D1*h)/(D1*h+(1-alpha[nn-1])*D*h1)#0
     bet1[nn-1]=(bet[nn-1]*D*h1)/(D1*h+(1-alpha[nn-1])*D*h1)#c00
     
     for i in range(nn-1,nn+j-1):
         F1=m1*c1[i,j-1]-(roa1[i,j]-roa1[i,j-1])
         bet1[i+1]=(F1+A1*bet1[i])/(B1-A1*alpha1[i])
         alpha1[i+1]=E1/(B1-A1*alpha1[i])
     #keyk qatlam uchun progonka
     c1[nn-1+j,j]=c00
     if j==0:
         c1[nn-1+j,j]=c00
     elif j==1:
             c1[nn-1+j,j]=c00
             c1[nn-1,j]=c00
     else:    
       for i in range(nn+j-2,nn-2,-1):
         c1[i,j]=alpha1[i+1]*c1[i+1,j]+bet1[i+1]
     
     if j==0:
         c[nn-1,j]=c00#(F+A*bet[nn-1])/(B-A*alpha[nn-1])
     elif j==1:
         c[nn-1,j]=c00
     else:
         c[nn-1,j]=c1[nn-1,j]
     for i in range(nn-2,-1,-1):
         c[i,j]=alpha[i+1]*c[i+1,j]+bet[i+1]
     #c[0,j]=c[1,j]

#grafiklar
c=c/c00
c1=c1/c00
x = np.linspace(0, h*(nn-1), num=nn, endpoint=True)
x1 = np.linspace((nn-1)*h,(nn-1)*h+h1*(tm-1), num=(tm), endpoint=True)
modul25.grafikX(x,x1,c,c1,roa,roa1,tm)
#c=c*c00

print(c1[nn-1,:])
print(c[nn-1,:])
print(np.max(c1))

print(np.max(alpha1))
print(np.max(bet1))
x = np.linspace(0, h1*(tm), num=(tm), endpoint=True)
#plt.plot(x,c1[(nn-1):tm+(nn-1):1,(tm-1)//3])
#plt.plot(x,c1[(nn-1):tm+(nn-1):1,2*(tm-1)//3])
#lt.plot(x,c1[(nn-1):tm+(nn-1):1,3*(tm-1)//3])

#plt.show
