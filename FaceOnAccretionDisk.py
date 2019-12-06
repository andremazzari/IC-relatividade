# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 14:06:47 2019

@author: Andre
"""
from numpy import sqrt,arctan,pi,linspace,zeros,exp,meshgrid,tile
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from time import time


M = 1 #Massa do buraco negro
d = 100*M #Distância do plano à origem

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y , Y[2]=phi
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*Y[0]**3, b*Y[0]**2]

#define as transfer functions 
def TF1(s, Y):
    return Y[2]-(pi/2)

def TF2(s, Y):
    return Y[2]-(3/2)*pi

def TF3(s, Y):
    return Y[2]-(5/2)*pi

#define o perfil de luminosidade
def SourceProfile(r):
    if r<6:
        x = 0
    else:
        x = 27.85*exp(-0.5545*r)
    return x

#define evento de passar o horizonte de eventos
def EventHorizon(s, Y):
    return 1/Y[0] - 2*M
EventHorizon.terminal = True #define este evento como terminal

#Lista com os parametros de impacto que serão utilizados
ListaParametros = linspace(0,15,300)
#Array onde será salvo a intensidade observada em cada parametro de impacto
ObsInt = zeros(len(ListaParametros))
EmittedInt = zeros(len(ListaParametros))

t0=time()
k=0
for b in ListaParametros:
    #calcula as condições inicias
    r0 = sqrt(b**2 + d**2)
    u0 = 1/r0
    phi0 = arctan(b/d)
    y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))
    
    #resolve o sistema
    sol = solve_ivp(F, [0, 150], [u0, y0, phi0], events=(EventHorizon, TF1, TF2, TF3), dense_output=True, max_step=0.01)
    
    #verifica as intersecções do raio de luz com o disco
    if len(sol.t_events[1])!=0:
        r=1/(sol.sol(sol.t_events[1][0])[0])
        if r<=2*M:
            print("m=1; r=",r,"; b=",b)
        else:
            ObsInt[k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
    if len(sol.t_events[2])!=0:
        r=1/(sol.sol(sol.t_events[2][0])[0])
        if r<=2*M:
            print("m=2; r=",r,"; b=",b)
        else:
            ObsInt[k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
    if len(sol.t_events[3])!=0:
        r=1/(sol.sol(sol.t_events[3][0])[0])
        if r<=2*M:
            print("m=3; r=",r,"; b=",b)
        else:
            ObsInt[k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
    
    EmittedInt[k] = SourceProfile(b)     
    k+=1

fig, ax = plt.subplots()
#PLota o perfil de emissão
ax.plot(ListaParametros, EmittedInt,'-b' , color='gray', label='Emitted profile')

#Plota o perfil de intensidade observado
ax.plot(ListaParametros, ObsInt,'-b' , color='blue', label='Observed profile')


ax.axis(xlim=(0,10),ylim=(0,10))
leg = ax.legend();
plt.show()

#Agora plota a imagem do buraco negro
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
azm = linspace(0, 2*pi, 100)
r, th = meshgrid(ListaParametros, azm)
z = tile(ObsInt, (r.shape[0], 1))

plt.pcolormesh(th, r, z)
plt.colorbar(label='Observed Intensity')


print(time()-t0)