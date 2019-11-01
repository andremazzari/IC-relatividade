# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 21:25:56 2019

@author: Andre
"""
from numpy import sqrt,arctan,pi,linspace,zeros,exp,meshgrid,tan,arccos,cos,concatenate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from time import time


M = 1 #Massa do buraco negro
d = 10*M #Distância do plano à origem
Inclination = 0.296705972 #Inclinação do disco de acreção, em radianos

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y , Y[2]=phi
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*Y[0]**3, b*Y[0]**2]

#define as transfer functions 

def TF1(s, Y):
    return Y[2] - (pi/2 - alpha)

def TF2(s, Y):
    return Y[2] - ((3/2)*pi - alpha)

def TF3(s, Y):
    return Y[2] - ((5/2)*pi - alpha)

def TF4(s, Y):
    return Y[2] + (pi/2 + alpha)

def TF5(s, Y):
    return Y[2] + ((3/2)*pi + alpha)

def TF6(s, Y):
    return Y[2] + ((5/2)*pi + alpha)

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

#Lista com os parametros de impacto que serão utilizados(também serão usados como os raios)
ListaParametros = linspace(0,15,200)
#Lista com os angulos
NumAngulos = 40
Angulos = concatenate((linspace(0,pi/2,NumAngulos), linspace(pi,(3/2)*pi,NumAngulos)))
#Array onde será salvo a intensidade observada em cada ponto de observação
ObsInt = zeros((4*NumAngulos-3,len(ListaParametros)))

t0=time()
k=0 #indice radial que sera usado na matriz ObsInt
for b in ListaParametros:
    i=0 #indice angular que será usada na matriz ObsInt
    for theta in Angulos:
        if theta>pi:
            b = -b
        
        #angulo alpha que será usado nas transfer functions
        alpha = arccos(sqrt((1+(tan(theta)**2)*(cos(Inclination)**2))/(1+tan(theta)**2)))
        
        #calcula as condições inicias
        r0 = sqrt(b**2 + d**2)
        u0 = 1/r0
        phi0 = arctan(b/d)
        y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))
        
        #resolve o sistema
        if b>=0:
            sol = solve_ivp(F, [0, 100], [u0, y0, phi0], events=(EventHorizon, TF1, TF2, TF3), dense_output=True, max_step=0.01)
        else:
            sol = solve_ivp(F, [0, 100], [u0, y0, phi0], events=(EventHorizon, TF4, TF5, TF6), dense_output=True, max_step=0.01)
        
        #Verifica se o raio de luz cruzou o horizonte de eventos
        if len(sol.t_events[0])!=0:
            EH = 1 #significa que o raio cruzou o horizonte de eventos
            ParametroEH = sol.t_events[0][0]
        else:
            EH = 0 #significa que o raio nao cruzou o horizonte de eventos
            ParametroEH = 0
        
        #Agora verifica as intersecções com o disco de acreção
        if len(sol.t_events[1])!=0:
            r=1/(sol.sol(sol.t_events[1][0])[0])
            if EH==0 or sol.t_events[1][0] < ParametroEH:
                ObsInt[i][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                if theta!=0 and theta!=pi/2 and theta<pi:
                    ObsInt[2*NumAngulos-i-2][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                elif theta>=pi and theta!=3/2*pi:
                    pos = i % (2*NumAngulos-2)
                    ObsInt[4*NumAngulos-4-pos][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
        if len(sol.t_events[2])!=0:
            r=1/(sol.sol(sol.t_events[2][0])[0])
            if EH==0 or sol.t_events[2][0] < ParametroEH:
                ObsInt[i][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                if theta!=0 and theta!=pi/2 and theta<pi:
                    ObsInt[2*NumAngulos-i-2][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                elif theta>=pi and theta!=3/2*pi:
                    pos = i % (2*NumAngulos-2)
                    ObsInt[4*NumAngulos-4-pos][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
        if len(sol.t_events[3])!=0:
            r=1/(sol.sol(sol.t_events[3][0])[0])
            if EH==0 or sol.t_events[3][0] < ParametroEH:
                ObsInt[i][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                if theta!=0 and theta!=pi/2 and theta<pi:
                    ObsInt[2*NumAngulos-i-2][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                elif theta>=pi and theta!=3/2*pi:
                    pos = i % (2*NumAngulos-2)
                    ObsInt[4*NumAngulos-4-pos][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
        if i!=(NumAngulos-1):        
            i+=1
        else:
            i=2*NumAngulos-2
    k+=1
    

Angulos= concatenate((linspace(0,pi/2,NumAngulos),linspace(pi/2,pi,NumAngulos)[1:NumAngulos-1],linspace(pi,3/2*pi,NumAngulos),linspace(3/2*pi,2*pi,NumAngulos)[1:NumAngulos]))

#Agora plota a imagem do buraco negro
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

r, th = meshgrid(ListaParametros, Angulos)

plt.pcolormesh(th, r, ObsInt)
plt.colorbar(label='Observed Intensity')


print(time()-t0)
