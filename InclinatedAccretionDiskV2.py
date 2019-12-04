# -*- coding: utf-8 -*-
"""
Com a substituição de variável dependente u=1/r, resolvemos o sistema:
    u''=(2u^3-3b^2u^5+7Mb^2u^6)
    phi'=lu^2
introduzindo y=u', podemos reescrever este sistema como :
    u'=y
    y'=(2-(3u^2-7Mu^3)b^2)u^3
    phi'=bu^2

O programa calcula as trajetorias de raios de luz nas direções [0,pi/2] e [3/2pi,2pi],
no plano da imagem do observador(plano xy), e usa a simetria do problema para obter os outros dois quadrantes
O 2 quadrante é simetrico ao primeiro, e o terceiro é simetrico ao quarto

Para os raios vistos em um certo angulo theta(angulo polar imagem), alpha é o angulo entre a reta do
disco de acreção que este raio atravessa e o plano xy(plano de observação). Este angulo alpha
é utlizado nas transfer functions

Para [0,pi/2]: b>0 e usamos as transfer function TF1, TF2, TF3
Para [3/2pi, 2pi]: b<0 e usamos as transfer functions TF4, TF5, TF6

As informações de intensidade sao salvas na matriz ObsInt. O primeiro indice indica o angulo
de observação, e o segundo o raio

Os parametros que devem ser alterados sao Inclination, ListaParametros, Numangulos e a funçao SourceProfile
"""
from numpy import sqrt,arctan,pi,linspace,zeros,exp,meshgrid,tan,arccos,cos,concatenate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from time import time


M = 1 #Massa do buraco negro
d = 10*M #Distância do plano à origem
Inclination = -pi/4 #Inclinação do disco de acreção, em radianos

#define a função que representa o lado direito do sistema de eq diferenciais
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

#define evento de atravessar o horizonte de eventos
def EventHorizon(s, Y):
    return 1/Y[0] - 2*M
EventHorizon.terminal = True #define este evento como terminal

#Lista com os parametros de impacto que serão utilizados(também serão usados como os raios)
ListaParametros = linspace(0,10,50)
#Lista com os angulos (o numero total de angulos utilizados na formaçao da imagem é 4*NumAngulos-3)
NumAngulos = 15
Angulos = concatenate((linspace(0,pi/2,NumAngulos), linspace(pi,(3/2)*pi,NumAngulos)))
#Array onde será salvo a intensidade observada em cada ponto de observação
ObsInt = zeros((4*NumAngulos-3,len(ListaParametros)))

t0=time()
k=0 #indice radial que sera usado na matriz ObsInt

if Inclination<0: #este fator f indica se o angulode inclinaçao é positivo ou negativo
    f = -1
else:
    f=1

for b in ListaParametros:
    i=0 #indice angular que será usada na matriz ObsInt
    for theta in Angulos:
        if theta>pi:
            b = -b
        
        #angulo alpha que será usado nas transfer functions
        alpha = f*arccos(sqrt((1+(tan(theta)**2)*(cos(Inclination)**2))/(1+tan(theta)**2)))
        
        
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
        if len(sol.t_events[1])!=0: #TF1 ou TF4
            r=1/(sol.sol(sol.t_events[1][0])[0])
            if EH==0 or sol.t_events[1][0] < ParametroEH:
                ObsInt[i][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                #Agora usa a simetria para salvar os valores do segundo e terceiro quadrante
                if theta!=0 and theta!=pi/2 and theta<pi:
                    ObsInt[2*NumAngulos-i-2][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                elif theta>=pi and theta!=3/2*pi:
                    pos = i % (2*NumAngulos-2)
                    ObsInt[4*NumAngulos-4-pos][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
        if len(sol.t_events[2])!=0:
            r=1/(sol.sol(sol.t_events[2][0])[0])
            if EH==0 or sol.t_events[2][0] < ParametroEH:
                ObsInt[i][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                #Agora usa a simetria para salvar os valores do segundo e terceiro quadrante
                if theta!=0 and theta!=pi/2 and theta<pi:
                    ObsInt[2*NumAngulos-i-2][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                elif theta>=pi and theta!=3/2*pi:
                    pos = i % (2*NumAngulos-2)
                    ObsInt[4*NumAngulos-4-pos][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
        if len(sol.t_events[3])!=0:
            r=1/(sol.sol(sol.t_events[3][0])[0])
            if EH==0 or sol.t_events[3][0] < ParametroEH:
                ObsInt[i][k] += (sqrt(1-(2*M)/r)**4)*SourceProfile(r)
                #Agora usa a simetria para salvar os valores do segundo e terceiro quadrante
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
