# -*- coding: utf-8 -*-
"""
Simula a imagem que observador no infinito ve de um buraco negro iluminado por um plano infinito atras dele.
Usamos a propriedade de que podemos inverter as trajetorias dos raios de luz,
e variando o parametro de impacto, calculamos se um dado raio de luz veio do plano de iluminação ou nao.
Usamos a simetria radial do problema para plotar a imagem da sombra do buraco negro que o observador ve.
Para uma lista com 1000 parametros b, o programa esta rodando em meia hora.
O Raio da sombra obtido foi de 6.087205523728037, mas segundo o artigo do Wald deveria ser 6.17.
O anel vermelhor interior tem raio de 5.2095154249579485, e o valor no artigo do Wald é de 5.20
"""

from numpy import sqrt,arctan,linspace,pi,cos,sin,zeros,meshgrid,tile
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from time import time

M = 1 #Massa do buraco negro
d = 100*M #Distância do plano à origem

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y , Y[2]=phi
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*Y[0]**3, b*Y[0]**2]

def HorizonteEventos(s, Y):
    return (1/Y[0])-2*M

#Lista com os parametros que serão calculados
ListaParametros = linspace(0,10*sqrt(2),1000)

t0=time()
theta = linspace(0,2*pi,100)
RaioSombra = 0
AnelFotons = []
ObsInt = zeros(len(ListaParametros))
k=0
for b in ListaParametros:
    #calcula as condições inicias
    r0 = sqrt(b**2 + d**2)
    u0 = 1/r0
    phi0 = arctan(b/d)
    y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))
    
    #resolve o sistema
    sol = solve_ivp(F, [0, 200], [u0, y0, phi0], events=HorizonteEventos,max_step=0.01)
    
    PhiFinal = sol.y[2][len(sol.y[2])-1]
    n=PhiFinal % (2*pi)
    #verifica se o raio de luz atravessou o horizonte de eventos e a variação angular
    if len(sol.t_events[0]) == 0 and n > (pi/2) and n < (3/2)*pi:
        #plot(b*cos(theta), b*sin(theta), color='red')
        ObsInt[k] = 1
        if b < 6.08:
            AnelFotons.append(b)
    else:
        #plot(b*cos(theta), b*sin(theta), color='black')
        RaioSombra = b
    
    if b<6.18 and b>6.07:
        print(b, " ; ",PhiFinal," ; ",n," ; ", len(sol.t_events[0]))
        
    k+=1

'''
axis("square")
xlim(-10,10)
ylim(-10,10)
show()
'''

#Agora plota a imagem do buraco negro
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
azm = linspace(0, 2*pi, 100)
r, th = meshgrid(ListaParametros, azm)
z = tile(ObsInt, (r.shape[0], 1))

plt.pcolormesh(th, r, z, cmap='Set1')
plt.colorbar(label='Observed Intensity')


print("raio da sombra:", RaioSombra)
print("anel de fotons:", AnelFotons)
print(time()-t0)