# -*- coding: utf-8 -*-
"""
Calcula e plota a trajetoria de raios de luz emitidos no mesmo ponto com diferentes parametros de impacto b.
O calculo é feito da mesma forma que nos outros programas, apenas alterando as condições iniciais.
"""

from numpy import sqrt,cos,sin
from scipy.integrate import solve_ivp
from matplotlib.pyplot import plot,show,xlim,ylim,axis,Circle,subplots,figure
from time import time

M = 1 #Massa do buraco negro
d = 10*M #Distância do plano à origem

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y , Y[2]=phi
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*Y[0]**3, b*Y[0]**2]

figure()
#plota um circulo que representa o buraco negro
circle1 = Circle((0, 0), 2*M, color='black')
fig ,ax = subplots()
#ax.set_aspect(1)
ax.add_patch(circle1)

#Lista com os parametros que serão calculados
ListaParametros = [1,2,3,4,5,6,7,7.348,3*sqrt(3)]
t0=time()

#condições iniciais fixas
r0 = 6*M #deve ser maior que 2*M
u0 = 1/r0
phi0 = 0
for b in ListaParametros:
    #Condição inicial faltante
    y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))
    
    #resolve o sistema
    sol = solve_ivp(F, [0, 100], [u0, y0, phi0], max_step=0.01)
    
    #plota a trajetória deste raio de luz
    plot((1/sol.y[0])*cos(sol.y[2]), (1/sol.y[0])*sin(sol.y[2]))
        
axis("square")
xlim(-10,10)
ylim(-10,10)

show()

print(time()-t0)