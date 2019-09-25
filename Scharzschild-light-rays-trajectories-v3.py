# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:14:23 2019

@author: Andre
"""
from numpy import sqrt,arctan,linspace,cos,sin,pi
from scipy.integrate import solve_ivp
from matplotlib.pyplot import plot,show,xlim,ylim,axis,Circle,subplots,figure


b = 3*sqrt(3) #Parametro de impacto
M = 1 #Massa do buraco negro
d = 10 #Distância do plano à origem

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y , Y[2]=phi
def F(theta, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*Y[0]**3, b*Y[0]**2]

#valores iniciais
r0 = sqrt(b**2 + d**2)
u0 = 1/r0
phi0 = arctan(b/d)
y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))

#angulos em que a função será calculada


sol = solve_ivp(F, [0, 100], [u0, y0, phi0], max_step=0.01)

figure()
#plota um circulo que representa o buraco negro
circle1 = Circle((0, 0), 2*M, color='black')
fig ,ax = subplots()
#ax.set_aspect(1)
ax.add_patch(circle1)

'''
#plota a trajetoria circular de fotons
theta=linspace(0,2*pi,100)
plot(3*M*cos(theta),3*M*sin(theta), color='red', linestyle='dashed')
'''

#Plota os pontos da trajetoria
plot((1/sol.y[0])*cos(sol.y[2]), (1/sol.y[0])*sin(sol.y[2]))


axis("square")
xlim(-(b+3),b+3)
ylim(-(b+3),b+3)

show()

DeltaPhi = sol.y[2][len(sol.y[2])-1] - sol.y[2][0]
print("Variação angular:", DeltaPhi," rad")
print("n = ", DeltaPhi/(2*pi))

