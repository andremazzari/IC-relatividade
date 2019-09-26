# -*- coding: utf-8 -*-
"""
Calcula e plota a trajetório de um raio de luz com parametro de impacto b no espaço-tempo de schwarzschild.
Os raios de luz iniciam de forma perpendicular a um plano no infinito.
A altura h que o raio inicia esta relacionada com o parametro de impacto por h=b.
O programa encontra as coordenadas r e phi parametrizadas por um parametro s.
Fazemos a substituição de variável dependente u=1/r, e precisamos resolver o seguinte sistema:
    u''=(2u^3-3b^2u^5+7Mb^2u^6)E^2
    phi'=lu^2
Sendo E a energia do foton, l o momento angular, e b=l/E o parametro de impacto. Escolhendo E=1, temos l=b,
e introduzindo y=u', podemos reescrever este sistema como um sistema de edo de primeira ordem:
    u'=y
    y'=(2-(3u^2-7Mu^3)b^2)u^3
    phi'=bu^2
"""
from numpy import sqrt,arctan,cos,sin,pi
from scipy.integrate import solve_ivp
from matplotlib.pyplot import plot,show,xlim,ylim,axis,Circle,subplots,figure


b = 6.17 #Parametro de impacto
M = 1 #Massa do buraco negro
d = 10*M #Distância do plano à origem

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y , Y[2]=phi
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*Y[0]**3, b*Y[0]**2]

#valores iniciais
r0 = sqrt(b**2 + d**2)
u0 = 1/r0
phi0 = arctan(b/d)
y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))

sol = solve_ivp(F, [0, 100], [u0, y0, phi0], max_step=0.01)

figure()
#plota um circulo que representa o buraco negro
circle1 = Circle((0, 0), 2*M, color='black')
fig ,ax = subplots()
#ax.set_aspect(1)
ax.add_patch(circle1)

'''
#plota a trajetoria circular dos fotons
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
print("Angulo final:", sol.y[2][len(sol.y[2])-1])
print("Variação angular:", DeltaPhi," rad")
print("n = ", DeltaPhi/(2*pi))

