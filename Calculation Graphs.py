# -*- coding: utf-8 -*-

#Calculation Graphs for SIMclean
import matplotlib.pyplot as plt
import math

g = 9.81
d = 5.6
h = 1.5
b = 2.37
a = 0.5

M_M = 5 + 1.5
M_S = 15
M_S_ex = 40

rho_B = 0.8
rho_R = 0.1

rope_E = 4e9
bar_E = 70e9

I_V = 14.77e-8
I_H = 1.64e-8

area_V = 384e-6
area_H = 282e-6
r = 0.01

F_M = g*M_M
F_S = g*M_S
F_B = g*d*rho_B
F_S_ex = g*M_S_ex

def alpha(theta):
    return math.atan2(d*math.cos(theta), h+d*math.sin(theta))

def phi(theta):
    return math.pi/2-alpha(theta)-theta

def F_T(x, theta):
    return math.cos(theta)/math.sin(phi(theta))*(x/d*F_M+F_B/2)

def delta_R(x, theta):
    return F_T(x, theta)/(rope_E*math.pi*r*r+F_T(x, theta))*(d*math.cos(theta)**2+(h+d*math.sin(theta))**2)**0.5

def F_X(x, theta):
    return F_T(x, theta)*math.cos(phi(theta)+theta)

def F_Y(x, theta):
    return F_M+F_B-F_T(x, theta)*math.sin(phi(theta)+theta)

def F_1(x, theta):
    return ((F_S+F_T(x, theta)*math.cos(alpha(theta))+F_Y(x, theta)+(h/a)*F_X(x, theta))/2)

def F_2(x, theta):
    return ((F_S+F_T(x, theta)*math.cos(alpha(theta))+F_Y(x, theta)-(h/a)*F_X(x, theta))/2) + F_S_ex

def Im(vel):
    return M_M*(vel/0.1)

def res(x, theta):
    if F_M > F_T(x, theta)*math.sin(alpha(theta)):
        return F_M-(F_T(x, theta)*math.sin(alpha(theta)))
    else:
        return 0
    
def flex_V(x, theta):
    return (F_T(x, theta)*math.cos(alpha(theta))*(h+b)**3)/(3*bar_E*I_V)

def flex_H(x, theta):
    return (res(x, theta)*d**3)/(3*bar_E*I_H)

def flex_Hm(x):
    return (F_M*x*(d**2-x**2)**1.5)/(9*(3**0.5)*d*bar_E*I_H)

def delta_V(x, theta):
    return (F_T(x, theta)*math.sin(alpha(theta))*(h+b))/(bar_E*area_V)

def delta_H(x, theta):
    return (F_T(x, theta)*math.cos(alpha(theta))*d)/(bar_E*area_H)

td = [t for t in range(-30, 61)]
ts = [t*(math.pi/180) for t in td]
vs = [v/100 for v in range(-100, 100)]
iz = [i/100 for i in range(0, 281)]

Ims = [Im(v) for v in vs]

FTs = [F_T(d, t) for t in ts]
FXs = [F_X(d, t) for t in ts]
FYs = [F_Y(d, t) for t in ts]
F1s = [F_1(d, t) for t in ts]
F2s = [F_2(d, t) for t in ts]

Dels = [delta_R(d, t)*1000 for t in ts]
Delv = [delta_V(d, t)*1000 for t in ts]
Delh = [delta_H(d, t)*1000 for t in ts]

DFv = [flex_V(d, t) for t in ts]
DFh = [flex_H(d, t) for t in ts]
DFhm = [flex_Hm(i) for i in iz]

def y_lb(args):
    if Ims in args:
        return 'Impulse (N)' 
    
    if DFv in args:
        return 'Length (m)'
    if DFh in args:
        return 'Length (m)'
    if DFhm in args:
        return 'Length (m)'
    
    if Dels in args:
        return 'Length (mm)'
    if Delv in args:
        return 'Length (mm)'
    if Delh in args:
        return 'Length (mm)'
    
    return 'Force (N)'

def x_lb(args):
    if td in args:
        return 'Angle of Bar to Horizontal (deg)'
    elif iz in args:
        return 'Distance of Module from End (m)'
    else:
        return 'Velocity (m/s)'
    
def save_plot(name, args):
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot()
    ax.axhline(linewidth=0.5, color='g')
    ax.axvline(linewidth=0.5, color='g')
    ax.plot(*args)
    plt.title(name, size=15)
    plt.ylabel(y_lb(args))
    plt.xlabel(x_lb(args))
    plt.savefig(name)

save_plot('Stopping Module Impulse', [vs, Ims])

save_plot('Rope Tension', [td, FTs])
save_plot('Horizontal Bar Reactions', [td, FXs, td, FYs])
save_plot('Base Support Reactions', [td, F1s, td, F2s])

save_plot('Rope Elongation', [td, Dels])
save_plot('Vertical Bar Elongation', [td, Delv])
save_plot('Horizontal Bar Elongation', [td, Delh])

save_plot('Vertical Bar End Deflection', [td, DFv])
save_plot('Horizontal Bar End Deflection', [td, DFh])
save_plot('Horizontal Bar Mid Deflection', [iz, DFhm])

x = d
theta = 30*(math.pi/180)

sol = [
    F_T(x, theta), delta_R(x, theta), F_X(x, theta), F_Y(x, theta), F_1(x, theta), F_2(x, theta), 
    flex_V(x, theta), flex_H(x, theta), flex_Hm(x), delta_V(x, theta), delta_H(x, theta)
]

print(sol)























