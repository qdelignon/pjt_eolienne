# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 13:39:04 2019

@author: 2017-0648
"""
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import animation


"""//////////////////////////////////////    Variables      ///////////////////////////////////////////////////"""
L = 3.2
S = 2.19e-3
rau = 1600
E = 21.3e9
J = 1.97e-6
"""///////////////////////////////////////////////////////////////////////////////////////////////////////////"""

L=3.2
rho=1600
S=2.19*(10**(-3))
E=21.3*(10**9)
I=1.97*(10**(-6))
Nk=100
Nt=50
omega = np.linspace(0,1000,Nk)


M=L*rho*S*np.array([[1/3,0,0],[0,1/5,1/7],[0,1/6,1/7]])
Kl=np.array([[E*S/L,0,0],[0,4*E*I/(L**3),6*E*I/(L**3)],[0,6*E*I/(L**3),12*E*I/(L**3)]])
M1=np.linalg.inv(M)

A=[]
for i in omega:
    N=i**2*L*rho*S*np.array([[1/3,0,0],[0,1/5,1/6],[0,1/6,1/7]])
    G=2*L*i*rho*S*np.array([[0,-1/4,-1/5],[1/4,0,0],[1/5,0,0]])
    Kp=i**2*rho*S*np.array([[L/3,0,0],[0,4*L/15,L/4],[0,L/4,9*L/35]])
    K=Kl+Kp-N
    Z=np.block([[np.zeros((3,3)),np.eye(3)],[-np.dot(M1,K),-np.dot(M1,G)]])
    X=np.linalg.eig(Z)
    A.append(X)

dbg = True

def campbell():
    lIm = np.zeros((Nk,6))
    for i in range(Nk):
        im = [A[i][0][k].imag for k in range(6)]
        im.sort()
        lIm[i] = im
    lIm = np.transpose(lIm)
    lIm=lIm[3:]
    for i in range(len(lIm)):
        plt.plot(omega, lIm[i])

"""///////////////////////////////////////   Matrices ////////////////////////////////////////////////////////"""
M = rau*S*L*np.array([[1/3,0,0],[0,1/5,1/6],[0,1/4,1/7]])
G0 = rau*S*L*np.array([[0,-1/2,-2/5],[1/2,0,0],[2/5,0,0]])
Kl = np.array([[E*S/L,0,0],[0,4*J*E/(L**3),6*E*J/(L**3)],[0,6*E*J/(L**3),12*E*J/(L**3)]])
Kp0 = rau*S*L*np.array([[1/3,0,0],[0,4/5,1/4],[0,1/4,9/35]])


def N0(x):
 return -rau*S*(L**2/2-x**2/2)

def G(omega):
    return G0*omega

def K(omega):
    return Kl+Kp0*omega**2-M*omega**2


I = np.identity(3)
O = np.zeros((3,3))
#A = np.block([[I,O],[O,M]])

def Z(omega):
    Ko = K(omega)
    Go = G(omega) 
    Mb = np.linalg.inv(M)
    return np.block([[O,I],[-Mb.dot(Ko),-Mb.dot(Go)]])



def B(omega):
    Ko = K(omega)
    Go = G(omega)
    return np.block([[O,I],[-Ko,-Go]])

def val(omega):
    a = np.linalg.inv(A)
    U = np.dot(a,B(omega))
    VAL = np.linalg.eig(Z(omega))[0]
    L=[]
    for k in VAL:
        if k.imag >=0:
            L.append(k.imag/2/np.pi*360)
    return L    
        





"""//////////////////////////////////////////////////////////CALCUL DES MODES//////////////////////////////////////////////////////"""
        


test =120

ValP = [A[k][0] for k in range(len(A))]
MU=[]
for k in range(len(A)):
    MU.append([A[k][0][i].real for i in range(6)])    
W= []
for k in range(len(A)):
    W.append([A[k][0][i].imag for i in range(6)])    
VecP = [A[k][1] for k in range(len(A))]
Yr = []
for k in range(len(A)):
    Yr.append([A[k][1][i].real for i in range(6)])
Yi = []
for k in range(len(A)):
    Yr.append([A[k][1][i].imag for i in range(6)])


def q(t,k,a,o):
    return  2*a*np.exp(MU[o][k])*(Yr[o][k]*np.cos(W[o][k]*t)-Yi[o][k]*np.sin(W[o][k]*t))

def phi(x):
    interp = np.array([[x/L,0,0],[0,(x/L)**2,(x/L)**2]])
    return interp

def u(x,t,k,a,o):
    return phi(x).dot(q(t,k,a,o))[0]

def v(x,t,k,a,o):
    return phi(x).dot(q(t,k,a,o))[1]

"""/////////////////////////////////////////////////////   ANIMATION DES MODES  /////////////////////////////////////// ////////// """
    
fig = plt.figure()
ax = plt.axes(xlim=(0,L+2),ylim=(-10,10))
line, = ax.plot([],[],lw=2)

global num_mode
global amp
num_mode = 0
amp = 2

def param(a,b):
    num_mode = a
    amp=b

#def init():
#    x = np.linspace(0,L,200)
#    y = [0 for i in x]
#    line.set_data([x],[y])
#    line.title("Animation du mode N°"+str(num_mode))
#    return line,
#
#def animate(i):
#    x = np.linspace(0,L,200)
#    V= [v(k,i,num_mode,amp) for k in x]
#    y= [u(k,i,num_mode,amp) for k in x]
#    x += V
#    line.set_data(x,y)
#    return line,ax
#
#
#anim = animation.FuncAnimation(fig, animate , frames=np.arange(1,200), interval=200)
#
#plt.show()


o = 2

#for i in np.linspace(0,100,50):
#   x = np.linspace(0,L,20)
#   V= [v(k,i,0,amp,o) for k in x]
#   y= [u(k,i,0,amp,o) for k in x]
#   x = x+V
#   plt.axes(xlim=(0,L),ylim=(-0.005,0.005))
#   plt.plot(x,y)
#plt.show()
#
#for i in np.linspace(0,100,50):
#   x = np.linspace(0,L,20)
#   V= [v(k,i,5,amp,o) for k in x]
#   y= [u(k,i,5,amp,o) for k in x]
#   x = x+V
#   plt.axes(xlim=(0,L),ylim=(-1,1))
#   plt.plot(x,y)
#plt.show()    
#
#for i in np.linspace(0,100,50):
#   x = np.linspace(0,L,20)
#   V= [v(k,i,2,amp,o) for k in x]
#   y= [u(k,i,2,amp,o) for k in x]
#   x = x+V
#   plt.axes(xlim=(0,L),ylim=(-1,1))
#   plt.plot(x,y)
#plt.show()    
    

  

    
"""/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"""
