# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:04:36 2019

@author: 2017-0648
"""

import numpy as np 
import matplotlib.pyplot as plt 


"""/////////////////////////////////////////// Données //////////////////////////////////////////////////////////////////////////"""
 
omega = np.arange(0,100,.1)
L = 3.2
S = 2.19e-3
rau = 1600
rho = 1600
E = 21.3e9
I = 1.97e-6
m = 380
mbeb = 25*10**-3


"""////////////////////////////////////////////  Matrice Globale de Masse ////////////////////////////////////////////////////////"""

M11 = rau * S * np.array([[L,0],[0,L]]) #ok
def phi(i):
    return i*2*np.pi/3
def M12(i): 
    return  rau * S * np.array([[L/2*np.cos(phi(i)),-L/3*np.sin(phi(i)),-L/4*np.sin(phi(i))],[L/2*np.sin(phi(i)),L/3*np.cos(phi(i)),L/4*np.cos(phi(i))]]) #ok
def M21(i):
    return M12(i).transpose()

M22 =L*rau*S*np.array([[1/3,0,0],[0,1/5,1/6],[0,1/6,1/7]]) #ok
Mm = np.array([[m,0],[0,m]]) #ok
O3 = np.zeros((3,3))

M = np.block([[3*M11+Mm,M12(0),M12(1),M12(2)],[M21(0),M22,O3,O3],[M21(1),O3,M22,O3],[M21(2),O3,O3,M22]])

"""/////////////////////////////////////// Matrice Globale de Dissipation (Ok) //////////////////////////////////////////////////////"""
def Gm(omega):
    return 2*m*np.array([[0,-omega],[omega,0]])
def G11(omega):
    return 2*rau * S * L * np.array([[0,-omega],[omega,0]])
def G22(omega):
    return 2*rau * S * omega *np.array([[0,-L/4,-L/5],[L/4,0,0],[L/5,0,0]])

def G12(i,omega):
    return 2*omega*rau * S * L*np.array([[-np.sin(phi(i))/2,-np.cos(phi(i))/3,-np.cos(phi(i))/4],[np.cos(phi(i))/2,-np.sin(phi(i))/3,-np.sin(phi(i))/4]])
def G21(i,omega):
    return -np.transpose(G12(i,omega))

def G(omega):
    return np.block([[3*G11(omega)+Gm(omega),G12(0,omega),G12(1,omega),G12(2,omega)],[G21(0,omega),G22(omega),O3,O3],[G21(1,omega),O3,G22(omega),O3],[G21(2,omega),O3,O3,G22(omega)]])


"""////////////////////////////////////////// Matrice Globale de Raideur /////////////////////////////////////////////////////////"""
k = 2.2*10**5

def Kg(omega):
    Km=[[k,0],[0,k]]
    Kl=np.array([[E*S/L,0,0],[0,4*E*I/(L*L*L),6*E*I/(L*L*L)],[0,6*E*I/(L*L*L),12*E*I/(L*L*L)]])
    Kp0=rho*L*S*np.array([[1/3,0,0],[0,4/15,1/4],[0,1/4,9/35]])
    Lk1=np.concatenate((Km,np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3))),axis=1)
    Lk2=np.concatenate((np.zeros((3,2)),Kl+omega**2*Kp0,np.zeros((3,3)),np.zeros((3,3))),axis=1)
    Lk3=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),Kl+omega**2*Kp0,np.zeros((3,3))),axis=1)
    Lk4=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),np.zeros((3,3)),Kl+omega**2*Kp0),axis=1)
    return np.concatenate((Lk1,Lk2,Lk3,Lk4),axis=0)


"""///////////////////////////////////////////Création de la matrice de l'équation /////////////////////////////////////////////"""
def om(omega): 
    Km=[[k,0],[0,k]]
    Kl=np.array([[E*S/L,0,0],[0,4*E*I/(L*L*L),6*E*I/(L*L*L)],[0,6*E*I/(L*L*L),12*E*I/(L*L*L)]])
    Kp0=rho*L*S*np.array([[1/3,0,0],[0,4/15,1/4],[0,1/4,9/35]])
    Lk1=np.concatenate((Km,np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3))),axis=1)
    Lk2=np.concatenate((np.zeros((3,2)),Kl+omega**2*Kp0,np.zeros((3,3)),np.zeros((3,3))),axis=1)
    Lk3=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),Kl+omega**2*Kp0,np.zeros((3,3))),axis=1)
    Lk4=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),np.zeros((3,3)),Kl+omega**2*Kp0),axis=1)
    Kglob=np.concatenate((Lk1,Lk2,Lk3,Lk4),axis=0)
    Nglob=-M*omega**2
    Z = np.block([[np.zeros((11,11)),np.identity(11)],[-np.dot(np.linalg.inv(M),Kglob+Nglob),-np.dot(np.linalg.inv(M),G(omega))]])
    EIG = np.linalg.eig(Z)
    valp = EIG[0]
    W=[]
    for i in valp:
        if i.imag>0:
            W.append(i.imag/2/np.pi)
        W.sort()
    return W




def amp(omega): 
    Kglob=Kg(omega)
    
    Nglob=-M*omega**2
    
    Fb = np.array([0 for k in range(11)])
    Fb[0] = mbeb * omega**2
    Fb.transpose()
    l = np.dot(np.linalg.inv(Kglob+Nglob),Fb)
    return  np.sqrt(l[0]**2+l[1]**2)





fig, axs = plt.subplots(2, 1)
axs[1].plot(omega,[amp(i)for i in omega])
axs[1].set_ylabel("Amplitude du mouvement")



Y= []
Y2 =[]
for i in omega:
    Y.append(om(i))
for k in range(0,5):
    axs[0].plot(omega,[Y[u][k] for u in range(len(Y))],label="mode n°"+str(k+1))
    axs[0].legend()
axs[0].set_title("Diagramme de Campbell des 8 premiers modes")
plt.xlabel("vitesse de rotation")
axs[0].set_ylabel("pulsation propre")
fig.tight_layout()
plt.show()

