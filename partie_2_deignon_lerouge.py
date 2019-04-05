# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
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



"""////////////////////////////////////////// Matrice Globale de Raideur /////////////////////////////////////////////////////////"""
k = 2.2*10**5
O22 = np.zeros((2,2))
O23  = np.zeros((2,3))
O32 = np.zeros((3,2))
Kl=np.array([[E*S/L,0,0],[0,4*E*I/(L*L*L),6*E*I/(L*L*L)],[0,6*E*I/(L*L*L),12*E*I/(L*L*L)]])

def testk(omega):
    Kp = omega*rau*S*L*np.array([[1/3,0,0],[0,4/15,1/4],[0,1/4,9/35]])
    Km = np.array([[k,0],[0,k]])
    Ki = Kl+Kp
    K = np.block([[Km,O23,O23,O23],[O32,Ki,O3,O3],[O32,O3,Ki,O3],[O32,O3,O3,Ki]])
    N = -omega**2*M
    K2 =K + N
    Km=[[k,0],[0,k]]
    Kp0=rho*L*S*np.array([[1/3,0,0],[0,4/15,1/4],[0,1/4,9/35]])
    Lk1=np.concatenate((Km,np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3))),axis=1)
    Lk2=np.concatenate((np.zeros((3,2)),Kl+omega**2*Kp0,np.zeros((3,3)),np.zeros((3,3))),axis=1)
    Lk3=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),Kl+omega**2*Kp0,np.zeros((3,3))),axis=1)
    Lk4=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),np.zeros((3,3)),Kl+omega**2*Kp0),axis=1)
    Kglob=np.concatenate((Lk1,Lk2,Lk3,Lk4),axis=0)
    return K2-Kglob
    
#Km=[[k,0],[0,k]]
#Kl=np.array([[E*S/L,0,0],[0,4*E*I/(L*L*L),6*E*I/(L*L*L)],[0,6*E*I/(L*L*L),12*E*I/(L*L*L)]])
#Kp0=rho*L*S*np.array([[1/3,0,0],[0,4/15,1/4],[0,1/4,9/35]])
#Lk1=np.concatenate((Km,np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3))),axis=1)
#Lk2=np.concatenate((np.zeros((3,2)),Kl+omega[i]**2*Kp0,np.zeros((3,3)),np.zeros((3,3))),axis=1)
#Lk3=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),Kl+omega[i]**2*Kp0,np.zeros((3,3))),axis=1)
#Lk4=np.concatenate((np.zeros((3,2)),np.zeros((3,3)),np.zeros((3,3)),Kl+omega[i]**2*Kp0),axis=1)
#Kglob=np.concatenate((Lk1,Lk2,Lk3,Lk4),axis=0)





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
    
def mu(omega):
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
            W.append(i.real)
        W.sort()
    return W


Y= []
Y2 =[]
for i in omega:
    Y.append(om(i))
    Y2.append(mu(i))
for k in range(0,8):
    plt.plot(omega,[Y[u][k] for u in range(len(Y))],label="mode n°"+str(k+1))
    plt.legend()
plt.title("Diagramme de Campbell des 8 premiers modes");plt.xlabel("vitesse de rotation"),plt.ylabel("pulsation propre")
plt.show()

for k in range(8,10):
    plt.plot(omega,[Y[u][k] for u in range(len(Y))],label="mode n°"+str(k+1))
    plt.legend()
plt.title("Diagramme de Campbell des 2 derniers modes");plt.xlabel("vitesse de rotation"),plt.ylabel("pulsation propre")
plt.show()

for k in range(0,10):
    plt.plot(omega,[abs(Y2[u][k]) for u in range(len(Y))],label="mode n°"+str(k+1))
    plt.title('Amortissement propre des modes en fonction de Oméga')
plt.xlabel("oméga");plt.ylabel("mu")
plt.show()




#paramètres
omega = 100
w=20
pi = np.pi
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


"""Représentation d'un epale"""

ValeursP=np.linalg.eig(Z)[0]
LVectP = [np.linalg.eig(Z)[1][i] for i in range(22)]
mu=np.real(ValeursP)
Yck=np.linalg.eig(Z)[1]
T=np.linspace(0,2*pi/w,10)
X=np.linspace(0,L,100)
Yk=np.zeros(((len(LVectP[:11]),len(T))))
u1=np.zeros((len(X),len(T)))
v1=np.zeros((len(X),len(T)))
u2=np.zeros((len(X),len(T)))
v2=np.zeros((len(X),len(T)))
u3=np.zeros((len(X),len(T)))
v3=np.zeros((len(X),len(T)))
P0=np.array([X,np.zeros(len(X))])
P1=np.zeros((100,2,100))
P2=np.zeros((100,2,100))
P3=np.zeros((100,2,100))


for j in range(len(T)):
    Yk[:,j]=1000*np.exp(mu[:11]*T[j])*(np.real(Yck[6][:11])*np.cos(w*T[j])-np.imag(Yck[6][:11])*np.sin(w*T[j]))
    for n in range(len(X)):
        u1[n][j]=(1/L)*X[n]*Yk[2][j]                                    #Calcul déplacement en 'x' de la pale 1
        v1[n][j]=(1/L**2)*X[n]**2*Yk[3][j]+(1/L**3)*X[n]**3*Yk[4][j]    #Calcul déplaceemnt en 'y' de la pale 1
        u2[n][j]=(1/L)*X[n]*Yk[5][j]
        v2[n][j]=(1/L**2)*X[n]**2*Yk[6][j]+(1/L**3)*X[n]**3*Yk[7][j]
        u3[n][j]=(1/L)*X[n]*Yk[8][j]
        v3[n][j]=(1/L**2)*X[n]**2*Yk[9][j]+(1/L**3)*X[n]**3*Yk[10][j]
        
        P1[j,0,n]+=P0[0,n]+u1[n][j]+Yk[0,j]                                     #Ajout de la position n sur la pale en 'x' au temps j avec prise en compte de la position décalée du centre de l'éolienne
        P1[j,1,n]+=P0[1,n]+v1[n][j]+Yk[1,j]                                     #Ajout de la position n sur la pale en 'y' au temps j avec prise en compte de la position décalée du centre de l'éolienne
    
        P2[j,0,n]+= P0[0,n]*np.cos(2*pi/3)+P0[1,n]*np.sin(2*pi/3) + u2[n][j] + Yk[0,j]
        P2[j,1,n]+= -P0[0,n]*np.sin(2*pi/3)+P0[1,n]*np.cos(2*pi/3) + v2[n][j] + Yk[1,j]
        
        P3[j,0,n]+= P0[0,n]*np.cos(4*pi/3)+P0[1,n]*np.sin(4*pi/3) + u3[n][j] + Yk[0,j]
        P3[j,1,n]+= -P0[0,n]*np.sin(4*pi/3)+P0[1,n]*np.cos(4*pi/3) + v3[n][j] + Yk[1,j]

    plt.plot(P1[j,0,:],P1[j,1,:])
    plt.plot(P2[j,0,:],P2[j,1,:])
    plt.plot(P3[j,0,:],P3[j,1,:])
    plt.grid()
    plt.axis('equal')
    plt.xlabel("x0 (m)")
    plt.ylabel("y0 (m)")
plt.grid()
plt.show()
