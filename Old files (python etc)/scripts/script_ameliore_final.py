import numpy as np
import matplotlib.pyplot as plt
import math

#Le code risque d'etre peu compréhensible sans avoir lu les slides


#Variables globales

#Constantes ******************************
alpha = 0.004 #0.004
beta = 0.02 #0.02
gamma = 0.015 #0.015
delta = 0.008 #0.008
K = 900 #900

#océan
alpha2 = 0.001 #0.001
omega = 0.0003 #0.0003
K2 = 40000 #40000

#exponentielle humains grâce à scipy
a = 0.022863115449118603
b = -21.769830005922582 

#Constantes ******************************

#Autres variables ************************
max_iter = 10000
eps = 1e-6
Tf = 2100 #temps final 


    #Conditions initiales
CA0 = 3306.296 #3306.296
CT0 = 370.23 #370.23
CS0 = 490.77 #490.77
CO0 = 40000 #40000

C0 = [CA0, CT0, CS0, CO0]
#Plus tard, le tableau de tableau C contiendra à chaque colonne n le tableau Cn
t0 = 2024 #temps initial
h = 0.1 #pas entre chaque temps

#Autres variables ************************

#fonctions********************************

def terme_source(t):
    if(t > 2020): 
        return np.array([a*2020 + b,0,0,0])
    else:
        return np.array([np.exp(t*a + b),0,0,0])/1000000000 #conversion de t en Gt

def f(t, Cn): 
    """
        fonction f définie dans la présentation.
        Args :
            Cn : numpy array ou python array (représente l'état du système au temps tn)
        Returns :
            numpy array représente f(Cn)
            
    """
    y1 = Cn[1] * (-alpha * (1 - (Cn[1])/K) + beta) + delta * Cn[2] - Cn[3] * (alpha2 * (1 - (Cn[3]/K2))) + omega * Cn[3]
    y2 = Cn[1] * (alpha * (1 - (Cn[1])/K) - beta - delta - gamma)
    y3 = Cn[1] * (gamma + delta) - delta * Cn[2]
    y4 = Cn[3] * (alpha2 * (1 - (Cn[3]/K2))) - omega * Cn[3]
    return np.array([y1,y2,y3,y4]) + terme_source(t)

def F2(t, Cnk_suiv, Cnk_prec):
    """
        fonction F2, même que F1 mais avec la méthode des trapèzes
        Args:
            Cnk_suiv : numpy array ou python array
            Cnk_prec : numpy array ou python array
        Returns:
            numpy array
    """
    Cnk_prec = np.array(Cnk_prec)
    Cnk_suiv = np.array(Cnk_suiv)
    return Cnk_prec + (h/2)*(f(t+h, Cnk_suiv) + f(t, Cnk_prec))

def f_newton(t, Cnk_suiv, Cnk_prec):
    """
        fonciton f_newton qui implémente celle décrite dans le rapport, on cherche à l'annuler
        Args:
            Cnk_suiv : numpy array ou python array
            Cnk_prec : numpy array ou python array
        Returns:
            numpy array 
    """
    return F2(t, Cnk_suiv, Cnk_prec) - Cnk_suiv

def df_newton(Cnk):
    """
        fonction df_newton qui représente la jacobienne de f_newton évaluée en le veccteur Cnk
        Args:
            Cnk : numpy array ou python array
        Returns:
            numpy array : jacobienne de f_newton évaluée en Cnk
    """
    A = np.zeros((4,4))
    A[:,1] = [-alpha + (2*alpha*Cnk[1])/K + beta, alpha - (2*alpha*Cnk[1])/K - beta - delta - gamma, gamma + delta, 0]
    A[:,2] = [delta, 0, -delta, 0]
    A[:,3] = [-alpha2 + (2*alpha2*Cnk[3])/K2 + omega, 0, 0, alpha2 - (2*alpha2*Cnk[3])/K - omega]
    return (h/2) * A - np.eye(4)

#fonctions********************************



#methodes numeriques**********************

def newton(t, X0, f, df, _eps=eps, _max_iter = max_iter):
    """
        fonction newton : résout f(X) = 0 par la méthode de newton
        Args:
            X0 : numpy array ou python array (premier terme de la suite de newton Xnk, (représente Cn-1))
            f : fonction de newton choisie 
            _eps : float (tolérance d'erreur)
            _max_iter : int (nombre maximum d'itérations de la méthode)
        Returns:
            numpy array : c'est le point où f s'annule (ici c'est Cn+1)
    """
    X0 = np.array(X0)
    Xk = X0
    for i in range(_max_iter):
        if(np.linalg.norm(f(t, Xk, X0)) < _eps):
            break
        Xk = Xk - np.linalg.inv(df(Xk)) @ f(t, Xk, X0)
    return Xk

def trapeze_newton(_C0):
    """
        applique la méthode de newton et des trapèzes combinés (comme décrit dans les slides) pour résoudre le problème
        Args:
            _C0 : numpy array ou python array (condition initiale de l'EDO)
        Returns :
            C : numpy array contenant à la colonne n le vecteur Cn
            T : python array qui contient tous les temps tn
    """
    t = t0
    C = np.zeros((4,1))
    T = [t0]
    C[:,0] = _C0
    k = 1
    while(t < Tf):
        t = t+h
        C = np.append(C, np.transpose( [newton(t, C[:,k-1],f_newton, df_newton)] ), axis = 1) #ajoute Cn au tableau C tel que à la colonne i on ait Ci
        T.append(t)
        k = k+1
    return C, T

#methodes numeriques**********************

C3, T3 = trapeze_newton(C0)

plt.plot(T3, C3[0], label='CA(t)')
plt.plot(T3, C3[1], label='CT(t)')
plt.plot(T3, C3[2], label='CS(t)')
#plt.plot(T3, C3[3], label='CO(t)')
plt.legend()
plt.xlabel("t")
plt.title("évolution des grandeurs")

plt.show()