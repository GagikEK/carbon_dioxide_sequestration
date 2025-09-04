import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#colonnes ['Entity', 'Code', 'Year', 'Annual CO₂ emissions']

print(np.exp(0))
data = pd.read_csv("annual-co2-emissions-per-country.csv") #emmissions en tonnes
tab_annees = []
tab_emmissions = []
for line in data.values:
    if line[0] == 'World': 
        tab_annees.append(line[2])
        tab_emmissions.append(line[3])

tab_annees = np.array(tab_annees)
tab_emmissions = np.array(tab_emmissions)
tab_emmissions = tab_emmissions

print(len(tab_annees))


def EQM(a,b):
    err = 0
    n = len(tab_annees)
    #tab_emmissions_normalise = (tab_emmissions - np.mean(tab_emmissions)) / np.std(tab_emmissions)
    
    tab_pred = np.exp(a * tab_annees + b)
    #tab_pred_normalise = (tab_pred - np.mean(tab_pred)) / np.std(tab_pred)

    tab_diff = tab_pred - tab_emmissions
    for i in range(n):
        err = err + (tab_diff[i])**2
    return err/n

def dEQM_da(a,b):
    somme = 0
    n = len(tab_annees)

    #tab_annees_normalise = (tab_annees - np.mean(tab_annees)) / np.std(tab_annees)
    #tab_emmissions_normalise = (tab_emmissions - np.mean(tab_emmissions)) / np.std(tab_emmissions)
    #tab_exp_normalise = (np.exp(a*tab_annees + b) - np.mean(np.exp(a*tab_annees + b)))/ np.std(np.exp(a*tab_annees + b))

    for i in range(n):
        ti = tab_annees[i]
        yi = tab_emmissions[i]

        somme = somme + np.exp(a*ti + b)*ti*(np.exp(a*ti + b) - yi)
    return 2*somme/n

def dEQM_db(a,b):
    somme = 0
    n = len(tab_annees)

    #tab_emmissions_normalise = (tab_emmissions - np.mean(tab_emmissions)) / np.std(tab_emmissions)
    #exp_normalise = (np.exp(a*tab_annees + b) - np.mean(np.exp(a*tab_annees + b)))/ np.std(np.exp(a*tab_annees + b))

    for i in range(n):
        yi = tab_emmissions[i]
        ti = tab_annees[i]

        somme = somme + np.exp(a*ti + b)*(np.exp(a*ti + b) - yi)
    return 2*somme/n

def descente_gradient(learning_rate_a, learning_rate_b):
    a = 0.01586
    b = -10
    print(a, b, EQM(a,b))
    plt.ion()
    fig, ax = plt.subplots()
    line1, = ax.plot(tab_annees, tab_emmissions/1000000000, label="Cible", color='blue')
    line2, = ax.plot(tab_annees, np.exp(a * tab_annees + b)/1000000000, label="Prédiction", color='red')
    ax.legend()
    for i in range(100):
        step_size_a = dEQM_da(a,b) * learning_rate_a
        a = a - step_size_a
        step_size_b = dEQM_db(a,b) * learning_rate_b
        b = b - step_size_b
        print(a, b, EQM(a,b))
        y_pred = np.exp(a * tab_annees + b)/1000000000
        line2.set_ydata(y_pred)
        ax.set_title(f"Iteration {i}")
        ax.relim()
        ax.autoscale_view()
        plt.pause(0.000001)
    plt.show()
    return a,b

#a2 = 0.0159521523623495
#b2 = -7.980498824947689

a,b = descente_gradient(1e-29, 1e-23)

a = 0.016763731933363734
b = -9.594328365363623


def f(x, a, b):
    return np.exp(a * x + b)

popt, _ = curve_fit(f, tab_annees, tab_emmissions, p0=(1e-2, -10))
a, b = popt
#print(f"Best fit: a={a}, b={b}")

#a=0.022863115449118603
#b=-21.769830005922582

plt.plot(tab_annees,tab_emmissions/1000000000)
plt.plot(tab_annees, np.exp(tab_annees*a +b)/1000000000)
plt.xlabel("Années")
plt.ylabel("C02 en Gt")
plt.title("Emissions mondiales par an en C02 dans l'atmosphère")
plt.show()
