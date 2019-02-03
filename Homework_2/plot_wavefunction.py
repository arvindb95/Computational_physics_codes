import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

a = 3.0
m = 1
V0 = 10
hbar = 1

delta = np.sqrt(2*m*V0/(hbar**2.0))

t = Table.read("Finite_well_solutions.tex",format="latex")

eta_even = t["Even solution"].data[0]
eta_odd = t["Odd solution"].data[0]
en_even = t["Even energy"].data[0]
en_odd = t["Odd energy"].data[0]


k_even = eta_even/a
k_odd = eta_odd/a
k1_even = np.sqrt((delta**2) + (eta_even**2))/a
k1_odd = np.sqrt((delta**2) + (eta_odd**2))/a

def even_wavefunction(x,k1,k,a):
    C = 1
    A = ((np.cos(k*a) - k*np.sin(k*a))*np.exp(k1*a))/(1-k1)
    print(A)
    if (np.all(x <= -a)):
        print("We got to die down at negative infinity!")
        return A*np.exp(k1*x)
    if (np.all(x > a)):
        return A*np.exp(-k1*x)
    if (-a < np.all(x) <= a):
        return C*np.cos(k*x)

def odd_wavefunction(x,k1,k,a):
    D = 1
    B = ((k*np.cos(k*a) - np.sin(k*a))*np.exp(k1*a))/(1-k1)
    if (np.all(x <= -a)):
        return -B*np.exp(k1*x)
    if (np.all(x > a)): 
        return B*np.exp(-k1*x)
    if (-a < np.all(x) <= a):
        return D*np.sin(k*x)




fig = plt.figure()
plt.title("Even and odd ground states of finite well")
#plt.ylim(-0.05,0.5)
plt.xlim(-10,10)
plt.axvline(x=0,color="k",linestyle="--",linewidth=0.5)
plt.hlines(y=[0,10,10],xmin=[-3,-5,3],xmax=[3,-3,5],color="k")
plt.vlines(x=[3,-3],ymin=0,ymax=10,color="k")

plt.axhline(y=en_even,color="b",linestyle="--",label="Even ground state")
plt.axhline(y=en_odd,color="y",linestyle="--",label="Odd ground state")

x1 = np.arange(-10,-a+0.01,0.01)
x2 = np.arange(-a+0.01,a+0.01,0.01)
x3 = np.arange(a+0.01,10+0.01,0.01)

plt.plot(x1,even_wavefunction(x1,k1_even,k_even,a)+en_even,"b-",label=r"Even $\psi$")
plt.plot(x2,even_wavefunction(x2,k1_even,k_even,a)+en_even,"b-")
plt.plot(x3,even_wavefunction(x3,k1_even,k_even,a)+en_even,"b-")
#plt.plot(x1,-A*np.exp(k1_odd*x1),"y-",label=r"Odd $\psi$")
#plt.plot(x2,np.sin(k_odd*x2),"y-")
plt.plot(x3,odd_wavefunction(x3,k1_odd,k_odd,a),"y-")



plt.legend()
plt.show()
