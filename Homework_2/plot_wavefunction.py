import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

a = 3.0
m = 1.0
V0 = 10.0
hbar = 1.0

delta = np.sqrt(2.0*m*V0*(a**2.0)/(hbar**2.0))

t = Table.read("Finite_well_solutions.tex",format="latex")

eta_even = t["Even solution"].data[0]
eta_odd = t["Odd solution"].data[0]
en_even = t["Even energy"].data[0]
en_odd = t["Odd energy"].data[0]

k_even = eta_even/a
k_odd = eta_odd/a
k1_even = np.sqrt((delta**2.0) - (eta_even**2.0))/a
k1_odd = np.sqrt((delta**2.0) - (eta_odd**2.0))/a

fig = plt.figure()
plt.title("Even and odd ground states of finite well")
plt.ylim(-1,11)
plt.xlim(-7.5,7.5)
plt.ylabel("Energy (eV)")
plt.xlabel(r"r ($\AA$)")
plt.axvline(x=0,color="k",linestyle="--",linewidth=0.5)
plt.hlines(y=[0,V0,V0],xmin=[-a,-7.5,a],xmax=[a,-a,7.5],color="k")
plt.vlines(x=[a,-a],ymin=0,ymax=V0,color="k")

plt.axhline(y=en_even,color="b",linestyle="--",label="Even ground state")
plt.axhline(y=en_odd,color="y",linestyle="--",label="Odd ground state")

### Defining the three regions

x1 = np.arange(-10,-a+0.01,0.01)
x2 = np.arange(-a,a+0.01,0.01)
x3 = np.arange(a,10,0.01)

### For even wavefunction 
C = 1.0
A  = ((np.cos(k_even*a) - k_even*np.sin(k_even*a))*np.exp(k1_even*a))/(1.0 - k1_even)
psi_1_even = A*np.exp(k1_even*x1)
psi_2_even = C*np.cos(k_even*x2)
psi_3_even = A*np.exp(-k1_even*x3)

### For odd wavefunction
D = 1.0
B = ((k_odd*np.cos(k_odd*a) - np.sin(k_odd*a))*np.exp(k1_odd*a))/(1.0 + k1_odd)
psi_1_odd = B*np.exp(k1_odd*x1)
psi_2_odd = D*np.sin(k_odd*x2)
psi_3_odd = -B*np.exp(-k1_odd*x3)

plt.plot(x1,psi_1_even+en_even,"b-",label=r"Even $\psi$")
plt.plot(x2,psi_2_even+en_even,"b-")
plt.plot(x3,psi_3_even+en_even,"b-")
plt.plot(x1,psi_1_odd+en_odd,"y-",label=r"Odd $\psi$")
plt.plot(x2,psi_2_odd+en_odd,"y-")
plt.plot(x3,psi_3_odd+en_odd,"y-")

plt.legend()
plt.show()
