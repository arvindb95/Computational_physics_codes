import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii

a = 3.0               # the half width of the well in angstroms 
m = 1.0               # mass of electron in 1 me units
V0 = 10.0             # height/depth of well in eV
hbar = 1.0            # hbar

delta = np.sqrt(2.0*m*V0*(a**2.0)/(hbar**2.0))     # Defined in the hand written solution
eta = np.arange(0.001,2*np.pi,0.01)  # Defined in the hand written solutiom

def f(f,eta): 
    """
    This is the lhs in the transcendental equations.
    """
    return np.sqrt(((delta/eta)**2.0)-1.0)

def f1(eta):
    """
    The equation whose zeros are to be found to obtain the even solution.
    """
    return f(eta) - np.tan(eta)

def f2(eta):
    """
    The equation whose zeros are to be found to obtain the odd solution.
    """
    return f(eta) + (1.0/np.tan(eta))

def firstDerivative_O4(function,x0,stepsize):
    """
    Returns first derivative of the function "function" at the point x0 by considering points, one and two steps on either side of x0.
    The Accuracy is of order (stepsize)^4.
    """
    return (function(x0 - 2.0*stepsize) - 8.0*function(x0 - stepsize) + 8.0*function(x0 + stepsize) - function(x0 + 2.0*stepsize))/(12.0*stepsize)

def Newton_Raphson(f,x0,stepsize):
    """
    Returns the zero of f. x0 is the guess. stepsize will be used for determining the first derivative of f.
    """
    fprime = firstDerivative_O4(f,x0,stepsize)
    x1 = x0 - (f(x0)/fprime)
    while((x0 - x1) >= stepsize**(8.0)):
        x0 = x1
        fprime = firstDerivative_O4(f,x0,stepsize)
        x1 = x0 - (f(x0)/fprime)
    return x1
    
def get_energy(eta):
    """
    Gets energy value for given eta value.
    """
    return (eta**2.0)*(hbar**2.0)/(2.0*m*(a**2.0))

##### Main portion begins #####

even_guess = np.pi/4.0   # Guess for even solution
odd_guess = 0.999*np.pi  # Guess for odd solution

even_sol = Newton_Raphson(f1,even_guess,10**(-3)) # Obtaining the solution using the Newton Raphson method
odd_sol = Newton_Raphson(f2,odd_guess,10**(-3))   # Obtaining the solution using the Newton Raphson method


print("eta for even sol : {etev:0.6f}".format(etev=even_sol))
print("eta for odd sol : {etod:0.6f}".format(etod=odd_sol))
print("Energy for even state is : {e2:0.6f}".format(e2=get_energy(even_sol)))
print("Energy for odd state is : {e1:0.6f}".format(e1=get_energy(odd_sol)))

## Plotting the graph

plt.plot(np.arange(0,(np.pi/2.0)-0.01,0.01),np.tan(np.arange(0,(np.pi/2.0)-0.01,0.01)),"b-",label=r"$tan(\eta)$")
plt.plot(np.arange((np.pi/2.0)+0.01,(3*np.pi/2.0)-0.01,0.01),np.tan(np.arange((np.pi/2.0)+0.01,(3*np.pi/2.0)-0.01,0.01)),"b-")
plt.plot(np.arange(0.01,np.pi-0.01,0.01),-1.0/np.tan(np.arange(0.01,np.pi-0.01,0.01)),"y-",label=r"$-cot(\eta)$")
plt.plot(np.arange(np.pi+0.01,2*np.pi-0.01,0.01),-1.0/np.tan(np.arange(np.pi+0.01,2*np.pi-0.01,0.01)),"y-")
plt.plot(eta,f(eta),"r-",label=r"$\sqrt{\frac{\Delta^{2}}{\eta^{2}}-1}$")
plt.ylim(-10,50)
plt.xticks(np.array([0,np.pi/4.0,np.pi/2.0,3*np.pi/4.0,np.pi,5*np.pi/4.0]),("0",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$",r"$\frac{5\pi}{4}$"))
plt.xlim(0,4)
plt.title("Even and odd solutions to finite well")
plt.xlabel(r"$\eta$")
plt.axhline(y=0,color="black",linewidth=0.5)
plt.axvline(x=even_guess,color="b",linestyle="--",linewidth=1,label="Guess even solution")
plt.axvline(x=even_sol,color="b",linewidth=1)
plt.axvline(x=odd_guess,color="y",linestyle="--",linewidth=1,label="Guess odd solution")
plt.axvline(x=odd_sol,color="y",linewidth=1)
plt.plot(even_sol,np.tan(even_sol),"bo",label="Even solution")
plt.plot(odd_sol,-1.0/np.tan(odd_sol),"yo",label="Odd solution")
plt.legend(loc='upper right',fontsize='xx-small')
plt.show()

t = Table([[even_sol],[odd_sol],[get_energy(even_sol)],[get_energy(odd_sol)]], names=["Even solution","Odd solution","Even energy","Odd energy"])

ascii.write(t,"Finite_well_solutions.tex",format="latex")

