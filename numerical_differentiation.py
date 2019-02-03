import numpy as np

def firstDerivative_O2(function,x0,stepsize):
    """
    Returns first derivative of the function "function" at the point x0 by considering points, one step on either side of x0.
    The Accuracy is of order (stepsize)^2.
    """
    return (function(x0+stepsize) - function(x0-stepsize))/(2*stepsize)

def firstDerivative_O4(function,x0,stepsize):
    """
    Returns first derivative of the function "function" at the point x0 by considering points, one and two steps on either side of x0.
    The Accuracy is of order (stepsize)^4.
    """
    return (function(x0 - 2*stepsize) - 8*function(x0 - stepsize) + 8*function(x0 + stepsize) - function(x0 + 2*stepsize))/(12*stepsize)

def secondDerivative_O4(function,x0,stepsize):
    """
    Returns second derivative of the function "function" at the point x0 by considering points, one step on either side of x0.
    The Accuracy is of order (stepsize)^4.
    """
    return (function(x0 - stepsize) - 2*function(x0) + function(x0 + stepsize))/(stepsize**2)

def secondDerivative_O6(function,x0,stepsize):
    """
    Returns first derivative of the function "function" at the point x0 by considering points, one and two steps on either side of x0.
    The Accuracy is of order (stepsize)^6.
    """
    return (-function(x0 - 2*stepsize) + 16*function(x0 - stepsize) - 30*function(x0) + 16*function(x0 + stepsize) - function(x0 + 2*stepsize))/(12*(stepsize**2))


### Define function for which derivative should be taken ###

def f(x):
    return x**3

### Define the point at which derivative should be taken and the stepsize to be used ###

x0 = 2
stepsize = 0.1

############################################################

print("The first derivative of the function at {a:0.2f} is {b:0.6f}, with accuracy of order 2 and is {c:0.6f}, with accuracy of order 4".format(a = x0, b = firstDerivative_O2(f,x0,stepsize), c = firstDerivative_O4(f,x0,stepsize)))
print("####################################################################################################")
print("The second derivative of the function at {a:0.2f} is {b:0.6f}, with accuracy of order 4 and is {c:0.6f}, with accuracy of order 6".format(a = x0, b = secondDerivative_O4(f,x0,stepsize), c = secondDerivative_O6(f,x0,stepsize)))

