#!/usr/bin/python

import numpy as np
from astropy.io import ascii
from astropy.table import Table

def trapezoidal(function,a,b,stepsize):
    """
    Function that takes a function, lower limit a, upper limit b and stepsize and calculates the integral using trapezoidal method.

    """
    interval = np.arange(a,b+stepsize,stepsize) ## array of the interval points
    integral = 0 ## variable to store the sum
    for i in range(1,len(interval)-1,2):
        integral = integral + (stepsize/2)*(function(interval[i-1]) + 2*function(interval[i]) + function(interval[i+1]))
    return integral

def simpsons(function,a,b,stepsize):
    """
    Function that takes a function, lower limit a, upper limit b and stepsize and calculates the integral using simpsons method.
    """
    interval = np.arange(a,b+stepsize,stepsize) ## array of interrval points
    integral = 0 ## variable to store the sum
    for i in range(1,len(interval)-1,2):
        integral = integral + (stepsize/3)*(function(interval[i-1]) + 4*function(interval[i]) + function(interval[i+1]))
    return integral                           

######### Define your function and limits of integration ############

def func_l(u):
    return 3*(1 - (u)**3)**(-1/3)

def func_r(s):
    return (3/2)*(1 - (s)**(3/2))**(-2/3)

st_l = np.zeros(4)
st_r = np.zeros(4)
nbin = np.zeros(4)
tmethod = np.zeros(4)
smethod = np.zeros(4)

filename = "Homework_1_out.tex"

theor = 2*np.pi/np.sqrt(3)

for n in range(1,5): 
    a_l = 0 
    b_l = (0.5)**(1/3)
    no_of_bins = 10**n
    stepsize_l = (b_l - a_l)/no_of_bins

    a_r = 0
    b_r = (0.5)**(2/3)
    stepsize_r = (b_r - a_r)/no_of_bins

    t_inte = trapezoidal(func_l,a_l,b_l,stepsize_l) + trapezoidal(func_r,a_r,b_r,stepsize_r)
    s_inte = simpsons(func_l,a_l,b_l,stepsize_l) + simpsons(func_r,a_r,b_r,stepsize_r)

    st_l[n-1] = stepsize_l
    st_r[n-1] = stepsize_r
    nbin[n-1] = no_of_bins
    tmethod[n-1] = t_inte
    smethod[n-1] = s_inte
    ##### Printing solution on screen #######################

    print("Step size for left half and right half = {step_l:0.6f} and {step_r:0.6f}".format(step_l = stepsize_l,step_r = stepsize_r))
    print("Trapezoidal method : {inte:0.6f}".format(inte = t_inte))
    print("Simpsons method : {inte:0.6f}".format(inte = s_inte))
    print("The value of the integral analytically : {inte:1.6f}".format(inte = theor))
    print("--------------------------------------------------------------------------------------")

table = Table([nbin, st_l, st_r, tmethod, tmethod/theor, smethod, smethod/theor],names=["No. of bins", "Step size (left)", "Step size (right)", "Trapezoidal", "Trapezoidal/Theoretical", "Simpsons", "Simpsons/Theoretical"])

ascii.write(table,filename,format="latex")
