import numpy as np
from astropy.table import Table
from astropy.io import ascii

## Comparing different algorithms to solve dy/dx = f(x,y)

def euler_method(function,x0,y0,x,stepsize):
    """
    function = dy/dx; x0,y0 are initial conditions; x is the point where y needs to be calculated.
    """
    x_array = np.arange(x0, x + stepsize, stepsize)
    y_array = np.zeros(len(x_array))
    y_array[0] = y0
    for i in range(1,len(x_array)):
        y_array[i] = y_array[i-1] + stepsize*function(x_array[i-1],y_array[i-1])
    
    return y_array[-1]

def partialx(function,x0,y0,stepsize):
    """
    Returns partial first derivative of the function "function" at the point x0,y0 with respect to x by considering points, one and two steps on either side of x0.
    The Accuracy is of order (stepsize)^4.
    """
    return (function(x0 - 2*stepsize, y0) - 8*function(x0 - stepsize, y0) + 8*function(x0 + stepsize, y0) - function(x0 + 2*stepsize, y0))/(12*stepsize)

def partialy(function,x0,y0,stepsize):
    """
    Returns partial first derivative of the function "function" at the point x0,y0 with respect to y by considering points, one and two steps on either side of y0.
    The Accuracy is of order (stepsize)^4.
    """
    return (function(x0, y0 - 2*stepsize) - 8*function(x0, y0 - stepsize) + 8*function(x0, y0 + stepsize) - function(x0, y0 + 2*stepsize))/(12*stepsize)


def taylor_method(function,x0,y0,x,stepsize):
    """
    function = dy/dx; x0,y0 are initial conditions; x is the point where y needs to be calculated.
    """
    x_array = np.arange(x0, x + stepsize, stepsize)
    y_array = np.zeros(len(x_array))
    y_array[0] = y0
    for i in range(1,len(x_array)):
        y_array[i] = y_array[i-1] + stepsize*function(x_array[i-1],y_array[i-1]) + (stepsize**2.0)*(partialx(function,x_array[i-1],y_array[i-1],10**(-7)) + function(x_array[i-1],y_array[i-1])*partialy(function,x_array[i-1],y_array[i-1],10**(-7)))
    
    return y_array[-1]

def adam_bashford_method(function,x0,y0,x,stepsize):
    """
    function = dy/dx; x0,y0 are initial conditions; x is the point where y needs to be calculated.
    """
    x_array = np.arange(x0, x + stepsize, stepsize)
    y_array = np.zeros(len(x_array))
    y_array[0] = y0
    y_array[1] = y0
    for i in range(2,len(x_array)):
        y_array[i] = y_array[i-1] + stepsize*(1.5*function(x_array[i-1],y_array[i-1]) - 0.5*function(x_array[i-2],y_array[i-2]))

    return y_array[-1]

def runge_kutta_method(function,x0,y0,x,stepsize):
    """
    function = dy/dx; x0,y0 are initial conditions; x is the point where y needs to be calculated.
    """
    x_array = np.arange(x0, x + stepsize, stepsize)
    y_array = np.zeros(len(x_array))
    y_array[0] = y0
    for i in range(1,len(x_array)):
        y_step = stepsize*function(x_array[i-1], y_array[i-1])
        y_array[i] = y_array[i-1] + stepsize*function(x_array[i-1] + (stepsize/2.0), y_array[i-1] + (y_step/2.0))
    
    return y_array[-1] 

def f(x,y):
    return -x*y

x0 = 0.0
y0 = 1.0
x1 = 1.0
x3 = 3.0

x = np.array([1.0,3.0])

print(euler_method(f,x0,y0,np.all(x),0.001))

stepsize_array = np.array([0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001])

#e1_array = np.zeros(len(stepsize_array))
#t1_array = np.zeros(len(stepsize_array))
#ab1_array = np.zeros(len(stepsize_array))
#rk1_array = np.zeros(len(stepsize_array))
#
#e3_array = np.zeros(len(stepsize_array))
#t3_array = np.zeros(len(stepsize_array))
#ab3_array = np.zeros(len(stepsize_array))
#rk3_array = np.zeros(len(stepsize_array))
#
#e1_diff = np.zeros(len(stepsize_array))
#t1_diff = np.zeros(len(stepsize_array))
#ab1_diff = np.zeros(len(stepsize_array))
#rk1_diff = np.zeros(len(stepsize_array))
#
#e3_diff = np.zeros(len(stepsize_array))
#t3_diff = np.zeros(len(stepsize_array))
#ab3_diff = np.zeros(len(stepsize_array))
#rk3_diff = np.zeros(len(stepsize_array))
#
#real_value_1 = np.exp(-1.0/2.0)
#real_value_3 = np.exp(-9.0/2.0)
#
#for i in range(len(stepsize_array)):
#    e1_array[i] = "{a:0.6f}".format(a=euler_method(f,x0,y0,x1,stepsize_array[i]))
#    e3_array[i] = "{a:0.6f}".format(a=euler_method(f,x0,y0,x3,stepsize_array[i]))
#    e1_diff[i] = "{b:0.3E}".format(b=abs(euler_method(f,x0,y0,x1,stepsize_array[i]) - real_value_1))
#    e3_diff[i] = "{b:0.3E}".format(b=abs(euler_method(f,x0,y0,x3,stepsize_array[i]) - real_value_3))
#
#    t1_array[i] = "{a:0.6f}".format(a=taylor_method(f,x0,y0,x1,stepsize_array[i]))
#    t3_array[i] = "{a:0.6f}".format(a=taylor_method(f,x0,y0,x3,stepsize_array[i]))
#    t1_diff[i] = "{b:0.3E}".format(b=abs(taylor_method(f,x0,y0,x1,stepsize_array[i]) - real_value_1))
#    t3_diff[i] = "{b:0.3E}".format(b=abs(taylor_method(f,x0,y0,x3,stepsize_array[i]) - real_value_3))
#
#    ab1_array[i] = "{a:0.6f}".format(a=adam_bashford_method(f,x0,y0,x1,stepsize_array[i]))
#    ab3_array[i] = "{a:0.6f}".format(a=adam_bashford_method(f,x0,y0,x3,stepsize_array[i]))
#    ab1_diff[i] = "{b:0.3E}".format(b=abs(adam_bashford_method(f,x0,y0,x1,stepsize_array[i]) - real_value_1))
#    ab3_diff[i] = "{b:0.3E}".format(b=abs(adam_bashford_method(f,x0,y0,x3,stepsize_array[i]) - real_value_3))
#
#    rk1_array[i] = "{a:0.6f}".format(a=runge_kutta_method(f,x0,y0,x1,stepsize_array[i]))
#    rk3_array[i] = "{a:0.6f}".format(a=runge_kutta_method(f,x0,y0,x3,stepsize_array[i]))
#    rk1_diff[i] = "{b:0.3E}".format(b=abs(runge_kutta_method(f,x0,y0,x1,stepsize_array[i]) - real_value_1))
#    rk3_diff[i] = "{b:0.3E}".format(b=abs(runge_kutta_method(f,x0,y0,x3,stepsize_array[i]) - real_value_3))
#
#
#t1 = Table([stepsize_array, e1_array, e1_diff, t1_array, t1_diff, ab1_array, ab1_diff, rk1_array, rk1_diff], names=["h", "e1", "e1_diff", "t1", "t1_diff", "ab1", "ab1_diff", "rk1", "rk1_diff"])
#t3 = Table([stepsize_array, e3_array, e3_diff, t3_array, t3_diff, ab3_array, ab3_diff, rk3_array, rk3_diff], names=["h", "e3", "e3_diff", "t3", "t3_diff", "ab3", "ab3_diff", "rk3", "rk3_diff"])
#
#ascii.write(t1,"final_table_x1.tex",format="latex")
#ascii.write(t3,"final_table_x3.tex", format="latex")

#print("Real value : ",np.exp(-1.0/2.0))
#print("Euler method : ",euler_method(f,0.0,1.0,1.0,0.001))
#print("Taylor method : ",taylor_method(f,0.0,1.0,1.0,0.001))
#print("Adam Bashford method : ",adam_bashford_method(f,0.0,1.0,1.0,0.001))
#print("Runge Kutta method : ",runge_kutta_method(f,0.0,1.0,1.0,0.001))
