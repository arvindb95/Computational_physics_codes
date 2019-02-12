import numpy as np
import matplotlib.pyplot as plt
import cmath as c
# Solving d^2x/dt^2 = -k*x^alpha
# We know that v = dx/dt and dv/dt = -k*x^alpha

def aofx(x,k,alpha):
    """
    Returns the function that equals d^2x/dt^2
    """
    return -k*(x**alpha)

def euler_cromer(x0, t0, t, delta_t, k, alpha):
    """
    Returns an array of x values from time t = t0 to t = t
    """    
    t_array = np.arange(t0, t + delta_t, delta_t)
    x_array = np.zeros(len(t_array))
    v_array = np.zeros(len(t_array))
    x_array[0] = x0

    for i in range(1, len(t_array)):
        v_array[i] = v_array[i-1] + delta_t*aofx(x_array[i-1],k,alpha)
        x_array[i] = x_array[i-1] + delta_t*v_array[i]

    return t_array, x_array, v_array

def FourierTransform(t_array, f_array, w_ini, w_fin):
    """

    """
    w_array = np.linspace(w_ini, w_fin, 1000)
    f_of_w = np.zeros(len(w_array))
    for wi in range(len(w_array)):
        integral = 0
        for ti in range(len(t_array)):
            f1 = f_array[ti]*(c.exp(-1j*w_array[wi]*t_array[ti]))
            integral = integral + f1/np.sqrt(2.0*np.pi)
        f_of_w[wi] = integral

    return w_array, f_of_w

x01 = 0.2
x02 = 0.5
x03 = 1.0
t0 = 0.0
delta_t = 0.1
t = 10.0
k = 1.0
alpha = 1.0

t_array1, x_array1, v_array1 = euler_cromer(x01, t0, t, delta_t, k, alpha)
t_array2, x_array2, v_array2 = euler_cromer(x02, t0, t, delta_t, k, alpha)
t_array3, x_array3, v_array2 = euler_cromer(x03, t0, t, delta_t, k, alpha)

wi = -2.0*np.pi/5
wf = 2.0*np.pi/5

#w_array1, f_of_w1 = FourierTransform(t_array1, x_array1, wi, wf)
#
#w_sel = w_array1[np.where(f_of_w1 == max(f_of_w1))[0][1]]
#T_period = 2.0*np.pi/w_sel
#
#fig1 = plt.figure()
#plt.title(r"Simple Harmonic Oscillator $\frac{d^{2}x}{dt^{2}}$ = - $x$")
#plt.plot(t_array1, x_array1, label="Amplitude = {a:0.2f} m, T = {t:0.2f} s".format(a=x01, t=T_period))
#plt.plot(t_array2, x_array2, label="Amplitude = {a:0.2f} m, T = {t:0.2f} s".format(a=x02, t=T_period))
#plt.plot(t_array3, x_array3, label="Amplitude = {a:0.2f} m, T = {t:0.2f} s".format(a=x03, t=T_period))
#plt.xlabel(r"$t (s)$")
#plt.ylabel(r"$x (m)$")
#plt.axhline(y=0,color="k", linewidth=0.7)
#plt.xlim(0,t)
#plt.legend(fontsize="xx-small")
#plt.savefig("harmonic_oscillator.pdf")

t = 100.0
alpha = 3.0

t_array1, x_array1, v_array1 = euler_cromer(x01, t0, t, delta_t, k, alpha)
t_array2, x_array2, v_array2 = euler_cromer(x02, t0, t, delta_t, k, alpha)
t_array3, x_array3, v_array2 = euler_cromer(x03, t0, t, delta_t, k, alpha)

w_array1, f_of_w1 = FourierTransform(t_array1, x_array1, wi, wf)
w_array2, f_of_w2 = FourierTransform(t_array2, x_array2, wi, wf)
w_array3, f_of_w3 = FourierTransform(t_array3, x_array3, wi, wf)

w_sel1 = w_array1[np.where(f_of_w1 == max(f_of_w1))[0][0]]
T_period1 = 2.0*np.pi/w_sel1
w_sel2 = w_array1[np.where(f_of_w2 == max(f_of_w2))[0][0]]
T_period2 = 2.0*np.pi/w_sel2
w_sel3 = w_array1[np.where(f_of_w3 == max(f_of_w3))[0][0]]
T_period3 = 2.0*np.pi/w_sel3

fig2 = plt.figure()
plt.title(r"Anharmonic Oscillator $\frac{d^{2}x}{dt^{2}}$ = - $x^{%0.1f}$"%alpha)
plt.plot(t_array1, x_array1, label="Amplitude = {a:0.2f} m, T = {t:0.2f} s".format(a=x01, t=abs(T_period1)))
plt.plot(t_array2, x_array2, label="Amplitude = {a:0.2f} m, T = {t:0.2f} s".format(a=x02, t=abs(T_period2)))
plt.plot(t_array3, x_array3, label="Amplitude = {a:0.2f} m, T = {t:0.2f} s".format(a=x03, t=abs(T_period3)))
plt.xlabel(r"$t$")
plt.ylabel(r"$x$")
plt.axhline(y=0,color="k", linewidth=0.7)
plt.xlim(0,t)
plt.legend(loc="upper left",fontsize="xx-small")
plt.savefig("anharmonic_oscillator.pdf")

