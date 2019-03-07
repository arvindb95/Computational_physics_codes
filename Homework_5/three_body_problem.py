import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
### Constants ###########

t_final = 10
dt = 10**(-2.0)

M_sun = 1.0
M_jupiter = 1.0 # 10**(-4)
M_earth = 10**(-6)
# Initial position
x0_sun = 0
y0_sun = 0
x0_jupiter = 5.2
y0_jupiter = 0
x0_earth = 1.0
y0_earth = 0
# Initial velocities
vx0_sun = 0
vy0_sun = 0
vx0_jupiter = 0
vy0_jupiter = 2*np.pi*13.07/30.0
vx0_earth = 0
vy0_earth = 2*np.pi

G = 4*(np.pi**2.0)/M_sun

##### Initializing arrays

time = np.arange(0, t_final + dt, dt)

def acceleration_1(Ms, Mj, pos1s, pos1j, r1s, r1j):

    return - (G*Ms*pos1s/(r1s**3.0)) - (G*Mj*pos1j/(r1j**3.0))

def three_body_problem(time_array):
    
    x_sun = np.zeros(len(time_array))
    x_sun[0] = x0_sun
    y_sun = np.zeros(len(time_array))
    y_sun[0] = y0_sun

    x_jupiter = np.zeros(len(time_array))
    x_jupiter[0] = x0_jupiter
    y_jupiter = np.zeros(len(time_array))
    y_jupiter[0] = y0_jupiter

    x_earth = np.zeros(len(time_array))
    x_earth[0] = x0_earth
    y_earth = np.zeros(len(time_array))
    y_earth[0] = y0_earth

    vx_sun = np.zeros(len(time_array))
    vx_sun[0] = vx0_sun
    vy_sun = np.zeros(len(time_array))
    vy_sun[0] = vy0_sun

    vx_jupiter = np.zeros(len(time_array))
    vx_jupiter[0] = vx0_jupiter
    vy_jupiter = np.zeros(len(time_array))
    vy_jupiter[0] = vy0_jupiter

    vx_earth = np.zeros(len(time_array))
    vx_earth[0] = vx0_earth
    vy_earth = np.zeros(len(time_array))
    vy_earth[0] = vy0_earth

    for i in range(len(time_array)-1):
        r_sj = np.sqrt(((x_sun[i] - x_jupiter[i])**2.0) + ((y_sun[i] - y_jupiter[i])**2.0))
        r_se = np.sqrt(((x_sun[i] - x_earth[i])**2.0) + ((y_sun[i] - y_earth[i])**2.0))
        r_je = np.sqrt(((x_earth[i] - x_jupiter[i])**2.0) + ((y_earth[i] - y_jupiter[i])**2.0))
        
        vx_sun[i+1] = vx_sun[i] + dt*acceleration_1(M_earth, M_jupiter, (x_sun[i] - x_earth[i]), (x_sun[i] - x_jupiter[i]), r_se, r_sj)
        vy_sun[i+1] = vy_sun[i] + dt*acceleration_1(M_earth, M_jupiter, (y_sun[i] - y_earth[i]), (y_sun[i] - y_jupiter[i]), r_se, r_sj)
        vx_jupiter[i+1] = vx_jupiter[i] + dt*acceleration_1(M_sun, M_earth, (x_jupiter[i] - x_sun[i]), (x_jupiter[i] - x_earth[i]), r_sj, r_je)
        vy_jupiter[i+1] = vy_jupiter[i] + dt*acceleration_1(M_sun, M_earth, (y_jupiter[i] - y_sun[i]), (y_jupiter[i] - y_earth[i]), r_sj, r_je)
        vx_earth[i+1] = vx_earth[i] + dt*acceleration_1(M_sun, M_jupiter, (x_earth[i] - x_sun[i]), (x_earth[i] - x_jupiter[i]), r_se, r_je)
        vy_earth[i+1] = vy_earth[i] + dt*acceleration_1(M_sun, M_jupiter, (y_earth[i] - y_sun[i]), (y_earth[i] - y_jupiter[i]), r_se, r_je)

        x_sun[i+1] = x_sun[i] + dt*vx_sun[i+1]
        y_sun[i+1] = y_sun[i] + dt*vy_sun[i+1]
        x_jupiter[i+1] = x_jupiter[i] + dt*vx_jupiter[i+1]
        y_jupiter[i+1] = y_jupiter[i] + dt*vy_jupiter[i+1]
        x_earth[i+1] = x_earth[i] + dt*vx_earth[i+1]
        y_earth[i+1] = y_earth[i] + dt*vy_earth[i+1]

    return x_sun, y_sun, x_jupiter, y_jupiter, x_earth, y_earth

x_sun, y_sun, x_jupiter, y_jupiter, x_earth, y_earth = three_body_problem(time)

fig1 = plt.figure(figsize=(6,6))

plt.title(r"Three body problem : M$_{Sun}$ = "+str(M_sun)+" M$_{Jupiter}$ = "+str(M_jupiter)+" M$_{Earth}$ = "+str(M_earth))
plt.plot(x_sun, y_sun, "y", label="Sun's path")
plt.plot(x_sun[0],y_sun[0],"yo",label="Sun start")
plt.plot(x_earth, y_earth, "b", label="Earth's path")
plt.plot(x_earth[0], y_earth[0], "bo", label="Earth start")
plt.plot(x_jupiter, y_jupiter, "r", label="Jupiter's path")
plt.plot(x_jupiter[0], y_jupiter[0], "ro", label="Jupiter start")
plt.xlabel("1AU")
plt.ylabel("1AU")
plt.legend(loc="upper right")

plt.savefig("Alpha_"+str(M_jupiter/M_sun)+".pdf")

fig, ax = plt.subplots()

axsun, = ax.plot([],[],"yo",label="Sun")
axjupiter, = ax.plot([],[],"ro", label="Jupiter")
axearth, = ax.plot([],[], "bo", label="Earth")
axsunpath, = ax.plot([],[],"y-",linewidth=1)
axjupiterpath, = ax.plot([],[],"r-",linewidth=1)
axearthpath, = ax.plot([],[], "b-", linewidth=1)


def init():
    ax.set_title(r"Three body problem : M$_{Sun}$ = "+str(M_sun)+" M$_{Jupiter}$ = "+str(M_jupiter)+" M$_{Earth}$ = "+str(M_earth))
    ax.set_xlim([min(min(x_sun), min(x_earth), min(x_jupiter)) - 2 ,max(max(x_sun), max(x_earth), max(x_jupiter)) + 2])
    ax.set_ylim([min(min(y_sun), min(y_earth), min(y_jupiter)) - 2 ,max(max(y_sun), max(y_earth), max(y_jupiter)) + 2])
    ax.set_xlabel("1AU")
    ax.set_ylabel("1AU")
    axsun.set_data([],[])
    axjupiter.set_data([],[])
    axearth.set_data([],[])
    return axsun, axjupiter, axearth
    

def animate(frame):
    print("time : ", frame)
    axsun.set_data([x_sun[frame]],[y_sun[frame]])
    axjupiter.set_data([x_jupiter[frame]], [y_jupiter[frame]])
    axearth.set_data([x_earth[frame]], [y_earth[frame]])
    axsunpath.set_alpha(1 - 0.7*frame/len(time))
    axsunpath.set_data([x_sun[:frame]], [y_sun[:frame]])
    axjupiterpath.set_alpha(1 - 0.7*frame/len(time))
    axjupiterpath.set_data([x_jupiter[:frame]],[y_jupiter[:frame]])
    axearthpath.set_alpha(1 - 0.7*frame/len(time))
    axearthpath.set_data([x_earth[:frame]], [y_earth[:frame]])
    return axsun, axjupiter, axearth, axsunpath, axjupiterpath, axearthpath


anim = animation.FuncAnimation(fig, animate, init_func = init, frames = len(time), interval=1, blit=True, repeat=False)
anim.save("Alpha_"+str(M_jupiter/M_sun)+".mp4",fps=25,dpi=180)
plt.show()

