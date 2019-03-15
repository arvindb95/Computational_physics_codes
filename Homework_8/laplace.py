import numpy as np

import matplotlib.pyplot as plt

h = 10**(-1)

x = np.arange(0,2+h,h)
y = np.arange(0,1+h,h)

V = np.zeros((len(x), len(y)))

# Boundary conditions

V[np.where(x==2),:] = 1
V[np.where(x==0),:] = -1
V[:,np.where(y==0)] = -1
V[:,np.where(y==1)] = 1

V_new = np.copy(V)

steps = 1000
s = 0

while(s <= steps):
    s += 1
    for i in range(1,(len(x)-1)):
        for j in range(1,(len(y)-1)):
            V_new[i,j] = (1.0/4.0)*(V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1])
            V = V_new

fig = plt.figure(figsize=(10,20))

plt.title("Solution to Laplace equation V(x,y)")
plt.grid(linestyle="--")
for i in range(len(x)):
    for j in range(len(y)):
        plt.scatter(x[i],y[j], c=V[i,j], cmap="viridis", vmin=-1, vmax =1)
plt.xlabel("x")
plt.ylabel("y")
plt.hlines(y=-0.02,xmin=-0.02,xmax=2.02,colors="k")
plt.hlines(y=1.02,xmin=-0.02,xmax=2.02,colors="k")
plt.vlines(x=-0.02,ymin=-0.02,ymax=1.02,colors="k")
plt.vlines(x=2.02,ymin=-0.02,ymax=1.02,colors="k")
plt.colorbar()        
plt.show()


