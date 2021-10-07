import numpy as np
from FluidsFFT import stablesolve
import matplotlib.pyplot as plt
plt.figure()
n = 128
p = 80
d = 0.05
u = np.zeros([n,n])
v = np.zeros([n,n])
visc = 10**-8
t = 1
dt = .01
for x in range(n):
    for y in range(n):
        if (y+.5)/n<=.5:
            u[x,y] = np.tanh(p*((y+.5)/n-.25))
        else:
            u[x,y] = np.tanh(p*(.75-(y+.5)/n))
        v[x,y] = d * np.sin(2*np.pi*((x+.5)/n+.25))
x,y = np.meshgrid(np.linspace(0,1,n),np.linspace(0,1,n))
plt.quiver(x,y,u,v,angles='xy')
plt.figure()
plt.plot(np.linspace(0,1,n),v[:,0])
for i in range(int(t/dt)):
    u,v = stablesolve(n, u, v, np.zeros_like(u),np.zeros_like(v),visc,dt)
    print(i)
dudy = np.zeros([n,n])
dvdx = np.zeros([n,n])
for x in range(n):
    for y in range(n):
        dudy[x,y] = (u[x,(y+1)%n] - u[x,(y-1)%n])/(2*1/n)
        dvdx[x,y] = (u[(x+1)%n,y] - u[(x-1)%n,y])/(2*1/n)
vorticity = dvdx-dudy
plt.figure()
plt.contour(np.linspace(0,1,n),np.linspace(0,1,n),vorticity,np.array(range(-7,8,1)))
