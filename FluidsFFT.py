# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:46:01 2021

@author: svenk
"""
import pyglet
from pyglet.window import mouse
import numpy as np
from indicator import Indicator
from scipy.fft import fft,ifft
from indicator import constrain
 
def stablesolve(n,u,v,u1,v1,visc,dt):
    u0 = u1.copy()
    v0 = v1.copy()
    for i in range(n**2):
            u[i] += dt*u0[i]
            u0[i] = u[i]
            v[i] += dt*v0[i]
            v0[i] = v[i]
    for x in range(n):
        for y in range(n):
                x0 = x-n*dt*u0[x+n*y]
                i0 = int(np.floor(x0))
                s = x0-i0
                i0 %= n
                i1 = (i0+1)%n
                
                y0 = y-n*dt*v0[x+n*y]
                j0 = int(np.floor(y0))
                t = y0-j0
                j0 %= n
                j1 = (j0+1)%n
                u[x+n*y] = (1-s)*((1-t)*u0[i0+n*j0]+t*u0[i0+n*j1])+\
                    s*((1-t)*u0[i1+n*j0]+t*u0[i1+n*j1])
                v[x+n*y] = (1-s)*((1-t)*v0[i0+n*j0]+t*v0[i0+n*j1])+\
                    s*((1-t)*v0[i1+n*j0]+t*v0[i1+n*j1])
    for i in range(n):
        for j in range(n):
            u0[i+(n+2)*j] = u[i+n*j] 
            v0[i+(n+2)*j] = v[i+n*j]
    u0 = fft(u0)
    v0 = fft(v0)
    U = [0,0]; V = [0,0]
    for i in range(0,n+1,2):
        x = 0.5*i
        for j in range(n):
            y = j if j<=n/2 else j-n
            r = x*x+y*y
            if r==0.0:
                continue
            f = np.exp(-r*dt*visc)
            U[0] = u0[i +(n+2)*j]; V[0] = v0[i +(n+2)*j]
            U[1] = u0[i+1+(n+2)*j]; V[1] = v0[i+1+(n+2)*j]
            u0[i +(n+2)*j] = f*( (1-x*x/r)*U[0]-x*y/r *V[0] )
            u0[i+1+(n+2)*j] = f*( (1-x*x/r)*U[1]-x*y/r *V[1] )
            v0[i+ (n+2)*j] = f*( -y*x/r *U[0] + (1-y*y/r)*V[0] )
            v0[i+1+(n+2)*j] = f*( -y*x/r *U[1] + (1-y*y/r)*V[1] )
                        
    u0 = ifft(u0).real
    v0 = ifft(v0).real
    f = 1.0/(n**2);
    for i in range(n):
        for j in range(n):
            u[i+n*j] = f*u0[i+(n+2)*j]
            v[i+n*j] = f*v0[i+(n+2)*j]
    return u,v
        

def updateGenerator(n,u0,v0,indicators,res):
    u = np.ones(n**2)*10
    v = np.zeros(n**2)
    def update(dt):
        nonlocal u,v,u0,v0
        u,v = stablesolve(n,u,v,u0,v0,0.001,dt)
        for ind in indicators:
            ind.update(u[int(ind.x/res)],v[int(ind.y/res)])
        u0/=2
        v0/=2
    return update
         
        


    
    
def main():
    n= 30
    u0 = np.zeros(n*(n+2))
    v0 = np.zeros(n*(n+2))
    size = 600
    res = size/n
    window = pyglet.window.Window(size, size)
    main_batch = pyglet.graphics.Batch()
    @window.event
    def on_draw():
        window.clear()
        main_batch.draw()
    
    indicators = [Indicator(n,i*res,j*res, i*res, j*res,batch=main_batch) for i in range(n) for j in range(n)]
    update = updateGenerator(n,u0,v0, indicators, res)
    pyglet.clock.schedule_interval(update, 1/4)

    @window.event
    def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
        nonlocal u0,v0
        if buttons & mouse.LEFT:
            x = constrain(x,0,size-1)
            y = constrain(y,0,size-1)
            xpos = int(x/res)
            ypos = int(y/res)
            u0[xpos+(n+2)*ypos] += dx*100
            v0[xpos+(n+2)*ypos] += dy*100
    @window.event
    def on_close():
        print("exiting")
        pyglet.clock.unschedule(update)
        pyglet.app.exit()
    pyglet.app.run()



if __name__ == '__main__':
    main()


     

    
