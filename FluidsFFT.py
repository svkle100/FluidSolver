# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:46:01 2021

@author: svenk
"""
import pyglet
from pyglet.window import mouse
import numpy as np
from indicator import Indicator
from scipy.fft import fft2,ifft2
from indicator import constrain
 
def stablesolve(n,u,v,uF,vF,visc,dt):

    #Add Force
    u+=dt*uF
    v+=dt*vF
    u0 = u.copy()
    v0 = v.copy()
    #Advect
    for x in range(n):
        for y in range(n):
                x0 = x-n*dt*u0[x,y]
                i0 = int(np.floor(x0))
                s = x0-i0
                i0 %= n
                i1 = (i0+1)%n
                
                y0 = y-n*dt*v0[x,y]
                j0 = int(np.floor(y0))
                t = y0-j0
                j0 %= n
                j1 = (j0+1)%n
                u[x,y] = (1-s)*((1-t)*u0[i0,j0]+t*u0[i0,j1])+\
                    s*((1-t)*u0[i1,j0]+t*u0[i1,j1])
                v[x,y] = (1-s)*((1-t)*v0[i0,j0]+t*v0[i0,j1])+\
                    s*((1-t)*v0[i1,j0]+t*v0[i1,j1])

    u = fft2(u)
    v = fft2(v)
    #calculate Wavenumbers
    freq = np.fft.fftfreq(n,1/n)
    kx = np.zeros(n*n).reshape([n,n])
    ky = np.zeros(n*n).reshape([n,n])
    for i in range(n):
        for j in range(n):
            kx[i,j] = 2*np.pi*freq[i]
            ky[i,j] = 2*np.pi*freq[j]
    #Diffuse
    for i in range(n):
        for j in range(n):
            x = kx[i,j]
            y = ky[i,j]
            r = x*x+y*y
            f = 1/(1+visc*dt*r)
            u[i,j] *= f
            v[i,j] *= f
    #Project
    for i in range(n):
        for j in range(n):
            x = kx[i,j]
            y = ky[i,j]
            r = x*x+y*y
            if r==0.0:
                continue
            u[i,j] = (1-x*x/r)*u[i,j]-x*y/r *v[i,j]
            v[i,j] = -y*x/r *u[i,j] + (1-y*y/r)*v[i,j]
                        
    u = ifft2(u).real
    v = ifft2(v).real
    return u.copy(),v.copy()
 
def substanceSolve(n,s,sS,u,v,diffRate,dissRate,dt):
    #Add Source 
    s+=dt*sS
    s0 = s.copy()
    #Advect
    for x in range(n):
        for y in range(n):
                x0 = x-n*dt*u[x,y]
                i0 = int(np.floor(x0))
                a = x0-i0
                i0 %= n
                i1 = (i0+1)%n
                
                y0 = y-n*dt*v[x,y]
                j0 = int(np.floor(y0))
                b = y0-j0
                j0 %= n
                j1 = (j0+1)%n
                s[x,y] = (1-a)*((1-b)*s0[i0,j0]+b*s0[i0,j1])+\
                    a*((1-b)*s0[i1,j0]+b*s0[i1,j1])
    s = fft2(s)
    #calculate Wavenumbers
    freq = np.fft.fftfreq(n,1/n)
    kx = np.zeros(n*n).reshape([n,n])
    ky = np.zeros(n*n).reshape([n,n])
    for i in range(n):
        for j in range(n):
            kx[i,j] = 2*np.pi*freq[i]
            ky[i,j] = 2*np.pi*freq[j]
    #Diffuse
    for i in range(n):
        for j in range(n):
            x = kx[i,j]
            y = ky[i,j]
            r = x*x+y*y
            f = 1/(1+diffRate*dt*r)
            s[i,j] *= f
    s = ifft2(s).real
    #Dissipation
    s/= (1+dt*dissRate)
    return s

def updateGenerator(n,u0,v0,s0,indicators,rectangles,res):
    u = np.zeros([n,n])
    v = np.zeros([n,n])
    s = np.zeros([n,n])
    def update(dt):
        nonlocal u,v,s,u0,v0,s0
        u,v = stablesolve(n,u,v,u0,v0,0.001,dt)
        s = substanceSolve(n,s,s0,u,v,0.00001,0.00001,dt)
        maxVal = max(np.amax(np.abs(u)),np.max(np.abs(v)))
        for i in range(n):
            for j in range(n):
                #indicators[i,j].update(u[i,j],v[i,j],maxVal)
                c = constrain(s[i,j],0,255)
                rectangles[i,j].color=(c,c,c)
        u0 /=2
        v0 /=2
        s0 /=2
    return update
         
        


    
    
def main():
    n= 40
    u0 = np.zeros([n,n])
    v0 = np.zeros([n,n])
    s0 = np.zeros([n,n])
    size = 600
    res = size/n
    main_batch = pyglet.graphics.Batch()
    window = pyglet.window.Window(size, size)
    rectangles = np.array([[pyglet.shapes.Rectangle(i*res,j*res,res,res,color=(0,0,0),batch=main_batch) for j in range(n)] 
                            for i in range(n)])
    
    @window.event
    def on_draw():
        window.clear()
        main_batch.draw()
    
    indicators = np.array([[Indicator(res,(i+.5)*res,(j+.5)*res, (i+.5)*res, (j+.5)*res,batch=main_batch) for j in range(n)] 
                            for i in range(n)])
    update = updateGenerator(n,u0,v0,s0, indicators,rectangles, res)
    pyglet.clock.schedule_interval(update, 1/120)

    @window.event
    def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
        nonlocal u0,v0,s0
        if buttons & mouse.LEFT:
            x = constrain(x,0,size-1)
            y = constrain(y,0,size-1)
            xpos = int(x/res)
            ypos = int(y/res)
            u0[xpos,ypos] += dx*3
            v0[xpos,ypos] += dy*3


        if buttons & mouse.RIGHT:
            x = constrain(x,0,size-1)
            y = constrain(y,0,size-1)
            xpos = int(x/res)
            ypos = int(y/res)
            s0[xpos,ypos] += 2000
    @window.event
    def on_close():
        print("exiting")
        pyglet.clock.unschedule(update)
        pyglet.app.exit()
    pyglet.app.run()



if __name__ == '__main__':
    main()


     

    
