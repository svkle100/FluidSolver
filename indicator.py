# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:49:48 2021

@author: svenk
"""
from scipy.interpolate import interp1d
import pyglet
def constrain(val, min_val, max_val):
    return min(max_val, max(min_val, val))

class Indicator(pyglet.shapes.Line):
    def __init__(self,n, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n = n
        self.m = interp1d([-.1,.1],[-1,1])
        self.factor = 500
        
    def update(self,vx,vy):
        vx = constrain(vx,-.1,.1)
        vy = constrain(vy,-.1,.1)

        self.x2 = self.x + constrain(self.m(vx)*self.factor,-self.n,self.n)
        self.y2 = self.y + constrain(self.m(vy)*self.factor,-20,self.n)
