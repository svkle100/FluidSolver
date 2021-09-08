# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:49:48 2021

@author: svenk
"""
import pyglet
import numpy as np

def constrain(val, min_val, max_val):
    return min(max_val, max(min_val, val))
def interpolate(value, a, b,newa,newb):
    if b-a==0:
        return 0
    return value/(b-a)*(newb-newa)+newb

class Indicator(pyglet.shapes.Line):
    def __init__(self,res, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.res = res
        self.min = 0.001
        
    def update(self,vx,vy,maxValue):
        vx *=(self.res/2.2)/maxValue**2
        vy *=(self.res/2.2)/maxValue**2
        norm = np.linalg.norm([vx,vy])
        self.x2 = self.x + vx/norm*self.res/2
        self.y2 = self.y + vy/norm*self.res/2
