# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:49:48 2021

@author: svenk
"""
import pyglet
import numpy as np
def constrain(val, min_val, max_val):
    return min(max_val, max(min_val, val))

class Indicator(pyglet.shapes.Line):
    def __init__(self,res, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.res = res
        self.factor = 500
        
    def update(self,vx,vy):
        norm =np.linalg.norm([vx,vy]) 
        if norm>self.res/2:
            vx = vx/norm*self.res/2
            vy = vy/norm*self.res/2
        self.x2 = self.x + vx
        self.y2 = self.y + vy
