# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 16:42:28 2021

@author: svenk
"""
import numpy as np
from numpy.fft import fft,ifft,fftfreq
import matplotlib.pyplot as plt
N = 100
L = 2*np.pi
x = np.arange(0,L,L/N)
u = np.sin(x)
u = fft(u)
k = 2*np.pi*fftfreq(N,L/N)
u = k*u *(1j)

u = ifft(u)
plt.plot(x,np.real(u))
plt.plot(x,np.cos(x))