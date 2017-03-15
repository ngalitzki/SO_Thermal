import sys
import numpy as np
import scipy.constants as PC
import os
import matplotlib.pyplot as plt
import pdb

pi=np.pi
h=6.626e-34
c=3e8
k_b=1.38e-23

def rho(f, T):
    B=pi*(2*h*f**3/c**2)/(np.exp(h*f/(k_b*T))-1)
    return B

def filter(f_in, T_in, f):
    i = f_in.tolist().index(f)
    T = T_in[i]
    f_in=np.array(f_in)    
    return T
