#! /usr/bin/python

# PlotSphericalAlbedo.py
# script to plot Spherical Albedo for 
# k=1.00,1.05,1.15,1.25

# based on material from,

# "Rough surfaces: Is the dark stuff just shadow?"
# Jeffrey N. Cuzzi, Lindsey B. Chambers, Amanda R. Hendrix
# Icarus (2016)

# calling sequence is,
# >python sphericalAlbedo.py
from scipy.integrate import tplquad
from scipy.integrate import quad
import numpy as np
import sys
from math import *
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# define some constants
pi = 4*atan(1.0)

# minimum and maximum values
Smin = 0.0
Smax = 2.0

kmin = 1.0
kmax = 1.5

curlyR = 1.0

# define values for shadowing parameters

shad_vals = np.linspace(Smin,Smax,21)
klam_vals = np.linspace(kmin,kmax,11)



curlyPi = []
cs = []

# define limits
mu0_min = 0.0
mu0_max = 1.0

def phi_min(mu0):
    return 0.0

def phi_max(mu0):
    return 2.0*pi

def mu_min(mu0,phi):
    return 0.0

def mu_max(mu0,phi):
    return 1.0
j=0    
for shad in shad_vals:
  curlyPi.append([])
  for klam in klam_vals:
    def integrand(mu0,phi,mu,args=(klam,shad)):
      alpha = acos(mu*mu0+cos(phi)*sqrt(1.0-mu**2)*sqrt(1.0-mu0**2))
#      return exp(-shad*alpha)*(mu**klam)*(mu0**klam)
      return exp(-shad*sqrt(tan(alpha/2.0)))*(mu**klam)*(mu0**klam)
   

    value, error = tplquad(integrand,mu0_min,mu0_max,phi_min,phi_max,
                           mu_min,mu_max,epsabs=1e-3,epsrel=1e-3)                          
    curlyPi[j].append(curlyR*value*(klam+1)/(pi))
#    print(klam,shad,curlyR*value*(klam+1)/(pi))
  j=j+1  
curlyPi = np.array(curlyPi)
shad_vals = np.array(shad_vals)
klam_vals = np.array(klam_vals)
print(np.min(curlyPi))
print(np.max(curlyPi))

print(curlyPi[:,0])

plt.plot(shad_vals,curlyPi[:,0],'b-',label='k=1.00')
plt.plot(shad_vals,curlyPi[:,1],'r-',label='k=1.05')
plt.plot(shad_vals,curlyPi[:,3],'g-',label='k=1.15')
plt.plot(shad_vals,curlyPi[:,5],'m-',label='k=1.25')
plt.xlabel('Shadowing Parameter')
plt.ylabel(r'$A_s$')
plt.show()

plt.plot(shad_vals,curlyPi[:,0]/curlyPi[0,0],'b-',label='k=1.00')
plt.plot(shad_vals,curlyPi[:,1]/curlyPi[0,1],'r-',label='k=1.05')
plt.plot(shad_vals,curlyPi[:,3]/curlyPi[0,3],'g-',label='k=1.15')
plt.plot(shad_vals,curlyPi[:,5]/curlyPi[0,5],'m-',label='k=1.25')
plt.xlabel('Shadowing Parameter')
plt.ylabel(r'$\frac{A_s}{A_{so}}$',rotation=90)
plt.show()

