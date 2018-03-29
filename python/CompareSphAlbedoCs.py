#! /usr/bin/python

# compareSphAlbedoCs.py
# script to compare shadowed lambert phase function to 
# spherical albedo for k=1.00,1.05,1.15,1.25

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

shad_vals = np.linspace(Smin,Smax,11)
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

#******compute Cs values

for shad in shad_vals:

     def phaseFn(alpha,args=(shad,)):
        pi = 4*atan(1.0)
        return (exp(-shad*sqrt(tan(alpha/2.0)))*(sin(alpha)+(pi-alpha)*cos(alpha))*sin(alpha))

     value, error = quad(phaseFn,0,pi)
     print(shad, value)
     cs.append(3*pi/(4*value))
 
cs=np.array(cs) 
plt.plot(shad_vals,1-cs*curlyPi[:,0])
plt.xlabel('Shadowing parameter')
plt.ylabel(r'1 - $c^{}_{SLH}A_s$')
#plt.ylim(0.99999,1.000001)
#plt.savefig('CompareCs_with_As.pdf',dpi=800,format='pdf')
plt.show()

print(cs)
