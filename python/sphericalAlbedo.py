#! /usr/bin/python

# sphericalAlbedo.py
# script to calculate and plot the spherical albedo
# based on the Mineart scattering law 

# based on material from,

# "Rough surfaces: Is the dark stuff just shadow?"
# Jeffrey N. Cuzzi, Lindsey B. Chambers, Amanda R. Hendrix
# Icarus (2016)

# calling sequence is,
# >python sphericalAlbedo.py
from scipy.integrate import tplquad
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

albedo_min = 0.0
albedo_max = 1.4

scale = (kmax-kmin)/(Smax-Smin)

curlyR = 1.0

# define values for shadowing parameters

shad_vals = np.linspace(Smin,Smax,21)
klam_vals = np.linspace(kmin,kmax,21)



curlyPi = []

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
    print(klam,shad,curlyR*value*(klam+1)/(pi))
  j=j+1  
curlyPi = np.array(curlyPi)
shad_vals = np.array(shad_vals)
klam_vals = np.array(klam_vals)
print(np.min(curlyPi))
print(np.max(curlyPi))
# print(curlyPi.shape)
# print(shad_vals.shape)
# print(klam_vals.shape)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(klam_vals, shad_vals, curlyPi,color='b')
# plt.show()

# plt.figure()
# CS = plt.contour(klam_vals, shad_vals, curlyPi, 10,
#                  linewidths=np.arange(.5, 4, .5),
#                  colors=('r', 'green', 'blue', (1, 1, 0), '#afeeee', '0.5')
#                  )
# plt.xlim(0.5,1.55)
# plt.ylim(0,2.0)
# plt.clabel(CS, inline=1, fontsize=9)
# plt.show()

if(True):
	fig=plt.figure(figsize=(8,8))
	im = plt.imshow(curlyPi, interpolation='bilinear', origin='lower',
					cmap=cm.gray, extent=(kmin,kmax, Smin,Smax),aspect=scale)
	levels = np.arange(albedo_min, albedo_max, 0.05)
	CS = plt.contour(curlyPi, levels,extent=(kmin, kmax, Smin, Smax))
	plt.clabel(CS,inline=1,fmt='%2.2f',fontsize=10)

	# make a colorbar for the contour lines
	CB = plt.colorbar(CS, shrink=1.0, extend='both')

	plt.title('Lines with colorbar')
	#plt.hot()  # Now change the colormap for the contour lines and colorbar
	plt.flag()

	# We can still add a colorbar for the image, too.
	CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)

	# This makes the original colorbar look a bit out of place,
	# so let's improve its position.

	l, b, w, h = plt.gca().get_position().bounds
	ll, bb, ww, hh = CB.ax.get_position().bounds
	CB.ax.set_position([ll, b + 0.1*h, ww, h*0.8])
	plt.xlabel('Minneart parameter')
	plt.ylabel('Shadowing parameter')
	plt.show()

# wireframe plot

klam_array = np.tile(klam_vals,(21,1))
shad_array = np.transpose(np.tile(shad_vals,(21,1)))

fig = plt.figure(figsize=(6,6), facecolor='w')
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(klam_array, shad_array, curlyPi, rstride=1, cstride=1)
ax.set_xlabel("Minneart parameter",linespacing=3.2)
ax.set_ylabel("Shadowing parameter",linespacing=3.2)
ax.set_zlabel('Spherical Albedo', linespacing=3.4)
ax.view_init(45, 145)
ax.dist = 10
plt.show()

# rotate viewpoint
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_wireframe(klam_array, shad_array, curlyPi, rstride=1, cstride=1)
# ax.set_xlabel("Minneart parameter")
# ax.set_ylabel("Shadowing parameter")
# for angle in range(0, 360):
#     ax.view_init(30, angle)
#     plt.draw()
#     plt.pause(.05)


fig = plt.figure(figsize=(6,6), facecolor='w')
ax = fig.add_subplot(111, projection='3d')
cset = ax.contour(klam_array, shad_array, curlyPi, levels=levels, cmap=cm.coolwarm)
ax.clabel(cset, fontsize=9,fmt='%2.2f', inline=1)
ax.set_xlabel("Minneart parameter")
ax.set_ylabel("Shadowing parameter")
CB = plt.colorbar(cset, shrink=0.5)
ax.view_init(45, 145)
ax.dist = 10
plt.show()