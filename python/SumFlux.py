#!/usr/bin/env python2
# SumFlux.py

import sys
import numpy as np
import math

flux1 = sys.argv[1]    # current Flux file  (flux.out)
flux2 = sys.argv[2]    # running flux sum   (flux_tot)
IoverF1 = sys.argv[3]  # current  i/f file  (if.out)
IoverF2 = sys.argv[4]  # running i/f sum (if_sum.out) 
iter = sys.argv[5]     # current iteration
maxiter = sys.argv[6]  # maximum number of iterations

# output filenames
flux3 = 'flux_new.out'
IoverF3 = 'if_new.out'

# define some vectors
order = []
flux_total = []
flux_solar = []
flux_satshine = []
fluxsum_total = []
fluxsum_solar = []
fluxsum_satshine = []
fluxsumsq_total = []
fluxsumsq_solar = []
fluxsumsq_satshine = []


albedo = []
IoverF_total = []
IoverF_solar = []
IoverF_satshine = []
IoverFsum_total = []
IoverFsum_solar = []
IoverFsum_satshine = []
IoverFsumsq_total = []
IoverFsumsq_solar = []
IoverFsumsq_satshine = []

# open files
f1 = open(flux1,'r')
f2 = open(flux2,'r')
f3 = open(flux3,'w')

g1 = open(IoverF1,'r')
g2 = open(IoverF2,'r')
g3 = open(IoverF3,'w')

# load columns of flux files

# save header

header1 = f1.readline() 
header2 = f1.readline() 
header3 = f1.readline() 
for line in f1:
    line = line.strip()
    columns = line.split()
    order.append(int(columns[0]))
    flux_solar.append(float(columns[1]))
    flux_satshine.append(float(columns[2]))    
    flux_total.append(float(columns[3]))

header1 = f2.readline() 
header2 = f2.readline() 
for line in f2:
    line = line.strip()
    columns = line.split()
    fluxsum_solar.append(float(columns[1]))
    fluxsum_satshine.append(float(columns[2]))    
    fluxsum_total.append(float(columns[3]))
    fluxsumsq_solar.append(float(columns[4]))
    fluxsumsq_satshine.append(float(columns[5]))    
    fluxsumsq_total.append(float(columns[6]))
    
    
fluxnew_total = np.array(fluxsum_total) + np.array(flux_total)
fluxnew_solar = np.array(fluxsum_solar) + np.array(flux_solar)
fluxnew_satshine = np.array(fluxsum_satshine) + np.array(flux_satshine)

fluxnewsq_total = np.array(fluxsumsq_total) + np.array(flux_total)**2
fluxnewsq_solar = np.array(fluxsumsq_solar) + np.array(flux_solar)**2
fluxnewsq_satshine = np.array(fluxsumsq_satshine) + np.array(flux_satshine)**2


if (iter == maxiter):
      fluxnew_total=fluxnew_total/float(maxiter)
      fluxnew_solar=fluxnew_solar/float(maxiter)
      fluxnew_satshine=fluxnew_satshine/float(maxiter)
      fluxnewsq_total=fluxnewsq_total/float(maxiter)
      fluxnewsq_solar=fluxnewsq_solar/float(maxiter)
      fluxnewsq_satshine=fluxnewsq_satshine/float(maxiter)

      
# write output to file
f3.write(header1)  # copy headers
f3.write(header2)

for i in range(len(fluxnew_total)):
    f3.write("%2.0d\t%11.6e\t%11.6e\t%11.6e\t%11.6e\t%11.6e\t%11.6e\n" % (order[i], \
               fluxnew_solar[i], \
               fluxnew_satshine[i], \
               fluxnew_total[i], \
               fluxnewsq_solar[i], \
               fluxnewsq_satshine[i], \
               fluxnewsq_total[i]))
#     str1 = str(order[i]) + '\t' + str(round(fluxnew_solar[i],15)) + '\t' + \
#     str(round(fluxnew_satshine[i],15)) + '\t' + str(round(fluxnew_total[i],15)) + \
#     '\t' + str(round(fluxnewsq_solar[i],15)) + '\t' + str(round(fluxnewsq_satshine[i],15)) + \
#     '\t' + str(round(fluxnewsq_total[i],15)) + '\n'
#     f3.write(str1) 

f1.close()
f2.close()
f3.close()

header1 = g1.readline() 
header2 = g1.readline()
for line in g1:
    line = line.strip()
    columns = line.split()
    albedo.append(float(columns[0]))
    IoverF_total.append(float(columns[1]))
    IoverF_solar.append(float(columns[2]))
    IoverF_satshine.append(float(columns[3]))

header1 = g2.readline() 
header2 = g2.readline()
for line in g2:
    line = line.strip()
    columns1 = line.split()
    IoverFsum_total.append(float(columns1[1]))
    IoverFsum_solar.append(float(columns1[2]))
    IoverFsum_satshine.append(float(columns1[3]))
    IoverFsumsq_total.append(float(columns1[4]))
    IoverFsumsq_solar.append(float(columns1[5]))
    IoverFsumsq_satshine.append(float(columns1[6]))
    
IoverFnew_total = np.array(IoverFsum_total) + np.array(IoverF_total)
IoverFnew_solar = np.array(IoverFsum_solar) + np.array(IoverF_solar)
IoverFnew_satshine = np.array(IoverFsum_satshine) + np.array(IoverF_satshine)

IoverFnewsq_total = np.array(IoverFsumsq_total) + np.array(IoverF_total)**2
IoverFnewsq_solar = np.array(IoverFsumsq_solar) + np.array(IoverF_solar)**2
IoverFnewsq_satshine = np.array(IoverFsumsq_satshine) + np.array(IoverF_satshine)**2

if (iter == maxiter):
      IoverFnew_total=IoverFnew_total/float(maxiter)
      IoverFnew_solar=IoverFnew_solar/float(maxiter)
      IoverFnew_satshine=IoverFnew_satshine/float(maxiter)
      IoverFnewsq_total=IoverFnewsq_total/float(maxiter)
      IoverFnewsq_solar=IoverFnewsq_solar/float(maxiter)
      IoverFnewsq_satshine=IoverFnewsq_satshine/float(maxiter)

g3.write(header1)
g3.write(header2)
for k in range(len(IoverFnew_total)):
    g3.write("%4.3f\t%11.6e\t%11.6e\t%11.6e\t%11.6e\t%11.6e\t%11.6e\n" % (albedo[k],IoverFnew_total[k], \
               IoverFnew_solar[k],IoverFnew_satshine[k],IoverFnewsq_total[k],IoverFnewsq_solar[k],IoverFnewsq_satshine[k]))

g1.close()
g2.close()
g3.close()
