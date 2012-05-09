#coding:utf-8

#===================
#LIBRARIES
#===================

import numpy as np
import scipy as sp
import pylab as pl
import math

import copy
import random
import array
import matplotlib.pyplot as plt

#===================
#OUTLINE
#===================

#<<Relaxation method with Ramsey Model>>
"""
1.derive final state
2.start loop and check defference
3.plot
"""

#===================
#VARIABLES
#===================

#parameters

alpha = 0.3
delta = 0.05
rho = 0.02
nPop = 0.01
theta = (delta + rho) / (alpha * (delta + nPop)) #it depends the hypothesis
kss = ((delta+rho)/alpha)**(1/(alpha-1)) #steady state of k
css = kss**alpha-(nPop+delta)*kss #steady state of c

print "steady state"
print "k =",kss
print "c =",css
print theta

#settings
m = 100 #mesh
n = 2 #number of variables
k1 = 3 #initial value of k
c1 = 0.5 #initial value of c
maxkount = 10000 #max kount
maxdif = 0.0000001 #max difference

#containers
y = np.zeros((m+1,n))
mm = range(0,m+1)
kk = np.zeros(m+1)
cc = np.zeros(m+1)
dd = np.zeros(m+1)

#===================
#FUNCTIONS
#===================



#===================
#MAIN CALCULATION
#===================

#--------1.derive final state--------

y[0][0] = k1
y[m][0] = kss
y[0][1] = c1
y[m][1] = css

for i in range(1,m):
    y[i][0] = (kss-k1) * i / m
    y[i][1] = (css-c1) * i / m


#--------2.start loop--------

kount = 1 #loop count
dif = 1000000 #check variable to converge

while kount<maxkount and dif>maxdif:
    dif = 0
    y_ch = copy.deepcopy(y) #for check
    #create converge curve from guess(y)
    for i in range(0,m):
        def_k = (y[i][0]**(alpha+1)) - y[i][0]*y[i][1] - (nPop + delta)*(y[i][0]**2)
        def_c = (y[i][1]**2) / theta * (alpha*y[i][0]**(alpha-1) - (delta+rho))
        y[i+1][0] = y[i][0] + def_k
        y[i+1][1] = y[i][1] + def_c
    #compute differ btw y and y_ch
    dif_m = y - y_ch
    """
    if kount >0:
        for i in range(0,m+1):
            kk[i] = y[i][0]
        plt.subplot(1,1,1)
        plt.plot(mm,kk,'x')
        plt.title('state of y')
        
        for i in range(0,m+1):
            cc[i] = y_ch[i][0]
        plt.subplot(1,1,1)
        plt.plot(mm,cc,'x')
        plt.title('y_ch')
        
        for i in range(0,m+1):
            dd[i] = 0.5 * (y[i][0] + y_ch[i][0])
        plt.subplot(1,1,1)
        plt.plot(mm,dd,'x')
        plt.title('state of c')
        
        plt.show()
    """

    for i in range(0,2):
        for j in range(0,m+1):
            dif = math.fabs(dif_m[j][i])
            y[j][i] = 0.5 * (y[j][i] + y_ch[j][i])
    
    print kount,"th loop is finished"
    print "error is",dif
    kount = kount + 1
    
    print y[50][0], y[50][1], def_k, def_c


#plotting graph

for i in range(0,m+1):
    kk[i] = y[i][0]
plt.subplot(2,1,1)
plt.plot(mm,kk,'x')
plt.title('state of k')

for i in range(0,m+1):
    cc[i] = y[i][1]
plt.subplot(2,1,2)
plt.plot(mm,cc,'x')
plt.title('state of c')
plt.show()

plt.subplot(1,1,1)
plt.plot(kk,cc,"x")
plt.title("k-c")
plt.show()
