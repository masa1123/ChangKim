#coding:utf-8

#===================
#LIBRARIES
#===================

import numpy as np
import scipy as sp
import pylab as pl
import math

from copy import deepcopy
import random
import array
import matplotlib.pyplot as plt
from pylab import *

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
#print theta

#settings
m = 1 #mesh
n = 2 #number of variables
k1 = 3 #initial value of k
c1 = 0.5 #initial value of c
maxkount = 10000 #max kount
maxdif = 0.0000001 #max difference

#settings of plotting
unit_plot_x = 0.2 #unit of x plotting
unit_plot_y = 0.04 #unit of y plotting
max_x = 20 #max size of x in graph
max_y = 2 #max size of y in grapf
min_x = unit_plot_x #min size of x in graph
min_y = unit_plot_y #min size of y in graph
dyn_arr = 0.1 #dynamics of arrow
max_arr = 0.2 #max length of arrow


#containers
y = np.zeros((m+1,n))
mm = range(0,m+1)
kk = np.zeros(m+1)
cc = np.zeros(m+1)
dd = np.zeros(m+1)
arr = np.zeros(2)

#===================
#FUNCTIONS
#===================

#derive angle
def vector_angle(p,q):
    inner_p = np.dot(p,q) #inner product
    p_len = np.linalg.norm(p)
    q_len = np.linalg.norm(q)
    cos = inner_p / (p_len * q_len)
    if p[0] > 0:
        d_angle = math.acos(cos)
    else:
        d_angle = math.acos(cos) + pi
    return d_angle
    

#derive change angle
def ram_arrow(k_ini,c_ini):
    y[0][0] = k_ini
    y[m][0] = kss
    y[0][1] = c_ini
    y[m][1] = css

    for i in range(1,m):
        y[i][0] = (kss-k1) * i / m
        y[i][1] = (css-c1) * i / m


    #--------2.start loop--------

    kount = 1 #loop count
    dif = 1000000 #check variable to converge

    while kount<maxkount and dif>maxdif:
        dif = 0
        y_ch = deepcopy(y) #for check
        #create converge curve from guess(y)
        for i in range(0,m):
            def_k = (y[i][0]**(alpha+1)) - y[i][0]*y[i][1] - (nPop + delta)*(y[i][0]**2)
            def_c = (y[i][1]**2) / theta * (alpha*y[i][0]**(alpha-1) - (delta+rho))
            y[i+1][0] = y[i][0] + def_k
            y[i+1][1] = y[i][1] + def_c

        
        #compute differ btw y and y_ch
        dif_m = y - y_ch

        for i in range(0,2):
            for j in range(0,m+1):
                dif = math.fabs(dif_m[j][i])
                y[j][i] = 0.5 * (y[j][i] + y_ch[j][i])
        
        #print kount,"th loop is finished"
        #print "error is",dif
        kount = kount + 1

    k_move = y[m][0] - y[0][0]
    c_move = y[m][1] - y[0][1] 

    return [k_move,c_move]

#===================
#MAIN CALCULATION
#===================
#plotting graph
title("ramsey model")
xlabel("k")
ylabel("c")
plot([min_x,max_x],[min_y,max_y],"x")

arr_x = arange(min_x,max_x + unit_plot_x,unit_plot_x)
arr_y = arange(min_y,max_y + unit_plot_y,unit_plot_y)

for i in arr_x:
    for j in arr_y:
        #create arrow
        arr = ram_arrow(i,j)
        arr_len = np.linalg.norm(arr) #length of arr
        arr /= arr_len #normalize
        arr = arr * max_arr * (1-1000**(-arr_len)) 
        dk = arr[0]
        dc = arr[1]

        #decide color
        v_color = vector_angle(arr,[1,0])
        v_color_10 =  int(math.ceil( v_color / (2 * pi) * 16777215 ) -1) #v_color_16
        col_r = v_color_10 % 256
        col_g = ((v_color-10 - col_r) % (256**2)) / 256
        col_b = (v_color_10 - col_r - col_g*256 ) / 256**2

        col_r = hex(int(col_r))
        col_g = hex(int(col_g))
        col_b = hex(int(col_b))


        col_r = col_r[col_r.rfind("x") +1: len(col_r)]
        while len(col_r) < 2:
            col_r = "0" + col_r
        col_g = col_g[col_g.rfind("x") +1: len(col_g)]
        while len(col_g) < 2:
            col_g = "0" + col_g
        col_b = col_b[col_b.rfind("x") +1: len(col_b)]
        while len(col_b) < 2:
            col_b = "0" + col_b

        col = "#" + col_r + col_g + col_b
        
        while len(col) < 3:
            col += "0"

        #draw arrow
        arrow(i,j,dk,dc,width = 0.003,shape="full", ec = col, fc = col)

show()

"""
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
"""
