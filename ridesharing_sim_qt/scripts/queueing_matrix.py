# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 17:19:12 2019

@author: ms119322
"""
import math
import scipy.optimize
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

cap = 5
x = 3.5
v = 1.0
lav = 6.5
q_calc = 7


def F(p):
    P_o = [[0 for i in range(cap+1)] for j in range(cap+1)]
    P_q = [[0 for i in range(q_calc+1)] for j in range(q_calc+1)]
    
    
    
    for i in range(cap+1):
        for j in range(cap + 1):
            for n in range(max(j-i,0)+1, q_calc+1):
                P_o[i][j] = P_o[i][j] + p[cap+1+n] * 1.0 / math.factorial(min(i+n, cap)-j) * (x*v/lav)**(min(i+n, cap)-j) * np.exp(-x*v/lav)
    
    for i in range(q_calc+1):
        for j in range(q_calc+1):
            if i <= cap:
                
                sum_occs = 0.0
                for l in range(cap - i + 1):
                    sum_occs = sum_occs + p[l]
                P_q[i][j] = P_q[i][j] + sum_occs * 1.0 / math.factorial(j) * (x*v/lav)**(j) * np.exp(-x*v/lav)
                
                for m in range(min(i,j) + 1):
                    P_q[i][j] = P_q[i][j] + p[cap-i+m] * 1.0 / math.factorial(j-m) * (x*v/lav)**(j-m) * np.exp(-x*v/lav)
            
            else:
                for l in range(min(j-i+cap,cap)):
                    P_q[i][j] = P_q[i][j] + p[l] * 1.0 / math.factorial(j-(i-(cap-l))) * (x*v/lav)**(j-(i-(cap-l))) * np.exp(-x*v/lav)
    
    rhs = np.zeros(2*(cap+1+q_calc+1)+2)
    total_occ = 0.0
    total_queue = 0.0
    for i in range(cap+1):
        total_occ = total_occ + p[i]
        prod = 0.0
        for k in range(cap+1):
            prod = prod + p[k] * P_o[k][i]
        rhs[i] = p[i] - prod
    for i in range(q_calc+1):
        total_queue = total_queue + p[cap+1+i]
        prod = 0.0
        for k in range(q_calc+1):
            prod = prod + p[cap+1+k] * P_q[k][i]
        rhs[i] = p[i] - prod
    
    for i in range(cap+1+q_calc+1, 2*(cap+1+q_calc+1)):
        if p[i-(cap+1+q_calc+1)] < 0.0:
            rhs[i] = p[i-(cap+1+q_calc+1)] * 100
        elif p[i-(cap+1+q_calc+1)] > 1.0:
            rhs[i] = (p[i-(cap+1+q_calc+1)] - 1.0) * 100

    
    rhs[2*(cap+1+q_calc+1)] = 1.0 - total_occ
    rhs[2*(cap+1+q_calc+1)+1] = 1.0 - total_queue    
    
    retvalue = 0.0
    for i in range(len(rhs)):
        retvalue = retvalue + abs(rhs[i])**2
    return retvalue
#    prod = np.zeros(cap+1)
#    sum1 = 0.0
#    sum2 = 0.0
#    for i in range(cap-1):
#        sum1 = sum1 + p[i]
#        sum2 = sum2 + i * p[i]
#    p[len(p)-1] = 1 + x - cap + (cap-1) * sum1 - sum2
#    p[len(p)-2] = 1 - p[len(p)-1] - sum1
    
    

initial_guess = np.zeros(2*(cap+1+q_calc+1)+2)
norm = 0.0

for i in range(cap+1):
    initial_guess[i] = np.exp(-1/2 * (i - x)**2 / x)
    norm = norm + initial_guess[i]
for i in range(cap+1):
    initial_guess[i] = initial_guess[i] / norm
    
initial_guess[cap+1] = 1
   
print(initial_guess)

s = scipy.optimize.minimize(F, initial_guess, method='Nelder-Mead', tol=1e-7)
sol = s.x
#for i in range(len(sol)):
    #print(sol[i])
    #print(F(sol)[i])
print('end')
sum1 = 0.0

for i in range(cap+1):
    sum1 = sum1 + sol[i]
for i in range(cap+1):
    sol[i] = sol[i] / sum1

sum2 = 0.0
for i in range(cap+1, cap+1+q_calc+1):
    sum2 = sum2 + sol[i]
for i in range(cap+1, cap+1+q_calc+1):
    sol[i] = sol[i] / sum2

for i in range(len(sol)):   
    print(sol[i])
#print(sol)
print(F(sol))


mean = 0.0
for i in range(cap+1):
    mean = mean + i * np.abs(sol[i])
    
print (mean)