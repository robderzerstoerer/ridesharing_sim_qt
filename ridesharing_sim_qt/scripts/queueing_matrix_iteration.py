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
import scipy.special

cap = 5
x = 3.5
v = 1.0
lav = 6.5
lmax = 12
q_calc = 20

p_2n = (1.0-1.0/12)/(5.5)


def F(p):
    P_o = [[0 for i in range(cap+1)] for j in range(cap+1)]
    P_q = [[0 for i in range(q_calc+1)] for j in range(q_calc+1)]
    
    for i in range(len(p)):
        p[i] = max(p[i],0)
    normo = 0.0
    p_norm = p
    for i in range(cap+1):
        normo = normo + p[i]
    for i in range(cap+1):
        p_norm[i] = p[i] / normo
    normq = 0.0
    for i in range(cap+1, cap+1+q_calc+1):
        normq = normq + p[i]
    for i in range(cap+1, cap+1+q_calc+1):
        p_norm[i] = p[i] / normq
    
    for i in range(cap+1):
        #P_o[i][0] = 1.0
        for j in range(cap + 1):
            for n in range(max(j-i,0), q_calc+1):
                #P_o[i][j] = P_o[i][j] + p_norm[cap+1+n] * scipy.special.binom(min(i+n,cap), j) * p_1n**(min(i+n,cap)-j) * (1-p_1n)**j
                for m in range(min(n,min(i+n,cap)-j)+1):
                    y = min(i+n,cap)-m-j
                    if y <= i:
                        P_o[i][j] = P_o[i][j] + scipy.special.binom(min(n,cap-i), m) * (1.0/lmax)**m * (1.0-1.0/lmax)**(min(n,cap-i)-m) * p_norm[cap+1+n] * scipy.special.binom(i, y) * p_2n**(y) * (1-p_2n)**(i-y)
                    
    for i in range(q_calc+1):
        for j in range(q_calc+1):
            if i <= cap:
                
                sum_occs = 0.0
                for l in range(cap - i + 1):
                    sum_occs = sum_occs + p_norm[l]
                P_q[i][j] = P_q[i][j] + sum_occs * 1.0 / math.factorial(j) * (x*v/lav)**(j) * np.exp(-x*v/lav)
                
                for m in range(1, min(i,j) + 1):
                    P_q[i][j] = P_q[i][j] + p_norm[cap-i+m] * 1.0 / math.factorial(j-m) * (x*v/lav)**(j-m) * np.exp(-x*v/lav)
            
            else:
                for l in range(min(j-i+cap,cap)+1):
                    P_q[i][j] = P_q[i][j] + p_norm[l] * 1.0 / math.factorial(j-(i-(cap-l))) * (x*v/lav)**(j-(i-(cap-l))) * np.exp(-x*v/lav)
    
    #print(p_norm)
    #print(P_o[cap])
    
    rhs = np.zeros(cap+1+q_calc+1)
    for i in range(cap+1):
        prod = 0.0
        for k in range(cap+1):
            prod = prod + p[k] * P_o[k][i]
        rhs[i] = prod
    for i in range(q_calc+1):
        prod = 0.0
        for k in range(q_calc+1):
            prod = prod + p[cap+1+k] * P_q[k][i]
        rhs[cap+1+i] = prod
    
    for i in range(len(rhs)):
        rhs[i] = max(rhs[i],0.0)
    
    return rhs
    
    

initial_guess = np.zeros(cap+1+q_calc+1)
norm = 0.0

#for i in range(cap+1):
#    initial_guess[i] = np.exp(-1/2 * (i - x)**2 / x)
#    norm = norm + initial_guess[i]
#for i in range(cap+1):
#    initial_guess[i] = initial_guess[i] / norm
    
initial_guess[cap-2] = 1.0
initial_guess[cap+2] = 1.0
   
print(initial_guess)
sol = F(initial_guess)
for it in range(100):
    sol = F(sol)
    
sum1 = 0.0


print("-----------------")
for i in range(cap+1):   
    print("p_o(%d) = %.6f" % (i, sol[i]))
print("-----------------")
for i in range(q_calc+1):
    print("p_q(%d) = %.6f" % (i, sol[cap+1+i]))
print("-----------------")
#print(sol)
#print(F(sol))

    
p_prime = np.zeros(cap+1)
sum_one = 0.0
sum_two = 0.0
for i in range(cap+1):
    #p_prime[i] = sol[i]
    for j in range(q_calc+1):
        p_prime[i+min(j, cap-i)] = p_prime[i+min(j, cap-i)] + sol[i] * sol[cap+1+j]
        sum_one = sum_one + (j - min(j, cap-i)) * sol[i] * sol[cap+1+j]
for j in range(q_calc+1):
    sum_two = sum_two + j * sol[cap+1+j]
    
mean = 0.0
total = 0.0
for i in range(cap+1):
    print("p\'_o(%d) = %.6f" % (i, p_prime[i]))
    mean = mean + i * p_prime[i]
    total = total + p_prime[i]

print("-----------------")
print ("mean: %.5f" % mean)
print("-----------------")
print ("p_full: %.5f" % (sum_one/sum_two))
print("-----------------")


print(total)

