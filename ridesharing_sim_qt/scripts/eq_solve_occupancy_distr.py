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
x = 3.0
v = 1.0
lav = 6.5
#tries = 10000


M = [[[0 for k in range(cap+1)] for j in range(cap+1)] for i in range(cap+1)]
for o in range(cap+1):
    for l in range(cap+1):
        for r in range(cap+1):
            for m in range(l+1):
                for n in range(o-(l-m)+1):
                    for q in range(r+1):
                        if (o - (l-m) - n) == 0:
                            for t in range(r-q, cap+1):
                                M[o][l][r] = M[o][l][r] + (x*v/lav)**(m+n-r+2*q+t) * np.exp(-4*x*v/lav) / (math.factorial(m) * math.factorial(n) * math.factorial(q) * math.factorial(-r+q+t))
                        else:
                            M[o][l][r] = M[o][l][r] + (x*v/lav)**(o-l+2*m-r+2*q+cap) * np.exp(-4*x*v/lav) / (math.factorial(m) * math.factorial(n) * math.factorial(q) * math.factorial(o-l+m-n-r+q+cap))
    #print(M[o])

def F(p):
    prod = np.zeros(cap+2)
    tot = 0.0
    for i in range(cap):
        tot = tot + p[i]
    p[cap] = 1-tot
    
    for o in range(cap+1):
        for l in range(cap+1):
            for r in range(cap+1):
                prod[o] = prod[o] + p[l] * M[o][l][r] * p[r]
    
    #for i in range(cap+1):
    #    if (p[i] < 0):
    #        return (np.abs(p) - prod) - 1000 * p[i]
    #    if (p[i] > 1):
    #        return (np.abs(p) - prod) + 1000 * (p[i]-1)
   
    mean = 0.0
    for i in range(cap+1):
        mean = mean + i * np.abs(p[i])
    
    p[cap+1] = (mean-x)**2
    return p - prod

def mini(initial_guess):
    sol = scipy.optimize.broyden1(F, initial_guess, verbose=2, f_rtol=1e-2, iter=100)
    minimum = F(sol.x)
    
#    for i in range(tries):
#        #print(i)
#        possible = False
#        while possible != True:
#            guess = sol.x * ((np.random.rand(cap+1) - 0.5) / 1.0 + 1.0)
#            total = 0.0
#            for j in range(cap):
#                total = total + guess[j]
#            if (total < 1.0):
#                possible = True
#                guess[len(guess) - 1] = 1.0 - total
#                #print(guess)
#                if (F(guess) < minimum):
#                    print("Here!")
#                    sol = scipy.optimize.minimize(F, guess, method='Nelder-Mead', tol=1e-7)
#                    minimum = F(sol.x)
    
    return sol.x

initial_guess = np.zeros(cap+2)
norm = 0.0

for i in range(cap+1):
    initial_guess[i] = np.exp(-1/2 * (i - x)**2 / x)
    norm = norm + initial_guess[i]
   

initial_guess = initial_guess / norm
initial_guess[cap+1] = 0
print(initial_guess)

sol = mini(initial_guess)
solnew = []
tot = 0.0
for i in range(cap+1):
    #solnew.append(sol[i])
    tot = tot + sol[i]
for i in range(cap+1):
    solnew.append(sol[i] / tot)
print(solnew)
#print(sol)
print(F(sol))


mean = 0.0
for i in range(cap+1):
    mean = mean + i * np.abs(solnew[i])
    
print (mean)