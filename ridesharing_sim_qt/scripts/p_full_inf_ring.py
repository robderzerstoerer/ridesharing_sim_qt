# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:20:02 2019

@author: ms119322
"""

sum1 = 0.0
p1 = 0.359055
lmax = 12
x = 6.0

for xe in range(12):
    prod1 = 1.0
    i = 0
    while i < xe + 1:
        prod1 = prod1 * (1 - p1 * (1 - i / lmax))
        i = i + 6.5 / (x-1)
    sum1 += prod1
        
print(1 - sum1 / lmax)