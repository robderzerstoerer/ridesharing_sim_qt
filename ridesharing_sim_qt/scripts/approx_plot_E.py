# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 11:42:23 2019

@author: ms119322
"""

import math
import scipy.optimize
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.special
import csv


topologies = [
    #"two_nodes_N_2",
    "torus_N_25",
    #"star_N_4",
    "simplified_city_N_16",
    "ring_N_25",
    "complete_graph_N_5",
    #"delaunay_random_torus_N_25",
    "3cayley_tree_N_46"
    ]
v = 1.0
q_calc = 25
num_iterations = 100


def fit_func_ln(lnB, lnBhalf):
    return lnBhalf - np.log(np.exp(lnBhalf) + np.exp(lnB))

def sort_data(xdata, ydata):
    if len(xdata) <= 1:
        return
    
    for i in range(len(xdata) - 1):
        smallest = 10000000
        smallest_index = -1
        for j in range(i, len(xdata)):
            if (xdata[j] < smallest):
                smallest = xdata[j]
                smallest_index = j
        
        if smallest_index != i:
            xdata[smallest_index] = xdata[i]
            xdata[i] = smallest
            tempy = ydata[smallest_index]
            ydata[smallest_index] = ydata[i]
            ydata[i] = tempy

def F(p, cap, x, lav, p_l1, p_2n):
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
        #if topo == "complete_graph_N_5":
        #    P_o[i][0] = 1.0
        #else:
        for j in range(cap + 1):
            for n in range(max(j-i,0), q_calc+1):
                #P_o[i][j] = P_o[i][j] + p_norm[cap+1+n] * scipy.special.binom(min(i+n,cap), j) * p_1n**(min(i+n,cap)-j) * (1-p_1n)**j
                for m in range(min(n,min(i+n,cap)-i)+1):
                    y = min(i+n,cap)-m-j
                    if y <= i and y >= 0:
                        #print("i = %d" % i)
                        #print("j = %d" % j)
                        #print("n = %d" % n)
                        #print("m = %d" % m)
                        #print("y = %d" % y)
                        P_o[i][j] = P_o[i][j] + scipy.special.binom(min(n,cap-i), m) * (p_l1)**m * (1.0-p_l1)**(min(n,cap-i)-m) * p_norm[cap+1+n] * scipy.special.binom(i, y) * p_2n**(y) * (1-p_2n)**(i-y)
                    
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
    
    
def p_full_approx(cap, x, lav, p_l1, p_2n, output=False):
    
    initial_guess = np.zeros(cap+1+q_calc+1)

    #norm = 0.0
    #for i in range(cap+1):
        #    initial_guess[i] = np.exp(-1/2 * (i - x)**2 / x)
        #    norm = norm + initial_guess[i]
        #for i in range(cap+1):
            #    initial_guess[i] = initial_guess[i] / norm
                              
    initial_guess[cap-1] = 1.0
    initial_guess[cap+1] = 1.0
   
    #print(initial_guess)
    sol = F(initial_guess, cap, x, lav, p_l1, p_2n)
    for it in range(num_iterations):
        sol = F(sol, cap, x, lav, p_l1, p_2n)

    if (output == True):
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
    
    p_full = sum_one/sum_two
    
    if (output == True):
        mean = 0.0
        total = 0.0
        for i in range(cap+1):
            print("p\'_o(%d) = %.6f" % (i, p_prime[i]))
            mean = mean + i * p_prime[i]
            total = total + p_prime[i]

        print("-----------------")
        print ("mean: %.5f" % mean)
        print("-----------------")
        print ("p_full: %.5f" % (p_full))
        print("-----------------")
        print(total)
    
    print("cap = %d, x = %.2f, p_full=%.4f" % (cap, x, p_full))
    
    return p_full

def main():
    
    norm_req_rates = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
    capacities = [1, 2, 3, 4, 5, 6, 7, 8]
    
    bhalves = []
    
    x_lnpoints = []
    y_lnpoints = []
    z_lnpoints = []  
    colors = []
    textlabels = []
    
    t = 0
    
    for topo in topologies:
        
        bhalves.append([])
        
        print(topo)
        
        # first read all the B_1/2
        for x in norm_req_rates:
            x_lndata = []
            y_lndata = []
            x_lndata_fit = []
            y_lndata_fit = []
        
            daten = open("allsim\\eff_" + topo + "_unl_x_" + ("%.1f" % x) + ".dateff.dat", newline = '')
            reader = csv.reader(daten, delimiter='\t')
            num_read = 0
            
            for punkt in reader:
                #print(punkt)
                x_lndata.append(np.log(float(punkt[0])))
                y_lndata.append(np.log(np.max([1.0 - float(punkt[1]), 0.001])))
                num_read = num_read + 1
                
                #if (1.0 - float(punkt[1])) < 0.15:
                    #    x_lndata_fit.append(np.log(float(punkt[0])))
                    #    y_lndata_fit.append(np.log(1.0 - float(punkt[1])))
            
            
            if num_read == 0:
                continue
            
            sort_data(x_lndata, y_lndata)
            # use 4 last values for fitting
            x_lndata_fit = x_lndata[(len(x_lndata) - 3):]
            y_lndata_fit = y_lndata[(len(y_lndata) - 3):]
            daten.close()

#            plt.figure(0, figsize=(12, 8))
#            plt.subplot(111, xscale="log", yscale="log")
#            plt.title("Fitting")
#            plt.xlabel("B")
#            plt.ylabel("Eff")
    
            #plt.plot(np.exp(x_lndata), np.exp(y_lndata), linestyle='-', 
            #         marker = 'o', markersize = 4, label=topo)
    
            popt, pcov = curve_fit(fit_func_ln, x_lndata_fit, y_lndata_fit)
            
            # standard deviation
            stddev = np.sqrt(np.diag(pcov))[0]
        
            bhalves[t].append(np.exp(popt[0]))
        
        
        if topo == "ring_N_25":
            lav = 6.5
            lmax = 12
            p_l1 = 1.0/lmax
            norm_const = 0.0
            for i in range(2, lmax+1):
                norm_const = norm_const + (1.0 - (i-1.0)/lmax)
            p_2n = (1.0-1.0/lmax)/norm_const
        elif topo == "torus_N_25":
            lav = 2.5
            p_l1 = 0.1666666666
            p_2n = 0.5555555555
        elif topo == "complete_graph_N_5":
            lav = 1
            p_l1 = 1
            p_2n = 0.1
        elif topo == "simplified_city_N_16":
            lav = 2.6
            p_l1 = 0.1666666666
            p_2n = 0.520833
        elif topo == "3cayley_tree_N_46":
            lav = 5.46957
            p_l1 = 0.043478
            p_2n = 0.214008
        
    
        print(p_full_approx(5,3.5, lav, p_l1, p_2n,True))
        
        p_measured = []
        p_approx = []
    
        for c in capacities:
        
            ix = 0
            for x in norm_req_rates:
        
                if x >= c:
                    break
                
                x_lndata = []
                y_lndata = []
                
                daten = open("allsim\\" + topo + "_cap_" + ("%d" % c) + "_x_" + ("%.1f" % x) + ".dateff.dat", newline = '')
                reader = csv.reader(daten, delimiter='\t')
                num_read = 0
                p_full2 = 1e10
            
                for punkt in reader:
                    x_lnpoints.append(np.log(float(punkt[0]) / bhalves[t][ix]))
                    z_lnpoints.append(np.log(float(punkt[1])))
                    p_full2 = float(punkt[9])
                    num_read = num_read + 1
                    #print(("%.4f" % punkt[9]) + ", " + ("%.4f" % pfull_approx(x, c)))
            
                # just use the p_full for (almost-)infinite B
                p_app = p_full_approx(c, x, lav, p_l1, p_2n)
                for read in range(num_read):
                    y_lnpoints.append(np.log(p_app))
                    colors.append(t + 1)
                    textlabels.append(topo)
                
                # for the 2D comparison plots
                if (p_full2 <= 1.0):
                    p_measured.append(p_full2)
                    p_approx.append(p_app)
                    
                ix = ix + 1
        
        #sort_data(p_measured, p_approx)
        plt.figure(0, figsize=(12, 8))
        plt.subplot(111)
        plt.title(topo)
        plt.xlabel("p measured")
        plt.ylabel("p approx")
        plt.autoscale(True)
        plt.plot(p_measured, p_approx, 'ro')    
        gerade = np.linspace(0.0, 1.0, 50)
        plt.autoscale(False)
        plt.plot(gerade, gerade, 'b--')
        plt.show()
    
        plt.figure(1, figsize=(12, 8))
        plt.subplot(111, xscale="log", yscale="log")
        plt.title(topo)
        plt.xlabel("p measured")
        plt.ylabel("p approx")
        plt.autoscale(True)
        plt.plot(p_measured, p_approx, 'ro')
        gerade = np.linspace(0.0, 1.0, 50)
        plt.autoscale(False)
        plt.plot(gerade, gerade, 'b--')
        plt.show()

        sum_sq = 0.0
        for i in range(len(p_measured)):
            sum_sq = sum_sq + ((p_measured[i] - p_approx[i])/p_measured[i])**2
    
        print("Relative deviation (mean_sq) = %.3f" % (math.sqrt(sum_sq/len(p_measured))))
    
        t = t + 1
        
    fig = {
        "data": [{"type": "scatter3d",
          "x": np.exp(x_lnpoints),
          "y": np.exp(y_lnpoints),
          "z": np.exp(z_lnpoints),
          "mode": "markers",
          "marker.color": colors,
          "marker.size": 2,
          "text": textlabels,
          "showlegend": True,
          }],
        "layout": {"title": {"text": topo}, 
                   "scene": {"xaxis": {"title": {"text": "B / B_{1/2}"}, "type": "log", "autorange": True, "dtick": 1}, 
                             "yaxis": {"title": {"text": "p_{full, B->inf}"}, "type": "log", "autorange": True, "dtick": 1},
                             "zaxis": {"title": {"text": "E"}, "type": "log", "range": (-1, 0), "dtick": 1}}}
    }
    
    plot(fig)
    
    
if __name__ == "__main__":
    main()