

from plotly.offline import init_notebook_mode, plot
import plotly.graph_objs as go
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
#matplotlib.rcParams['text.usetex'] = True
init_notebook_mode()


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
            
def pfull_approx0(x, cap):
    sum1 = 0.0
    
    for i in range(0, 100):
        sum1 = sum1 + min(i, cap) * x**i * np.exp(-x) / math.factorial(i)
        
    sum_total = 0.0
    for i in range(0, 100):
        sum_total = sum_total + i * x**i * np.exp(-x) / math.factorial(i)
        
    return 1.0 - sum1/sum_total

def pfull_approx1(x, cap):
    sum1 = 0.0
    
    for i in range(0, 100):
        p_i = 0.0
        for k in range(0, cap + 1):
            p_i = p_i + x**(i+k) * np.exp(-2*x) / (math.factorial(i) * math.factorial(k))
        for j in range(1, 100):
            p_i = p_i + x**(max(0,i-j)) * x**(cap+j) * np.exp(-2*x) / (math.factorial(max(0,i-j)) * math.factorial(cap+j))
        
        sum1 = sum1 + min(i,cap) * p_i
        
    sum_total = 0.0
    for i in range(0, 100):
        sum_total = sum_total + i * x**i * np.exp(-x) / math.factorial(i)
        
    return 1.0 - sum1/sum_total

def main():

    
    
    topologies = [
    #"two_nodes_N_2",
    #"torus_N_25",
    #"star_N_4",
    #"simplified_city_N_16",
    #"ring_N_25",
    "complete_graph_N_5",
    #"delaunay_random_torus_N_25",
    #"3cayley_tree_N_46"
    ]
    
    norm_req_rates = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
    capacities = [1, 2, 3, 4, 5, 6, 7, 8]
    
    bhalves = []
    
    
    t = 0
    
    for topo in topologies:
    
        print(topo)
        
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
                p_full2 = -1.0
            
                for punkt in reader:
                    p_full2 = float(punkt[9])
                    #print(("%.4f" % punkt[9]) + ", " + ("%.4f" % pfull_approx(x, c)))
                
                p_measured.append(p_full2)
                p_approx.append(pfull_approx1(x, c))
        
        #sort_data(p_measured, p_approx)
        plt.figure(0, figsize=(12, 8))
        plt.subplot(111)
        plt.title(topo)
        plt.xlabel("p measured")
        plt.ylabel("p approx")
        plt.plot(p_measured, p_approx, 'ro')
        
        gerade = np.linspace(0.0, p_measured[len(p_measured)-1], 50)
        plt.plot(gerade, gerade, 'b--')
        
        plt.show()
        
        
        
        t = t + 1

if __name__ == "__main__":
    main()
