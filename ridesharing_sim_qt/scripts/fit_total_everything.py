"""

Curve fitting for different topologies

"""

from plotly.offline import init_notebook_mode, plot
import plotly.graph_objs as go
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
#matplotlib.rcParams['text.usetex'] = True
init_notebook_mode()

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
            

def main():

    print (__doc__)
    
    
    topologies = [
    #"two_nodes_N_2",
    "torus_N_25",
    #"star_N_4",
    "simplified_city_N_16",
    "ring_N_25",
    "complete_graph_N_5",
    "delaunay_random_torus_N_25",
    "3cayley_tree_N_46"
    ]
    
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
        
        for x in norm_req_rates:
            x_lndata = []
            y_lndata = []
            x_lndata_fit = []
            y_lndata_fit = []
        
            # first read all the B_1/2
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
    
            fit_x_ln = np.linspace(x_lndata[0], x_lndata[len(x_lndata)-1], 50)
            fit_y_ln = fit_func_ln(fit_x_ln, popt[0])
    
           # plt.plot(np.exp(fit_x_ln), np.exp(fit_y_ln), 'r-',
           #          label= topo + ' B_{1/2}=' + str(np.exp(popt[0])) + ' +/- ' + str(np.exp(popt[0])*stddev))
        
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
                p_full2 = -1.0
            
                for punkt in reader:
                    #print(punkt)
                    if topo == "two_nodes_N_2":
                        x_lnpoints.append(np.log(float(punkt[0])))
                    else:
                        x_lnpoints.append(np.log(float(punkt[0]) / bhalves[t][ix]))
                    z_lnpoints.append(np.log(float(punkt[1])))
                    
                    p_full2 = float(punkt[9])
                    num_read = num_read + 1
                
                    #if (1.0 - float(punkt[1])) < 0.15:
                        #    x_lndata_fit.append(np.log(float(punkt[0])))
                        #    y_lndata_fit.append(np.log(1.0 - float(punkt[1])))
            
                # just use the p_full for (almost-)infinite B
                for read in range(num_read):
                    y_lnpoints.append(np.log(p_full2))
                    colors.append(t + 1)
                    textlabels.append(topo)
                
                ix = ix + 1
             
        print(topo)
        print(bhalves[t])
        print(len(x_lnpoints))
        print(len(y_lnpoints))
        
        
        t = t + 1

    #fig = go.Figure(data=[go.Scatter3d(x=x_lnpoints, y=y_lnpoints, z=z_lnpoints, mode='markers')])
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
