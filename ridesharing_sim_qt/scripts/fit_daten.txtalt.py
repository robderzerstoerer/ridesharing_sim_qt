"""

Numerische Integration

In diesem Programm werden Mittelpunkts-, Trapez- und Simpsonregel zur
numerischen Berechnung verschiedener Integrale benutzt. Angezeigt wird der
relative Fehler des Integrals ueber der genutzten Intervallbreite h (doppelt
logarithmisch) fuer die Funktion f(x) = cosh(2x), integriert von -pi/2 bis pi/3.
Weiterhin ist das erwartete Skalierungsverhalten fuer die verschiedenen
Methoden eingezeichnet.

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
#matplotlib.rcParams['text.usetex'] = True

def fit_func_ln(lnB, m):
    return lnB + np.log(m)


def main():

    print (__doc__)
    
    x_lndata = []
    y_lndata = []
    x_lndata_fit = []
    y_lndata_fit = []
    
    
    
    daten = open('occdata.txt', newline = '')
        
    reader = csv.reader(daten, delimiter='\t')
    for punkt in reader:
        print(punkt)
        x_lndata.append(np.log(float(punkt[0])))
        y_lndata.append(np.log(float(punkt[1])))
        
        x_lndata_fit.append(np.log(float(punkt[0])))
        y_lndata_fit.append(np.log(float(punkt[1])))
          
    daten.close()

    plt.figure(0, figsize=(12, 8))
    plt.subplot(111)
    plt.title("Fitting")
    plt.xlabel("B")
    plt.ylabel("Eff")
    
    plt.plot(np.exp(x_lndata), np.exp(y_lndata), linestyle='-', 
             marker = 'o', markersize = 4, label='Ring 100')
    
    popt, pcov = curve_fit(fit_func_ln, x_lndata_fit, y_lndata_fit)
    
    # standard deviation
    stddev = np.sqrt(np.diag(pcov))[0]
    
    fit_x_ln = np.linspace(x_lndata[0], x_lndata[len(x_lndata)-1], 50)
    fit_y_ln = fit_func_ln(fit_x_ln, popt[0])
    
    plt.plot(np.exp(fit_x_ln), np.exp(fit_y_ln), 'r-',
             label='fit: B_{1/2}=' + str(np.exp(popt[0])) + ' +/- ' + str(np.exp(popt[0])*stddev))
    
    test_x_ln = np.linspace(x_lndata_fit[0], x_lndata_fit[len(x_lndata_fit)-1], 50)
    test_y_ln = np.log(1/(np.exp(fit_x_ln) + 1))
    
    plt.plot(np.exp(test_x_ln), np.exp(test_y_ln), 'g-',
             label='test')

    

    plt.autoscale(True)    # Das Skalierungsverhalten der Methoden
                                   # soll die Plotgrenzen nicht beeinflussen.

    plt.legend(loc = 2)  # Legende oben links
    plt.show()

if __name__ == "__main__":
    main()
