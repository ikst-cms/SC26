import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
from scipy import optimize


def func(x, a, b, c):
    return a * np.exp( -((x - b) ** 2.0) / (2.0 * c ** 2) )


def plot_histogram(indata, labels):
    """
    This method is used to plot histogram with Normal fit
    :param indata:
    :param labels:
    :return: None.
    """
    color = iter( cm.summer( np.linspace( 0, 1, 3 ) ) )
    NBins = 30
    rcol = next( color )
    x0, x1 = min( indata ), max( indata )
    x0 = x0 - (x1 - x0) * 0.2
    x1 = x1 + (x1 - x0) * 0.2
    yv, bins, Patches = plt.hist( indata, NBins, range=[x0, x1], label=labels, facecolor='grey',
                                  alpha=0.5 )
    xv = [0.5 * (bins[i] + bins[i + 1]) for i in range( len( bins ) - 1 )]
    if sum( _ > 0.0 for _ in yv ) < 2 or abs( xv[0] - xv[1] ) < 1.0e-4:
        print(
          "\nPlot_histogram :: Only one bin with frequency > 0, skipping histogram plot.\n"
        )
        return
    norm_x, norm_y = np.average( xv ), np.average( yv )
    xv = [_ / norm_x for _ in xv]
    yv = [_ / norm_y for _ in yv]
    popt, pcov = optimize.curve_fit( func, xv, yv, p0=[1.0, np.average( xv ), 1.0] )
    delx = (xv[1] - xv[0])
    x_fit = np.linspace( xv[0] - 2 * delx, xv[-1] + 2 * delx, num=70 )
    y_fit = func( x_fit, popt[0], popt[1], popt[2] )
    x_fit = [_ * norm_x for _ in x_fit]
    y_fit = [_ * norm_y for _ in y_fit]
    plt.plot( x_fit, y_fit, lw=2, color='blue' )
    plt.xlabel( 'CE-calculated energies' )
    plt.ylabel( 'Normalized Frequency' )
    plt.title( 'Histogram of Energies' )
    plt.grid( True )
    plt.legend( loc='upper right', numpoints=1 )
    nB = int( 0.6 * len( x_fit ) )
    nBt = int( 0.7 * len( x_fit ) )
    plt.annotate( r'$a e^{-(x-\mu)^2/2\sigma^2}$', xy=(x_fit[nB], y_fit[nB]), xycoords='data',
                  xytext=(x_fit[nBt], y_fit[nB]),
                  textcoords='data', arrowprops=dict( arrowstyle="->" ) )
    plt.show()

    return
