import numpy as np
import matplotlib.pyplot as plt
import glob
from numpy import trapz
from scipy.optimize import curve_fit
from matplotlib import cm
from astropy.io import ascii
from pylab import *


def main(path_to_cat, xmin, xmax, magbins, ymin, ymax, zbins, cat, cs, ds):
    """
    Generates plots for C(m), S(z,m), and S(z,m)C(m)
    Args:
        path_to_cat (string) = Path to the folder with the science images.
                               Given in the parameters file.
        xmin (float) = Minimum input magnitude assigned by the user.
        xmax (float) = Maximum input magnitude assigned by the user.
        magbins (int) = Number of magnitude bins dessignated by user.
        ymin (float) = Minimum input redshift assigned by the user.
        ymax (float) = Maximum input redshift assigned by the user.
        zbins (int) = Number of redshift bins dessignated by user.
        cat (string) = Name of the field for which the simulation is run.
        cs (int array) = Array with fraction of recovered sources for each
                         magnitude and redshift.
        ds (int array) = Array with fraction of dropouts for each magnitude
                         and redshift.
    """

    x = linspace(xmin, xmax, magbins)
    y = linspace(ymin, ymax, zbins)

    spacing_x = round(5*(xmax - xmin)/(magbins), 1)
    spacing_y = round(5*(ymax - ymin)/(zbins), 1)

    with np.errstate(divide='ignore', invalid='ignore'):
        csm = np.true_divide(ds, cs)  # cover in case of division by 0.

    plt.imshow(cs, extent=[min(y), max(y), max(x), min(x)], cmap='RdPu',
               origin="upper")
    plt.xlabel('$z$', fontsize=16)
    plt.ylabel('Detection Band (AB mag)', fontsize=16)
    plt.yticks(np.arange(xmin, xmax, spacing_x), fontsize=16)
    plt.xticks(np.arange(ymin, ymax, spacing_y), fontsize=16)
    plt.colorbar().set_label(label='$C(m)$', size=16)
    plt.savefig(path_to_cat+'Results/Plots/Completeness_Field'+cat+'.pdf')
    plt.close()

    # Enter the IF below any dropout is different than 0 ()
    if np.sum(ds) != 0.0:
        plt.imshow(ds, extent=[min(y), max(y), max(x), min(x)],
                   cmap='GnBu', origin="upper")
        plt.xlabel('$z$', fontsize=16)
        plt.ylabel('Detection Band (AB mag)', fontsize=16)
        plt.yticks(np.arange(xmin, xmax, spacing_x), fontsize=16)
        plt.xticks(np.arange(ymin, ymax, spacing_y), fontsize=16)
        plt.colorbar().set_label(label='$C(m)S(z,m)$', size=16)
        plt.savefig(path_to_cat+'Results/Plots/Dropouts_Field'+cat+'.pdf')
        plt.close()

        plt.imshow(csm, extent=[min(y), max(y), max(x), min(x)],
                   cmap='YlGn', origin="upper")
        plt.xlabel('$z$', fontsize=16)
        plt.ylabel('Detection Band (AB mag)', fontsize=16)
        plt.yticks(np.arange(xmin, xmax, spacing_x), fontsize=16)
        plt.xticks(np.arange(ymin, ymax, spacing_y), fontsize=16)
        plt.colorbar().set_label(label='$S(z,m)$', size=16)
        plt.savefig(path_to_cat+'Results/Plots/CS_Field'+cat+'.pdf')
        plt.close()
