"""
Contains the main equations and calculations needed for the program.
"""

import numpy as np
import random
import os
import pysysp
import scipy.special
import scipy.optimize
import scipy.integrate
from astropy.io import fits
from astropy.table import Table


def gammainc(s, x):
    """
    Define and return the value of the gamma incomplete function.
    """
    def integrand(t, s): return t**(s-1) * np.exp(-1*t)
    gi = scipy.integrate.quad(integrand, 0, x, args=s)[0]
    return gi


def gamma(s):
    """
    Define and return the value of the gamma function.
    """
    def integrand(t, s): return t**(s-1) * np.exp(-1*t)
    gi = scipy.integrate.quad(integrand, 0, np.inf, args=s)[0]
    return gi


def get_bn(n0):
    """
    Calculates the parameter bn from the Sersic profile.

    Args:
        n0 (int) = Sersic index.

    Returns:
        bn (float) = Value of bn.
    """
    def errfunc(bn, n):
        return abs(scipy.special.gamma(2*n0) -
                   2*scipy.special.gammainc(2*n0, bn) *
                   scipy.special.gamma(2*n0))
    bn = scipy.optimize.fmin(errfunc, 1., args=(n0,), disp=False)[0]
    return bn


def get_Ie_n1(flux0, minor_axis, major_axis):
    """
    Calculates de parameter Ie for the Sersic profile with n = 1,
    which corresponds to the intensity at the radius that encloses
    half of the total light of the galaxy,  the effective radius.

    Args:
        flux0 (float) = Total flux for the galaxy.
        minor_axis (float) = Minor axis of the ellipse in pixels.
        major_axis (float) = Major axis of the ellipse in pixels.

    Returns:
        ie (float) = Value of Ie
    """
    ie = 1 / (2 * minor_axis * major_axis * np.pi)
    return ie


def get_Ie(bn0, n0, flux0, re):
    """
    Calculates de parameter Ie for the Sersic profile with n = 4,
    which corresponds to the intensity at the radius that encloses
    half of the total light of the galaxy, the effective radius.

    Args:
        bn0 (float) = bn parameter.
        flux0 (float) = Total flux for the galaxy.
        re (float) = effective radius in pixels.
    Returns:
        ie (float) = Value of Ie.
    """
    ie = ((bn0)**(2 * n0)) / (re**2 * 2 * np.pi * n0 * gamma(2*n0))
    return ie


def calculate_distance(x1, y1, x2, y2):
    """
    Calculates the distance between two pairs of coordinates.
    Args:

    Returns:
       dis (float) = distance in the same units as the coordinates.
    """
    dis = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    return dis


def makeSersic(n0, bn0, re, ell, inc_angle, size_galaxy):
    """
    Calculates the flux for each pixel following a Sersic profile.

    Args:
        n0 (int) = Sersic index.
        bn0 (float) = bn parameter.
        re (float) = Effective radius in pixels.
        ell (float) = Eccentricity. Varies between 0 and 1.
        inc_angle (float) = Inclination angle in radians. Varies
                            between 0 and Pi/2.
        size_galaxy (int) = Diameter of the galaxy stamp in pixels.
    Returns:
        fl (float) = Flux for each pixel.
    """
    stamp = np.zeros((size_galaxy,size_galaxy))
    s2 = size_galaxy / 2
    major_axis = re
    minor_axis = re * (1.-ell)
    I_e = ((bn0)**(2*n0)) / (2*np.pi*n0*major_axis*minor_axis*gamma(2*n0))
    def f(x, y):
        x_aux = (x-s2)*np.cos(inc_angle) + (y-s2)*np.sin(inc_angle)
        y_aux = -(x-s2)*np.sin(inc_angle) + (y-s2)*np.cos(inc_angle)
        radius = np.sqrt((x_aux/major_axis)**2 + (y_aux/minor_axis)**2)
        return I_e * np.exp(-bn0*((radius)**(1./n0)))

    for i in xrange(size_galaxy):
        def g(x):
            return i - 1./2.
        def h(x):
            return i + 1./2.
        for j in xrange(size_galaxy):
            fl = scipy.integrate.dblquad(f, j-1./2., j+1./2., g, h, 
                                 epsabs=1.49e-08, epsrel=1.49e-08)[0]
            stamp[i,j] = fl
    return stamp


def galaxies_positions(image_data, nsources, size, re):
    """
    Provides the position in the image for the simulated galaxies.

    Args:
        image_data (float array) = Science image. It corresponds to
                                   an array with the counts value for
                                   each pixel.
        nsources (int) = Number of simulated galaxies per iteration.
        size (int) = Diameter of the galaxy stamp in pixels.
        re (float) = Effective radius (in pixels).
    Returns:
        xpos (int) = Position of the simulated galaxy in the x axis.
        ypos (int) = Position of the simulated galaxy in the y axis.
    """
    s2 = size / 2
    xpos, ypos = np.zeros(nsources), np.zeros(nsources)
    for i in range(nsources):
        xr, yr = 1141, 1156
        while ((xr == 0 and yr == 0) or image_data[xr, yr] == 0):
            xr = random.randrange(s2 + 1, image_data.shape[0] - s2 - 1, 1)
            yr = random.randrange(s2 + 1, image_data.shape[1] - s2 - 1, 1)
            if (image_data[xr, yr] != 0):
                xpos[i] = int(xr)
                ypos[i] = int(yr)
    for j in range(nsources):
        d2 = calculate_distance(xpos[j], ypos[j], xpos, ypos)
        w1 = np.where(np.logical_and(np.array(d2) <= re,
                      np.array(d2) > 0.0))[0]
        while ((xpos[j] == 0 and ypos[j] == 0) or
               (image_data[int(xpos[j]), int(ypos[j])] == 0) or
               (len(w1) >= 1)):
            xr = random.randrange(s2 + 1, image_data.shape[0] - s2 - 1, 1)
            yr = random.randrange(s2 + 1, image_data.shape[1] - s2 - 1, 1)
            xpos[j] = int(xr)
            ypos[j] = int(yr)
            d2 = calculate_distance(xpos[j], ypos[j], xpos, ypos)
            w1 = np.where(np.logical_and(np.array(d2) <= 5 * re,
                          np.array(d2) > 0.0))[0]
    xpos[0] = int(1138)
    ypos[0] = int(1152)
    return xpos, ypos


def spectrum(a, x, b, redshift):
    """
    Creates spectrum following a Lyman break galaxy model.
    Args:
        a (float) = Normalisation factor for the spectrum.
        x (float array) = Array with points for the wavelength axis.
        b (float) = Value of the UV spectral slope.
        redshift (float) = Input redshift of the artificial source.
    Returns:
        f (float array )= flux for each wavelength for the spextrum
    """
    x = x / (1.+redshift)
    f = np.piecewise(x, [x < 1216, x >= 1216], [0, lambda x: a*(x**b)/
        (1.+redshift)])
    return f


def write_spectrum(lambda_detection, mag, beta, redshift):
    """
    Saves the spectrum (spec.fits) of the simulated Lyman break galaxy
    so the magnitudes expected in each filter can be calculated.
    Args:
        lambda_detection (float) = Central wavelength of the detection band 
                                   in Angstroms.  
        mag (float) = Input mangitude of the artificial source in the 
                      detection band.
        beta (float) = Value of the UV spectral slope.
        redshift (float) = Input redshift of the artificial source.
    """
    wavelength = np.arange(0, 30000, 1.)
    c = 2.99792458e18  # In angstrom
    flux_detection = (c * 10**(-mag/2.5)) / (lambda_detection**2)
    norm = flux_detection / (lambda_detection**beta * (1.+redshift))
    spec = spectrum(norm, wavelength, beta, redshift)
    data = [[wavelength],[spec]]
    hdu = fits.PrimaryHDU()
    hdu.header['CTYPE1'] = 'Wavelength'
    hdu.header['CTYPE2'] = 'Flux'
    t = Table(data, names=('Wavelength', 'flux'))
    data = np.array(t)
    fits.writeto('Files/spec.fits', data, clobber=True)


def mag_band(name_band, zp):
    """
    Args
        name_band (string) = name of the band for which the magnitude of the
                             simulated galaxy is calculated.
    Returns
        mag (float): Expected magnitude in "name_band".
    """
    synth_spec = pysysp.StarSpectrum('Files/spec.fits')
    flux = pysysp.BandPass('Files/' + name_band + '.tab')
    mag = synth_spec.apmag(flux, mag='AB') + 48.6
    return mag
