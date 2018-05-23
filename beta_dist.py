from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import random
import os
import time
import scipy.special, scipy.optimize, scipy.integrate
import pysysp

def spectrum(a,x,b,z):
	x = x/(1.+z)
	f = np.piecewise(x, [x < 1216, x >= 1216], [0, lambda x: a*(x**b)/(1.+z)])
	return f


def main(magnitude, redshift, beta):

	f160w = pysysp.BandPass('f160w.tab')
	f098m = pysysp.BandPass('f098m.tab')
	f125w = pysysp.BandPass('f125w.tab')
	f600lp = pysysp.BandPass('f600lp.tab')
	f606w = pysysp.BandPass('f606w.tab')
	wavelength = np.arange(0, 30000, 1.);
	lambda_f160w=15369.14
	c=3e18 #In angstrom
	flux_f160w=(c*10**(-magnitude/2.5))/(lambda_f160w**2)
    
	norm = flux_f160w/(lambda_f160w**beta*(1.+redshift))
	spec=spectrum(norm,wavelength,beta,redshift)
	data=[[wavelength],[spec]]

	hdu = fits.PrimaryHDU()
	hdu.header['CTYPE1'] = 'Wavelength'
	hdu.header['CTYPE2'] = 'Flux'
	t=Table(data, names=('Wavelength', 'flux'))
	data=np.array(t)

	fits.writeto('spec.fits',data, clobber=True)
	synth_spec = pysysp.StarSpectrum('spec.fits')

	m_f160w = synth_spec.apmag(f160w,mag='AB')+48.6
	m_f098m = synth_spec.apmag(f098m,mag='AB')+48.6
	m_f125w = synth_spec.apmag(f125w,mag='AB')+48.6
	m_f600lp = synth_spec.apmag(f600lp,mag='AB')+48.6
	m_f606w = synth_spec.apmag(f606w,mag='AB')+48.6

	#m_f098m2 = 1.0674516*m_f160w+6.6643267
	#m_f125w2 = 1.0951267*m_f160w-2.3757875
	#m_f606w2 = 1.1528192*m_f160w+20.214228
	#m_f600lp2 = 1.152062*m_f160w+7.04734
	
	#print m_f160w,m_f098m,m_f125w,m_f606w,m_f600lp
	return m_f160w,m_f098m,m_f125w,m_f606w,m_f600lp      

if __name__ == "__main__":
     main()
