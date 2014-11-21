#!/usr/bin/env python 
'''
	University of Southampton -- JM -- November 2014

				py_plot_output.py

Synopsis:
	various plotting routines for making standard plots
	from Python outputs 

Usage:
    Either import as a module in a python session e.g.
    import py_plot_output as p 

    or run from the command line e.g. 

    py_plot_output root mode


	
Arguments:
    root 
        root filename to analyse

    mode 
        mode of plotting 
        wind        plot of common wind quantites
        ions        plot of common ions 
        spec        spectrum for different viewing angles
        spec_comps  components contributing to total spectrum e.g. disk, wind
        all         make all the above plots
'''

#import pylab as p 
import pylab as p
import py_read_output as r 
import numpy as np 
import os, sys
import py_plot_util as util

has_astropy = True 
try:
    import astropy
except ImportError:
    has_astropy = False



def make_spec_plot(s, fname, smooth_factor = 10, angles = True, components = False):

	'''
    make a spectrum plot from astropy.table.table.Table object 

    Parameters
    ----------
    s: astropy.table.table.Table
    	table containing spectrum data outputted from Python 

    fname: str 
    	filename to save as e.g. sv

    smooth_factor: int 
    	factor you would like to smooth by, default 10

    angles: Bool 
    	Would you like to plot the viewing angle spectra? 

    components: Bool 
    	would you like to plot the individual components e.g. Disk Wind 
    
    Returns
    ----------
    Success returns 0
    Failure returns 1

    Saves output as "spectrum_%s.png" % (fname)
    '''

	if type(s) != astropy.table.table.Table:
		raise TypeError("make_spec_plot takes astropy.table.table.Table object as first arg")
		return 1

	if s.colnames[8] != "Scattered":
		print "Warning- colnames are not in expected order! %s != Scattered" % (s.colnames[8])

	ncomponents = 9

	if angles:

		# first make viewing angle plot
		p.figure(figsize=(8,12))

		nspecs = len(s.dtype.names) - ncomponents

		nx = 1
		ny = nspecs
		if nspecs > 4:
			nx = 2
			ny = (1 + nspecs) / nx 


		print "Making a %i by %i plot, %i spectra" % (nx, ny, nspecs)

		for i in range(nspecs):

			p.subplot(ny, nx, i+1)

			p.plot(s["Lambda"], util.smooth(s[s.dtype.names[ncomponents + i]], window_len = smooth_factor))

			p.title(s.dtype.names[ncomponents + i])
			p.xlabel("Wavelength")
			p.ylabel("Flux")

		p.savefig("spectrum_%s.png" % (fname))
		p.clf()

	if components:

		p.figure(figsize=(8,12))
		p.subplot(211)
		p.plot(s["Lambda"], util.smooth(s["Created"], window_len = smooth_factor), label="Created")
		p.plot(s["Lambda"], util.smooth(s["Emitted"], window_len = smooth_factor), label="Emitted")

		p.subplot(212)

		for i in range(4,9):
			p.plot(s["Lambda"], util.smooth(s[s.dtype.names[i]], window_len = smooth_factor), label=s.dtype.names[i])

		p.xlabel("Wavelength")
		p.ylabel("Flux")
		p.legend()

		p.savefig("spec_components_%s.png" % (fname))
		p.clf()
		
	return 0


def make_spec_plot_from_class(s, fname, smooth_factor = 10, angles = True, components = False):

	'''
    make a spectrum plot from py_classes.specclass object

    Parameters
    ----------
    s: specclass object
    	table containing spectrum data outputted from Python 

    fname: str 
    	filename to save as e.g. sv

    smooth_factor: int 
    	factor you would like to smooth by, default 10

    angles: Bool 
    	Would you like to plot the viewing angle spectra? 

    components: Bool 
    	would you like to plot the individual components e.g. Disk Wind 
    
    Returns
    ----------
    Success returns 0
    Failure returns 1

    Saves output as "spectrum_%s.png" % (fname)
    '''


	if angles:

		# first make viewing angle plot
		p.figure(figsize=(8,12))
		nspecs = len(s.spec)

		nx = 1
		ny = nspecs
		if nspecs > 4:
			nx = 2
			ny = (1 + nspecs) / nx 


		print "Making a %i by %i plot, %i spectra" % (nx, ny, nspecs)

		for i in range(nspecs):


			p.subplot(ny, nx, i+1)

			p.plot(s.wavelength, util.smooth(s.spec[i], window_len = smooth_factor))

			p.xlabel("Wavelength")
			p.ylabel("Flux")

		p.savefig("spectrum_%s.png" % (fname))
		p.clf()

	if components:

		p.figure(figsize=(8,12))
		p.subplot(211)
		p.plot(s.wavelength, util.smooth(s.created, window_len = smooth_factor), label="Created")
		p.plot(s.wavelength, util.smooth(s.emitted, window_len = smooth_factor), label="Emitted")

		p.subplot(212)
		p.plot(s.wavelength, util.smooth(s.censrc, window_len = smooth_factor), label="CenSrc")
		p.plot(s.wavelength, util.smooth(s.disk, window_len = smooth_factor), label="Disk")
		p.plot(s.wavelength, util.smooth(s.wind, window_len = smooth_factor), label="Wind")
		p.plot(s.wavelength, util.smooth(s.hitsurf, window_len = smooth_factor), label="HitSurf")
		p.plot(s.wavelength, util.smooth(s.scattered, window_len = smooth_factor), label="Scattered")
		p.xlabel("Wavelength")
		p.ylabel("Flux")
		p.legend()

		p.savefig("spec_components_%s.png" % (fname))
		p.clf()

	return 0




def make_wind_plot(d, fname, var=None, shape=(4,2), axes="log", den_or_frac=0, fname_prefix="wind", lims=None):
    '''
    make a wind plot from astropy.table.table.Table object 

    Parameters
    ----------
    d: astropy.table.table.Table
    	table containing wind data outputted from Python 
    	if == None then this routine will get the data for you

    fname: str 
    	filename to save as e.g. sv

    var: array type 
    	array of string colnames to plot 

    angles: Bool 
    	Would you like to plot the viewing angle spectra? 

    components: Bool 
    	would you like to plot the individual components e.g. Disk Wind 

    axes: str 
        lin or log axes

    den_or_frac: int
        0 calculate ion densities
        1 calculate ion fractions

    lims: array-like
        limits of plot, specified as ((xmin,xmax), (ymin, tmax))
        can be array or tuple. Default is Nonetype.

    Returns
    ----------
    Success returns 0
    Failure returns 1

    Saves output as "spectrum_%s.png" % (fname)
    '''
    
    if d == None:
        util.get_pywind_summary(fname, den_or_frac=den_or_frac)
        d = r.read_pywind(fname)

    if var == None:
        var = ["ne", "te", "tr", "IP", "nphot", "v", "w", "ionC4"]

    if axes != "lin" and axes != "log":
        print "Error: didn't understand mode %s, defaulting to log" % axes
        axes = "log"

    nplots = len(var)

    # check shape is ok for variables required
    if shape[0] * shape[1] < nplots:
        print "Error: shape is less than length of var array"
        return 1

    if shape[0] * shape[1] > nplots:
        print "Warning: shape is more than length of var array"

    p.figure(figsize=(8,12))

    # cycle over variables and make plot
    for i in range(nplots):

        p.subplot(shape[0], shape[1], i+1)

        value_string = var[i]

        x,z,v = util.wind_to_masked(d, value_string)

        p.contourf(z,x,np.log10(v), extend="both")
        p.colorbar()
        p.title("Log(%s)" % value_string)

        if lims != None:
            p.xlim(lims[0][0], lims[0][1])
            p.ylim(lims[1][0], lims[1][1])

        if axes == "log":
            # log axes
            p.semilogy()
            p.semilogx()

    p.savefig("%s_%s.png" % (fname_prefix, fname))

    return 0



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys

    ishow='yes'
    if len(sys.argv)>2:
        fname = sys.argv[1]
        mode = sys.argv[2]
    else:
        print __doc__
        sys.exit(1)

    if mode == "spec":
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname)

    elif mode == "wind":
        make_wind_plot(None, fname)

    elif mode == "spec_comps":
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname, angles = False, components = True)

    elif mode == "ions":
        make_wind_plot(None, fname, 
                       var = ["ionh1", "ionh2", "ionC3", "ionC4", "ionC5", "ionSi4", "ionN5", "ionO6"],
                       fname_prefix="ions", den_or_frac = 1)

    elif mode == "all":
        s = r.read_spectrum(fname)
        make_spec_plot(s, fname, components = True)
        make_wind_plot(None, fname)
        make_wind_plot(None, fname, 
                       var = ["ionh1", "ionh2", "ionC3", "ionC4", "ionC5", "ionSi4", "ionN5", "ionO6"],
                       fname_prefix="ions", den_or_frac = 1)

    else:
        print "didn't understand mode %s" % mode
        print __doc__










