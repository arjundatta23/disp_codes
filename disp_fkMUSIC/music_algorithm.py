#!/usr/bin/python

import numpy as np
import cmath as cm
import matplotlib.pyplot as plt

##########################################################################################################################

""" Variable names inspired from original paper -
    Multiple Emitter Location and Signal Parameter Estimation, IEEE transactions on Antennas and Propagation,
    Vol. AP-34, No. 3, March 1986

    Recorded data vector: xvect
    Spatial correlation matrix: rmat
"""

def do_single_freq(xvect,rec_locs,fhz,nsignals):

	delta_x = rec_locs[1] - rec_locs[0]
	k_nyquist = 3/(2*delta_x)
	# k_nyquist controls the range of the trial kvalues so make it larger (change the numerator) if your result looks truncated.
	aperture = rec_locs[-1] - rec_locs[0]
	#kstep = 1/(10*aperture)
	kstep=1/(600*delta_x)
	kvalues = 2*np.pi*np.arange(-k_nyquist,k_nyquist+kstep,kstep)
	if (len(kvalues)%2)==0:
		""" Because of the way arange works, sometimes you may get an extra element at the end
		    This is to ensure that you always have an odd number of elements in kvalues
		    (equal number on either side of 0, plus 0) """
		kvalues=kvalues[:-1]
	power_mus = np.zeros(len(kvalues))
	xvalues = [p-rec_locs[0] for p in rec_locs]
	num_sensors = xvect.shape[0]
	freq_samples = xvect.shape[1]
	print("Number of receivers is: ", num_sensors)
	print("Number of frequency samples is: ", freq_samples)
	if freq_samples>1:
		xvect = xvect - (xvect.mean(1) * np.ones((1,freq_samples)) )
	rmat = (xvect*(xvect.H)) / freq_samples
	umat,s,vpython = np.linalg.svd(rmat)
	vmat = vpython.H
	# print("SHAPE OF VMAT: ", vmat.shape)
	# print("NO. OF RECEIVERS: ", num_sensors)
	nss_vmat = vmat[:,nsignals:num_sensors]
	""" Remember numpy's svd returns the conjugate transpose of the V matrix rather than the V matrix itself
	    ( svd(A) = USV.H  and python returns V.H)
	    nss_vmat below stands for noise subspace v matrix """
	print("Shape of noise-subspace eigenvector matrix is: ", nss_vmat.shape)
	approx_rinv = nss_vmat*(nss_vmat.H)
	print("Shape of approximate R-inverse matrix is: ", approx_rinv.shape)
	print("Length of kvalues:", kvalues.size)
	for j,k in enumerate(kvalues):
		steer = np.matrix([ cm.rect(1,-k*x) for x in xvalues ]).T
		denom = steer.H * approx_rinv * steer
		quad_form = denom[0,0]
		# print(final.shape, final, final[0,0], type(final))
		power_mus[j] = abs(1 / quad_form)
	# print("Min and max of power_mus: ", min(power_mus), max(power_mus))
	power_mus=power_mus/(max(power_mus))
	return kvalues, power_mus
