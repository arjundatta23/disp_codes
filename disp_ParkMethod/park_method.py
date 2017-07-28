#!/usr/bin/python

import cmath as cm
import numpy as np

######################################################################################

def dosinglefreq(absdist,true_phase,ktrials,tempprint=False):

	rel_dist=absdist-absdist[0]
	kxmat=np.outer(ktrials,rel_dist)
	apsm=np.exp(+1j*kxmat)
	# apsm stands for applied_phase_shift_matrix
	tpmat=true_phase.reshape(len(true_phase),1)
	#tpmat=np.flipud(true_phase.reshape(len(true_phase),1))
	tpcm=np.exp(1j*tpmat)
	# tpcm stands for true_phase_column_matrix
	stacked_result=np.dot(apsm,tpcm)
	mod_stack=np.abs(stacked_result)
	return mod_stack

######################################################################################
# The below was used for debugging
#	if tempprint:
#		print "Distances are ", rel_dist
#		print "Row of kxmat: ", kxmat[100,:]
#		print "Row of apsm: ", apsm[100,:]
#		print "tpmat: ", tpmat
#		print "tpcm: ", tpcm
#		print "Individual results of matrix multiplication: "
#		summ=0
#		for x in range(apsm.shape[1]):
#			prod= apsm[100,x]*tpcm[x]
#			print prod
#			summ+=prod
#		print "Sum is ", summ
