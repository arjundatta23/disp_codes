#!/usr/bin/python

# this is a script to read (any number of) SECOND-ORDER ucd pickles.
# ucd pickles are of types - first-order and second-order

# Standard Python modules

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

import read_surf96_io as rs96

######################################### Main Program #######################################################
fig_allev=plt.figure()
ax_allev=fig_allev.add_subplot(111)
pickledir=sys.argv[1]
jarlist=[n for n in os.listdir(pickledir) if n.endswith('.pckl')]
num_pickles=len(jarlist)
print "Found %d pickle files" %(num_pickles)
fsam_allp=range(num_pickles)
mtot=6
# total number of modes in the pickle files needs to be known before hand
usrc=raw_input("See theoretical dispersion ? (y/n): ")
if usrc=='y':
	thdpfile=[raw_input('File containing theoretical dispersion: ')]

val_mode_allp=[[] for i in range(mtot)]

#print val_mode_allp

for mode in range(mtot):
	for j,jarname in enumerate(jarlist):
		jarfile=os.path.join(pickledir,jarname)
		jar=open(jarfile)
		cookie1 = pickle.load(jar)
		cookie2 = pickle.load(jar)
		#print cookie1.shape, cookie2.shape

		fstp = cookie1
		# fstp = freq_samples_this_pickle
		vtp = cookie2[:,mode]
		# vtp = values_this_pickle
		#print vtp
		
		val_mode_allp[mode].append(vtp)
			
		if mode==0:
			fsam_allp[j]=fstp
			# for a given pickle, all modes have the same frequency samples (many of which are not picked)
		
		jar.close()
	
#print val_mode_allp

################################################################################################################

mks=['o','s','^','d','p','v']
wm='0.9'
wm2='0.8'

shown_already=[False for mode in range(mtot)]
for mn in range(mtot):
	curve_name="Mode %d" %(mn)
	allf_thismode=[]
	fv_rel=range(num_pickles)
	for pkl,fsp in enumerate(fsam_allp):
		# fsp is frequency_samples_pickle
		#print val_mode_allp[mn][pkl] #, len(val_mode_allp[mn]), len(fsam_allp)
		if shown_already[mn] or np.all(val_mode_allp[mn][pkl]==-1):
			ax_allev.plot(fsp,val_mode_allp[mn][pkl],mks[mn],zorder=-100,color=wm,markeredgecolor=wm2)
		else:
			ax_allev.plot(fsp,val_mode_allp[mn][pkl],mks[mn],zorder=-100,color=wm,markeredgecolor=wm2,label=curve_name)
		if not shown_already[mn] and not np.all(val_mode_allp[mn][pkl]==-1):
			shown_already[mn]=True

		fv=zip(fsp,val_mode_allp[mn][pkl])
		fpicked=[x for x,y in fv if y>0]
		vpicked=[y for x,y in fv if y>0]
		allf_thismode = np.union1d(allf_thismode,fpicked)
		fv_rel[pkl]=zip(fpicked,vpicked)
	#print "Picked frequency samples for mode ", mn
	#print allf_thismode

	# have to loop over pickles a 2nd time, after the first time is completed
	# this is because, only after the 1st loop has completed do we know all the frequency samples picked for a mode.

	if usrc=='y':
		rs96obj = rs96.read_disp(thdpfile)
        	theor_cdisp = rs96obj.disp[0]
		solidcurve = theor_cdisp
		print "Theoretical curve mode number: ", mn
		try:
			f = [x for x,y,z in solidcurve[mn]]
			v = [y for x,y,z in solidcurve[mn]]
		except ValueError:
			f = [x for x,y in solidcurve[mn]]
			v = [y for x,y in solidcurve[mn]]
		curve_name="Mode %d" %mn
		ax_allev.plot(f,v,'.-') #,label=curve_name)
	elif len(allf_thismode)>0:
		mode_pick_matrix=np.zeros((num_pickles,len(allf_thismode)))
		#mode_pick_matrix[:]=np.nan
		#print mode_pick_matrix.shape
		for pkl in range(num_pickles):
			pfs_tpm=np.array([x for x,y in fv_rel[pkl]])
			# pfs_tpm is picked_freq_sample_this_pickle_&_mode
			loc_thisp=np.searchsorted(allf_thismode,pfs_tpm)
			#print pfs_tpm, loc_thisp
			mode_pick_matrix[pkl,loc_thisp]=[y for x,y in fv_rel[pkl]]

		# ideally, I would fill all the "unpicked" elements of "mode_pick_matrix" with nans. And then use np.nanmean
		# and np.nanstd on each column of the matrix
		# but it turns out np.nanmean and np.nanstd are not available on the old version of numpy.
		# so I have to fill the "unpicked" elements of "mode_pick_matrix" with zeros and then loop through each
		# column of the matrix, extract the nonzero elements and calculate their mean. No way to ignore zeros
		# when calculating mean etc.
		# Unfortunately, this means AN EXTRA LOOP!!

		mode_mean=[]
		mode_std=[]
		for col in range(mode_pick_matrix.shape[1]):
			samples_thisfreq=mode_pick_matrix[:,col]
			#nstf=np.count_nonzero(samples_thisfreq)
			# nstf is num_samples_this_freq
			stf=samples_thisfreq[np.nonzero(samples_thisfreq)]
			# stf is samples_this_freq
			#if mn==3:
			#	print allf_thismode[col], np.mean(stf), np.std(stf)
			mode_mean.append(np.mean(stf))
			mode_std.append(np.std(stf))
		
		for i,freqout in enumerate(allf_thismode[::-1]):
			print "%s %s %s %s %3d %11.4f %10.4f %12.6f" %('SURF96','R','C','T',mn,1/freqout,mode_mean[::-1][i],mode_std[::-1][i])
		ax_allev.errorbar(allf_thismode,mode_mean,yerr=mode_std,color='k',ecolor='0.4')

ax_allev.legend(loc=1)
ax_allev.set_xlabel('Frequency [Hz]')
ax_allev.set_ylabel('Phase Velocity [km/s]')
ax_allev.set_xlim(0.004,0.064)
ax_allev.set_ylim(4,7)
plt.show()
