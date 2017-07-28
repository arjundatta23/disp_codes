#!/usr/bin/python

# Standard Python modules

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Modules written by me
import read_surf96_io as rs96

pickledir=sys.argv[1]
jarlist=[n for n in os.listdir(pickledir) if n.endswith('.pckl')]
#print "Found %d pickle files" %(len(jarlist))
numf=1
for j,jarname in enumerate(jarlist):
	jarfile=os.path.join(pickledir,jarname)
	jar=open(jarfile)
	#print "Reading ", jarfile
	cookie = pickle.load(jar)
	jar.close()
	if len(cookie)>numf:
		numf=len(cookie)
	f=[x for x,y,e in cookie]
	c=[y for x,y,e in cookie]
	cerr=[2*e for x,y,e in cookie]		# NB: MULTIPLYING THE ERROR BY 2
	modenum=0 #int(jarname[1])
	#print "Mode number ", modenum
	revind=range(len(cookie))[::-1]
	# NB: don't want to output every sample as that gets too much for surf96. Want to have finer sampling at the 
	# short periods and coarser sampling at long periods
	sample_num=0
	for i in revind:
		if f[i]>0.05:
			coarseness=6
		elif f[i] <= 0.05 and f[i] > 0.02:
			coarseness=4
		elif f[i] <= 0.02 and f[i] > 0.01:
			coarseness=2
		else:
			coarseness=1
		#coarseness=1
		if sample_num%coarseness == 0:
			print "%s %s %s %s  %d	%9.5f	%.2f	%.5f" %('SURF96','R','C','T',modenum,1/(f[i]),c[i],cerr[i])
		sample_num+=1

	plt.xlim(0.005,0.06)
	plt.ylim(3,7)
	#plt.errorbar(f,c,yerr=cerr)
	plt.plot(f,c,'o',color='0.9',markeredgecolor='0.6')
#print "Found total of %d frequency samples" %(numf)

usrc=raw_input("See theoretical dispersion too ? (y/n): ")
if usrc=='y':
	thdpfile=[raw_input('File containing theoretical dispersion: ')]
	rs96obj = rs96.read_disp(thdpfile)
	theor_cdisp = rs96obj.disp[0]
	for k in range(len(theor_cdisp)):
		print "mode number: ", k
		try:
			f_solid = [x for x,y,z in theor_cdisp[k]]
			v = [y for x,y,z in theor_cdisp[k]]
		except ValueError:
			f_solid = [x for x,y in theor_cdisp[k]]
			v = [y for x,y in theor_cdisp[k]]
		curve_name="Mode %d" %k
		plt.plot(f_solid,v,'-',linewidth=2) #,label=curve_name)

plt.xlim(0,0.06)
plt.ylim(3,5)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase Velocity [km/s]')
plt.title('Results from %d profiles' %(len(jarlist)))
plt.grid(True,color='0.5')
plt.show()
