#!/usr/bin/python

# Standard Python modules

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Modules written by me

import music_supplement as ms

#############################################################################################################

def show_plot(withpicks):
	spec=plt.figure()
	axspec=spec.add_subplot(111)
	X,Y = np.meshgrid(freqs,pvels)
	cax=axspec.pcolor(X,Y,stack)
	if withpicks:
		axspec.plot(fwinpick,pv,'ko')
		axspec.plot(fwinpick,pvlb,color='k',lw=1.5,) #,'*')
		axspec.plot(fwinpick,pvub,color='k',lw=1.5) #,'*')
	axspec.set_xlim(f1,f2)
	axspec.set_ylim(3,7)
	axspec.set_xlabel('Frequency [Hz]')
	axspec.set_ylabel('Phase Velocity [km/s]')
	axspec.set_title('Stack of %d events' %(len(jarlist)))
	spec.colorbar(cax)
	return axspec

##############################################################################################################

def store_result():

	jarname=raw_input("Name of pickle file (must end with .pckl) : ")
	jarfile=os.path.join(os.getcwd(),jarname)
	jar=open(jarfile,'w')
	pickle.dump(spobj.finalans,jar)
	jar.close()

######################################### Main Program #######################################################

pickledir=sys.argv[1]
jarlist=[n for n in os.listdir(pickledir) if n.endswith('.pckl')]
print "Found %d pickle files" %(len(jarlist))
for j,jarname in enumerate(jarlist):
	jarfile=os.path.join(pickledir,jarname)
	jar=open(jarfile)
	print "Reading ", jarfile
	cookie1=pickle.load(jar)
	cookie2=pickle.load(jar)
	cookie3 = pickle.load(jar)
	jar.close()
	if j==0:
		f1=cookie1[0] #0.0065
		f2=cookie1[1] #0.0615
		cmin=cookie2[0] #3
		cmax=cookie2[1]
		nrows=cookie3.shape[0]
		ncols=cookie3.shape[1]
		freqs=np.linspace(f1,f2,ncols) 
		pvels=np.linspace(cmin,cmax,nrows)
		stack=np.zeros((nrows,ncols))
		setv=[cookie1,cookie2]
		print "Numrows is %d" %(nrows)
	else:
		print cookie1, setv[0]
		print cookie2, setv[1]
		if (cookie1!=setv[0] or cookie2!=setv[1]):
			sys.exit('Stored pickles are inconsistent !')
	try:
		stack+=cookie3
	except ValueError:
		sys.exit("Problem with shape of matrix read in from %s " %(jarname))
		
stack=stack/(j+1)
figax=show_plot(False)
usrc_extra=raw_input("Plot ucd pickle too ? (y/n): ")
if usrc_extra=='y':
	ucdps=raw_input('UCD pickle file/directory: ')
	if os.path.isdir(ucdps):
		jarlist=[n for n in os.listdir(ucdps) if n.endswith('c.pckl')]
		nump=len(jarlist)
		print "Found %d ucd pickle files" %(nump)
		for j,jarname in enumerate(jarlist):
			jarfile=os.path.join(ucdps,jarname)
			jar=open(jarfile)
			cookie = pickle.load(jar)
			vel=np.array([y for x,y in cookie])
			if j==0:
				fsamples=np.array([float("%.4f" %(x)) for x,y in cookie])
				pfae=np.zeros((len(jarlist),len(fsamples)))         		# pfae is for picks_from_all_events
			npicks=vel.shape[1]
			for pick in range(npicks):
				cpick=vel[:,pick]
				figax.plot(fsamples,cpick,'-o')
			jar.close()
	else:
		jar=open(ucdps)
		print "Reading ", ucdps
		cookie = pickle.load(jar)
		fsam=np.array([float("%.4f" %(x)) for x,y in cookie])
		vel=np.array([y for x,y in cookie])
		npicks=vel.shape[1]
		for pick in range(npicks):
			cpick=vel[:,pick]
			figax.plot(fsam,cpick,'-o')
plt.show()
#print "Shape of original spectrum is ", stack.shape
usrc=raw_input("Pick curve (with errors) on spectrum ? (y/n): ")
if usrc=='y':
	spobj=ms.ypick(stack,freqs,pvels,0.95)
	fwinpick=spobj.frel
	pv=spobj.ypicks
	pvub=spobj.pickub
	pvlb=spobj.picklb
	show_plot(True)
	usrc2=raw_input("Do you want to pickle the result ? (y/n): ")
	if usrc2=='y':
		store_result()
