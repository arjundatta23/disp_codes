#!/usr/bin/python

#########################
#
# Arjun Datta, March 2014
#
#########################


# Available Python modules
import os
import sys
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('../modules_common')

# Modules written by me
import seisarray_data as sad
import park_method as pm
import read_earth_io as reo
#import read_surf96_io as rs96

####################################################################################################################

def plot_spectrum(showimage,scaleshow,showth,pp):

	mcol=['b','g','r','c','m','y','k','b','g','r']
        X,Y = np.meshgrid(tdd.freqpoints, ctrials)
        pmfig=plt.figure()
        ax=pmfig.add_subplot(111)
        cax=ax.pcolor(X, Y, showimage)#, cmap=plt.cm.jet_r)
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Phase Velocity [km/s]')
        if scaleshow:
        	#cb = pmfig.colorbar(cax,orientation='horizontal',fraction=0.05,pad=0.15)
               cb = pmfig.colorbar(cax)
	if showth:
		dfile=raw_input('File containing theoretical dispersion: ')
		#dispfile=os.path.join(script_dir,dfile)
		mnums=raw_input("Start and end mode numbers to plot: ")
		#minc=int(raw_input("Incident mode to highlight: "))
		ml=int(mnums.split()[0])
		mh=int(mnums.split()[1])
		theor_disp = reo.read_disp([dfile],ml,mh)
		# if dispfile happens to be from surf96 instead of earthsr, that is taken care of by reo
		solidcurves = theor_disp.modcdisp[0]
		for i,mode in enumerate(theor_disp.rel_modes):
		        try:
				f = [x for x,y in solidcurves[i]]
			        v = [y for x,y in solidcurves[i]]
			except ValueError:
				f = [x for x,y,z in solidcurves[i]]
			        v = [y for x,y,z in solidcurves[i]]
		        curve_name="Mode %d" %mode
			#if mode==minc:
			ax.plot(f,v,'-',linewidth=2,color=mcol[mode],label=curve_name)
			#else:
		        #ax.plot(f,v,'--',linewidth=2,color=mcol[mode],label=curve_name)
		ax.legend(loc='best')
        if pp==1:
		ucdpfname=raw_input('UCD pickle file: ')
		ucdpfile=os.path.normpath(os.path.join(script_dir,ucdpfname))
		try:
			jar=open(ucdpfile)
			print "Reading ", ucdpfile
			cookie3 = pickle.load(jar)
			fsam=np.array([float("%.4f" %(x)) for x,y in cookie3])
			vel=np.array([y for x,y in cookie3])
			npicks=vel.shape[1]
			for pick in range(npicks):
				cpick=vel[:,pick]
				ax.plot(fsam,cpick,'wo')
		except ValueError:
			print "Unable to read ucd pickle file"
			pass
        df=tdd.freqpoints[1]-tdd.freqpoints[0]
        ax.set_xlim(lf,hf)
        defxticks=ax.get_xticks()
        ax.set_xticks(defxticks[:-1])
        print "defxticks are: ", defxticks
	ax.set_xlim(tdd.freqpoints[0],tdd.freqpoints[-1])
        ax.set_ylim(3.0,8.0)
        ax.set_title("Aperture %d km" %(tdd.epdist[-1]-tdd.epdist[0]))
        plt.show()
	#plt.savefig('pm_result.png')

##########################################################################################################

def store_result():

	#pickledir=raw_input("Name of pickle directory: ")
	pickledir='.'
	resmatrix=specnorm	# resmatrix is result matrix
	freqs_final=tdd.freqpoints
	frange=[freqs_final[0],freqs_final[-1]]
	crange=[cmin,cmax]
	jarname="pm_fc.pckl"
	jarfile=os.path.join(workdir,pickledir,jarname)
	jar=open(jarfile,'w')
	pickle.dump(frange,jar)	
	pickle.dump(crange,jar)
	pickle.dump(resmatrix,jar)
	jar.close()
	print "Stored file ", jarfile

##########################################################################################################

class get_freqdomain_info():

	def __init__(self,dir_data):
		
		dlistsad=[]
		dlistsad.append(dir_data)
		#sadobj=sad.read_data(dlistsad,fid,'win',ftype)
		sadobj=sad.read_data(dlistsad,fid,'whole',ftype)
		usrc_st=raw_input("Use all stations (y/n) ?: ")
		#usrc_st='y'
                if usrc_st=='y':
			try:
				tdodata=sadobj.windata
				print "using windata"
			except AttributeError:
				print "using fulldata"
				tdodata=sadobj.fulldata
                        self.epdist=sadobj.ed
                else:
                        st_range=raw_input("Start and end station numbers: ")
                        s1=int(st_range.split()[0])
                        s2=int(st_range.split()[1])
                        try:
				tdodata=sadobj.windata[s1-1:s2]
			except AttributeError:
				tdodata=sadobj.fulldata[s1-1:s2]
                        self.epdist=sadobj.ed[s1-1:s2]
		anobj=sad.add_noise(tdodata,sadobj.si,1)
		tdudata=anobj.newsig
		tdudata=tdodata
		# above line to be commented in order to run with noise added to data
		si=sadobj.si
		#td_data=np.flipud(tdudata)
		#"""  The above line rearranges the td_data matrix so as to have data from the farthest station
		#     in the first row and that from the closest (to the source) station in the last row
		#"""
		td_data=tdudata
		num_ts=td_data.shape[1]
		td_data = td_data - td_data.mean()
		""" Remove data mean before Fourier transforming"""
		fs=np.fft.fftfreq(num_ts,si)
		fd_data=np.fft.rfft(td_data)
		fs_positive=fs[:fd_data.shape[1]]
		if num_ts%2==0:
			fs_positive[-1]=-1*fs_positive[-1]
		# when number of samples is even, the last value returned by rfft is the Nyquist frequency term
		# whereas the corresponding value returned by fftfreq is the negative Nyquist frequency
		#else:
		#	sys.exit('Odd number of time samples: %d' %(num_ts))
		rel_indices=np.intersect1d(np.where(fs_positive>=lf)[0],np.where(fs_positive<=hf)[0])
		print "length of rel_indices is ", len(rel_indices)
		self.freqpoints = fs_positive[rel_indices]
		self.afdm = np.matrix(fd_data[:,rel_indices])
		print self.afdm.shape 
		#print "Relevant FFT frequency samples are: ", self.freqpoints
		self.cumphase=np.unwrap(np.angle(self.afdm))	
		#self.cumphase=np.angle(self.afdm)	
#		print "Shape of cumphase: ", self.cumphase.shape
#		tempfig=plt.figure()
#		axtemp=tempfig.add_subplot(111)
#		for i in range(self.cumphase.shape[0]):
#			axtemp.plot(self.freqpoints,self.cumphase[i,:])

################################### Main Program #############################################

shot=False
cmin=3.0
cmax=8.0
step=0.01
ctrials=np.arange(cmin,cmax+step,step)
script_dir=os.getcwd()
datadir=sys.argv[1]
workdir=os.path.normpath(os.path.join(script_dir,datadir))
fid=sys.argv[2]
try:
	ftype=sys.argv[3]
except IndexError:
	ftype=None
frange = raw_input("Enter frequency range of interest (Hz): ")
lf = float(frange.split()[0])
hf = float(frange.split()[1])
tdd=get_freqdomain_info(datadir)
specraw=np.zeros((ctrials.size,len(tdd.freqpoints)))
#usrc2=raw_input("Normalize v-f spectrum ? (y/n): ")
usrc2='y'
for j in range(len(tdd.freqpoints)):
	ktry=(2*np.pi*tdd.freqpoints[j])/ctrials
	print "Working on frequency %f" %(tdd.freqpoints[j])
	stack_thisfreq=pm.dosinglefreq(tdd.epdist,tdd.cumphase[:,j],ktry,False)
	specraw[:,j]=stack_thisfreq[:,0]
#	if j==11:
#		ptitle="Freq %f Hz " %(tdd.freqpoints[j])
#		tempfig=plt.figure()
#		axtemp=tempfig.add_subplot(111)
#		axtemp.plot(ctrials,specraw[:,j])
#		axtemp.set_title(ptitle)
#		axtemp.set_xlabel('c [km/s]')
#		axtemp.set_ylabel('Stack amplitude')
		#axtemp.set_ylim(0,5)
print "Plotting..."
usrc3='n'
usrc='y'
if usrc2=='y' or usrc2=='Y':
	specnorm=specraw/np.max(specraw,0)
	plot_spectrum(specnorm,False,False,0)
	#if fid.endswith('SAC') or fid.endswith('sac'):
	if '.SAC' in fid or '.sac' in fid:
		while usrc != 'n':
			usrc=raw_input('Plot theoretical dispersion ? (y/n): ')
			if usrc=='y' or usrc=='Y':
				shot=True
				plot_spectrum(specnorm,False,shot,0)
	while usrc3 != 'n':
		usrc3=raw_input("Plot a pickled file on top (y/n): ")
		if usrc3=='y' or usrc3=='Y':
			plot_spectrum(specnorm,False,shot,1)
else:
	plot_spectrum(specraw,True,shot,0)
	while usrc3 != 'n':
		usrc3=raw_input("Plot a pickled file on top (y/n): ")
		if usrc3=='y' or usrc3=='Y':
			plot_spectrum(specraw,True,shot,1)
usrc_store=raw_input("Do you want to save the results ? (y/n) ")
if usrc_store=='y':
	store_result()
