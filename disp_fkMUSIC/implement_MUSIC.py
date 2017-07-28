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

# Modules written by me
import seisarray_data as sad
import music_algorithm as ma
import pairwise_interpolate as pi
import read_earth_io as reo
import read_surf96_io as rs96

##########################################################################################################

def plot_spectrum(showimage,imtype,seeth):

	aperture=int(locs[-1]-locs[0])
	fig=plt.figure()
	ax=fig.add_subplot(111)
	if imtype==1:	
		X,Y = np.meshgrid(kfinal,freqs_final)
		ax.set_ylabel('Frequency [Hz]')
		ax.set_xlabel('Wavenumber [Radian/km]')
		ax.set_ylim(max(freqs_final),min(freqs_final))
		ax.set_xlim(min(kfinal),max(kfinal))
		#ax.set_title('Aperture %d km' %(aperture))
		print "kfinal min and max values are: ", kfinal[0], kfinal[-1]
	elif imtype==2:
		legcols=0
		X,Y = np.meshgrid(freqs_final,cfinal)
		ax.set_xlabel('Frequency [Hz]')
		ax.set_ylabel('Phase Velocity [km/s]')
		ax.set_ylim(cmin,cmax)
			#usrc=raw_input("See theoretical dispersion too ? (y/n): ")
			#if usrc=='y':
		if seeth:
			mcol=['b','g','r','c','m','y','k','b','g','r']
			os.listdir(os.getcwd())
			thdpfile=[raw_input('File containing theoretical dispersion: ')]
			mnums=raw_input("Start and end mode numbers to plot: ")
			#minc=int(raw_input("Incident mode to highlight: "))
			ml=int(mnums.split()[0])
			mh=int(mnums.split()[1])
			theor_disp = reo.read_disp(thdpfile,ml,mh)
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
				#	ax.plot(f,v,'--',color=mcol[mode],label=curve_name)
			if i>3 and i<=7:
				legcols=2
				legsp=0.5
			elif i>7:
				legcols=3
				legsp=0.25
			else:
				legcols=1
				legsp=1
			ax.legend(loc=1,ncol=legcols,labelspacing=legsp)
		ax.set_xlim(freqs_final[0],freqs_final[-1])
		#ax.set_xlim(0.02,freqs_final[-1])
		#ax.set_title('%.2f km, param = %d' %(aperture,nmodes))
		ax.set_title('Aperture %d km' %(aperture))
		if legcols>1:
			leg = plt.gca().get_legend()
        	        ltext  = leg.get_texts()
                	#plt.setp(ltext, fontsize='small')
	print "ARJUN ", len(X), len(Y), showimage.shape
	cax=ax.pcolor(X,Y,showimage)
	#fig.colorbar(cax)
	return ax

##########################################################################################################

def get_windows(ts,wl,nsig):

	""" Function to determine the number of frequency windows, which of course depends on the relationship 
    	    between number of relevant frequency samples and length of (no. of samples in) each window.
    	    Results from MUSIC analysis of each window will be assigned to the centre of the window """

# ts stands for total_samples, wl for window_length
	counter=0
	for s in range(0,ts,nsig):
		if s+wl-1<ts:
			counter+=1
		else:
			break
	return counter

##########################################################################################################

def shift_negative(inspec,shift):
	""" Shifts negative ktrial values over to the positive side thereby altering the range
	    of trial wavenumbers """

	inspec=np.array(inspec)
	if shift==True:
		ncols=inspec.shape[1]
		pos_ks=((ncols-1)/2)+1
		outspec=np.zeros(inspec.shape)
		outspec[:,:(pos_ks-1)]=inspec[:,pos_ks:]
		outspec[:,(pos_ks-1):]=inspec[:,:pos_ks]
		ks=ktrials[1]-ktrials[0]
		numk=len(ktrials)
		kfirstnonneg=ktrials[ktrials>=0][0]
		""" Extra precaution because Python behaves funnily sometimes - replaces 0 with a very
		    small positive number """
		powten=math.floor(math.log(ks,10))
		if kfirstnonneg<(10**(powten)):
			kfirstnonneg=0.0
		print "kfirstnonneg is: ", kfirstnonneg
		kmax=kfirstnonneg+(numk-1)*ks
		knew=np.arange(kfirstnonneg,kmax,ks)
	else:
		outspec=inspec
		knew=ktrials
	return outspec, knew

##########################################################################################################

def phasevel_spectrum(fks):

	""" Function to convert the fk-spectrum into a phase velocity spectrum. Discards negative k values """
	cfixed=np.arange(cmin,cmax+cstep,cstep)
	rel_portion=np.where(kfinal>0)
	kpos=kfinal[rel_portion[0]]
	print "Number of positive k values is: ", len(rel_portion[0])
	print "Length of kpos is: ", len(kpos)
	pvs=np.zeros((len(cfixed),len(freqs_final)))
	for j in range(fks.shape[0]):
		omega=2*np.pi*freqs_final[j]
		cthisf=np.array([ omega/k for k in kpos ])
		fks_thisf=fks[j,rel_portion[0]]
		try:
			pvs[:,j]=pi.lagrange_linear(cthisf[::-1],fks_thisf[::-1],cfixed)
			#print "Interpolated values are: ", pvs[:,j]
			if np.any(np.isnan(pvs[:,j])):
				sys.exit('Problem with interpolation: x and y inputs to function must have the same length')
		except IndexError:
			sys.exit('Problem interpolating at frequency sample number %d, %f Hz' %(j,freqs_final[j]))
	return pvs, cfixed

##########################################################################################################

def store_result():

	#pickledir=raw_input("Name of pickle directory: ")
	pickledir='.'
	resmatrix=fc_spectrum	# resmatrix is result matrix
	frange=[freqs_final[0],freqs_final[-1]]
	crange=[cmin,cmax]
	jarname="music_fc.pckl"
	jarfile=os.path.join(workdir,pickledir,jarname)
	jar=open(jarfile,'w')
	pickle.dump(frange,jar)	
	pickle.dump(crange,jar)
	pickle.dump(resmatrix,jar)
	jar.close()
	print "Stored file ", jarfile

##########################################################################################################

class time_domain_data():

	def __init__(self,dir_data):
		
		dlistsad=[]
		dlistsad.append(dir_data)
		sadobj=sad.read_data(dlistsad,fid,'win',ftype)
		#sadobj=sad.read_data(dlistsad,fid,'whole',ftype)
		usrc_st=raw_input("Use all stations (y/n) ?: ")
		if usrc_st=='y':
			try:
				tdodata=sadobj.windata
			except AttributeError:
				print "Using full data..."
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
		# above line to be commented out to run MUSIC with random moise added to data
		si=sadobj.si
		#td_data=np.flipud(tdudata)
		#"""  Rearrange the td_data matrix so as to have data from the farthest station
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
		self.freqpoints = fs_positive[rel_indices]
		self.afdm = np.matrix(fd_data[:,rel_indices])
		print self.afdm.shape 
		#print "Relevant FFT frequency samples are: ", self.freqpoints
		

################################### Main Program #############################################

script_dir=os.getcwd()
datadir=sys.argv[1]
workdir=os.path.normpath(os.path.join(script_dir,datadir))
evdir=os.path.basename(workdir)
outdir=os.path.dirname(workdir)
fid=sys.argv[2]
try:
	ftype=sys.argv[3]
except IndexError:
	ftype=None
nmodes=int(raw_input("Number of signals (modes) in the data: "))
if nmodes%2==0:
	efs=1
# efs stands for extra_freq_samples
else:
	efs=2
frange = raw_input("Enter frequency range of interest (Hz): ")
lf = float(frange.split()[0])
hf = float(frange.split()[1])
crange = raw_input("Phase velocity range: ")
cmin = float(crange.split()[0])
cmax=float(crange.split()[1])
cstep=float(raw_input("Phase velocity step: "))
tdd_obj = time_domain_data(datadir)
full_datamat = tdd_obj.afdm
freq_values = tdd_obj.freqpoints
locs = tdd_obj.epdist
nfreqs = len(freq_values)
fk_spectrum=[]
freqs_final=[]

fsew=nmodes+efs
# fsew stands for freq_samples_each_window
nw=get_windows(nfreqs,fsew,nmodes)
wstart=0
for w in range(nw):
	ws=wstart
	we=ws+fsew
	wc=(ws+we)/2
	wstart+=nmodes
	print "Working on frequency %f, with columns %d to %d" %(freq_values[wc],ws,we-1)
	musmat = full_datamat[:,ws:we]
	this_freq = freq_values[wc]
	ktrials, peaks_thisfreq = ma.do_single_freq(musmat,locs,this_freq,nmodes)
	fk_spectrum.append(peaks_thisfreq)
	freqs_final.append(this_freq)	
fk_spectrum, kfinal = shift_negative(fk_spectrum,True)
print "Number of windows is: ", nw, nfreqs
print "Shape of fk_spectrum is: ", fk_spectrum.shape
freqs_final=[float("%.4f" %x) for x in freqs_final]
print "Frequencies are: ", freqs_final
#print "Used station numbers %d to %d" %(s1,s2)
print "Finished MUSIC calculations, computing phase velocity spectrum..."
fc_spectrum, cfinal = phasevel_spectrum(fk_spectrum)
#plot_spectrum(fk_spectrum,1,False)
fcspec = plot_spectrum(fc_spectrum,2,False)
#cpicks = cfinal[np.argmax(fc_spectrum,axis=0)]
#cpicks=[float("%.5f" %x) for x in cpicks]
#fcspec.plot(freqs_final,cpicks,'wo')
plt.show()
plt.savefig('music_result.png')
usrc='n'
if fid.endswith('SAC') or fid.endswith('sac'):
	# usually means you're looking at synthetic data
	usrc=raw_input('Plot theoretical dispersion too ? (y/n): ')
	#usrc='n'
	if usrc=='y' or usrc=='Y':
		fcspec = plot_spectrum(fc_spectrum,2,True)
		plt.show()
		#plt.savefig('music_fc.png')
usrc2=raw_input("Do you want to save the results ? (y/n) ")
if usrc2=='y':
	store_result()
