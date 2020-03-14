#!/usr/bin/python

##################################################################################################
# Module meant specifically for analysis of seismic array data
#
# Arjun Datta, February 2014
##################################################################################################

import os
import sys
import math
import obspy.core
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as clt

# Non-standard modules

import window_functions as wf

###################################################################################################

class read_data():

	def __init__(self,dirlist,fe,howmuch,ftype=None):
	
		if not (howmuch=='whole' or howmuch=='win'):
			sys.exit("It is unclear what portion of the time series is to be read. Please try again")
		self.filext=fe
		self.mnames=dirlist
		self.arraydata=range(len(dirlist))
		for dn,dl in enumerate(dirlist):
			self.arraydata[dn]=self.do_single_dir(dl,howmuch)
		print "Read %d directory" %(dn+1)

	def do_single_dir(self,filesdir,hm):
		
		""" Function to create a matrix containing all seismograms from the station array
	    	with each row of the matrix corresponding to a single station.

	    	Arguments: a directory containing sac files, a string indicating filename extension, 
	    	and optionally, a string to identify a particular type of file (eg. a certain component)
	    	to selectively extract """

		
		names_list=os.listdir(filesdir)
		evfile=[n for n in names_list if n.endswith('.txt')]
		try:
			evfile=evfile[0]
			self.evname='Event %s: Vertical Component' %(evfile[:-4])
		except IndexError:
			self.evname='Synthetic Data: Vertical Component'
		print "filext is ", self.filext
		#print names_list
		names_list=[ n for n in names_list if n.endswith(self.filext) ]
		num_stations=len(names_list)
		if num_stations==0:
			sys.exit("Terminating program - No relevant seismograms in directory")
		file_list = []
		epdist_unsorted=[]
		self.stnames=[]
		""" Sort seismogram list by epicentral distance of station """
		for i,x in enumerate(names_list):
			file_list.append(os.path.join(filesdir,x))
			try:
				tr = obspy.core.read(file_list[-1])[0]
				epdist_unsorted.append(float(tr.stats.sac.dist))
				self.stnames.append(tr.stats.station)
			except TypeError:
				sys.exit("Python cannot read %s" %(file_list[-1]))
		combined = zip(epdist_unsorted,file_list,self.stnames)
		ep_dist = [x for x,y,z in sorted(combined)]
		records = [y for x,y,z in sorted(combined)]
		self.stnames = [z for x,y,z in sorted(combined)]
		""" Time series parameters obtained from the first station in the array
	    	Ideally recording parameters are the same for all stations """
		example_station = records[0]
		tr = obspy.core.read(example_station)[0]
		samples = tr.stats.npts
		self.si = tr.stats.delta
		#self.tstart=0
		""" Start time determined from nearest station and end time from farthest station """
		self.tstart = tr.stats.sac.b
		self.tend = obspy.core.read(records[-1])[0].stats.sac.e
		print "start and end times are: ", self.tstart, self.tend
		tlen = self.tend + self.si - self.tstart
		#self.fulldata = np.matrix(np.zeros((len(records),samples)))
		self.fulldata = np.matrix(np.zeros((len(records),int(math.ceil(tlen)))))
		print "Shape of fulldata is ", self.fulldata.shape
		if __name__ == '__main__':
			print "\nMin ep dist of array is for station %s: %.3f" %(self.stnames[0],ep_dist[0])
			print "Max ep dist of array is for station %s: %.3f" %(self.stnames[-1],ep_dist[-1])
			print "Station spacing for array: %.3f " %(ep_dist[1]-ep_dist[0])
		for i,x in enumerate(records):
			tr=obspy.core.read(x)[0]
			try:
				self.fulldata[i,:] = tr.data
			except ValueError:
				""" This means that
					1. Either all records in the array do not have the same number of samples as the first one
					   OR
					2. Start and end times are different for each record
				"""
				if abs(self.fulldata.shape[1]-len(tr.data))<0.01*samples:
					""" If the discrepancy is less than 1% of the expected value it is accepted with
						truncation or zero padding as the case may be """
					if self.fulldata.shape[1]<len(tr.data):
						#print "tlen and samples are: ", tlen, samples
						try:
							self.fulldata[i,:] = tr.data[:samples]
						except ValueError:
							self.fulldata[i,:] = tr.data[:samples-1]
					elif self.fulldata.shape[1]>len(tr.data):
						self.fulldata[i,:len(tr.data)] = tr.data
						self.fulldata[i,len(tr.data):] = 0.0
				elif self.tend > (self.tstart + samples*self.si):
					# this means different records have different start and end times
					trst=tr.stats.sac.b
					tret=tr.stats.sac.e
					#print i, tr.stats.station, trst, tret, tret-trst
					self.fulldata[i,(trst-self.tstart):(tret+self.si-self.tstart)] = tr.data
				else:
					# if it is a case of difference in no. of samples and the discrepancy is greater than 1%, the program quits """
					print "NPTS in %s is %d but number of samples in %s is %d" %(os.path.basename(records[0]),samples,os.path.basename(x),len(tr.data))
					sys.exit("Quitting - cannot read data")
		self.ed=np.array(ep_dist)
		#return self.fulldata, np.array(ep_dist), si
		if hm=='win':
			#self.win_wholesec(2,8)
			self.win_eachrec(4,8)
			try:
				return self.windata
			except AttributeError:
				return self.fulldata
		else:
			return self.fulldata

	def win_wholesec(self,umin,umax):
	
		""" Function that trims the matrix generated by read_sac_data down to just
	    	the surface wave arrivals (using a group velocity window) across the array """
	
		t1=int(math.floor(self.ed[0]/umax))
		t2=int(math.ceil(self.ed[-1]/umin))
		#if (t2-t1)>=1500 and (t2-t1)<=2000:
		#	t2=t1+2000
		#else:
		#	sys.exit("%d samples in surface wave window" %(t2-t1))
		
		ntfull=self.fulldata.shape[1]	# total number of samples
		print "ntfull is", ntfull, self.tstart
		tfull=np.arange(self.tstart,ntfull*self.si,self.si)
		twin=np.intersect1d(np.where(tfull>=t1)[0],np.where(tfull<=t2)[0])
		samples_win=t2-t1 # number of samples in the window
		if samples_win==len(twin):
			samples_win-=1
		print "samples_win is: ", samples_win
		print "length of twin is: ", len(twin)
		boxcardata=np.array(self.fulldata[:,twin[0]:twin[-1]])
		window=wf.tukeywin(samples_win,0.5)
		window=window.reshape(1,len(window))
		print "Shape of 'window' is ", window.shape
		wtmat=np.dot(np.ones((self.fulldata.shape[0],1)),window)
		try:
			#print "shapes are ", boxcardata.shape, wtmat.shape, self.fulldata.shape, window.shape
			self.windata=np.matrix(boxcardata*wtmat)
		except ValueError:
			print "Time window is %d to %d so number of time samples is %d" %(t1,t2,(t2-t1)/self.si)
			print "But length of twin is: ", len(twin)
			print "First and last elements of twin are: ", twin[0],twin[-1]
			if __name__=='__main__':
				sys.exit('The data are not long enough for the specified minimum group velocity (%.1f km/s)' %(umin))
		self.tstart=tfull[twin[0]]

	def win_eachrec(self,umin,umax):
	
		""" Also a windowing function like win_wholesec but with the difference that the group velocity window 
		    is now applied to each record individually rather than to the entire record section as a whole
		
		    Ultimately we still get the same length of time series from each record by zero padding at the front
		    or at the end (or both) as appropriate.
		"""
	    	
		master_t1=int(math.floor(self.ed[0]/umax))
		master_t2=int(math.ceil(self.ed[-1]/umin))
				
		ntfull=self.fulldata.shape[1]	# total number of samples
		tfull=np.arange(self.tstart,ntfull*self.si,self.si)

		master_twin=np.intersect1d(np.where(tfull>=master_t1)[0],np.where(tfull<=master_t2)[0])
		boxcardata=np.array(self.fulldata[:,master_twin[0]:master_twin[-1]])
		
		wtmat=np.zeros((self.fulldata.shape[0],len(master_twin)))
		print "Shape of wtmat is ", wtmat.shape
		# NB: desired length of time series is known before-hand, so we start with a matrix of 0s
		# and fill the matrix with the window function for each record (each row of matrix)
		# which means zero padding (on either end) is automatic, it does not have to be done explicitly.

		# check for minor dicrepancies due to rounding of floats
		if boxcardata.shape[1]<wtmat.shape[1]:
			boxcardata=np.array(self.fulldata[:,master_twin[0]:master_twin[-1]+1])
		print "Shape of boxcardata is ", boxcardata.shape	
		
		# Loop through each record in the section
		for j in range(self.fulldata.shape[0]):
			rec_t1=int(math.floor(self.ed[j]/umax))
			rec_t2=int(math.ceil(self.ed[j]/umin))
			samples_win = rec_t2 - rec_t1 # number of samples in the window for this record
			#if samples_win==len(twin):
			#	samples_win-=1
			#print "Arrival window for record no. %d is %d seconds" %(j,samples_win)
				
			window=wf.tukeywin(samples_win,0.5)
			winstart = rec_t1 - master_t1
			winend = winstart + len(window)
			#print "..and win start and end are: ", winstart, winend
			try:
				wtmat[j,winstart:winend]=window
			except ValueError:
				print "Time window is %d to %d so number of time samples is %d" %(master_t1,master_t2,(master_t2-master_t1)/self.si)
				print "But length of twin is: ", len(master_twin)
				print "First and last elements of twin are: ", master_twin[0],master_twin[-1]
				if __name__=='__main__':
					sys.exit('The data are not long enough for the specified minimum group velocity (%.1f km/s)' %(umin))

		self.windata=np.matrix(boxcardata*wtmat)
		self.tstart=tfull[master_twin[0]]
	
##################################################################################################################

class add_noise():

	def __init__(self,puresig,si,nt):

		print "Shape of original data is ", puresig.shape
		pcentnoise=20 # strength of noise relative to pure signal at the frequency of signal peak
		self.pcent=pcentnoise/100.0
		nsam=puresig.shape[1]
		fs=np.fft.fftfreq(nsam,si)
		self.sspec=np.fft.rfft(puresig)
		self.fspos=fs[:self.sspec.shape[1]]
		print "Length of freq samples is ", len(self.fspos)
		if nsam%2==0:
			self.fspos[-1]=-1*self.fspos[-1]
		shapefd=(puresig.shape[0],len(self.fspos))
		noisyspec=self.compute_noise(nt,shapefd,puresig.shape)
		self.newsig = np.fft.irfft(noisyspec)
		self.nsyampspec=np.abs(np.fft.rfft(self.newsig))

	def compute_noise(self,ntype,shapef,shapet):

		sampspec=np.abs(self.sspec)
		ntd=np.random.normal(0,1,shapet)
		nft=np.fft.rfft(ntd)
		self.nphspec=np.angle(nft)
		if ntype==1:
			# noise at Earth microseism frequencies
			#t=np.arange(shapereq[1])
			#tmat=np.zeros(shapereq)
			#for i in range(tmat.shape[0]):
			#	tmat[i,:]=t
			#noisemat=np.sin(2*np.pi*tmat/14)+2*np.sin(2*np.pi*tmat/7)
			#noisemat=noisemat+0.8*np.random.normal(0,1,shapereq)
			self.nampspec=np.ones(shapef)
			indeachst=np.argmax(sampspec,1)
			sd=0.02
			for i in range(self.nampspec.shape[0]):
				sumgauss=0.5*(np.exp(-(self.fspos)**2/(2*(sd)**2)))+np.exp(-(self.fspos-0.07)**2/(2*(sd)**2))+2*np.exp(-(self.fspos-0.14)**2/(2*(1.5*sd)**2))
				self.nampspec[i,:]+=sumgauss
				nasp=self.nampspec[i,indeachst[i]] # nasp stands for noise_at_signal_peak
				sp=sampspec[i,indeachst[i]]
				self.nampspec[i,:]=(self.pcent*sp/nasp)*self.nampspec[i,:]
		elif ntype==2:
			self.nampspec=np.abs(nft)
			for i in range(self.nampspec.shape[0]):
				sigmax=np.amax(sampspec[i,:])
				nmax=np.amax(self.nampspec[i,:])
				self.nampspec[i,:]=(self.pcent*sigmax/nmax)*self.nampspec[i,:]
		#newsigampspec = self.nampspec + sampspec
		#newsigspec = newsigampspec*(np.exp(1j*np.angle(self.sspec)))
		newsigspec = self.sspec + self.nampspec*(np.exp(1j*self.nphspec))
		return newsigspec

##################################################################################################################

class plot_data():
	
	""" Plot a record section, i.e. data from a linear seismic array in a way that honours station spacing """

	def __init__(self,master,sdist,snames,deltat,tstart,figtitle,pn):
	
		""" Generates a collection object compatible with matplotlib.collections """
	
		egdatamat=master[0]
		self.dataset_list=range(len(master))
		self.pn=pn
		max_amp=egdatamat.max()
		print "max_amp is: ", max_amp
		""" Normalize w.r.t. maximum amplitude across array """
		nmlzr=(-1/max_amp)*np.matrix(np.identity(egdatamat.shape[1]))
		# the minus sign here is so that the plot comes out correct- the plot in make_plot has its y-axis flipped
		# in the sense that largest epicentral distance comes at the bottom. So what is plotted is an upside down
		# flipped image of the collection object. This minus sign corrects for that.
		tsamples=egdatamat.shape[1]
		self.epdist=sdist
		self.time = np.arange(tstart,tstart+tsamples*deltat,deltat)
		#print "Number of samples and sampling interval are: ", tsamples,deltat
		#print "Time window for plotting: ", self.time[0], self.time[-1]
		seis_list = list()
		dlist = list()
		self.st_plotted = list()
		dscolors=[[0,0,0],[1,0,0],[0,0,1]]
		for ed,d in enumerate(self.dataset_list):
			try:
				#print "Performing multiplication: ", d, master[d].shape,  nmlzr.shape
				normdata=master[d]*nmlzr
			except ValueError:
				print "Shapes are: ", master[d].shape,nmlzr.shape
				sys.exit()
			for i in range(0,master[d].shape[0]):
				seismgm = np.array(normdata[i,:])[0]
				seismgm = 50*seismgm                   # a final amplification for plotting
				single_trace = zip(self.time,seismgm)  # for sideways plotting
				seis_list.append(single_trace)
				dlist.append((0,self.epdist[i]))
				self.st_plotted.append(snames[i])
			print "RGB colors are: ", dscolors[ed]
			col = clt.LineCollection(seis_list,offsets=dlist,color=dscolors[ed])
			self.dataset_list[d] = col
		self.make_plot(figtitle)

	def make_plot(self,caption):

		""" Plot the collection object(s) """
		plopt=1
		try:
			minspacing=min([self.epdist[xi+1]-xel for xi,xel in enumerate(self.epdist) if xi+1<len(self.epdist)])
		except ValueError:
			minspacing=50
		#print self.epdist
		ytop=self.epdist[0]-minspacing
		ybot=self.epdist[-1]+minspacing
		#print ytop,ybot
		st_locs=self.dataset_list[0].get_offsets()[:,1]
		#fig = plt.figure(figsize=(12,8))
		fig=plt.figure()
		ax = fig.add_subplot(111)
		#self.dataset_list[1].set_color('r')
		#self.dataset_list[0].set_color('b')
		for c,col in enumerate(self.dataset_list):
			ax.add_collection(col)
		#ax.yaxis.tick_right()
		ax.set_xlim([self.time[0],self.time[-1]]) #1100
		#ax.set_xlim(400,900) 
		ax.set_ylim([ybot,ytop])
		ax.set_xlabel("Time [s]")
		ax.set_ylabel("Distance from source [km]")#, labelpad=30)
		#ax.set_ylabel("Distance in model [km]")
		#ax.set_title(caption)
		#axr=ax.twinx()
		y1,y2=ax.get_ylim()
		x1,x2=ax.get_xlim()
		#axr.set_ylim(y1,y2)
		if plopt==1:
			try:
				namenum=["%5s %2d" %(stpl,(self.pn*nspp)+i+1) for i,stpl in enumerate(self.st_plotted)]
			except NameError:
				namenum=["%5s %2d" %(stpl,i+1) for i,stpl in enumerate(self.st_plotted)]
			#axr.set_yticks(st_locs)
			#axr.set_yticklabels(namenum,size=10)
		else:
			mev=5
			mst=np.array([(mev*i)-1 for i in range(1,(len(st_locs)+mev)/mev)])
			try:
				mstlocs=st_locs[mst]
			except IndexError:
				mstlocs=[st_locs[-1]]
			#print "mstlocs is ", mstlocs
			#axr.set_yticks(mstlocs)
			#axr.set_yticklabels(mst+1,size=10)
		plt.show()

###################################################################################################################

def fk_transform(dmat,epd,sit):
	nseis = dmat.shape[0]	
	lent=dmat.shape[1]
	fkt=np.fft.fftshift(np.fft.fft2(dmat))
	dx=epd[1]-epd[0]
	""" Check that spatial sampling is uniform """
	epd_dummy=np.arange(epd[0],epd[-1]+dx,dx)
	if not np.all(epd_dummy)==np.all(epd):
		print "Cannot take 2D Fourier transform as spatial sampling is uneven"
		return None
	kvals=np.fft.fftshift(2*np.pi*np.fft.fftfreq(nseis,dx))
	fvals=np.fft.fftshift(np.fft.fftfreq(lent,sit))
	X,Y=np.meshgrid(kvals,fvals)
	#print X.shape, Y.shape, fkt.shape
	fkt=np.transpose(fkt)
	plt.pcolor(X,Y,fkt)
	plt.ylim(0.1,0)
	plt.xlabel('Wavenumber')
	plt.ylabel('Frequency [Hz]')
	plt.show()

##################################################################################################################

def usage():
	print "What do you want to do ?"
	print "Choose from the following options"
	print "		 1 - View data"
	print "		 2 - Add synthetic noise to data"
	print "          3 - Exit"

##################################################################################################################	

def see_subset(desd,rdo):

	fd=desd[0]
	nspp=int(raw_input("Number of seismograms per plot: "))
	nstot=fd.shape[0]
	nplots=nstot/nspp
	if not nstot%nspp==0:
		nplots+=1
	#print "total, per plot and nplots: ", nstot, nspp, nplots
	for plotn in range(nplots):
		s1=nspp*plotn
		s2=min(nstot,s1+nspp)
		print "Plotting station numbers %d to %d, %d plots remaining" %(s1+1,s2,nplots-plotn-1)
		showdata=range(len(desd))
		for dd in range(len(desd)):
			showdata[dd]=desd[dd][s1:s2]
		showed=rdo.ed[s1:s2]
		shownames=rdo.stnames[s1:s2]
		pdobj=plot_data(showdata,showed,shownames,rdo.si,rdo.tstart,rdo.evname,plotn)

##################################################################################################################	

if __name__ == '__main__':
	ndir=len(sys.argv)-1
	dirnums=range(1,ndir)
	xnames=[]
	for dn in dirnums:
		xnames.append(sys.argv[dn])
	#datadir=sys.argv[1]
	try:
		fext=sys.argv[-1]
	except IndexError:
		fext=raw_input("Extension for sac files in directory (Enter eg: .sac, .SAC, .LHZ ): ")
	desired_output = 0
	while desired_output != 3:
		usage()
		desired_output = int(raw_input("Enter your choice (number) here: "))
		if desired_output>=3:
			sys.exit('Thank you')
		elif desired_output==1:
			usrc_ww=int(raw_input("See entire data (1) or a windowed portion (2) ? "))
			if usrc_ww==1:
				rdobj=read_data(xnames,fext,'whole')
			elif usrc_ww==2:
				rdobj=read_data(xnames,fext,'win')
			#reldata=getattr(rdobj,'fulldata') if usrc_ww==1 else getattr(rdobj,'windata')
			desireddata=rdobj.arraydata
			print "No. of directories is: ", len(desireddata)
			print "No. of seismograms in first directory is ", desireddata[0].shape[0]
			usrc_hm=raw_input("See only few stations at a time ? (y/n): ")
			if usrc_hm=='y':
				see_subset(desireddata,rdobj)
			elif usrc_hm=='n':		
				pdobj=plot_data(desireddata,rdobj.ed,rdobj.stnames,rdobj.si,rdobj.tstart,rdobj.evname,0)
			else:
				sys.exit('Thank you\n')
		elif desired_output==2:
			rdobj=read_data(xnames[0],fext,'whole')
			# addition of noise will only be done on one directory (set of seismograms) at a time
			usrcn=int(raw_input("Enter 1 for Earth microseism noise or 2 for random Gaussian noise: " ))
			anobj=add_noise(rdobj.arraydata[0],rdobj.si,usrcn)
			fig = plt.figure()
			#figtemp = plt.figure()
			#axtemp = figtemp.add_subplot(111)
			axn = fig.add_subplot(111)
			showsspec = np.abs(anobj.sspec[0,:])
			shownspec = anobj.nampspec[0,:]
			shownsigspec = anobj.nsyampspec[0,:]
			fsnspec = anobj.fspos
			#axtemp.plot(fsnspec,np.unwrap(anobj.nphspec[0,:]))
			#axn.set_title("Amplitude spectrum of noise added")
			axn.plot(fsnspec,showsspec,label='Original signal')
			axn.plot(fsnspec,shownspec,color='g',label='Noise added')
			axn.plot(fsnspec,shownsigspec,label='Noisy signal')
			axn.set_xlabel('Frequency [Hz]')
			axn.set_ylabel('Amplitude')
			axn.set_xlim(0.005,0.2)
			axn.legend()
			#showdata=range(2)
			#showdata[0]=rdobj.arraydata[0]
			#showdata[1]=anobj.newsig
			# use above 3 lines if you want to plot original data and the artificially noisy data on the same plot
			# if you want to plot only the artificially noisy data, use following line
			showdata=[anobj.newsig]
			usrc_hm=raw_input("See only few stations at a time ? (y/n): ")
			if usrc_hm=='y':
				see_subset(showdata,rdobj)
			elif usrc_hm=='n':		
				pdobj=plot_data(showdata,rdobj.ed,rdobj.stnames,rdobj.si,rdobj.tstart,rdobj.evname,0)
			#fk_transform(rdobj.fulldata,rdobj.ed,rdobj.si)
