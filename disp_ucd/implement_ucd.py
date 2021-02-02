#!/usr/bin/python

###################################################################################
# Program to perform UC-diagram analysis on data from a seismic array
#
# Input to program: Directory containing sac files from a single station-array
#
# Arjun Datta, March 2014
#
# Algorithm and notation based on Duputel et al. 2010:
#
# Improving the analysis and inversion of multimode Rayleigh-wave dispersion by using
# group-delay time information observed on arrays of high-frequency sensors
# Geophysics, Vol. 75, No. 2, 2010
###################################################################################

# Available Python modules
import os
import sys
import math
import obspy.core
import numpy as np
import scipy.signal as ss
import matplotlib.pyplot as plt

# sys.path.append('../modules_common')
sys.path.append(os.path.expanduser('~/code_general/modules.python'))
# path to extraneous modules

# Modules written by me
import ucd_supplement as ucs
import array_proc.seisarray_data as sad
import SW1D_earthsr.read_earthsr_io as reo
# import read_surf96_io as rs96

# Command line arguments
script_dir=os.getcwd()
data_dir=sys.argv[1]
fid=sys.argv[2]
try:
	ftype=sys.argv[3]
except IndexError:
	ftype=None

################################################################################

if not os.path.exists(data_dir):
	sys.exit("Terminating program - No such directory -%s" %(data_dir))

################################################################################
# User input
################################################################################
finput=input("Enter frequency at which or frequency range in which to perform analysis: ")
#finput='0.01 0.08'
numf=len(finput.split())
if numf==1:
	freqs=[]
	freqs.append(float(finput))
elif numf==2:
	flow=float(finput.split()[0])
	fhigh=float(finput.split()[1])
	fstep=float(input("Enter frequency increment: "))
	freqs=np.arange(flow,fhigh+fstep,fstep)
	#freqs=np.insert(freqs,1,0.005)        # TEMPORARY CHANGE FOR TA DATA ANALYSIS (JAN 2016)
	print("Frequency values for analysis: ", freqs)
crange = input("Phase velocity range: ")
urange = input("Group velocity range: ")
totpick = int(input("Total number of modes to pick ? "))
cmin = float(crange.split()[0])
cmax=float(crange.split()[1])
umin=float(urange.split()[0])
umax=float(urange.split()[1])
if fid.endswith('SAC') or fid.endswith('sac'):
	# means the program is dealing with synthetic data
	#usrc=input("See theoretical dispersion too ? (y/n): ")
	usrc='y'
	ml=0
	mh=totpick-1
	if usrc=='y':
		#showmnums=input("Start and end mode numbers to plot: ")
		showmnums="0 5"
		ml=int(showmnums.split()[0])
		mh=int(showmnums.split()[1])
		#minc=int(input("Incident mode to highlight: "))
		#thdpfile=[input('File containing theoretical dispersion: ')]
		thdpfile=["../disp.bg.ray.gz"]
		#try:
		theor_disp = reo.read_disp(thdpfile,ml,mh)
		theor_cdisp = theor_disp.modcdisp[0]
		theor_udisp = theor_disp.modudisp[0]
		#except IndexError:
		#	rs96obj = rs96.read_disp('cdisp_surf96.ray')
		#	theor_cdisp = rs96obj.disp[0]
		#	rs96obj = rs96.read_disp('udisp_surf96.ray')
		#	theor_udisp = rs96obj.disp[0]
	else:
		theor_cdisp=0
		theor_udisp=0
else:
	usrc=None
synth_case=True if usrc=='y' else False
#usrc_st=input("Use all stations (y/n) ?: ")
usrc_st='n'

################################################################################
# Reading the data
################################################################################

data_dir=os.path.normpath(os.path.join(script_dir,data_dir))
dlistsad=[]
dlistsad.append(data_dir)
#sadobj=sad.read_data(dlistsad,fid,'win',ftype)
sadobj=sad.read_data(dlistsad,fid,'whole',ftype)
if usrc_st=='y':
	try:
		tdodata=sadobj.windata
	except AttributeError:
		tdodata=sadobj.fulldata
	epdist=sadobj.ed
else:
	#st_range=input("Start and end station numbers: ")
	st_range="17 58"
	s1=int(st_range.split()[0])
	s2=int(st_range.split()[1])
	try:
		tdodata=sadobj.windata[s1-1:s2]
	except AttributeError:
		tdodata=sadobj.fulldata[s1-1:s2]
	epdist=sadobj.ed[s1-1:s2]
anobj=sad.add_noise(tdodata,sadobj.si,1)
td_data=anobj.newsig
# comment the below line in order to run the code with artificially added noise
td_data=tdodata
si=sadobj.si
#begint=sadobj.tstart
begint=int(sadobj.tstart)

################################################################################
# Assignment of main parameters
################################################################################
nrec=td_data.shape[0]					# number of stations in the network
nsamp=td_data.shape[1]					# number of time samples in the records
delta_x = epdist[1] - epdist[0]				# interstation spacing
aperture=epdist[-1] - epdist[0]
k_nyquist = 1/(2*delta_x)		# *** Remember this is an approximation as station spacing is usually not uniform ***
kstep=1/(nrec*delta_x)
kvalues = 2*np.pi*np.arange(-k_nyquist,k_nyquist,kstep) # wavenumber
t=np.arange(begint,begint+(nsamp*si),si)		# time samples
fs=np.fft.fftshift(np.fft.fftfreq(nsamp,si))		# FFT frequency samples (Hz)
omega=2*np.pi*fs					# angular frequency (rad/s)
xbar=np.mean(epdist)					# Mean epicentral distance of the network
print("Number of frequency samples is: ", len(fs))

################################################################################
# Preliminary processing
################################################################################
# ***** Generate weight function for the stations to attenuate secondary lobes
def gauss(x,mu,sigma):
	return np.exp(-0.5*((x-mu)/sigma)**2)
wn=gauss(epdist,(epdist[-1]+epdist[0])/2,aperture/4)
#wn=np.ones((nrec))
wn_uss=ss.gaussian(nrec,nrec/4)                         # can be used only in case of uniformly spaced stations
#***** Remove mean
xtdata=ss.detrend(td_data)
#***** Normalize each trace with respect to greatest dynamic range in entire data
xtdata=np.matrix(xtdata)
maxet=np.matrix.max(xtdata,axis=1)
minet=np.matrix.min(xtdata,axis=1)
ret=maxet-minet
divisor=ret*(np.matrix(np.ones((1,nsamp))))
xtdata=xtdata/divisor
#***** Fourier transform with respect to time
xwdata=np.fft.fftshift(np.fft.fft(xtdata),axes=1)

#####################################################################################################

class do_single_frequency():

	""" Implements the UC-diagram processing technique for a single frequency """

	def __init__(self,fhz,prevuc=None):

		print("Working on frequency %f Hz" %(fhz))
		self.omega0=2*np.pi*fhz
		self.ktrials = np.linspace(self.omega0/cmax,self.omega0/cmin,100)# trial wavenumbers for this frequency
		par=0.9
		alpha = (4*(math.log(10**1.5)))/((self.omega0*par)**2) 		# for a filter whose total bandwidth at 30dB is equal to its central frequency
		self.gf=np.exp(-alpha*(omega-self.omega0)**2) 			# Gaussian filter
		heaviside=lambda x: 0 if x < 0 else 1
		self.H=np.array(list(map(heaviside,omega)))
		self.H=self.H.reshape(1,len(self.H))
		self.gf=self.gf.reshape(1,len(self.gf))
		#********** Mute the negative frequencies and apply the Gaussian filter on each trace ********#
		muting=np.dot(np.ones((nrec,1)),self.H)
		filtering=np.dot(np.ones((nrec,1)),self.gf)
		ansig=xwdata*muting*filtering					# analytical signal
		self.uc_est=self.taup(ansig)
		if fhz<0.02:			# NB: TEMPORARY CHANGE FOR TA DATA ANALYSIS
			self.uc_est=self.taup(ansig)
		elif fhz<0.03:
			self.uc_est=4.5
		else:
			self.uc_est=4.25
		if prevuc != None and abs(self.uc_est-prevuc)>0.4*prevuc:
			self.uc_est=prevuc
		print("Using Uc = %f units" %(self.uc_est))
		self.cgat=xbar/self.uc_est
		self.main_processing()
		self.build_ucd()

	def taup(self,anasig):

		""" estimates the central group velocity of the multi-mode wave packet """

		xt_modif=np.matrix(abs(np.fft.ifft(np.fft.ifftshift(anasig,axes=1))))
		div=np.matrix.max(xt_modif,axis=1)*np.matrix(np.ones((1,nsamp)))
		xtmodif=xt_modif/div
		if begint>0:
			bz=np.zeros((nrec,begint))
			xtmodif=np.concatenate((bz,xtmodif),axis=1)
		od=math.log(umin,10)
		du=10**(math.floor(od)-1)
		uc=np.arange(umin,umax+du,du)
		tau=np.arange(0,(nsamp+1)*si,si)
		#tau=np.arange(begint,begint+((nsamp+1)*si),si)
		tau=tau.reshape(1,len(tau))
		taumat=np.dot(np.ones((nrec,1)),tau)
		tp=np.zeros((len(uc),tau.shape[1]))
		#taumat=taumat-taumat[0][0]
		# print("taumat is ", taumat[10,:10])
		offsets=epdist.reshape(len(epdist),1)
		for l in range(len(uc)):
			tt = taumat + np.dot(offsets/uc[l],np.ones((1,tau.shape[1])))
			tti=np.around(tt/si)
			tti=tti.astype(int)
			#if l==0:
			#	print("tt is: ", tt)
			counter=0
			for m in range(nrec):
				xindices=np.where(tti[m,:]<=nsamp)[0]
				q=tau.shape[1]-len(xindices)
				z=np.zeros((1,q))
				col=tti[m,xindices]-1
				col=col.reshape(1,len(xindices))
				tp[l,:]=tp[l,:] + np.concatenate((xtmodif[m,col],z),axis=1)
				#if l==0 and m==0:
					# print(tp.shape)
				if q>0:
					counter+=1
			if counter>0:
				tp[l,:]=tp[l,:]/counter
		locofmax=np.where(tp==tp.max())
		return uc[locofmax[0][0]]

	def main_processing(self):

		""" calculates the G(k,w) function: denoted here by wkdata """

		ltte=(omega-self.omega0)/self.uc_est	# Linear Term in the Taylor Expansion of k(w)
		ltte.reshape(1,len(ltte))
		apst=(omega*self.cgat)
		# apst stands for Additional Phase Shift Term - this term does not appear in the
		# formula for G(k,w) in the paper. It is this term that produces a time shift of
		# xbar/uc so that the time samples in the inverse F.T. (self.tkdata) are shifted
		# from the original time samples in the data by the amount equal to xbar/uc
		apst=apst.reshape(1,len(apst))
		self.wkdata=np.zeros((len(self.ktrials),len(ltte)))
		for rec in range(nrec):
			knot_term=np.exp(1j*self.ktrials*epdist[rec])
			knot_term=knot_term.reshape(len(knot_term),1)
			rest=xwdata[rec,:]*self.gf*self.H*np.exp(1j*(ltte*epdist[rec]-apst))*wn[rec]
			self.wkdata = self.wkdata + np.dot(knot_term,rest)
		# print(np.where(self.wkdata[0,:]>0))

	def build_ucd(self):

		""" generates the final content of the UC diagram """

		self.tkdata=np.abs(np.fft.ifft(np.fft.ifftshift(self.wkdata,axes=1)))
		mini=self.tkdata.min()
		le=(self.tkdata-mini).max()
		self.tkdb=20*np.log10(self.tkdata/le)

#####################################################################################################

def make_ucd_plot(frad,wavnum,toplot,ucest,showucd):

	# print("Time samples are: ", t)
	# ***** Calculate array response function
	theorex=True
	wn_mup=wn.reshape(1,len(wn))
	k_mup=wavnum.reshape(len(wavnum),1)
	k_mup=k_mup-np.mean(k_mup)
	x_mup=epdist.reshape(1,len(epdist))
	rk=np.sum(np.dot(np.ones((len(wavnum),1)),wn_mup)*np.exp(1j*np.dot(k_mup,x_mup)),axis=1)
	rkdb=20*np.log10(abs(rk)/nrec)
	# print("Shape of response function ", rk.shape)
	#ax_ucd=plt.subplot(111)
	# for array response function plotted on right use following:
	ax_rk=plt.subplot2grid((1,5),(0,4))
	ax_ucd=plt.subplot2grid((1,5),(0,0),colspan=4)
	# for array response function plotted on top use following:
	#ax_rk=plt.subplot2grid((6,1),(0,0),rowspan=1)
	#ax_ucd=plt.subplot2grid((6,1),(1,0),rowspan=5)
	# for not plotting the array response function at all use the following:
	#fig=plt.figure()
	#ax_ucd=fig.add_subplot(111)
	# Uncomment below to plot a legend for the contours
	#fig.colorbar(cs)
	thisf=float('%.4f' %(frad/(2*np.pi)) )
	cs=ax_ucd.contour(t,wavnum,toplot,levels=np.arange(contourmin,0))
	pv=np.linspace(cmin,cmax,10)
	gv=np.linspace(umin,umax,6)
	X=xbar/gv
	Y=frad/pv
	ax_ucd.set_xlabel('U: Group Velocity')
	ax_ucd.set_ylabel('C: Phase Velocity')
	ax_ucd.set_xticks(X)
	ax_ucd.set_xticklabels(map(lambda z: "%.2f" %z, gv))
	ax_ucd.set_yticks(Y)
	ax_ucd.set_yticklabels(map(lambda z: "%.2f" %z, pv))
	ax_ucd.set_xlim(xbar/umin,xbar/umax)
	ax_ucd.set_ylim(frad/cmin,frad/cmax)
	#ax_ucd.set_title('%.2f km/s' %(ucest))
	ax_ucd.set_title('%.3f Hz, Uc=%.2f km/s' %(thisf,ucest))
	#ax_ucd.grid(True)
	if synth_case:
		for tm in range(len(theor_cdisp)):
			try:
				ctm = [b for a,b,c in theor_cdisp[tm] if a==thisf][0]
				utm = [q for p,q,r in theor_udisp[tm] if p==thisf][0]
			except ValueError:
				try:
					ctm = [b for a,b in theor_cdisp[tm] if a==thisf][0]
					utm = [q for p,q in theor_udisp[tm] if p==thisf][0]
				except IndexError:
					theorex=False
					pass
			except IndexError:
				theorex=False
				print("Could not find theoretical values for frequency: ", thisf)
				pass
			# print("Values for mode %d are %.3f and %.3f" %(tm+ml,ctm,utm))
			if theorex:
				utmx = xbar/utm
				ctmy = frad/ctm
				# following marks the location of the theoretical mode c/u values
				ax_ucd.plot(utmx,ctmy,'ko',ms=5,mew=2.0)
				dc=ax_ucd.transData.transform((utmx,ctmy))
				inv=ax_ucd.transAxes.inverted()
				ac=inv.transform((dc[0],dc[1]+4))
				tms="%d" %(tm)
				# following controls where exactly the mode numbers are marked
				ax_ucd.text(ac[0],ac[1],tms,transform=ax_ucd.transAxes,fontsize=14)
	# following 8 lines need to be uncommented in order to plot the response function
	ax_rk.plot(rkdb,wavnum)
	ax_rk.set_ylim(frad/cmin,frad/cmax)
	ax_rk.yaxis.tick_right()
	ax_rk.xaxis.tick_top()
	ax_rk.set_xlabel('r(k)')
	rk_xtauto=ax_rk.get_xticks()
	rk_xtshow=(rk_xtauto[1],rk_xtauto[-2])
	ax_rk.set_xticks(rk_xtshow)
	if showucd:
		figpath=None
		stlocs_im=np.linspace(epdist[0],epdist[-1],nrec)
		figwts=plt.figure()
		axwts=figwts.add_subplot(111)
		axwts.plot(stlocs_im,wn_uss,'-o',label='Ideally spaced stations')
		axwts.plot(epdist,wn_uss,'-o',label='Scipy signal gaussian')
		axwts.plot(epdist,wn,'-o',label='Actually used')
		axwts.set_ylabel('Station Weight')
		axwts.set_xlabel('Epicentral Distance')
		plt.legend(loc=8)
	else:
		figname='ucd_%.3f_hz.png' %(thisf)
		figpath=data_dir+'/'+figname
		#plt.close(fig)
	return ax_ucd,figpath

#####################################################################################################
# Apply UC-diagram analysis to each frequency and pick dispersion
#####################################################################################################

contourmin=-29
pv_allf=[]
gv_allf=[]
for fnum,f0 in enumerate(freqs):
	# f0 is the current frequency (centre of gaussian filter)
	if len(freqs)>1:
		if fnum<=1:
			dsfobj=do_single_frequency(f0)
			ucused=dsfobj.uc_est
		else:
			dsfobj=do_single_frequency(f0,ucused)
			ucused=dsfobj.uc_est
		#currax,savenm=make_ucd_plot(dsfobj.omega0,dsfobj.ktrials,dsfobj.tkdb,dsfobj.uc_est,False)
		print("Done. Picking modes on the UC diagram...")
		peaks=ucs.get_peaks(dsfobj.tkdb,umin,umax,xbar,t,contourmin,totpick,False)
		pv_thisf=[None for i in range(totpick)]
		gv_thisf=[None for i in range(totpick)]
		for e,m in enumerate(peaks):
			mc=dsfobj.omega0/dsfobj.ktrials[m[0]]
			mu=xbar/t[m[1]]
			pv_thisf[e]=mc
			gv_thisf[e]=mu
		#	currax.plot(t[m[1]],dsfobj.ktrials[m[0]],'k+',ms=25,mew=2.0)
		#plt.savefig(savenm)
		# print(pv_thisf)
		# print(gv_thisf)
		pv_allf.append(pv_thisf)
		gv_allf.append(gv_thisf)
	elif len(freqs)==1:
		dsfobj=do_single_frequency(f0)
		currax,savenm=make_ucd_plot(dsfobj.omega0,dsfobj.ktrials,dsfobj.tkdb,dsfobj.uc_est,True)
		peaks=ucs.get_peaks(dsfobj.tkdb,umin,umax,xbar,t,contourmin,totpick,True)
		for pk in peaks:
			pvel=dsfobj.omega0/dsfobj.ktrials[pk[0]]
			gvel=xbar/t[pk[1]]
			print("Picked phase and group velocity: ", pvel, gvel)
			# this is where you mark the picked values with '+' signs
			currax.plot(t[pk[1]],dsfobj.ktrials[pk[0]],'k+',ms=25,mew=2.0)
		plt.show()
if len(freqs)>1:
	for xfig in range(2):
		# making 2 figures - one each for phase and group velocity
		dispfig=plt.figure()
		axdisp=dispfig.add_subplot(111)
		#axdisp.grid(True)
		if xfig==0:
			figname='ucd_cdisp.eps'
			axdisp.set_ylabel('Phase Velocity [km/s]')
			#axdisp.set_ylim(cmin,cmax+0.5)
			axdisp.set_ylim(cmin,cmax)
			whichvel='pvd'
			pw=pv_allf
			lloc=1
		else:
			figname='ucd_udisp.eps'
			axdisp.set_ylabel('Group Velocity [km/s]')
			axdisp.set_ylim(umin,umax)
			whichvel='gvd'
			pw=gv_allf
			lloc=4
		if synth_case: #and xfig==0: # TEMPORARY. PLEASE CHANGE.
			mcol=['b','g','r','c','m','y','k','b','g','r']
			print("Length of theor_cdisp is: ", len(theor_cdisp))
			solidcurve = theor_cdisp if xfig==0 else theor_udisp
			print("Length of solidcurve is: ", len(solidcurve))
			for k in range(len(solidcurve)):
				print("k is ", k)
				try:
					f = [x for x,y,z in solidcurve[k]]
					v = [y for x,y,z in solidcurve[k]]
				except ValueError:
					f = [x for x,y in solidcurve[k]]
					v = [y for x,y in solidcurve[k]]
				amn=theor_disp.rel_modes[k]
				curve_name="Mode %d" %(amn)
				#if amn==minc:
				axdisp.plot(f,v,'-',color=mcol[amn],label=curve_name)
				#else:
				#	axdisp.plot(f,v,'--',color=mcol[amn],label=curve_name)
		for mdet in range(totpick):
			pq=[pvf[mdet] for pvf in pw]
			#if mdet<(totpick-1):
			curve_name="UCD pick %d" %(mdet+1)
			axdisp.plot(freqs,pq,'*',label=curve_name)
		try:
			pl=len(solidcurve) #+totpick
		except NameError:
			pl=totpick
		if pl<=5:
			legcols=1
			legsp=1
		elif pl<=10:
			legcols=2
			legsp=0.25
		else:
			legcols=3
			legsp=0.25
		axdisp.legend(loc=lloc,labelspacing=legsp,ncol=legcols)
		axdisp.set_xlim(freqs[0],freqs[-1])
		#axdisp.set_title('Aperture %d km' %(int(aperture)))
		axdisp.set_xlabel('Frequency [Hz]')
		if legcols > 1:
			leg = plt.gca().get_legend()
			ltext  = leg.get_texts()
			#plt.setp(ltext, fontsize='small')
		plt.savefig(data_dir+'/'+figname)
if len(freqs)>1:
	#usrc2=input("Do you want to pickle the results ? (y/n) ")
	usrc2='y'
	if usrc2=='y':
		ucs.make_pickle(data_dir,freqs,pv_allf,gv_allf)
