#!/usr/bin/python

import os
import sys
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Modules written by me
import read_earth_io as reo
import read_surf96_io as rs96

#####################################################################################################

def get_peaks(inmat,gvmin,gvmax,dist,ts,thresh,npick,verbose):
	
	""" First you need to determine the (very small) portion of the input matrix that constitutes 
	the UC diagram - that small part of the entire time series over which group velocity values
	are meaningful """
	t1=dist/gvmax
	t2=dist/gvmin
	relind=np.intersect1d(np.where(ts>t1)[0],np.where(ts<t2)[0])
	ucdim=inmat[:,relind[0]:relind[-1]]
	wwidth=max(7,int(math.ceil(0.001*ucdim.shape[1])))
	wheight=int(math.ceil(0.05*ucdim.shape[0]))
	if wwidth%2==0:
		wwidth+=1
	if wheight%2==0:
		wheight+=1
	# ******** Define the edges of the sub image over which the window will move ******
	le=relind[0] #(wwidth-1)/2
	re=relind[-1] #inmat.shape[1]-le-1
	te=(wheight-1)/2
	be=ucdim.shape[0]-te-1
	relim=inmat[te:be,le:re]
	locgm=np.where(relim==relim.max())
	rowgm=locgm[0][0]+te
	colgm=locgm[1][0]+le#+relind[0]
	if verbose:
		print relind[0],relind[-1]
		print "Shape of original matrix: ", inmat.shape
		print "Location of global maxima of UC diagram: ", rowgm, colgm
		print "Shape of relevant matrix: ", ucdim.shape
		print "Window size is ", wwidth, wheight
		print "Edges wrt original matix are ",le,re,te,be
	r=te
	prosp=0
	peaklocs=[]
	peakv=[]
	while r<=be:
		c=le
		while c<=re:
			wmid=float("%.5f" %(inmat[r,c]))
			win=inmat[r-te:r+te+1,c-(wwidth-1)/2:c+(wwidth-1)/2+1]
			wmax=float("%.5f" %(np.amax(win)))
			if wmid>(0.8*thresh) and wmax==wmid:
				if verbose:
					print "Found peak at %d,%d with value %f" %(r,c,wmid)
				prosp+=1
				""" Don't consider those peaks which correspond to higher phase velocity but lower group velocity than the global maxima """
				if not (r<rowgm and c>colgm):
					peakv.append(wmid)
					peaklocs.append([r,c])
				c+=(wwidth-1)/2
			elif wmax==wmid: 
				# it is a peak but of small amplitude
				c+=(wwidth-1)/2
			else:
				c+=1
		r+=1
	""" Sort the peaks in ascending order of magnitude """
	peaklocs=[b for a,b in sorted(zip(peakv,peaklocs))]
	
	""" Filter out those peaks that occur very close to the global maxima """
	# the way I do this is very ad-hoc, must replace with something more theoretically sound
	neargmc=range(colgm-wwidth,colgm+wwidth)
	neargmr=range(rowgm-wwidth,rowgm+wwidth)
	relpeaks=[[i,j] for i,j in peaklocs if (not i in neargmr) and (not j in neargmc)]
	if verbose:
		print "Sorted peak locations: ", peaklocs
		print "Location of prospective peaks other than global maxima: ", relpeaks
	
	""" Ensure you pick only npick number of modes """
	if len(relpeaks)>(npick-1):
	# drop the lowest amplitude picks that are superfluous
		extra=len(relpeaks)-(npick-1)
		relpeaks=relpeaks[extra:]
	
	""" Make sure the global maxima is included """
	#relpeaks.insert(0,[rowgm,colgm])
	# modif. ARJUN Aug 2018: global maxima needs to go to the END of the list
	# as it is in ASCENDING order.
	relpeaks.append([rowgm,colgm])
	
	""" Discard any peak that is extremely close, in phase velocity, to another one of greater amplitude """
	pn=0
	while pn<len(relpeaks)-1:
		tp=relpeaks[pn]
		rtp=tp[0]
		ctp=tp[1]
		rpks=relpeaks[pn+1:]
		print "Comparing ", tp, "with ", rpks
		for i,op in enumerate(rpks):
			rdiff=abs(rtp-op[0])
			if rdiff<3:
				if verbose:
					print "Discarding peak at %d,%d" %(rtp,ctp)
				del relpeaks[pn]
				# tp will be lower in magnitude between the two because the peaks have already been sorted in ascending order of magnitude
				break
		pn+=1
	
	# modif. ARJUN August 2018: finally, sort the selected peaks in descending order of magnitude
	relpeaks=relpeaks[::-1]
	if verbose:
		print "Location of finally selected peaks: ", relpeaks
	return relpeaks
	
#####################################################################################################
# Function to pickle (store) UCD dispersion results
#####################################################################################################

def make_pickle(wdir,fsamples,cdisp,udisp):

	""" Function to pickle UCD dispersion results. Offers the option of storing entire dispersion curves or only parts of them as specified by the user"""

	parentdir=os.path.dirname(wdir)
	evdir=os.path.basename(wdir)
	#pickledir=raw_input("Name of pickle directory: ")
	pickledir='.'
	storethis=zip(fsamples,cdisp)
	jarname=evdir+"_c"+".pckl"
	jarfile=os.path.join(wdir,pickledir,jarname)
	print "pickledir is: ", pickledir
	print "Name of pickle file: ", jarfile
	jarc=open(jarfile,'w')
	pickle.dump(storethis,jarc)
	jarc.close()
	storethis=zip(fsamples,udisp)
	jarname=evdir+"_u"+".pckl"
	jarfile=os.path.join(wdir,pickledir,jarname)
	jaru=open(jarfile,'w')
	pickle.dump(storethis,jaru)
	jaru.close()

#####################################################################################################
# Class to perform various operations with previously stored results (pickles)
#####################################################################################################

class workon_pickle():

	def __init__(self,pfilelist,todo):
		
		if len(pfilelist)==1:
			jar=open(pfilelist[0])
			print "Reading ", pfilelist[0]
			cookie = pickle.load(jar)
			self.fsam=np.array([float("%.4f" %(x)) for x,y in cookie])
			self.vel=np.array([y for x,y in cookie])
			self.npicks=self.vel.shape[1]
			jar.close()
		self.totp=len(pfilelist)
		self.plist=pfilelist
		if todo==1:
			self.view_pickle()
		elif todo==2:
			fsam, pickmatrix = self.pick_modes()
			newpfile='p'+pfilelist[0]
			jarnew=open(newpfile,'w')
			pickle.dump(fsam,jarnew)
			pickle.dump(pickmatrix,jarnew)
			jarnew.close()
			print "New pickle file created"
		elif todo==3:
			self.write_pickle_surf96()
		#jar.close()

	def view_pickle(self):
	
		#usrc=raw_input("See theoretical dispersion too ? (y/n): ")
		usrc='y'
		mcol=['b','g','r','c','m','y','k','b','g','r']
		ncols=1
		r=self.totp/ncols if self.totp%ncols==0 else (self.totp/ncols)+1
		usr_t=raw_input("Enter title of plot: ")
		#usr_t="whatever"
		if r==1:
			sizy=4.75
		elif r==2:
			sizy=9
		elif r==3:
			sizy=12
		fig=plt.figure()
		#fig=plt.figure(figsize=(12,sizy))
		#fig.suptitle(usr_t)
		for i,pckl in enumerate(self.plist):
			jar=open(pckl)
			cookie = pickle.load(jar)
			fsam=np.array([float("%.4f" %(x)) for x,y in cookie])
			vel=np.array([y for x,y in cookie])
			npicks=vel.shape[1]
			ax=fig.add_subplot(r,ncols,i+1)
			if usrc=='y' and len(thdpfile)>0:
				#for k in range(len(solidcurve)):
				for k,mode in enumerate(reoobj.rel_modes):
					print "mode number: ", k
					try:
						f = [x for x,y,z in solidcurve[k]]
						v = [y for x,y,z in solidcurve[k]]
					except ValueError:
						f = [x for x,y in solidcurve[k]]
						v = [y for x,y in solidcurve[k]]
					curve_name="Mode %d" %mode
					ax.plot(f,v,'-',color=mcol[mode],label=curve_name)
			for pick in range(npicks):
				cpick=vel[:,pick]
				curve_name="UCD pick %d" %(pick+1)
				ax.plot(fsam,cpick,'*',ms=8,label=curve_name)
			try:
				pl=npicks+len(solidcurve)
			except NameError:
				pl=npicks
			#pl=npicks
			if pl<=4:
				legcols=1
				legsp=1
			elif pl<=10:
				legcols=2
				legsp=0.25
			else:
				legcols=3
				legsp=0.25
			ax.set_xlabel('Frequency [Hz]')
			#ax.set_ylabel('Group Velocity [km/s]')
			#ax.set_ylim(2,6)
			ax.set_ylabel('Phase Velocity [km/s]')
			ax.set_ylim(3,7.5)
			ax.set_xlim(0.005,0.08)
			ax.set_title(usr_t)
			#plt.grid(True)
			if i==0:
				ax.legend(loc=1,labelspacing=legsp,ncol=legcols,prop={'size':12})
			jar.close()
		plt.show()
		#plt.savefig('ucd_pickle_plot.eps')

	def pick_modes(self):

		print "Original cookie: "
		for i in range(len(self.fsam)):
			print self.fsam[i], self.vel[i]
		usrc_pick='y'
		fsamples_modes=[]
		picks_modes=[]
		picked_mn=[]
		while usrc_pick=='y':
			pmn=int(raw_input("Mode number to pick: "))
			mfr=raw_input("Frequency range for this pick: ")
			mf1=float(mfr.split()[0])
			mf2=float(mfr.split()[1])
			find1=np.where(self.fsam==mf1)[0][0]
			find2=np.where(self.fsam==mf2)[0][0]
			frel=self.fsam[find1:find2+1]
			#print frel
			mode_picks=range(len(frel))
			try:
                        	ftheor = [x for x,y,z in solidcurve[pmn] if x in frel]
                                vtheor = [y for x,y,z in solidcurve[pmn] if x in frel]
                        except ValueError:
                        	ftheor = [x for x,y in solidcurve[pmn] if x in frel]
                                vtheor = [y for x,y in solidcurve[pmn] if x in frel]
			if frel[0]==ftheor[-1]:
				# this means your theoretical file has a frequency ordering opposite to that of your pickle file
				# (one is in ascending, the other in descending order of frequency)
				vtheor=vtheor[::-1]
			#print "Theoretical values are: ", zip(frel,vtheor)
			for j,freq in enumerate(frel):
				#print freq, self.vel[j+find1], vtheor[j]
				if np.all(self.vel[j+find1]>0):
					discrep=np.abs(self.vel[j+find1]-vtheor[j])
				else:
					dtmp=range(len(self.vel[j+find1]))
					for i in range(len(self.vel[j+find1])):
						try:
							dtmp[i]=np.abs(self.vel[j+find1][i]-vtheor[j])
						except TypeError:
							# means you've encountered a pickle value = None
							dtmp[i]=100
					discrep=np.array(dtmp)
				#print "Frequency, differences: ", freq, discrep
				mode_picks[j] = self.vel[j+find1][np.argmin(discrep)]
			plt.plot(frel,mode_picks,'*')
			plt.xlim(0.004,0.062)
			plt.ylim(3,7)
			plt.show()
			fsamples_modes.append(frel)
			picks_modes.append(mode_picks)
			picked_mn.append(pmn)
			usrc_pick=raw_input("Pick another mode ? (y/n): ")
		print "Picked %d modes " %(len(picked_mn))
		allf=[]
		for fsam_mode in fsamples_modes:
			allf=np.union1d(allf,fsam_mode)

		#print "allf: ", allf

		# all picks made, now build final matrix of picks
		# this is a matrix with columns = modes
		#                   and rows = frequency samples
		all_picks=np.zeros((len(allf),6))
		all_picks[:]=-1
		for k,col in enumerate(picked_mn):
			rrtc=np.searchsorted(allf,fsamples_modes[k])
			# rrtc = rel_rows_this_col
			#print "rel indices for mode %d " %(col), rrtc
			all_picks[rrtc,col]=picks_modes[k]
			
		print "Final answer: "
		print all_picks

		return allf, all_picks
		
		#for j in range(self.npicks):
		#	rf=raw_input("Frequency range for pick %d: " %(j+1))
		#	f1=float(rf.split()[0])
		#	f2=float(rf.split()[1])
		#	for k in range(len(self.fsam)):
		#		if self.fsam[k]<f1 or self.fsam[k]>f2:
		#			self.vel[k][j]=None
		#print "Modified cookie: "
		#for i in range(len(self.fsam)):
		#	print self.fsam[i], self.vel[i]
		#new_cookie=zip(self.fsam,self.vel)
		#return new_cookie

	def write_pickle_surf96(self):

		for i in range(self.npicks):
			for j in range(len(self.fsam))[::-1]:
				if self.vel[j,i] != None:
					print "%s %s %s %s  %d	%9.5f	%.2f	%.5f" %('SURF96','R','C','T',i,1/(self.fsam[j]),self.vel[j,i],0.0)
				else:
					continue

#####################################################################################################
# Function to modify (basically screening) pickled values
#####################################################################################################

#def modify_pickle(pfile):
#
#	jar=open(pfile)
#	print "Reading ", pfile
#	cookie = pickle.load(jar)
#	fsam=np.array([float("%.4f" %(x)) for x,y in cookie])
#	vel=np.array([y for x,y in cookie])
#	npicks=vel.shape[1]
#	print "Original cookie: "
#	for i in range(len(fsam)):
#		print fsam[i], vel[i]
#	for j in range(npicks):
#		rf=raw_input("Frequency range for pick %d: " %(j+1))
#		f1=float(rf.split()[0])
#		f2=float(rf.split()[1])
#		for k in range(len(fsam)):
#			if fsam[k]<f1 or fsam[k]>f2:
#				vel[k][j]=None
#	print "Modified cookie: "
#	for i in range(len(fsam)):
#		print fsam[i], vel[i]
#	jar.close()
#	new_cookie=zip(fsam,vel)
#	jarnew=open(pfile,'w')
#	pickle.dump(new_cookie,jarnew)
#	jarnew.close()
#	print "Pickle file has been modified"

####################################################################################################

def usage():
	print "What do you want to do ?"
	print "Choose from the following options"
	print "		 1 - View pickle"
	print "		 2 - Pick modes"
	print "          3 - Write pickle out in surf96 format"
	print "		 4 - Exit"

#####################################################################################################
# Main program - for using this module by itself
#####################################################################################################

if __name__=='__main__':
	nfiles=len(sys.argv)
	filnums=range(1,nfiles)
	flist=[]
	for fn in filnums:
		fname=sys.argv[fn]
		flist.append(fname)
	thdpfile=raw_input('File containing theoretical dispersion (just hit ENTER if N/A): ')
	if len(thdpfile)>0:
		mnums=raw_input("Start and end mode numbers to plot: ")
		#minc=int(raw_input("Incident mode to highlight: "))
		ml=int(mnums.split()[0])
		mh=int(mnums.split()[1])
		#try:
		reoobj = reo.read_disp([thdpfile],ml,mh)
		theor_cdisp = reoobj.modcdisp[0]
		#theor_udisp = reoobj.modudisp[0]
		#except IndexError:
		#rs96obj = rs96.read_disp([thdpfile])
		#theor_cdisp = rs96obj.disp[0]
		#rs96obj = rs96.read_disp('udisp_surf96.ray')
		#theor_udisp = rs96obj.disp[0]
		solidcurve = theor_cdisp
	task = 0
	while task != 4:
		usage()
		task = int(raw_input("Enter your choice (number) here: "))
		if task>=4:
			sys.exit('Thank you')
		elif task==1 or task==2 or task==3:
			wpobj=workon_pickle(flist,task)
