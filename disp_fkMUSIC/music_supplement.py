#!/usr/bin/python

import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Modules written by me
import read_earth_io as reo
import read_surf96_io as rs96

###################################################################################################

def show_single(inmat,withpicks):

	spec=plt.figure()
	axspec=spec.add_subplot(111)
	X,Y = np.meshgrid(freqs,pvels)
	cax=axspec.pcolor(X,Y,inmat)
	if withpicks:
		axspec.plot(fwinpick,pv,'ko')
		axspec.plot(fwinpick,pvlb,color='k',lw=1.5,) #,'*')
		axspec.plot(fwinpick,pvub,color='k',lw=1.5) #,'*')
	axspec.set_xlim(f1,f2)
	axspec.set_ylim(3,7)
	axspec.set_xlabel('Frequency [Hz]')
	axspec.set_ylabel('Phase Velocity [km/s]')
	#axspec.set_title('Stack of %d events' %(len(jarlist)))
	#spec.colorbar(cax)
	plt.show()

##############################################################################################################

def store_result():

	fname=flist[0]
	jarname=fname.replace("music","mpic")
	jarfile=os.path.join(os.getcwd(),jarname)
	jar=open(jarfile,'w')
	pickle.dump(ypobj.finalans,jar)
	jar.close()

#####################################################################################################

class ypick():

	def __init__(self,inspec,xvalues,yvalues,threshold):

		# inputs are f-c spectrum, freqs, pvels, th<1

		xdom=raw_input("frequency bounds (or just lower bound) for picking: ")
		ydom=raw_input("phase velocity bounds (or just lower bound) for picking: ")
		try:
			x1=float(xdom.split()[0])
			x2=float(xdom.split()[1])
			relcols=np.intersect1d(np.where(xvalues>=x1)[0],np.where(xvalues<=x2)[0])
		except IndexError:
			x2=None
			relcols=np.where(xvalues>x1)[0]
		try:
			y1=float(ydom.split()[0])
			y2=float(ydom.split()[1])
			relrows=np.intersect1d(np.where(yvalues>y1)[0],np.where(yvalues<=y2)[0])
		except IndexError:
			y2=None
			relrows=np.where(yvalues>=y1)[0]
		relspec=inspec[relrows[0]:relrows[-1]+1,relcols[0]:relcols[-1]+1]
		crel=yvalues[relrows[0]:relrows[-1]+1]
		self.frel=xvalues[relcols[0]:relcols[-1]+1]
		#print "Shape of selected portion is ", relspec.shape
		#print "Relevant c-values: ", crel
		#print "Relevant frequencies: ", self.frel
		indmax_relspec=np.argmax(relspec,axis=0)
		indmax_orig=[ ir+relrows[0] for ir in indmax_relspec]
		self.ypicks=[float("%.5f" %x) for x in crel[indmax_relspec]]
		self.picklb=[]
		self.pickub=[]
		for pn,ptf in enumerate(self.ypicks):
			col_orig=pn+relcols[0]						# pn is for pick number
			mvtf=inspec[indmax_orig[pn],col_orig]				# mvtf is for maximum_value_this_frequency
			ostf=inspec[:,col_orig]						# ostf is for original_spectrum_this_frequency
			bw=np.where(ostf>=threshold*mvtf)[0] 				# bw is for bandwidth
			belowmax=ostf[bw[0]:indmax_orig[pn]][::-1]
			abovemax=ostf[indmax_orig[pn]:bw[-1]]
			print xvalues[col_orig], mvtf, ptf, yvalues[bw[0]], yvalues[bw[-1]]
			lb=yvalues[bw[0]]
			ub=yvalues[bw[-1]]
			if not all(belowmax[i]>=belowmax[i+1] for i in range(len(belowmax)-1)):
				print "Intelligent lower bound picking at frequency ", xvalues[col_orig]
				hwbelow=bw[bw<indmax_orig[pn]][::-1] # hw is for half-width
				cbelow=yvalues[bw[0]:indmax_orig[pn]][::-1]
				""" First, check if there's any gap at all in hwbelow """
				mi=[j for j in range(len(hwbelow)-1) if hwbelow[j+1]<(hwbelow[j]-1)]
				if len(mi)>0:
					cutoff=mi[0]
					lb=cbelow[cutoff]
					belowmax=belowmax[:cutoff]
				else:
					cutoff=None
				sh=[ j for j in range(len(belowmax)-1) if belowmax[j+1]>belowmax[j] ]
				if len(sh)>0:
					lb=cbelow[sh[0]]
				#if float("%.4f" %(xvalues[col_orig]))==0.0405:
				#	print "cutoff is ", cutoff
				#	print lb
				#	print belowmax
				#	print "sh is ", sh
			if not all(abovemax[i]>=abovemax[i+1] for i in range(len(abovemax)-1)):
				print "Intelligent upper bound picking at frequency ", xvalues[col_orig]
				hwabove=bw[bw>=indmax_orig[pn]] # hw is for half-width
				cabove=yvalues[indmax_orig[pn]:bw[-1]]
				""" First, check if there's any gap at all in hwabove """
				mi=[j for j in range(len(hwabove)-1) if hwabove[j+1]>(hwabove[j]+1)]
				if len(mi)>0:
					cutoff=mi[0]
					ub=cabove[cutoff]
					abovemax=abovemax[:cutoff]
				sh=[ j for j in range(len(abovemax)-1) if abovemax[j+1]>abovemax[j] ]
				if len(sh)>0:
					ub=cabove[sh[0]]
			print "Bounds are ", lb, ub
			self.picklb.append(lb)
			self.pickub.append(ub)
		self.calculate_error()
	
	def calculate_error(self):
		
		print "Getting error estimate... "
		#errtop=[i-j for i,j in zip(self.pickub,self.ypicks)]
		#errbot=[i-j for i,j in zip(self.ypicks,self.picklb)]
		#errav=[(i+j)/2 for i,j in zip(errbot,errtop)]
		errpick=range(len(self.frel))
		ystep=0.01
		for f,fr in enumerate(self.frel):
			ydist=np.arange(self.picklb[f],self.pickub[f],ystep)
			sd=np.std(ydist)
			errpick[f]=sd
		self.finalans=zip(self.frel,self.ypicks,errpick)
		#print self.finalans

##########################################################################################################

class view_pickle():

	def __init__(self,pfilelist):
	
		self.totp=len(pfilelist)
		self.plist=pfilelist

		self.freqs=range(self.totp)
		self.pvels=range(self.totp)
		self.resmat=range(self.totp)

                self.make_plot()

	def make_plot(self):

		usrc='n'
		usrc=raw_input("See theoretical dispersion too ? (y/n): ")
		if usrc=='y':
			thdpfile=[raw_input('File containing theoretical dispersion: ')]
			#try:
			#reoobj = reo.read_disp(thdpfile,0,5)
			#theor_cdisp = reoobj.modcdisp[0]
			#theor_udisp = reoobj.modudisp[0]
			rs96obj = rs96.read_disp(thdpfile)
			theor_cdisp = rs96obj.disp[0]
			
			solidcurve = theor_cdisp
		else:
			theor_cdisp=0
			theor_udisp=0

		if self.totp==1:
			pncols=1
			r=1
			spec=plt.figure()
		else:
			pncols=2
			r=self.totp/pncols if self.totp%pncols==0 else (self.totp/pncols)+1
			if r==1:
				sizy=4.75
			elif r==2:
				sizy=9
			elif r==3:
				sizy=12
			spec=plt.figure(figsize=(12,sizy))

		usr_t=raw_input("Enter title of plot: ")
		spec.suptitle('Event '+usr_t)

		for i,pfile in enumerate(self.plist):
			jar=open(pfile)
			cookie1=pickle.load(jar)
			f1=cookie1[0] #0.0065
			f2=cookie1[1] #0.0615
			cookie2=pickle.load(jar)
			cmin=cookie2[0] #3
			cmax=cookie2[1] #8
			cookie3 = pickle.load(jar)
			nrows=cookie3.shape[0]
			ncols=cookie3.shape[1]
			print "No. of rows and columns: ", nrows, ncols
			self.freqs[i]=np.linspace(f1,f2,ncols)
			self.pvels[i]=np.linspace(cmin,cmax,nrows)
			print "Length of pvels and freqs: ", len(self.pvels[i]), len(self.freqs[i])
			self.resmat[i]=np.zeros((nrows,ncols))
			try:
				self.resmat[i]=cookie3
			except IndexError:
				sys.exit('Please check the dimensions of the stored matrix and try again')
			jar.close()
		
			axspec=spec.add_subplot(r,pncols,i+1)
			#spec=plt.figure()
			#axspec=spec.add_subplot(111)
			X,Y = np.meshgrid(self.freqs[i],self.pvels[i])
			#print "Shape of matrix is: ", self.resmat[i].shape
			cax=axspec.pcolor(X,Y,self.resmat[i])
			axspec.set_xlim(self.freqs[i][0],self.freqs[i][-1])
			axspec.set_ylim(3,7)
			axspec.set_xlabel('Frequency [Hz]')
			axspec.set_ylabel('Phase Velocity [km/s]')
			#axspec.set_title('self.resmat of %d events' %(len(jarlist)))
			#spec.colorbar(cax)

			if usrc=='y' and i==1:
				for k in range(len(solidcurve)):
					print "mode number: ", k
					try:
						f = [x for x,y,z in solidcurve[k]]
						v = [y for x,y,z in solidcurve[k]]
					except ValueError:
						f = [x for x,y in solidcurve[k]]
						v = [y for x,y in solidcurve[k]]
					curve_name="Mode %d" %k
					axspec.plot(f,v,'-',linewidth=2) #,label=curve_name)
		
		#usrc_extra=raw_input("Plot ucd pickle too ? (y/n): ")
		#usrc_extra='n'
		#if usrc_extra=='y':
		#	ucdpfname=raw_input('UCD pickle file: ')
		#	ucdpfile=os.path.normpath(os.path.join(script_dir,ucdpfname))
			#ucdpfname=[n for n in os.listdir(ucdpdir) if self.pname+'_c' in n][0]
		#	try:
		#		jar=open(ucdpfile)
		#		print "Reading ", ucdpfile
		#		cookie3 = pickle.load(jar)
		#		fsam=np.array([float("%.4f" %(x)) for x,y in cookie3])
		#		vel=np.array([y for x,y in cookie3])
		#		npicks=vel.shape[1]
		#		for pick in range(npicks):
		#			cpick=vel[:,pick]
		#			axspec.plot(fsam,cpick,'wo')
		#	except ValueError:
		#		print "Unable to read ucd pickle file"
		#		pass
		plt.show()
		#plt.savefig('music_pickle_plot.png')

#####################################################################################################
# Main program - for using this module by itself
#####################################################################################################

if __name__=='__main__':

	script_dir=os.getcwd()
	nfiles=len(sys.argv)
	filnums=range(1,nfiles)
	flist=[]
	for fn in filnums:
		fname=sys.argv[fn]
		flist.append(fname)
	#picklefile=sys.argv[1]		obsolete; for use when there was necessarily only one input
	#vpobj=view_pickle(picklefile)
	if len(flist)==1:
		usrc=raw_input("Pick curve (with errors) on spectrum ? (y/n): ")
		if usrc=='y':
			jar=open(flist[0])
			cookie1=pickle.load(jar)
			f1=cookie1[0]
			f2=cookie1[1]
			cookie2=pickle.load(jar)
			cmin=cookie2[0]
			cmax=cookie2[1]
			cookie3 = pickle.load(jar)
			nrows=cookie3.shape[0]
			ncols=cookie3.shape[1]
			print "No. of rows and columns: ", nrows, ncols
			freqs=np.linspace(f1,f2,ncols)
			pvels=np.linspace(cmin,cmax,nrows)
			print "Length of pvels and freqs: ", len(pvels), len(freqs)
			resmat=np.zeros((nrows,ncols))
			try:
				resmat=cookie3
			except IndexError:
				sys.exit('Please check the dimensions of the stored matrix and try again')
			jar.close()
			show_single(resmat,False)
			ypobj=ypick(resmat,freqs,pvels,0.95)
			fwinpick=ypobj.frel
			pv=ypobj.ypicks
			pvub=ypobj.pickub
			pvlb=ypobj.picklb
			show_single(resmat,True)
			store_result()
		else:
			vpobj=view_pickle(flist)
	else:
		vpobj=view_pickle(flist)
