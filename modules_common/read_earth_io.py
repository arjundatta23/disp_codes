#!/usr/bin/python

##########################
# Arjun Datta, August 2013
##########################

import os
import re
import gzip
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# Modules written by me

import read_surf96_io as rs96

##################################################################################################################

def plot_mod(plotthese):

	mplcol=['k','r','b','y']
	mnames=['Regular','Perturbed','Other']
	fig=plt.figure()
	ax=fig.add_subplot(111)
	print "Plotting %d model(s)..." %(len(plotthese))
	for pn,pt in enumerate(plotthese):
		dep=[x for x,y in pt]
		val=[y for x,y in pt]
		#ax.plot(val,dep,label=mnames[pn])
		#ax.plot(val,dep,'-',color=mplcol[pn],label=mnames[pn])
		try:
			ax.plot(val,dep,'-',color=mplcol[pn],label='Medium %d' %(pn+1))
		except IndexError:
			ax.step(val,dep,'-',label='Medium %d' %(pn+1))
	ax.set_xlim(3.5,7)
	ax.set_ylim(700,0)
	ax.set_xlabel('Vs [km/s]')
	ax.set_ylabel('Depth [km]')
	plt.legend(loc='best',prop={'size':20})
	#plt.grid(True)
	plt.show()
	#plt.savefig('reo_mod.eps',orientation='portrait')

##################################################################################################################

def plot_pdegn(mode,period,plotcontent,xname,dothis,xmin=0,xmax=0):

	""" Function to plot partial derivatives or eigenfunctions, used only when this module is run
	    as a script by itself
	"""

	fig=plt.figure()
	ax=fig.add_subplot(111)
	cnames=['Background Model','Perturbed Region']
	for pn,pt in enumerate(plotcontent):
		try:
			dep=[x for x,y in pt]
			val=[y for x,y in pt]
			#ax.plot(val,dep,label='File %d' %(pn+1))
			ax.plot(val,dep,label=cnames[pn])
		except TypeError:
			sys.exit('Requested component does not exist')
		#print "Max. depth is: ", dep[-1]
	ax.set_ylim([dep[-1],0])
	#ax.set_ylim(2000,0)
	figinfo = "Mode %d Period %.1f s" %(mode,period)
	ax.set_ylabel('Depth [km]')
	ax.grid(True,color='0.6')
	ax.legend(loc=3)
	if dothis=='show':
		ax.set_xlabel(xname)
		ax.set_title(figinfo)
		plt.show()
	elif dothis=='save':
		locinfo = " %d km" %(location)
		ax.set_xlabel(locinfo)
		ax.set_xlim([xmin,xmax])
		ax.set_title(xname)
		ax.text(0.1,0.1,figinfo,transform=ax.transAxes,fontsize=18)
		plt.savefig('localeigenfn.jpg')

###################################################################################################################

def plot_disp(modes,firstfig,secondfig=None):
	
	""" Function to plot phase or group velocity dispersion, used only when this module is run
	    as a script by itself
	"""
	def singlefig(thisax,plot_content,plot_name=None):
		mplcol=['k','r','b','y']
		if wvel=='p':
			yappropriate='Phase Velocity [km/s]'
		else:
			yappropriate='Group Velocity [km/s]'
		for pn,pt in enumerate(plot_content):
			for i,thismode in enumerate(pt):
				try:
					f = [x for x,y in thismode]
					v = [y for x,y in thismode]
				except ValueError:
					f = [x for x,y,z in thismode]
					v = [y for x,y,x in thismode]
				if len(plot_content)==1:
					curve_name="Mode %d" %modes[i]
					ax.plot(f,v,'-',label=curve_name)
				else:
					if i==0:
						ax.plot(f,v,color=mplcol[pn],label='Medium %d' %(pn+1))
					else:
						ax.plot(f,v,color=mplcol[pn])
		thisax.set_xlabel('Frequency [Hz]')	
		thisax.set_ylabel(yappropriate)
		thisax.set_ylim([3,7])
		#thisax.set_xlim([0,0.08])
		#thisax.set_title('Love')
		if not plot_name==None:
			thisax.set_title(plot_name)
		if plot_name=="Group Velocity" or plot_name==None:
			thisax.legend(loc=1)
		thisax.grid(True,color='0.6')
	if secondfig==None:
		ax = plt.subplot(111)
		singlefig(ax,firstfig)
	else:
		for f in range(1,3):
			ax = plt.subplot(1,2,f)
			if f==1:
				singlefig(ax,firstfig,'Phase Velocity')
			if f==2:
				singlefig(ax,secondfig,'Group Velocity')
		plt.suptitle("Dispersion in Model", fontsize=16)
	plt.show()

###################################################################################################################

def plot_excitation(x,y1,y2,y3):

	""" Function to plot source excitation for all modes, used only when this module is run
	    as a script by itself
	"""

	plt.plot(x,y1,'--o',label='py1')
	plt.plot(x,y2,'--o',label='py2')
	plt.plot(x,y3,'--o',label='py3')
	plt.legend(loc=1)
	plt.show()

###################################################################################################################

class read_modfile:

	def __init__(self,filelist):
		self.mnames=filelist
		self.struc=range(len(filelist))
		for fn,fl in enumerate(filelist):
			try:
				self.struc[fn]=self.read_single_file(fl)
				print "Reading file as normal: ",fl
			except ValueError:
				listforearth=[fl]
				print "Trying to read file %s using rs96" %(fl)
				rs96obj=rs96.read_models(listforearth)
				self.struc[fn]=rs96obj.vs_struc[0]
		print "reo finished reading %d file(s)" %(fn+1)

	def read_single_file(self,modfile):

		try:
			fobj = open(modfile,'r')
		except IOError:
			sys.exit('Cannot read file %s' %(modfile))
		nlayers=int((fobj.readline()).split()[0])
		dep=[]
		dep.append(0.0)
		vp=range(nlayers)
		vs=range(nlayers)
		rho=range(nlayers)
		fcontent=fobj.readlines()
		for ln,line in enumerate(fcontent):
			if ln==nlayers:
				break
			else:
				thk=float(line.split()[0])
				vp[ln]=float(line.split()[1])
				vs[ln]=float(line.split()[2])
				rho[ln]=float(line.split()[3])
				dep.append(dep[-1]+thk)
	
		result=zip(dep,vs)
		fobj.close()	
		del fobj
		return result
	
###################################################################################################################

class read_egnfile:
	
	""" Reads the eigenfunction for a specified mode and period from one or more eigenfunction files
	    Inputs to class: list of eigenfunction files, mode number, period

	    NB: If the exact specified period is not present in the file, this class will read the next higher
		period that is in the file
	"""

	def __init__(self,flist,m_concerned,p_concerned):
		self.mpp = m_concerned
		self.ppp = p_concerned
		self.uz=range(len(flist))
		self.ur=range(len(flist))
		self.tz=range(len(flist))
		self.tr=range(len(flist))
		self.ut=range(len(flist))
		self.tt=range(len(flist))
		for fn,fl in enumerate(flist):
			allcomp=self.read_single_file(fl)
			try:
				self.uz[fn]=allcomp[0]
				self.ur[fn]=allcomp[1]
				self.tz[fn]=allcomp[2]
				self.tr[fn]=allcomp[3]
				self.ut[fn]=allcomp[4]
				self.tt[fn]=allcomp[5]
			except TypeError:
				# period too long for mode in question
				pass

	def read_single_file(self,inyifile):
		moddep=self.parse_file(inyifile)
		x=self.pick_right_slice()
		return x
	def parse_file(self,egn_file):

		""" Goes through the eigen file & reads it into memory. Returns a list containing 
	   	    the entire file below the lines containing the model at the top

		    In the eigen file each mode has a header of two lines (specifying number of periods) and each
		    period has a header of 1 line
		"""

		if egn_file.endswith('.gz'):
		            egn_contents = gzip.GzipFile(egn_file,'r')
		else:
	        	    egn_contents = open(egn_file,'r')
	
		dep=[]
		ncol_ph=7
		# First line in the file is a string
		egn_contents.readline()
		# Second line specifies number of layers in model
		model_lines = int(egn_contents.readline().split()[0])

		# Read the file into memory but get rid of the initial model-specifying part
		egn_entire=egn_contents.readlines()
		for i,line in enumerate(egn_entire):
			if i < model_lines: 
				values_line = line.split()
				dep.append(float(values_line[0]))
			else: break
		#del egn_entire[:model_lines]
		self.whole_egn = egn_entire[model_lines+2:]
		# Get the period samples for each mode
		self.nps=[]
		self.ps=[]
		ps_mode=[]
		for line in self.whole_egn:
			if "samples" in line.split():
				self.nps.append(int(line.split()[1]))
				if __name__=='__main__':
					print "no. of periods is ", self.nps[-1]
			if (len(line.split())==ncol_ph) and (not line.split()[2].isalpha()):
				per=float(line.split()[1])
				per_lines=int(line.split()[-2])
				ps_mode.append([per,per_lines])
			#	print "length of ps_mode is ", len(ps_mode)
			if len(ps_mode)==self.nps[-1]:
				if __name__=='__main__':
					print "Found mode %d" %(len(self.ps))
				self.ps.append(ps_mode)
				ps_mode=[]
		#print self.nps, self.ps
		self.modes_lines=[]
		for ps_m in self.ps:
			m_lines = sum([y for x,y in ps_m]) + (len(ps_m)) + 2
			self.modes_lines.append(m_lines)
		egn_contents.close()
		del egn_contents
		return np.array(dep)

	def pick_right_slice(self):

		result=range(6)
		islov = False
		isray = False		
		modes = range(len(self.nps))
		prev_mlines = sum([ self.modes_lines[x] for x,y in enumerate(modes) if y<self.mpp ])
		if __name__=='__main__':
			print "mpp & ps are ", self.mpp
		prev_ps = [y for x,y in self.ps[self.mpp] if x<self.ppp]
		prev_plines = sum(prev_ps) + len(prev_ps)
		prev_lines = int(prev_mlines + prev_plines)
		try:
			self.ppp = float(self.whole_egn[prev_lines+1].split()[1])
			print "Extracting eigenfunctions for mode %d, period %f" %(self.mpp,self.ppp)
		except (ValueError, IndexError):
			if __name__ == '__main__':
				sys.exit("Mode number %d does not exist at this period !!" %(self.mpp))
			else:
				return 1
		rel_lines = self.ps[self.mpp][len(prev_ps)][1]
		rel_slice = self.whole_egn[(prev_lines+1):(prev_lines+1+rel_lines+1)]
		nfac=1 #float(rel_slice[0].split()[-4])
		print "From reo: nfac is ", nfac
		d=[]
		y1=[]
		y2=[]
		y3=[]
		y4=[]
		for line in rel_slice[1:]:
			values_line = [ float(i) for i in line.split() ]
			d.append(values_line[0])
			y1.append(nfac*values_line[1])
			y2.append(nfac*values_line[2])
			try:
				y3.append(nfac*values_line[3])
				y4.append(nfac*values_line[4])
				isray = True
			except IndexError:
				islov = True
		if isray:
			result[0]=zip(d,y1)
			result[1]=zip(d,y2)
			result[2]=zip(d,y3)
			result[3]=zip(d,y4)
			result[4]=None
			result[5]=None
		elif islov:
			result[4]=zip(d,y1)
			result[5]=zip(d,y2)
			result[0]=None
			result[1]=None
			result[2]=None
			result[3]=None
		return result
		
####################################################################################################################

class read_egnfile_per:

	""" Class to read a single eigenfunction file and extract the eigenfunction FOR ALL MODES AT A SPECIFIED
	    PERIOD

	NB: The difference between this class and the read_egnfile class in terms of operation is that this
	    class will only work if the EXACT period (exact upto a certain precision) specified as input is 
	    present in the eigenfunction file. This class is NOT capable of improvising to read the closest period
	    to the one requested """
	
	def __init__(self,infile,p_concerned):
		if infile.endswith('.gz'):
		            egn_contents = gzip.GzipFile(infile,'r')
		else:
	        	    egn_contents = open(infile,'r')
		ncol_ph=7
		# First line in the file is a string
		egn_contents.readline()
		# Second line contains number of layers
		model_deps = int(egn_contents.readline().split()[0])
		self.dep=np.arange(model_deps,dtype=float)
		self.mu=np.arange(model_deps,dtype=float)
		self.lamda=np.arange(model_deps,dtype=float)
		self.rho=np.arange(model_deps,dtype=float)
		alpha=range(model_deps)
		beta=range(model_deps)
		#self.rho=range(model_deps)
		i=0
		j=0
		k=0
		startreading=False
		for line in egn_contents:
			if i < model_deps: 
				values_line = line.split()
				self.dep[i]=float(values_line[0])
				alpha[i]=float(values_line[3])
				beta[i]=float(values_line[1])
				self.rho[i]=float(values_line[2])
				self.mu[i]=self.rho[i]*(beta[i]**2)
				self.lamda[i]=self.rho[i]*(alpha[i]**2)-(2*self.mu[i])
				i+=1
			if "modes listed" in line:
				self.totm=int(line.split()[0])
				self.wavnum=np.zeros(self.totm)
				self.uzmat=np.zeros((model_deps,self.totm))
				self.urmat=np.zeros((model_deps,self.totm))
				self.tzmat=np.zeros((model_deps,self.totm))
				self.trmat=np.zeros((model_deps,self.totm))
				self.utmat=np.zeros((model_deps,self.totm))
				self.ttmat=np.zeros((model_deps,self.totm))
				
			if startreading:
				values_line = [ float(i) for i in line.split() ]
				y1[j]=values_line[1]
				y2[j]=values_line[2]
				try:
					y3[j]=values_line[3]
					y4[j]=values_line[4]
					isray=True
					islov=False
				except IndexError:
					islov=True
					isray=False
				j+=1
				if j==lyrsthism:
					startreading=False
					if isray:
						self.uzmat[:len(y1),k]=y1
						self.urmat[:len(y2),k]=y2
						self.tzmat[:len(y3),k]=y3
						self.trmat[:len(y4),k]=y4
						self.utmat=None
						self.ttmat=None
					elif islov:
						self.utmat[:len(y1),k]=y1
						self.ttmat[:len(y2),k]=y2
						self.uzmat=None
						self.urmat=None
						self.tzmat=None
						self.trmat=None
					k+=1
			#if "mode=" in line.split()[0]:
			if (len(line.split())==ncol_ph) and (not line.split()[0].isalpha()):
				perto5th=round(float(line.split()[1]),5)
				if perto5th==round(p_concerned,5):
					if __name__=='__main__':
						print "Extracting eigenfunction for mode %d and period %.6f" %(int(line.split()[0]),float(line.split()[1]))
					j=0
					self.wavnum[k]=2*np.pi/(p_concerned*float(line.split()[2]))
					lyrsthism=int(line.split()[5])
					startreading=True
					y1=np.arange(lyrsthism,dtype=float)
					y2=np.arange(lyrsthism,dtype=float)
					y3=np.arange(lyrsthism,dtype=float)
					y4=np.arange(lyrsthism,dtype=float)

		# At this stage, i.e. after the entire file has been read, the value of k will be equal to the number
		# of modes that exist (are present in the file) for the given period. The matrices containing the 
		# eigenfunctions (each column of matrix = separate mode) should therefore be trimmed down to k
		# columns. This is done by the function final_result
		
		self.final_result(k,isray)
		
	def final_result(self,relm,ray):
		
		extracols=np.arange(relm,self.totm)
		self.wavnum=np.delete(self.wavnum,extracols)
		if ray:
			self.uzmat=np.delete(self.uzmat,extracols,1)
			self.urmat=np.delete(self.urmat,extracols,1)
			self.tzmat=np.delete(self.tzmat,extracols,1)
			self.trmat=np.delete(self.trmat,extracols,1)
		else:
			self.utmat=np.delete(self.utmat,extracols,1)
			self.ttmat=np.delete(self.ttmat,extracols,1)
	
####################################################################################################################

class read_disp:

	""" Class to read the dispersion file (Love or Rayleigh) produced by earthsr, and extract either the:
		
		1. The partial derivatives of phase velocity w.r.t. model parameters for a specified mode & period
			OR		
		2. Phase and Group velocity dispersion curves for all modes in between 2 specified mode numbers
			OR
		3. The flattened model which is only present in the disp file
	
	    Which of these 3 operations is performed depends on the input to the class
	    If the last argument (2nd integer) is a period (in seconds), the class will perform operation 1,
            if it is a higher mode number, the class will perform operation 2.
	    If the last 2 arguments are both None (i.e. user runs script with disp file(s) as only argument), the class will perform operation 2

	    Input is: list of dispersion files, first integer, second integer (the 2 'integers' may be = None)
		      When they are indeed integers, the first one will always be interpreted as a mode number,
			2nd integer will be interpreted as higher mode number if it is <=30 or as period if it is > 30 """
	
	def __init__(self,flist,int1,int2):

		self.getpd=False
		self.getdisp=False
		self.getfmod=False

		if int1==None and int2==None:
			# operation 3 - read flattened model
			self.getfmod=True
			self.fstruc=range(len(flist))
			for fn,fl in enumerate(flist):
                                self.fstruc[fn]=self.read_single_file(fl)
                                print "Finished reading file: ", fl
		elif int2>30:
			# operation 1 - reading partial derivatives
			self.getpd=True
			self.mnum=int1
			self.psec=int2
			self.modpd=range(len(flist))
                        for fn,fl in enumerate(flist):
                                self.modpd[fn]=self.read_single_file(fl)
                                print "Finished reading file: ", fl
		else:
			# operation 2 - reading dispersion
			self.getdisp=True
			self.mode_l = int1
			self.mode_h = int2
			self.modcdisp=range(len(flist))
			self.modudisp=range(len(flist))
			for fn,fl in enumerate(flist):
				try:
					[self.modcdisp[fn],self.modudisp[fn]]=self.read_single_file(fl)
				except TypeError:
					listforrs96=[fl]
					rs96obj=rs96.read_disp(listforrs96)
					self.modcdisp[fn]=rs96obj.disp[0]
					self.modudisp[fn]=self.modcdisp[fn] # because with surf96, a particular file will contain either phase or group velocity, not both
					self.rel_modes=range(len(self.modcdisp[fn]))
	def read_single_file(self,indfile):
		
		self.ncol_ph=8
		if self.getpd:
			self.parse_filepd(indfile)
			return zip(self.deporig,self.reqdpd)
		elif self.getdisp:
			self.tlm=[]	# tlm stands for Total_Lines_Mode - total lines in file before the end of any mode
			self.parse_filedisp(indfile)
		        #print "tlm is ", self.tlm
			if len(self.tlm)>0:
				[cdisp,udisp]=self.pick_right_slice()
				return cdisp,udisp
			else:
				return None
			#	if __name__=='__main__':
			#		self.extra_analysis(cdisp)
		elif self.getfmod:
			self.parse_filepd(indfile,True)
			return zip(self.depf,self.vsf)

	def parse_filepd(self,disp_file,modonly=False):

		""" Reads the file to extract partial derivatives for the required mode & period.
		    Stops reading once it has found the right mode-period combination and read in
		    all its partial derivatives
		"""

		if disp_file.endswith('.gz'):
		            disp_contents = gzip.GzipFile(disp_file,'r')
		else:
	        	    disp_contents = open(disp_file,'r')
		
		pd_one_par=[]
         	pd_par=[]
		npar=int(disp_contents.readline().split()[0])
		nlo=int(disp_contents.readline().split()[0])
		# Second line contains num. of layers in original model
		
		# Number of columns in file for the period headers, and number of model parameters
		# w.r.t which partial derivatives are present in the file, need to be known before-hand
		
		self.deporig=[0.0]
		self.depf=[0.0] 
		self.vsf=[]
		cl=0
		pardone=0
		mainreading=False
		foundtarget=False
		for line in disp_contents:
			if not mainreading:
				if cl<nlo:
					thk=float(line.split()[0])
					self.deporig.append(self.deporig[-1]+thk)
					cl+=1		
				else:
					try:
						thk_fltnd=float(line.split()[0])
						vs_fltnd=float(line.split()[-1])
						self.vsf.append(vs_fltnd)
						try:
							#self.depf[cl-nlo]=self.depf[cl-nlo-1]+thk_fltnd
							self.depf.append(self.depf[-1]+thk_fltnd)
							cl+=1
						except IndexError:
							pass
					except (IndexError, ValueError):
						pass
					if "modes listed" in line:
						mainreading=True
						if modonly:
							return
							print "returning now!"
			else:
				if foundtarget:
					if pardone==npar:
						break
					else:
						if (len(line.split())==self.ncol_ph and line.split()[0].isdigit()):
							sys.exit('Partial derivatives for period %d not listed in file' %(self.psec))
						values_on_line = [ float(i) for i in line.split() ]
       	               				pd_one_par.extend(values_on_line)
				                if len(pd_one_par)==depsthismp:
       	                       				pd_par.append(np.array(pd_one_par))
							pd_one_par=[]
       	                        			pardone+=1
                               	else:
					if (len(line.split())==self.ncol_ph): #and ("mode" not in line.split()):
						# on a period header line
						try:
							perto5th=round(float(line.split()[1]),5)
						except ValueError:
							print "problem with line: ", line
							sys.exit()
						if (int(line.split()[0])==self.mnum) and (perto5th==self.psec):
							foundtarget=True
							depsthismp=int(line.split()[-2])
							pardone=0
						else:
							continue
					else:
						continue
		# out of for loop
		if len(pd_par)==0:
			sys.exit('Period %d not listed in file' %(self.psec))
		print "Extracted partial derviatives for %d parameters" %(len(pd_par))
		pd_all=np.zeros((nlo,len(pd_par)))
		for col,par in enumerate(pd_par):
			pd_all[:len(par),col]=par
		#print len(self.deporig), len(self.depf), self.vsf
		#print pd_par, pd_all
		self.reqdpd=pd_all[:,2]
		disp_contents.close()
		del disp_contents

	def parse_filedisp(self,disp_file):

		""" Parses the file by first reading all of it into memory 
		"""

		if disp_file.endswith('.gz'):
		            disp_contents = gzip.GzipFile(disp_file,'r')
		else:
	        	    disp_contents = open(disp_file,'r')
		
		# Read the file into memory but get rid of the initial model-specifying part
		disp_entire=disp_contents.readlines()
		if "SURF96" in disp_entire[0]:
			# this file is produced by surf96
			return
		for i,line in enumerate(disp_entire):
			if "modes listed" in line:
				total_modes = int(line.split()[0])
				if total_modes < (self.mode_h+1):
					self.mode_h = total_modes-1
				break
		self.whole_disp = disp_entire[i+1:]	
		
		# Get the period samples and number of lines for each mode
		per_linescounter=0
		#ps_mode=[]
		for line in self.whole_disp:
			per_linescounter+=1
			if "****" in line:
				self.tlm.append(per_linescounter)
				per_counter=0	
		disp_contents.close()
		del disp_contents

	def pick_right_slice(self):
		
		all_modes = range(len(self.tlm))
		try:
			prev_lines = [ self.tlm[x] for x,y in enumerate(all_modes) if y<self.mode_l ][-1]
		except IndexError:
			prev_lines = 0
		last_rel_line = [ self.tlm[x] for x,y in enumerate(all_modes) if y>=self.mode_l and y <= self.mode_h ][-1]
		self.rel_modes = range(self.mode_l,self.mode_h+1)
		rel_slice = self.whole_disp[prev_lines:last_rel_line]
		#print self.tlm
		print "Reading between lines %d and %d " %(prev_lines, last_rel_line)
		pvd=range(len(self.rel_modes)) # pvd stands for phase velocity dispersion
		gvd=range(len(self.rel_modes)) # gvd stands for group velocity dispersion
		pv=[]
		gv=[]
		freq=[]
		mdone=0
		for line in rel_slice:
			if "****" in line:
				try:
					pvd[mdone]=zip(freq,pv)
				except IndexError:
					break
				gvd[mdone]=zip(freq,gv)
				print "Finished reading dispersion of mode ", self.rel_modes[mdone]
				mdone+=1
				pv=[]
				gv=[]
				freq=[]
			else:
				if (len(line.split())==self.ncol_ph) and (line.split()[0].isdigit()) and (line.split()[-1].isdigit()):
					values_line = [ i for i in line.split() ]
					#print values_line
					freq_temp = '%.4f' %(1./float(values_line[1]))
					freq.append(float(freq_temp))
					pv.append(float(values_line[2]))
					gv.append(float(values_line[3]))
		return pvd,gvd
			
####################################################################################################################

class read_excitation:
	
	""" Class to read the Rayleigh wave excitation file produced by earthsr, and
	    extract excitation amplitudes of all modes (all those considered by earth) at a particular frequency
	    Input to class: excitation file, period """
	
	def __init__(self,exfile,freq_interest):
		self.omega = 2*np.pi*freq_interest
		self.rel_modes = []
		self.a=[]
		self.b=[]
		self.c=[]
		if exfile.endswith('.gz'):
		            ex_contents = gzip.GzipFile(exfile,'r')
		else:
	        	    ex_contents = open(exfile,'r')
		mode_exists=0
		foundline=False
		for line in ex_contents:
			if foundline:
				print line
				self.a.append(float(line.split()[0]))
				self.b.append(float(line.split()[1]))
				self.c.append(float(line.split()[2]))
				foundline=False
			col1 = float(line.split()[0])
			if len(line.split())==1 and not float(line.split()[0])==-1:
				self.rel_modes.append(int(float(line.split()[0])))
			if abs(col1-self.omega)<0.000001:
				foundline=True
				mode_exists+=1
		self.rel_modes=self.rel_modes[1:mode_exists+1]
		
####################################################################################################################

""" If running this module as a script by itself, it produces a plot of the eigenfunctions of a particular mode 
    at a particular period or a plot of phase or group velocity dispersion for a set of specified modes, or of the
    excitation amplitudes of all modes at a particular frequency

    This module when run by itself also outputs some auxillary information that may be useful for certain applications
"""

if __name__ == '__main__':
	nfiles=len(sys.argv)
	filenums=range(1,nfiles)
	xnames=[]
	for fn in filenums:
		xnames.append(sys.argv[fn])
	egfile=xnames[0]
	if os.path.basename(egfile).startswith('eigen') or os.path.basename(egfile).startswith('yi.'):
		egn = {'Ux': 'ur', 'Uz': 'uz', 'Uy': 'ut', 'szz': 'tz', 'sxz': 'tr', 'syz': 'tt'}
		egn_name = {'Ux': 'Radial displacement', 'Uz': 'Vertical displacement', 'Uy': 'Transverse displacement', 'szz': 'Normal stress', 'sxz': 'Shear stress radial', 'syz': 'Shear stress transverse'}
		try:
			location=int([ k for k in egfile.split('.') if k.isdigit() ][0])
			m = int(sys.argv[-6])
			p = float(sys.argv[-5])
			component = sys.argv[-4]
			dowhat=sys.argv[-3]
			minx=float(sys.argv[-2])
			maxx=float(sys.argv[-1])
			fklist=[egfile]
			refobj=read_egnfile(fklist,m,p)
		except (IndexError,ValueError):
			dowhat="show"
			try:
				m = int(sys.argv[-3])
				p = float(sys.argv[-2])
				component = sys.argv[-1]
				xnames=xnames[:-3]
				if not (component in egn):
					print "Eigenfunction component must be one of the following:"
					for d in egn:
						print d
					print "Please try again"
					sys.exit()
				refobj=read_egnfile(xnames,m,p)
			except ValueError:
				sys.exit("You must specify mode,period & which component of the eigenfunction you want to look at !! Please try again")
		plot_quantity = getattr(refobj,egn[component])
		print "Plotting the: ", egn_name[component]
		#print "Number of points in eigenfunction array: ", plot_quantity.size
		#rel_depth = refobj.dep[:plot_quantity.size]
		if dowhat=="save":
			plot_pdegn(refobj.mpp,refobj.ppp,plot_quantity,egn_name[component],dowhat,minx,maxx)
		else:
			plot_pdegn(refobj.mpp,refobj.ppp,plot_quantity,egn_name[component],dowhat)
	elif os.path.basename(egfile).startswith('disp') or os.path.basename(egfile).startswith('dp.'):
		try:
			m1 = int(sys.argv[-2])
			m2 = int(sys.argv[-1])
			xnames=xnames[:-2]
		except ValueError:
			# no mode numbers specified means we read the flattened model from the disp file
			m1=None
			m2=None
		rdobj = read_disp(xnames,m1,m2)
		if hasattr(rdobj,'modpd'):
			plot_pdegn(rdobj.mnum,rdobj.psec,rdobj.modpd,'dc/dVs','show')
		elif hasattr(rdobj,'modcdisp'):
			wvel = raw_input("Do you want to see dispersion of group (g) or phase (p) velocity or both (b) ? ")
			if wvel=='p':
				plot_disp(rdobj.rel_modes,rdobj.modcdisp)
			elif wvel=='g':
				plot_disp(rdobj.rel_modes,rdobj.modudisp)
			elif wvel=='b':
				plot_disp(rdobj.rel_modes,rdobj.modcdisp,rdobj.modudisp)
		else:
			print "plot flattened model"
			plot_mod(rdobj.fstruc)
		#	for mn,md in enumerate(rdobj.pvd):
		#		v = [y for x,y in md if x==1./p ]
		#		print rdobj.rel_modes[mn],p, v
	elif os.path.basename(egfile).startswith('excitation'):
		try:
			see_freq = float(sys.argv[-1])
		except ValueError:
			print "You must specify a frequency at which excitation amplitudes are to be extracted !!"
			print "Please try again"
			sys.exit()
		source_ex = read_excitation(infile,see_freq)
		plot_excitation(source_ex.rel_modes,source_ex.a,source_ex.b,source_ex.c)
	elif os.path.basename(egfile).startswith('mod'):
		rmfobj=read_modfile(xnames)
		plot_mod(rmfobj.struc)
	else:
		print "This script needs a disp, eigen or excitation file produced by earthsr"
		print "Please supply one and try again"
