#!/usr/bin/python

""" Module that takes as input an array of x-values, an array of y-values, and an array of 'new' x-values; 
    produces an output of y-values interpolated to the new x-values """

import sys
import numpy as np

def lagrange_linear(xdata,ydata,xnew,precision=6):

	""" The only condition on the inputs to this function is that the x-coordinates of the data
	    be sorted in ascending order. The 4th argument is optional and it defines the precision
	    at which two float values are deemed to be equal """
	if not len(xdata)==len(ydata):
		print "x-input has length %d, y-input has length !" %(len(xdata),len(ydata))
		return None
	else:
		ynew=np.zeros(xnew.shape)
		pos_xnew=np.searchsorted(xdata,xnew)
		repstring="%."+str(precision)+"f"
		for j,ind in enumerate(pos_xnew):
			try:
				xd_actual=float(repstring %(xdata[ind]))
				xnew_actual=float(repstring %(xnew[j]))
				#print "xnew is ", xnew_actual
			except IndexError:
				#print pos_xnew
				sys.exit('Either the input x data is not sorted in ascending order or you are trying to extrapolate rather than interpolate')
			if xd_actual==xnew_actual:
				#print "this x-point is common: ", xnew[j]
				ynew[j]=ydata[ind]
			else:
				#print "this x-point is new: ", xnew[j]
				x_upper=xdata[ind]
				x_lower=xdata[ind-1]
				y_lower=ydata[ind-1]
				try:
					y_upper=ydata[ind]
				except IndexError:
					print "Trying to access element number %d out of %d elements in the data" %(ind,len(ydata))
					#print pos_xnew
				ynew[j]=y_lower*((x_upper-xnew[j])/(x_upper-x_lower)) + y_upper*((xnew[j]-x_lower)/(x_upper-x_lower))
				#print "Bounds are:", x_lower, x_upper
		return ynew
	


