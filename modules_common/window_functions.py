#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
 
def tukeywin(window_length, alpha=0.5):
	'''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    	that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    	at \alpha = 0 it becomes a Hann window.
 
    	We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    	output
 
   	Reference
    	---------
 
	http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
 
    	'''
	# Special cases
    	if alpha <= 0:
    	    return np.ones(window_length) #rectangular window
   	elif alpha >= 1:
     	    return np.hanning(window_length)
 
	# Normal case
	x = np.linspace(0, 1, window_length)
	w = np.ones(x.shape)
 
   	# first condition 0 <= x < alpha/2
    	first_condition = x<alpha/2
   	w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
 
    	# second condition already taken care of (central region where window is 1)
 
    	# third condition 1 - alpha / 2 <= x <= 1
    	third_condition = x>=(1 - alpha/2)
    	w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))
	return w

###############################################

if __name__=='__main__':
	import sys
	param=float(sys.argv[1])
	x=range(2000)
	y=tukeywin(len(x),param)
	plt.plot(x,y)
	plt.show()
