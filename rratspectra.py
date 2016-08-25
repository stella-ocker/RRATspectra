#!/usr/bin/env python

'''
rratspectra.py

Plots RRAT spectra and calculates their slope for single or multiple .spd files produced by make_spd.py. Optionally stacks spectra from multiple .spd files to increase the S/N. Default set to analyze data at zero DM. Background spectra plotted in grey.

Usage on the command line:
>>python rratspectra.py -infile INFILE(.spd file(s)) [OPTIONS]

Stella Koch Ocker (socker@oberlin.edu) - Aug 22, 2016 

'''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from scipy import stats
import math
import argparse
from ConfigParser import ConfigParser
import sys

#defining command line arguments
parser = argparse.ArgumentParser('python rratspectra.py')
parser.add_argument('-infile',help='List of .spd or .npz file(s) to analyze',nargs='*',action='store',dest='infile',required=True)
parser.add_argument('-stack',help='If flag is set and multiple spd files are given, spectra will be stacked. Default is False.',action='store_true',dest='stack')
parser.add_argument('-m_start',help='Number of subbands to ignore at the beginning of the data array when calculating spectral slope. Default: 1.',action='store',dest='m_start')
parser.add_argument('-m_end',help='Number of subbands to ignore at the end of the data array when calculating spectral slope. Default: 1.',action='store',dest='m_end')
parser.set_defaults(stack=False,m_start=1,m_end=1)
args = parser.parse_args()
if len(args.infile)==1:
	spd_list = args.infile
if len(args.infile)>1:
	spd_list = args.infile
stack = args.stack
m_start = int(args.m_start) * -1
m_end = int(args.m_end)

def main(spd_list):
	
	#specifying figure settings depending on amount of input files and value of -stack
	if len(spd_list) > 1:
		fig, axs = plt.subplots(len(spd_list),1,sharex=True,sharey=True)
		fig.subplots_adjust(hspace=0,wspace=0)
		axs=axs.ravel()
	
	if len(spd_list) == 1:
		fig = plt.figure()
		ax = fig.add_subplot(111)
	
	if stack == True:
		fig = plt.figure()
		ax = fig.add_subplot(111)

	if stack == True:
		npz_fn = spd_list[0]
		spectrum = np.load(npz_fn)
		dat = spectrum['Data_dedisp_nozerodm'].astype(np.float32)
		onstack = np.zeros(np.size(dat,0))
		offstack = np.zeros(np.size(dat,0))
	
	mstack = np.zeros(1)

	for i in spd_list:
		#loading data
		npz_fn = i
		spectrum = np.load(npz_fn)
		dat = spectrum['Data_dedisp_zerodm'].astype(np.float32) #Data_dedisp_nozerodm for no 0 DM filtering
		freqs = spectrum['freqs_nozerodm']
		text_array = spectrum['text_array']

		#pulse peak index, sampling time, pulsewidth, no. samples
		indpeak = (np.argmax(np.sum(dat,axis=0))).astype(np.float32)
		sampling_time = float(text_array[14])
		pulsewidth = float(text_array[13])
		indfin = float(text_array[7])

		#calculating minimum and maximum index of pulse
		indmin = indpeak - int(math.ceil(pulsewidth/sampling_time))
		indmax = indpeak + int(math.ceil(pulsewidth/sampling_time))

		#subtracting the off-pulse mean from the total data array and dividing by the off-pulse standard dev.
		sm1 = np.mean(dat[:,0:indmin])
		sm2 = np.mean(dat[:,indmax:indfin])
		sm = (sm1+sm2)/2
		strd1 = np.std(dat[:,0:indmin])
		strd2 = np.std(dat[:,indmax:indfin])
		strd = (strd1+strd2)/2
		dat = (dat-sm)/strd

		#Creating normalized off- and on-pulse spectra. Spectra need to be reversed due to format of .spd files.
		offpulse = np.sum(dat[:,0:indmin],axis=1)/float(indmin)
		offpulse = offpulse[::-1]
		onpulse = np.sum(dat[:,indmin:indmax],axis=1)/float(indmax-indmin)
		onpulse = onpulse[::-1]

		#Calculating linear fit to spectrum (ignoring first and last frequency channels).
		sum1 = 1+abs(min(onpulse))
		sum2 = onpulse + sum1
		x = np.log(freqs[m_end:m_start])
		y = np.log(sum2[m_end:m_start])
		m, intercept, r_value, p_value, std_error = stats.linregress(x,y)
		m = float('%.2f'%(m))
		std_error = float('%.2f'%(std_error))
		slope = 'Slope:' + ' ' + str(m) + ' ' + '$\pm$' + ' ' + str(std_error)
		
		#Adding calculated spectra and slope to designated figure.
		if len(spd_list) == 1:
			ax.plot(freqs[m_end:m_start],onpulse[m_end:m_start],'black',label=slope)
			ax.plot(freqs[m_end:m_start],offpulse[m_end:m_start],'grey')
			ax.legend(frameon=False,fontsize=12)

		if len(spd_list) > 1:
			axs[spd_list.index(i)].plot(freqs[m_end:m_start],onpulse[m_end:m_start],c='black',label=slope)
			axs[spd_list.index(i)].plot(freqs[m_end:m_start],offpulse[m_end:m_start],c='grey')
			axs[spd_list.index(i)].legend(frameon=False,fontsize=10)

		#Stacks spectra if -stack is True.
		if stack == True:
			onstack = np.vstack((onstack,onpulse))
			offstack = np.vstack((offstack,offpulse))
			mstack = np.vstack((mstack,m))

		npz_fn = None
		spectrum =  None
		dat = None
		text_array = None
		indfin = None
		sampling_time = None
		pulsewidth = None
		sm1 = None
		sm2 = None
		sm = None
		s = None
		strd1 = None
		strd2 = None
		strd = None
		indpeak = None
		indmin = None
		indmax = None
		onpulse = None
		offpulse = None

	plt.rc('text',usetex=True)
	plt.rc('font',family='serif')

	if stack == True:
		onstack1 = (np.sum(onstack,axis=0))/float(len(spd_list))
		offstack1 = (np.sum(offstack,axis=0))/float(len(spd_list))
		
		x = np.log(freqs[m_end:m_start])
		y = np.log(onstack1[m_end:m_start])
		m, intercept, r_value, p_value, std_error = stats.linregress(x,y)
		m = float('%.2f'%(m))
		std_error = float('%.2f'%(std_error))
		slope = 'Slope:' + ' ' + str(m) + ' ' + '$\pm$' + ' ' + str(std_error)	

		ax.plot(freqs[m_end:m_start],onstack1[m_end:m_start],'black',label=slope)
		ax.plot(freqs[m_end:m_start],offstack1[m_end:m_start],'grey')
		ax.set_xlabel(r'Frequency (MHz)')
		ax.legend(frameon=False,fontsize=12)
	if len(spd_list)>1:
		axs[-1].set_xlabel(r'Frequency (MHz)')
	if len(spd_list)==1:
		ax.set_xlabel(r'Frequency (MHz)')

	plt.show()

main(spd_list)





	







