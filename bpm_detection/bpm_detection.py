# Copyright 2012 Free Software Foundation, Inc.
#
# This file is part of The BPM Detector Python
#
# The BPM Detector Python is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# The BPM Detector Python is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with The BPM Detector Python; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

import math
import soundfile as sf

import matplotlib.pyplot as plt
import numpy
import pywt
from scipy import signal
from tkinter import filedialog
from tkinter import *


def read_wav():
	root = Tk()
	root.withdraw()

	file_selected = filedialog.askopenfilename()

	#file_selected = 'C:/Users/Ivan/Desktop/a.wav'
	#file_selected = 'C:/Users/Ivan/Desktop/Simple-drum-loop-116-bpm.wav'

	if file_selected == '':
		exit("File not selected. Exiting.")

	data, fs = sf.read(file_selected)

	# if stereo wav file
	if data.ndim == 2:
		# avrage from both channels
		data = data.T[:][0] + data.T[:][1]
		data /= 2

	return data, fs


# simple peak detection
def peak_detect(data):
	max_val = numpy.amax(abs(data))
	peak_ndx = numpy.where(data == max_val)
	if len(peak_ndx[0]) == 0:  # if nothing found then the max must be negative
		peak_ndx = numpy.where(data == -max_val)

	return peak_ndx


def bpm_detector(data, fs):
	LEVELS = 4
	max_decimation = 2 ** (LEVELS - 1)
	min_ndx = math.floor(60.0 / 220 * (fs / max_decimation))
	max_ndx = math.floor(60.0 / 40 * (fs / max_decimation))
	
	cA, cD = pywt.dwt(data, "db4")
	cD_minlen = math.floor(len(cD) / max_decimation + 1)
	cD_sum = numpy.zeros(cD_minlen)

	for i in range(1, LEVELS):
		# 0) DWT
		cA, cD = pywt.dwt(cA, "db4")

		# 1) Filter
		cD = signal.lfilter([0.01], [1, -0.99], cD)

		# 2) FWR
		cD = abs(cD)

		# 3) Downsampling
		cD = abs(cD[:: (2 ** (LEVELS - i - 1))])

		# 4) Norm
		cD -= numpy.mean(cD)
		cA -= numpy.mean(cA)

		# Recombine the signal before ACF
		cD_sum += cD[0 : cD_minlen]

	# 5) ACRL
	correl = numpy.correlate(cD_sum, cD_sum, "full")

	"""What you need to do is take the last half of your correlation 
	result, and that should be the autocorrelation you are looking for.
	See: https://stackoverflow.com/a/676302
	"""
	midpoint = len(correl) // 2
	correl_midpoint_tmp = correl[midpoint:]

	peak_ndx = peak_detect(correl_midpoint_tmp[min_ndx:max_ndx])

	peak_ndx_adjusted = peak_ndx[0] + min_ndx
	bpm = 60.0 / peak_ndx_adjusted * (fs / max_decimation)
	bpm = round(numpy.median(bpm), 1)

	return bpm


def plot_signal(samps, fs, bpm):
	fig, ax1 = plt.subplots()

	time = numpy.linspace(0, len(samps)/fs, num=len(samps))
	samps_norm = samps / numpy.max(samps)

	ax1.plot(time, samps_norm, linewidth=1)

	ax1.set_xlabel('time (s)')
	ax1.set_title(f'Normalized sound amplitude\nEstimated Beats Per Minute: {bpm}')
	ax1.set_xlim(0, time[-1])
	plt.subplots_adjust(hspace=1)
	plt.show()


if __name__ == "__main__":
	samps, fs = read_wav()

	bpm = bpm_detector(samps, fs)

	print("Estimated Beats Per Minute:", bpm)
	plot_signal(samps, fs, bpm)
