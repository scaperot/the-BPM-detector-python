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
	LEVELS = 5
	max_decimation = 2 ** (LEVELS - 1)
	min_ndx = math.floor(60.0 / 220 * (fs / max_decimation))
	max_ndx = math.floor(60.0 / 40 * (fs / max_decimation))
	
	cA, cD = pywt.dwt(data, "db4")
	cD_minlen = math.floor(len(cD) / max_decimation + 1)
	cD_sum = numpy.zeros(cD_minlen)

	for loop in range(1, LEVELS):
		# 1) DWT
		cA, cD = pywt.dwt(cA, "db4")

		# 2) Filter
		a = 0.99
		cD = signal.lfilter([1-a], [a], cD)

		# 5) Decimate for reconstruction later.
		cD = abs(cD[:: (2 ** (LEVELS - loop - 1))])

		# 4) Subtract out the mean.
		cD -= numpy.mean(cD)
		cA -= numpy.mean(cA)

		# 6) Recombine the signal before ACF
		cD_sum += cD[0 : cD_minlen]

	# ACF
	correl = numpy.correlate(cD_sum, cD_sum, "full")

	midpoint = math.floor(len(correl) / 2)
	correl_midpoint_tmp = correl[midpoint:]

	peak_ndx = peak_detect(correl_midpoint_tmp[min_ndx:max_ndx])

	peak_ndx_adjusted = peak_ndx[0] + min_ndx
	bpm = 60.0 / peak_ndx_adjusted * (fs / max_decimation)

	return bpm


if __name__ == "__main__":
	WINDOWS = 3

	samps, fs = read_wav()

	nsamps = len(samps)
	window_samps = int(WINDOWS * fs)
	samps_ndx = 0  # First sample in window_ndx
	max_window_ndx = math.floor(nsamps / window_samps)
	bpms = numpy.zeros(max_window_ndx)

	# Iterate through all windows
	for window_ndx in range(0, max_window_ndx):

		# Get a new set of samples
		data = samps[samps_ndx : samps_ndx + window_samps]
		if not ((len(data) % window_samps) == 0):
			raise AssertionError(str(len(data)))

		bpm = bpm_detector(data, fs)

		bpms[window_ndx] = round(int(bpm))

		# Iterate at the end of the loop
		samps_ndx += window_samps

	bpm_median = int(numpy.median(bpms))
	print(sorted(bpms))
	print("Estimated Beats Per Minute:", bpm_median)

	fig, ax1 = plt.subplots()
	time = numpy.linspace(0, nsamps/fs, num=nsamps)
	ax1.plot(time, samps/numpy.max(samps), linewidth=1)
	ax1.set_xlabel('time (s)')
	ax1.set_title(f'Normalized sound amplitude\nEstimated Beats Per Minute: {bpm_median}')

	plt.subplots_adjust(hspace=1)
	plt.show()
