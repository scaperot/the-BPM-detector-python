from __future__ import division, print_function
import wave, array, math, time, argparse, sys
import numpy, pywt
import wavio
from scipy import signal
import pdb
import matplotlib.pyplot as plt
import seaborn
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.integrate import cumtrapz


def smooth(data):
    """Smooth data with 3rd order Butterworth filter."""
    b, a = signal.butter(3, 7 / len(data))
    data = signal.filtfilt(b, a, data)
    return data


def read_wav(filename):
    wf = wavio.read(filename)
    samps = wf.data
    fs = wf.rate
    return samps.flatten(), fs


# print an error when no data can be found
def no_audio_data():
    print("No audio data for sample, skipping...")
    return None, None


# simple peak detection
def peak_detect(data):
    max_val = numpy.amax(abs(data))
    peak_ndx = numpy.where(data==max_val)
    if len(peak_ndx[0]) == 0: #if nothing found then the max must be negative
        peak_ndx = numpy.where(data == -max_val)
    return peak_ndx


def bpm_detector(data, fs):
    cA = []
    cD = []
    correl = []
    cD_sum = []
    levels = 4
    max_decimation = 2**(levels - 1)
    min_ndx = 60./ 220 * (fs/max_decimation)
    max_ndx = 60./ 40 * (fs/max_decimation)

    for loop in range(0, levels):
        cD = []
        # 1) DWT
        if loop == 0:
            [cA,cD] = pywt.dwt(data,'db4')
            cD_minlen = len(cD) // max_decimation + 1
            cD_sum = numpy.zeros(cD_minlen)
        else:
            [cA,cD] = pywt.dwt(cA, 'db4')
        # 2) Filter
        cD = signal.lfilter([0.01], [1 -0.99], cD)

        # 4) Subtractargs.filename out the mean.

        # 5) Decimate for reconstruction later.
        cD = abs(cD[::(2**(levels - loop - 1))])
        cD = cD - numpy.mean(cD)
        # 6) Recombine the signal before ACF
        #    essentially, each level I concatenate
        #    the detail coefs (i.e. the HPF values)
        #    to the beginning of the array
        cD_sum = cD[0:cD_minlen] + cD_sum

    if [b for b in cA if b != 0.0] == []:
        return no_audio_data()
    # adding in the approximate data as well...
    cA = signal.lfilter([0.01], [1 -0.99], cA)
    cA = abs(cA)
    cA = cA - numpy.mean(cA)
    cD_sum = cA[0:cD_minlen] + cD_sum

    # ACF
    correl = numpy.correlate(cD_sum, cD_sum, 'full')

    midpoint = len(correl) // 2
    correl_midpoint_tmp = correl[midpoint:]
    peak_ndx = peak_detect(correl_midpoint_tmp[int(min_ndx):int(max_ndx)])
    if len(peak_ndx) > 1:
        return no_audio_data()

    peak_ndx_adjusted = peak_ndx[0] + min_ndx
    bpm = 60./ peak_ndx_adjusted * (fs/max_decimation)
    return bpm, correl


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process .wav file to '
                                     'determine the Beats Per Minute.')
    parser.add_argument('filename', help='.wav file for processing')
    parser.add_argument('--window', type=float, default=3,
                        help='size of the the window (seconds) that will be '
                        'scanned to determine the bpm. Typically less than 10 '
                        'seconds. [3]')
    parser.add_argument('--plot', '-p', action='store_true', default=False,
                        help='Plot tempo with matplotlib')
    parser.add_argument('--write-midi', '-w', help='Write MIDI to file')
    parser.add_argument('--midi-note', '-n', help='MIDI note number for click '
                        'track output', type=int, default=67)

    args = parser.parse_args()
    samps, fs = read_wav(args.filename)

    data = []
    correl = []
    bpm = 0
    n = 0
    nsamps = len(samps)
    seconds = numpy.arange(len(samps)) / fs
    window_samps = int(args.window*fs)
    samps_ndx = 0  #first sample in window_ndx
    max_window_ndx = nsamps // window_samps
    bpms = numpy.zeros(max_window_ndx)
    seconds_mid = numpy.zeros(max_window_ndx)

    #iterate through all windows
    for window_ndx in range(0, max_window_ndx):

        #get a new set of samples
        #print n,":",len(bpms),":",max_window_ndx,":",fs,":",nsamps,":",samps_ndx
        data = samps[samps_ndx:samps_ndx + window_samps]
        seconds_mid[window_ndx] = \
                seconds[samps_ndx:samps_ndx + window_samps].mean()
        if not ((len(data) % window_samps) == 0):
            raise AssertionError(str(len(data)))

        bpm, correl_temp = bpm_detector(data, fs)
        if bpm is None:
            continue
        bpms[window_ndx] = bpm
        correl = correl_temp

        #iterate at the end of the loop
        samps_ndx = samps_ndx + window_samps
        n = n + 1 #counter for debug...

    # Smoothing from http://stackoverflow.com/questions/28536191\
    # /how-to-filter-smooth-with-scipy-numpy
    filtered = lowess(bpms, seconds_mid, is_sorted=True, frac=0.3, it=0)[:, 1]
    # Create beats array
    bps = filtered / 60
    beats = cumtrapz(bps, seconds_mid, initial=0)

    bpm = numpy.median(bpms)
    print('Completed. Estimated Beats Per Minute:', bpm)

    if args.write_midi is not None:
        from midiutil.MidiFile import MIDIFile
        pitch = ags.midi_note
        track = 0
        channel = 0
        duration = 0.125 # In beats
        volume = 127     # 0-127, as per the MIDI standard
        mf = MIDIFile(1) # One track, defaults to format 1 (tempo track
                         # automatically created)
        # Set time signature to 1/4
        # mf.addTimeSignature(track, 0, 1, 2, 24, notes_per_quarter=8)
        # Create new beats array that's evenly spaced and map to this
        beats_even = numpy.arange(0, beats.max()*1.1, 0.5)
        bpm_mapped = numpy.interp(beats_even, beats, filtered)
        for tempo, beat in zip(bpm_mapped, beats_even):
            mf.addTempo(track, beat, tempo)
            mf.addNote(track, channel, pitch, beat, duration, volume)
        with open(args.write_midi, "wb") as output_file:
            mf.writeFile(output_file)

    if args.plot:
        plt.plot(beats, bpms)
        plt.plot(beats, filtered)
        plt.xlabel("Beat")
        plt.ylabel("BPM")
        plt.show()
