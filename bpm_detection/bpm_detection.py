import wave, array, math, time, argparse, sys
import numpy, pywt
from scipy import signal
import pdb
import matplotlib.pyplot as plt


#class BpmDetection:
#    def __init__(self,filename='/dev/rand'):
    

def read_wav(filename):

    #open file, get metadata for audio
    try:
        wf = wave.open(filename,'rb')
    except IOError, e:
        print e
        return

    # typ = choose_type( wf.getsampwidth() ) #TODO: implement choose_type
    nsamps = wf.getnframes();
    assert(nsamps > 0);

    fs = wf.getframerate()
    assert(fs > 0)

    # read entire file and make into an array
    samps = list(array.array('i',wf.readframes(nsamps)))
    #print 'Read', nsamps,'samples from', filename
    try:
        assert(nsamps == len(samps))
    except AssertionError, e:
        print  nsamps, "not equal to", len(samps)
    
    return samps, fs
    
# simple peak detection
def peak_detect(data):
    max_val = numpy.amax(abs(data)) 
    peak_ndx = numpy.where(data==max_val)
    if len(peak_ndx[0]) == 0: #if nothing found then the max must be negative
        peak_ndx = numpy.where(data==-max_val)
    return peak_ndx
    
def bpm_detector(data):
    cA = [] 
    cD = []
    correl = []
    cD_sum = []
    levels = 4
    max_decimation = 2**(levels-1);
    min_ndx = 60./ 220 * (fs/max_decimation)
    max_ndx = 60./ 40 * (fs/max_decimation)
    
    for loop in range(0,levels):
        cD = []
        # 1) DWT
        if loop == 0:
            [cA,cD] = pywt.dwt(data,'db4');
            cD_minlen = len(cD)/max_decimation+1;
            cD_sum = numpy.zeros(cD_minlen);
        else:
            [cA,cD] = pywt.dwt(cA,'db4');
        # 2) Filter
        cD = signal.lfilter([0.01],[1 -0.99],cD);

        # 4) Subtractargs.filename out the mean.

        # 5) Decimate for reconstruction later.
        cD = abs(cD[::(2**(levels-loop-1))]);
        
        cD = cD - numpy.mean(cD);
        # 6) Recombine the signal before ACF
        #    essentially, each level I concatenate 
        #    the detail coefs (i.e. the HPF values)
        #    to the beginning of the array
        cD_sum = cD[0:cD_minlen] + cD_sum;

    # adding in the approximate data as well...    
    cA = signal.lfilter([0.01],[1 -0.99],cA);
    cA = abs(cA);
    cA = cA - numpy.mean(cA);
    cD_sum = cA[0:cD_minlen] + cD_sum;
    
    # ACF
    correl = numpy.correlate(cD_sum,cD_sum,'full') 
    
    midpoint = len(correl) / 2
    correl_midpoint_tmp = correl[midpoint:]
    peak_ndx = peak_detect(correl_midpoint_tmp[min_ndx:max_ndx]);
    peak_ndx_adjusted = peak_ndx[0]+min_ndx;
    bpm = 60./ peak_ndx_adjusted * (fs/max_decimation)
    print bpm
    return bpm,correl
    
    
    
def bpm_detection(samps,fs,window=3):
    bpm = 0
    nsamps = len(samps)
    window_samps = int(window*fs)         
    window_ndx = int(1); #current window we are processing
    samps_ndx = 0;  #first sample in window_ndx 
    max_window_ndx = nsamps / window_samps;

    
    data = []
    correl=[]
    bpms = numpy.zeros(max_window_ndx)

    #iterate through all windows
    while window_ndx < max_window_ndx:


        #get a new set of samples
        data = samps[samps_ndx:samps_ndx+window_samps]
        if not ((len(data) % window_samps) == 0):
            raise AssertionError( str(len(data) ) ) 
        #print 'Read', len(data), 'samples. From',str(samps_ndx+1),'to',str(samps_ndx+window_samps)
        
        bpms[window_ndx],correl = bpm_detector(data)

        # integration...
        #if window_ndx == 1:
        #    accum_correl = numpy.zeros(len(correl))
        #accum_correl = (correl+accum_correl) / len(correl)
        
        #iterate at the end of the loop
        window_ndx = window_ndx + 1;
        samps_ndx = samps_ndx+window_samps;

    # Peak detection - just involves only picking a peak within the window that we are 
    #accum_zero = len(correl) / 2
    #correl_final = accum_correl[accum_zero:]
    
    #val = numpy.amax(abs(correl_final[min_ndx:max_ndx]))    
    #k = numpy.where(correl_final[min_ndx:max_ndx]==val)
    #if len(k[0]) == 0:
    #    k = numpy.where(correl_final[min_ndx:max_ndx]==-val)
    #j = k[0]+min_ndx;
    #bpm_final = 60./ j * (fs/max_decimation)
    
    bpm = numpy.median(bpms)
    return bpm, correl

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process .wav file to determine the Beats Per Minute.')
    parser.add_argument('--filename', required=True,
                   help='.wav file for processing')
    parser.add_argument('--window', type=float, default=3,
                   help='size of the the window (seconds) that will be scanned to determine the bpm.  Typically less than 10 seconds. [3]')

    args = parser.parse_args()
    samps,fs = read_wav(args.filename)
    bpm,plot_data = bpm_detection(samps,fs, args.window);
    print 'Completed.  Estimated Beats Per Minute:', bpm
    
    n = range(0,len(plot_data))
    plt.plot(n,abs(plot_data)); 
    plt.show(False); #plot non-blocking
    time.sleep(10);
    plt.close();
