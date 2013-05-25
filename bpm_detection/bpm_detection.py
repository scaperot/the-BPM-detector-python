import wave
import array
import math
import numpy
import pywt
from scipy import signal
import pdb
import time
import matplotlib.pyplot as plt

#Type code 	C Type 	Python 	Type 	Minimum size in bytes
#'c' 	  char 	character 	1
#'b' 	  signed char 	int 	1
#'B' 	  unsigned char 	int 	1
#'u' 	  Py_UNICODE 	Unicode character 	2 (see note)
#'h' 	  sig      limitless     ned short 	int 	2
#'H' 	  unsigned short 	int 	2
#'i' 	  signed int 	int 	2
#'I' 	  unsigned int 	long 	2
#'l' 	  signed long 	int 	4
#'L' 	  unsigned long 	long 	4
#'f' 	  float 	float 	4
#'d' 	  double 	float 	
def choose_type(nbytes):
    return; 


def bpm_detection(filename,window=0):
    bpm = 0
    levels = 4; #TODO: should be arg
    overlap = 0.5; #TODO: should be arg

    #open file, get metadata for audio
    try:
        wf = wave.open(filename,'rb')
    except IOError, e:
        print e
        return 0,0,0,0
    
    # typ = choose_type( wf.getsampwidth() ) #TODO: implement choose_type
    nsamps = wf.getnframes();
    assert(nsamps > 0);
    
    
    fs = wf.getframerate()
    assert(fs > 0)

    # read entire file and make into an array
    samps = list(array.array('i',wf.readframes(nsamps)))
    print 'Read', nsamps,'samples from', filename
    try:
        assert(nsamps == len(samps))
    except AssertionError, e:
        print  nsamps, "not equal to", len(samps)
        #return 0,0,0,0
        
    #iterate through samples based on window size.  overlapping by 50% each time.
    window_samps = window*fs; #TODO fix me...
    if (window_samps % 2) == 0:
        window_step_size = int(window_samps * overlap);
    else:
        window_step_size = int((window_samps-1) * overlap);

    window_ndx = 1; #current window we are processing
    max_window_ndx = nsamps / window_samps;
    samps_ndx = 0;  #first sample in window_ndx 
    accum_correl = []
    #iterate through all windows
    while window_ndx < max_window_ndx:
        print 'Window #:',window_ndx;
        cA = [] 
        cD = []
        final_signal = []#numpy.zeros((levels+1)*window_samps).reshape((levels+1),window_samps);
        data = []
        correl = []

        #get a new set of samples
        data = samps[samps_ndx:samps_ndx+window_samps]
        if not ((len(data) % window_samps) == 0):
            raise AssertionError( str(len(data) ) ) 
        print 'Read', len(data), 'samples. From',str(samps_ndx+1),'to',str(samps_ndx+window_samps)
        
        
        #4 level loop - #TODO: make N level loop...
        cD_sum = []
        max_decimation = 2**(levels-1);
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

            # 4) Subtract out the mean.

            # 5) Decimate for reconstruction later.
            cD = abs(cD[::(2**(levels-loop-1))]);
            
            cD = cD - numpy.mean(cD);
            # 6) Recombine the signal before ACF
            #    essentially, each level I concatenate 
            #    the detail coefs (i.e. the HPF values)
            #    to the beginning of the array
            cD_sum = cD[0:cD_minlen] + cD_sum;


        cA = signal.lfilter([0.01],[1 -0.99],cA);
        cA = abs(cA);
        cA = cA - numpy.mean(cA);
        cD_sum = cA[0:cD_minlen] + cD_sum;
        
        # ACF
        correl = numpy.correlate(cD_sum,cD_sum,'full')
        # integration...
        if window_ndx == 1:
            accum_correl = numpy.zeros(len(correl))
        accum_correl = (correl+accum_correl) / len(correl) 
        
        #iterate at the end of the loop
        window_ndx = window_ndx + 1;
        samps_ndx = samps_ndx+window_step_size;

    # Peak detection
    accum_zero = len(correl) / 2
    correl_final = accum_correl[accum_zero:]
    
    min_ndx = 60./ 220 * (fs/max_decimation)
    max_ndx = 60./ 40 * (fs/max_decimation)

    val = numpy.amax(abs(correl_final[min_ndx:max_ndx]))
    
    k = numpy.where(correl_final[min_ndx:max_ndx]==val)
    if len(k[0]) == 0:
        k = numpy.where(correl_final[min_ndx:max_ndx]==-val)
    j = k[0]+min_ndx;
    bpm = 60./ j * (fs/max_decimation)
    print 'Completed.  Estimated BPM =', bpm
 
    n = range(0,len(correl_final))
    plt.plot(n,abs(correl_final)); 
    #plt.xlim(min_ndx,max_ndx)
    plt.show();
    #plt.show(False);
    #time.sleep(10);
    #plot.close();
        
     
    return 0, 0, bpm, 0

if __name__ == '__main__':
    window_secs = 1; #time window for wavelets over the whole file
    signal,corr,bpm,cd = bpm_detection('GreatIam.wav', window_secs); 
