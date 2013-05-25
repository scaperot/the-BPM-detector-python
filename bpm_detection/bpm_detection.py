import wave, array, math, time, argparse, sys
import numpy, pywt
from scipy import signal
import pdb
import matplotlib.pyplot as plt


def bpm_detection(filename,window=3):
    bpm = 0
    levels = 4; 
    overlap = 0.5; 

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
    #print 'Read', nsamps,'samples from', filename
    try:
        assert(nsamps == len(samps))
    except AssertionError, e:
        print  nsamps, "not equal to", len(samps)
                
    #iterate through samples based on window size.  overlapping by 50% each time.
    window_samps = window*fs; 
    if (window_samps % 2) == 0:
        window_step_size = int(window_samps * overlap);
    else:
        window_step_size = int((window_samps-1) * overlap);

    window_ndx = 1; #current window we are processing
    samps_ndx = 0;  #first sample in window_ndx 
    accum_correl = []
    max_window_ndx = nsamps / window_samps;

    
    #iterate through all windows
    while window_ndx < max_window_ndx:
        sys.stdout.write('=');
        cA = [] 
        cD = []
        data = []
        correl = []
        cD_sum = []

        #get a new set of samples
        data = samps[samps_ndx:samps_ndx+window_samps]
        if not ((len(data) % window_samps) == 0):
            raise AssertionError( str(len(data) ) ) 
        #print 'Read', len(data), 'samples. From',str(samps_ndx+1),'to',str(samps_ndx+window_samps)
        
        
        #4 level loop
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

        # adding in the approximate data as well...    
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

    sys.stdout.write('\n');

    # Peak detection - just involves only picking a peak within the window that we are 
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
     
    return bpm, correl_final

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process .wav file to determine the Beats Per Minute.')
    parser.add_argument('--filename', required=True,
                   help='.wav file for processing')
    parser.add_argument('--window', type=int, default=3,
                   help='size of the the window (seconds) that will be scanned to determine the bpm.  Typically less than 10 seconds. [3]')

    args = parser.parse_args()
    
    bpm,plot_data = bpm_detection(args.filename, args.window);
    print 'Completed.  Estimated Beats Per Minute:', bpm
    
    n = range(0,len(plot_data))
    plt.plot(n,abs(plot_data)); 
    plt.show(False); #plot non-blocking
    time.sleep(10);
    plt.close();
