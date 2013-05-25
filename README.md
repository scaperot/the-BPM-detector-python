the-BPM-detector-python
=======================
Implementation of a Beats Per Minute (BPM) detection algorithm, as presented in the paper of G. Tzanetakis, G. Essl and P. Cook titled: "Audio Analysis using the Discrete Wavelet Transform".

You can find it here: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.63.5712

Based on the work done in the MATLAB code located at github.com/panagiop/the-BPM-detector-python.

Process .wav file to determine the Beats Per Minute.

optional arguments:
  -h, --help           show this help message and exit
  --filename FILENAME  .wav file for processing
  --window WINDOW      size of the the window (seconds) that will be scanned
                       to determine the bpm. Typically less than 10 seconds.
                       [3]

