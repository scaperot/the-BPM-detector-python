BPM Detector in Python
=======================


## Installation

### Activate virtual environment
``` bash
source ./bin/activate
```

### Dependencies
``` bash
pip install --upgrade numpy
pip install --upgrade PyWavelets
pip install --upgrade scipy
pip install --upgrade pdb
pip install --upgrade matplotlib
```

## Test
``` bash
python ./bpm_detection/bpm_detection.py --filename ./data/*.wav 
```

## Original Readme
Implementation of a Beats Per Minute (BPM) detection algorithm, as presented in the paper of G. Tzanetakis, G. Essl and P. Cook titled: "Audio Analysis using the Discrete Wavelet Transform".

You can find it here: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.63.5712

Based on the work done in the MATLAB code located at github.com/panagiop/the-BPM-detector-python.

Process .wav file to determine the Beats Per Minute.

Dependencies: scipy, numpy, pywt, matplotlib
