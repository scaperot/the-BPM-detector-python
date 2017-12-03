BPM Detector in Python
=======================


## Installation

### Activate virtual environment
``` bash
source ./bin/activate
```

### Dependencies
- numpy
- PyWavelets
- scipy
- pdb
- matplotlib
``` bash
pip install -r requirements.txt
```

## Test
``` bash
python ./bpm_detection/bpm_detection.py --filename ./data/*.wav 
```

## Help

### Converting mp3
``` bash
mpg123 -w [path/to/original.mp3] [./data/wavfile.wav]
```

## Original Readme
Implementation of a Beats Per Minute (BPM) detection algorithm, as presented in the paper of G. Tzanetakis, G. Essl and P. Cook titled: "Audio Analysis using the Discrete Wavelet Transform".

You can find it here: http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.63.5712

Based on the work done in the MATLAB code located at github.com/panagiop/the-BPM-detector-python.

Process .wav file to determine the Beats Per Minute.

Dependencies: scipy, numpy, pywt, matplotlib
