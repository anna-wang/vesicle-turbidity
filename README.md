# vesicle-turbidity
### A. Wang
### awang@molbio.mgh.harvard.edu
### anna.wang@unsw.edu.au

This code is for calculating the turbidity of a sample of vesicles. This turbidity is the
optical depth of the sample, and is related to the 'absorbance' as measured on a 
spectrophotometer. For further details, please see the related paper

Wang, Chan Miller, and Szostak, Biophysical Journal, 2019
https://www.sciencedirect.com/science/article/pii/S0006349519300220

How to use this code:
The relevant constants are located in 'constants.py', and the functions are in 'turbidity.py'.
To calculate the turbidity, we use the scattering cross sections calculated from HoloPy.

1. Please first install HoloPy:
   https://holopy.readthedocs.io/en/latest/
2. Once HoloPy is installed and running, open the example_calculations iPython notebook and try
   to calculate the turbidty spectra for different samples of vesicles.
