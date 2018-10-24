import scipy
import numpy as np

########################
## required constants ##
########################

## length unit: micrometres

# Avogadro's number
Na = 6.022e23

# area per lipid um^2
a = 0.627e-6 # N. Kucerka et al. / Biochimica et Biophysica Acta 1808 (2011) 2761–2771 

# half the thickness of bilayer, in micrometres (half the steric bilayer thickness)
# Leftin A, Molugu TR, Job C, Beyer K, Brown MF. Biophysical Journal 107(10):2274–86 (2014) 
# Huber T, Rajamoorthi K, Kurze VF, Beyer K, Brown MF. J Am Chem Soc. 124(2):298–309 (2002) 

l = 5e-3/2 #in micrometres

# default centre-to-centre separation between bilayers
s = 4*l

# lipid index: wavelength-dependent refractive index of egg PC
# eq 10 from B.N Khlebtsov, L.A Kovler, V.A Bogatyrev, N.G Khlebtsov, S.Yu Shchyogolev,
# JQSRT, 79–80, 825-838, (2003)
def li(wl):
    return 1.4713+1.31/(wl*1e3)+4309/(wl*1e3)**2

# medium index: wavelength-dependent refractive index of water
# eq 10 from Amelia G. Van Engen, Scott A. Diddams, and Tracy S. Clement
# Appl. Opt. 37, 5679-5686 (1998)
def mi(wl):
    return 1.313242+15.7834/(wl*1e3)-4382/(wl*1e3)**2+1.1455e6/(wl*1e3)**3

# set up optics for HoloPy
illum_wavelen = 0.4 # set default wavelength to 400 nm
illum_polarization = (0,1)

# set up lognormal probability distribution
def lognormS(mn,sd):
    return np.log(mn/np.sqrt(1+(sd/mn)**2)), np.sqrt(np.log(1+(sd/mn)**2))

def p(x,mu,sigma):
    return np.exp(-(np.log(x)-mu)**2/(2*sigma**2))/(x*sigma*np.sqrt(2*np.pi))