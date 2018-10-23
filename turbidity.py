import scipy
import pylab
import numpy as np
import holopy as hp
from holopy.scattering import calc_holo, Sphere
from holopy.scattering import calc_cross_sections
from constants import *

#############################
## total sample scattering ##
#############################

# length unit: micrometres

# set up oligolamellar vesicles Sphere object
# s is the separation between bilayers, centre-to-centre
# r_olv should go [r-2s-l,r-2s+l,r-s-l,r-s+l,r-l,r+l]
# i.e. [r-ns-l,r-ns+l] with n = [0,nlayers]

def n_olv(layers, medium_index, lipid_index):
    return [medium_index, lipid_index]*layers

def r_olv(layers, r, s, t=l):
    r_temp = np.array([r-t,r+t]*layers)
    s_temp = np.array([x for j in zip(np.arange(layers-1,-1,-1),np.arange(layers-1,-1,-1)) for x in j])*s
    return r_temp-s_temp

# surface areas of vesicles: calculate radius from middle of each bilayer
def area_guv(Rs):
    return 2*4*np.pi*Rs**2

def area_olv(layers,s,Rs):
    areas=[]
    for r in Rs:
        list_r = r_olv(layers, r, s)
        list_r = (np.array(list_r[::2])+np.array(list_r[1::2]))/2 # take middle of each bilayer
        total_area_for_r = (2*4*np.pi*list_r**2).sum()
        areas.append(total_area_for_r)
    return areas

def extinction_olv_wld(layers, wls, rad, s, t, flag=0):
    extinction_o = []
    for illum_wavelen in wls:
        coated_sphere = Sphere(n = n_olv(layers, mi(illum_wavelen), li(illum_wavelen)), r = r_olv(layers,rad,s,t))
        x_secs = calc_cross_sections(coated_sphere, mi(illum_wavelen), illum_wavelen, illum_polarization)
        extinction_o.append(x_secs[flag])
    return np.array(extinction_o)

def extinction_guv_wld(wls,rad,t=l,flag=0):
    extinction_v = []
    for illum_wavelen in wls:
        coated_sphere = Sphere(n = (mi(illum_wavelen),li(illum_wavelen)), r = (rad-t, rad+t))
        x_secs = calc_cross_sections(coated_sphere, mi(illum_wavelen), illum_wavelen, illum_polarization)
        extinction_v.append(x_secs[flag])
    return np.array(extinction_v)

def sample_extinction(extinction,area,concentration=5e-3):
    return 10**3*0.01*concentration*Na*a*1e-12*extinction/area

# we calculate the sample extinction for each radius, weigh it by the size distribution, and add it all together

def weighted_turb(t, c, m, s, R, wls):
    '''
    t is half the thickness of the membrane (in micrometres)
    c is the concentration of the amphiphile, in M (moles/L)
    mn, sd, are the lognormal 
    '''
    step = R[1]-R[0]
# initialise the cumulative extinction, weighted by the size distribution
    sdist = 0
    psum = 0
    for r in R:
        sdist+=(step*p(r, m, s)*sample_extinction(extinction_guv_wld(wls=wls,t=t,rad=r),area_guv(r),concentration=c))
        psum+=step*p(r, m, s)
    return sdist,psum