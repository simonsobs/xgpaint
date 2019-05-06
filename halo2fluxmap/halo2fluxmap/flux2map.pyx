import healpy as hp
import numpy as np
cimport numpy as np

from numpy.random import *
from .hod   import *
from utils import *

def makemap(np.ndarray x,np.ndarray y, np.ndarray z,np.ndarray F, nside):
    
    pixind = hp.vec2pix(nside,x,y,z) #to get same orientation as pks2map
    hmap = np.zeros(hp.nside2npix(nside))
    for i in range(np.shape(F)[0]):
        hmap[pixind[i]] += F[i]
        
    return hmap
    
def makemapflat(np.ndarray x,np.ndarray y, np.ndarray z,np.ndarray F, nside,fov):
    fov = fov*np.pi/180           
    dp  = fov/nside  
    #+fov/2 to ensure numbers in pixind are positive
    
    dm = (x>0)
    z = z[dm]
    y = y[dm]         
    x = x[dm] 
    F = F[dm]
    
    r = np.sqrt(z**2+x**2+y**2)
    theta = y / r
    phi   = z / r
    theta = np.arcsin(theta) 
    phi   = np.arcsin(phi)
    
    dm = (abs(theta)<fov/2)
    theta = theta[dm]         
    phi   = phi[dm] 
    F = F[dm]
    dm = ((abs(phi)<fov/2))
    theta = theta[dm] 
    phi   = phi[dm]
    F     = F[dm]
    
    pixind = np.floor((theta+fov/2)/dp) + nside*np.floor((phi+fov/2)/dp) 
    pixind = pixind.astype(int)
    
    hmap = np.zeros(nside**2)
    for i in range(np.shape(F)[0]):
        hmap[pixind[i]]+=F[i] 

    return hmap
    
def cen2sat(np.ndarray cen, np.ndarray n):
    
    '''
    Parameters
    cen[:,0]:  array of central masses
    n: 	       array of central number of satellites for each central
    
    Returns
    sat[:,0]: array containing parent halo mass for each satellite
    '''
    
    cdef int N_sat = np.sum(n)
    cdef int N_cen = np.shape(cen)[0]
    cdef int N_prp = np.shape(cen)[1]
    
    sat = np.zeros((N_sat,N_prp),dtype='float32')
    
    cdef int count = 0
    cdef int i
    for i in range(N_cen):
        sat[count:count+n[i],:] = cen[i,:]
        count += n[i]
        
    return sat    

def cen2sat_masses(np.ndarray cen, np.ndarray n, np.ndarray nmean):
    
    '''
    Parameters
    cen[:,0]: array of central masses
    n: 	      array of number of satellites for each central
    
    Returns
    msat[:]: array containing subhalo mass for each satellite
    '''
    
    cdef int N_sat = np.sum(n)
    cdef int N_cen = np.shape(cen)[0]
    cdef int N_prp = np.shape(cen)[1]
    
    msat = np.zeros(N_sat,dtype='float32')
    
    # Make function of mass fraction as a function of number of satellites
    muofn = make_muofn()
    
    cdef int count = 0
    cdef int i
    for icen in range(N_cen):
        N_sat    = n[icen]
        N_satbar = nmean[icen]
        M_cen    = cen[icen,0]
        Rank = uniform(0.0,N_satbar,N_sat)
        mu   = muofn(Rank)
        
        for isat in range(N_sat):
            msat[count+isat] = mu[isat] * M_cen
        count += n[icen]
        
    return msat

