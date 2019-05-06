#from    mods import * 
#from   utils import *
#from halocat import *
#from globals import *
import mpi4py 
import sys

import healpy as hp

from   utils   import *
from   halocat import *

def report_local(string,i,Ncen):
    if params.rank == 0:
        sys.stdout.write('\n '+params.justify+str(i+1)+' of '+str(params.Nreads)+
                         ' chunks '+string+' Ncen = '+str(Ncen))
        sys.stdout.flush()
        
def read_catalog():
    for i in range(params.Nreads):  

        ceni = get_catalog(i)
        report_local('read',i,np.shape(ceni)[0])
        ceni = cull_catalog(ceni)
        report_local('culled',i,np.shape(ceni)[0])
        ceni = shuffle_catalog(ceni)

        if i==0:
            cen = distribute_catalog(ceni)
        else:
            cen  = np.concatenate((cen,distribute_catalog(ceni)))
            
        report_local('read, shuffled, and culled',i,np.shape(cen)[0])

    if params.rank == 0: print ''
 
    del ceni
    return cen

def get_catalog(i):
    '''
    Gets the file. Checks what the format given in the parameters file 
    is and returns x,y,z,M values as column stack data accordingly.
    '''
    extension  = params.format 
    if extension=='pksc':
        peakdata  = ReadPkscLightCone(i)
        data      = peakdata.data            
    else:
        print" extension "+extension+" not recognized"

    return data

def get_catalog_Non():
    Non  = ReadLightConeHeader()
    
    return Non

def writemap(base,intensity):

    if(params.rank != 0): return
        
    if params.flat==0:
        filename = base+'.fits' 
        hp.write_map(filename,intensity)
    elif params.flat==1:
        filename = base+'.map' 
        intensity = intensity.astype('float32')

        fov      = np.array([params.fov*np.pi/180]).astype('float32')
        npix     = np.array([params.nside]).astype('int32')

        outfile = open(filename,'wb') 
        npix.tofile(outfile)
        npix.tofile(outfile)
        fov.tofile(outfile)
        fov.tofile(outfile)
        intensity.tofile(outfile)

    report('Wrote Map '+filename,2)
