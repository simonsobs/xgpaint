import os
import datetime
import h5py
import numpy              as np
import pylab              as pl
import scipy              as sp

from astropy.io import fits

class HaloLightCone():
    """
    @brief Class describing a halo light cone
    """
    def __init__(self, **kwargs):
        pass

    def copy(self):
        """
        @brief Creates a copy of the halo light cone
        """
        return copy.copy(self)

def ReadPkscLightCone(filename, LightIO=True, omegam=0.31, h=0.68):
    """
    @brief Reads peakpatch output in binary format into an object.
    """    
    halos = HaloLightCone()

    Nreads   = 1
    i = 0

    pkfile   = open(filename,"rb")

    halos.Non        = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]
    halos.RTHmax     = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]
    halos.redshiftin = np.fromfile(pkfile, dtype=np.float32, count=1).astype('float')[0]

    start = int(np.ceil(float(halos.Non)/Nreads)*i)
    end   = int(min( halos.Non , np.ceil(float(halos.Non)/Nreads)*(i+1) ))

    Nread = end-start

    outnum  = 10 #7 floats per halo
    massindex = 6
    if LightIO: 
        outnum = 4
        massindex = 3

    npkdata = Nread*outnum

    pkfile.seek( (3+start*outnum)*4 ,os.SEEK_SET)  #12 byte header plus array
    halos.data = np.fromfile(pkfile, dtype=np.float32, count=npkdata)
    halos.data = np.reshape(halos.data,(Nread,outnum)).astype('float32')
    pkfile.close()
    
    omegam = 0.25
    h = 0.7
    
    print(omegam, h)
    
    rho       = 2.775e11*omegam*h**2
    halos.data[:,massindex] = (4*(np.pi)/3)*(halos.data[:,massindex]**3)*rho # convert 6 from RTH to M

    halos.data = np.column_stack((
        halos.data[:,0], # x in comoving MPC
        halos.data[:,1], # y in comoving MPC
        halos.data[:,2], # z in comoving MPC
        halos.data[:,massindex], # mass in Solar Mass (no h anyhere)!
    ))            

    return halos

def ReadLightConeHeader():

    filename = params.folder+params.inputfile

    pkfile = open(filename,"rb")
    Non    = np.fromfile(pkfile, dtype=np.int32, count=1).astype('int')[0]
    
    return Non

import sys
file = sys.argv[1]

halos = ReadPkscLightCone(file)
f = h5py.File(file.replace('.pksc', '.hdf5'), "w")
f.create_dataset('halos', data=halos.data)
f.close()

