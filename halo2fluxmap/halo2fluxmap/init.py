import mpi4py.rc
import datetime

from   pdict     import *
from   cosmolopy import fidcosmo

def getparameters(filename):
    
    import params
    import globals

    dict=pdict()
    dict.read_from_file(filename)

    for item in dir(params):
        if not item.startswith("__"):
            if item in dict: 
                exec("params." + item + " = dict['" + item + "']")

    # derived parameters
    params.norm_freq = params.freq_list[0]
    params.omegam    = params.omegab+params.omegac
    params.h         = float(params.H0)/100
    
    # for cosmolopy 
    params.cosmo                   = fidcosmo
    params.cosmo['h']              = params.h
    params.cosmo['n']              = params.scalar_index
    params.cosmo['omega_M_0']      = params.omegam
    params.cosmo['omega_b_0']      = params.omegab
    params.cosmo['omega_lambda_0'] = params.omegal
    params.cosmo['sigma_8']        = params.sigma8

    # convert model parameters to fiducial units
    params.shang_Mpeak /= globals.Mf

    return params

def initialize(parameterfile):

    params = getparameters(parameterfile)

    if params.serial: mpi4py.rc.initialize = False
    from mpi4py import MPI
    
    if MPI.Is_initialized():
        params.comm     = MPI.COMM_WORLD
        params.rank     = params.comm.Get_rank()
        params.size     = params.comm.Get_size()
        params.parallel = True
    else:
        params.rank     = 0
        params.size     = 1
        params.parallel = False
    
    fmt       = '%H:%M:%S on %m/%d/%Y'
    timestamp = datetime.datetime.now().strftime(fmt)

    if(params.rank==0):
        print ''
        bar = 72*'-'
        print bar
        print 'Running on',params.size,'processor(s)'
        print 'Time:      '+timestamp
        print 'Directory: '+os.getcwd()
        print bar
        print ''
