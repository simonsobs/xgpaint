import mpi4py.rc
import datetime

from   pdict     import *
from   cosmolopy import fidcosmo

def getparameters(filename):
    
    import params
    import globals

    dict=pdict()
    dict.read_from_file(filename)

    if 'verbose' in dict: params.verbose = dict['verbose']

    if 'serial'    in dict: params.serial    = dict['serial']
    if 'Nreads'    in dict: params.Nreads    = dict['Nreads']
    if 'format'    in dict: params.format    = dict['format'] 
    if 'folder'    in dict: params.folder    = dict['folder']
    if 'inputfile' in dict: params.inputfile = dict['inputfile']
    if 'outbase'   in dict: params.outbase   = dict['outbase']
    if 'LightIO'   in dict: params.LightIO   = dict['LightIO']

    if 'H0'           in dict: params.H0           = dict['H0']
    if 'omegab'       in dict: params.omegab       = dict['omegab']
    if 'omegac'       in dict: params.omegac       = dict['omegac']
    if 'omegal'       in dict: params.omegal       = dict['omegal']
    if 'scalar_index' in dict: params.scalar_index = dict['scalar_index']
    if 'sigma8'       in dict: params.sigma8       = dict['sigma8']
    if 'min_redshift' in dict: params.min_redshift = dict['min_redshift']
    if 'max_redshift' in dict: params.max_redshift = dict['max_redshift']
    if 'num_redshift' in dict: params.num_redshift = dict['num_redshift']
    if 'octx'         in dict: params.octx         = dict['octx']
    if 'octy'         in dict: params.octy         = dict['octy']
    if 'octz'         in dict: params.octz         = dict['octz']

    if 'box_size' in dict: params.box_size = dict['box_size']
    if 'fov'      in dict: params.fov      = dict['fov']
    if 'flat'     in dict: params.flat     = dict['flat']
    if 'numdens'  in dict: params.numdens  = dict['numdens']
    if 'iwantcat' in dict: params.iwantcat = dict['iwantcat']

    if 'hod' in dict: params.hod = dict['hod']
    if 'LM'  in dict: params.LM  = dict['LM']

    if 'shang_Td'     in dict: params.shang_Td     = dict['shang_Td']
    if 'beta'         in dict: params.shang_beta   = dict['shang_beta']
    if 'shang_eta'    in dict: params.shang_eta    = dict['shang_eta']
    if 'shang_I0'     in dict: params.shang_I0     = dict['shang_I0']
    if 'shang_Mmin'   in dict: params.shang_Mmin   = dict['shang_Mmin']
    if 'shang_Msmin'  in dict: params.shang_Msmin  = dict['shang_Msmin']
    if 'shang_Mpeak'  in dict: params.shang_Mpeak  = dict['shang_Mpeak']
    if 'shang_sigmaM' in dict: params.shang_sigmaM = dict['shang_sigmaM']
    if 'shang_nc2ns'  in dict: params.shang_nc2ns  = dict['shang_nc2ns']
    
    if 'freq_list' in dict: params.freq_list = dict['freq_list'] 

    if 'Inu_norm' in dict: params.Inu_norm = dict['Inu_norm']
    if 'nside'    in dict: params.nside    = dict['nside']

    if 'min_mass'   in dict: params.min_mass   = dict['min_mass']
    if 'M1'         in dict: params.M1         = dict['M1']
    if 'M0'         in dict: params.M0         = dict['M0']
    if 'sigma_logM' in dict: params.sigma_logM = dict['sigma_logM']
    if 'logM_min'   in dict: params.logM_min   = dict['logM_min']
    if 'alpha'      in dict: params.alpha      = dict['alpha']
    if 'mass_type'  in dict: params.mass_type  = dict['mass_type']

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
