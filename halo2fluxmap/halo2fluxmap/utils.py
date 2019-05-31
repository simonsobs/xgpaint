import datetime
import params

import numpy              as np
import scipy              as sp
import cosmolopy.distance as cd

from   scipy.integrate    import *
from   scipy.interpolate  import *
from   cosmolopy          import *

# hard coded (rho_{crit} / h^2)
# should replace in params
RHO_CRIT_OVER_H2 = 2.78e11

def r2m(radius):
    """Compute mass of a sphere of of critical density, given radius.
    
    Parameters
    ----------
    radius : float
        radius in Mpc
        
    Returns
    -------
    float
        mass in Msun
    """
    # rho is physical density in Msun / Mpc^3
    rho = RHO_CRIT_OVER_H2 * params.omegam * (params.h**2)
    return 4./3.*np.pi*rho*radius**3

def m2r(mass):
    """Compute radius of a sphere of critical density, given mass.
    
    Parameters
    ----------
    mass : float
        mass in Msun
        
    Returns
    -------
    float
        radius in Mpc
    """
    rho = RHO_CRIT_OVER_H2 * params.omegam * (params.h**2)
    return (3.*mass/4./np.pi/rho)**(1./3.)

def random_phi(N):
    """Sample uniformly in (0,2pi).
    
    This is used for phi in spherical coordinates, the 
    azimuthal angle in physics convention. 
    
    Parameters
    ----------
    N : int
        number of samples to generate
        
    Returns
    -------
    ndarray of float32 of length N
        random samples
    """
    return (np.random.uniform(size=N)*2*np.pi).astype('float32')

def random_theta(N):
    """Sample theta uniformly over the sphere.
    
    This is used for theta in spherical coordinates, the 
    polar angle in physics convention. 
    
    Parameters
    ----------
    N : int
        number of samples to generate
        
    Returns
    -------
    ndarray of float32 of length N
        random samples
    """
    return (np.arccos(2*np.random.uniform(size=N)-1)).astype('float32')

def mz2c(mass, z):
    """Get halo concentration from mass and redshift.
    
    These are from fits to the WMAP9 cosmology for Duffy et al. 
    2008 eq. 4, and describing the full (F) population for redshift
    0-2, described in Duffy et al. 2008 Table 1.
    
    Parameters
    ----------
    mass : float
        halo mass in Msun
    z : float
        redshift
    
    Returns
    -------
    float
        halo concentration parameter
    """
    return 7.85 * (m / (2e+12/params.h))**-0.081 / (1+z)**0.71

def r2z(r):
    """Compute comoving distance from redshift.
    
    This currently computes the distance by generating an array of
    distances and redshifts from cosmolopy.
    
    Parameters
    ----------
    r : float
        Radius (Mpc)
        
    Returns
    -------
    float
        comoving distance
    """
    zrange   = np.linspace(0,6,1000)
    r_to_z   = sp.interpolate.interp1d(
        cd.comoving_distance_transverse(
            zrange, **params.cosmo), zrange)
    return r_to_z(r).astype('float32')

def report(description,verbosity):
    """Print a string only if on MPI rank 0 and sufficient verbosity.
    
    Parameters
    ----------
    description : string
        string to be printed by the rank 0 process
    verbosity : int
        will be compared to params.verbosity
    """
    if(params.rank==0 and params.verbose>=verbosity): 
        print(params.justify,description)

def resetoverhead():
    params.overhead = params.proc.get_memory_info().rss        

def check_memory(description,N):
  
    memt = params.proc.get_memory_info().rss
    mem = memt - params.overhead
    mem_per_object = mem/(N*1.0)
    flt_per_object = mem_per_object / 4

    mem /= 1024.**3
    memt /= 1024.**3

    params.maxmem=max(params.maxmem,memt)

    print('                     ',description,memt,'GB total')
    print('                         ',mem_per_object,'bytes per array )element'
    print('                         ',flt_per_object,'floats per array element')
  
    return

def distribute_catalog(data):

#    report('Distributing...',2)

    N = np.shape(data)[0]

    Nextra = N % params.size

    Nlocal = (N-Nextra) / (params.size)
        
    start = Nlocal*(params.rank  )
    end   = Nlocal*(params.rank+1)

    if(params.rank==params.size-1): end=N

    Nlocal = end - start

    return data[start:end,:]


def shuffle_catalog(data):

    np.random.seed(13579) #added seed or each process does its own permutation
                          # this was an error before and halos were double counted
    p = np.random.permutation(np.shape(data)[0])

    return data[p,:]

def cull_catalog(data):
    
    r = np.sqrt(data[:,0]**2+data[:,1]**2+data[:,2]**2)
    redshift = r2z(r)
    
    # filtering halos in the sphere of r<box/2 and z<max_redshift
    dm = (
            (redshift  > params.min_redshift ) & 
            (redshift  < params.max_redshift ) & 
            (  abs(r)  < (params.box_size)/2 ) & 
            (data[:,3] > params.min_mass     )
            )    

    data = data[dm]    

    if params.flat == 1:
        thetaxc = np.abs(np.arctan(data[:,1]/data[:,0]))*2
        thetayc = np.abs(np.arctan(data[:,2]/data[:,0]))*2	
        dm = ((thetaxc < np.radians(params.fov)) & (thetayc < np.radians(params.fov))
              & (data[:,0]>0))
        data = data[dm]
    else:        
        xcmin=-1e10; xcmax=1e10; ycmin=-1e10; ycmax=1e10; zcmin=-1e10; zcmax=1e10;

        if params.octx == 0: xcmin=-1e10; xcmax=0
        if params.octy == 0: ycmin=-1e10; ycmax=0
        if params.octz == 0: zcmin=-1e10; zcmax=0

        if params.octx == 1: xcmin=0; xcmax=1e10
        if params.octy == 1: ycmin=0; ycmax=1e10
        if params.octz == 1: zcmin=0; zcmax=1e10

        dm = ( (data[:,0] > xcmin) & (data[:,0] < xcmax) & 
               (data[:,1] > ycmin) & (data[:,1] < ycmax) & 
               (data[:,2] > zcmin) & (data[:,2] < zcmax) )
        data = data[dm]

    return data

def jiang_shmf(m,M_halo):
    gamma_1    = 0.13
    alpha_1    = -0.83
    gamma_2    = 1.33
    alpha_2    = -0.02
    beta_2     = 5.67
    zeta       = 1.19
    
    dndm = (((gamma_1*((m/M_halo)**alpha_1))+
             (gamma_2*((m/M_halo)**alpha_2)))*
             (np.exp(-(beta_2)*((m/M_halo)**zeta))))
    
    return dndm

def write_time(string_in, rank):
    if rank==0:
        fmt       = '%H:%M:%S on %m/%d/%Y'
        timestamp = datetime.datetime.now().strftime(fmt)
        bar = 72*'-'
        print ''
        print bar
        print string_in
        print 'Time:      '+timestamp
        print bar
        print ''

    return

