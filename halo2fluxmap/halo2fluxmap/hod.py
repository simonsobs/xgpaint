import time
import scipy

from utils   import *

def hod(M): 

    if (params.hod == "manera"):
        return hod_manera(M)
    elif (params.hod == "shang"):
        return hod_shang(M)
    else:
        print "hod " + params.hod + " not found, exiting ..."
        exit()

def hod_manera(M):

    if params.mass_type=='log':
        M = 10**M
    f=(M-params.logM_min)/(params.sigma_logM)  #function inside erf
    N_cen = 0.5*(1+special.erf(f))      #<N_cen>
    N_sat_list=[]

    for i in range(len(M)):              #sets N_sat to zero if M halo<M0
        if M[i]<params.M0:
            N_sat=0
        else:
            N_sat = N_cen[i] * (((M[i]-params.M0)/(params.M1))**params.alpha)
        N_sat_list.append(N_sat)

    N_sat_list = np.array(N_sat_list)	

    N_gal = N_cen + N_sat_list          #weights

    return N_sat_list,N_cen,N_gal

def hod_shang(M):

    if params.mass_type=='log':
        M = 10**M

    N_cen = M * 0 + 1

    #<n_sat> interpolation values
    def integrand_m(lm,lM_halo):
        m=np.exp(lm)
        M_halo=np.exp(lM_halo)

        dns_dm = jiang_shmf(m,M_halo)

        return dns_dm

    x_m = np.linspace(np.log(M.min()),np.log(M.max()),1000)
    N_sat_i = np.zeros(len(x_m))
    for i in range(len(x_m)):
        N_sat_i[i], foo = scipy.integrate.quad(integrand_m,np.log(params.shang_Msmin),
                       x_m[i],args=x_m[i])
        N_sat_i[i] = max(0,N_sat_i[i])

    f_m = interpolate.interp1d(x_m,N_sat_i)
    N_sat = f_m(np.log(M)) 

    N_gal = N_cen + N_sat

    return N_sat,N_cen,N_gal

def fnfw(x):
    fnfw = np.log(1+x)-x/(1+x)
    return fnfw


def g(x,c,m):
    '''
    g-f(x)=0
    '''
    g = m - (fnfw(x)/fnfw(c))
    return g


def root_finder(c,lnm):
    '''
    finding root of a given function
    '''
    m=np.exp(lnm)
    root = scipy.optimize.fsolve(g,c*0+1,args=(c,m))/c
    return root

def cen2ns(cen):
    nsmean = hod(cen[:,3])[0]
    return np.random.poisson(nsmean).astype('int32'),nsmean

def clnm2r(c,lnm):

    mmin = -7.1
    mmax = 0
    cmin = np.min(c*0.999)
    cmax = np.max(c*1.001)

    ctab  = np.linspace(cmin,cmax,30)
    mtab  = np.linspace(mmin,mmax,30)
    cc,mm = np.meshgrid(ctab,mtab)

    rr    = np.zeros((len(ctab),len(mtab)))
    rr    = np.reshape(rr,(len(ctab)*len(mtab)))
    cc    = np.reshape(cc,(len(ctab)*len(mtab)))
    mm    = np.reshape(mm,(len(ctab)*len(mtab)))
    rr    = root_finder(cc,mm)

    return scipy.interpolate.griddata((cc,mm),rr,(c,lnm))

def populate_centrals(cen,ns):

    from mpi4py import MPI
    import flux2map

    t1=time.time()

    N   = np.array(0,dtype='int64')
    Nt  = np.array(0,dtype='int64')
    N  += np.sum(ns) # N = total number of satellites

    if params.parallel:
        params.comm.Reduce([N,MPI.INT],[Nt,MPI.INT],op = MPI.SUM,root=0)
    else:
        Nt = N

    report('Number of satellites: '+str(Nt),-1)

    report('Populating...',2)

    if(params.rank==0 and params.verbose>0): 
        print ''
        print params.justify,'Centrals:  ',np.shape(cen)[0]
        print params.justify,'Satellites:',N

    sat = flux2map.cen2sat(cen,ns)
    mhsat = sat[:,3].copy() # host halo mass of sattelite

    # satarr_1 = concentration of satellite
    # satarr_2 = log satellite inner mass in units of host mass

    satarr_1  = mz2c(sat[:,3],
                     r2z((sat[:,0]**2+sat[:,1]**2+sat[:,2]**2)**0.5))
    satarr_2  = np.log(np.random.uniform(size=N).astype('float32'))

    satarr_2[satarr_2<-7] = -7

    sat[:,3] = m2r(sat[:,3]) * clnm2r(satarr_1,satarr_2) 
    sat[:,3] = sat[:,3] * 200**(-1./3) #Lagrangian to Eulerian

    # satarr_1 = phi of satellite
    # satarr_2 = theta of satellite

    satarr_1  = random_phi(N)
    satarr_2  = random_theta(N)

    sat[:,0] += sat[:,3]*(np.sin(satarr_2))*(np.cos(satarr_1))
    sat[:,1] += sat[:,3]*(np.sin(satarr_2))*(np.sin(satarr_1))
    sat[:,2] += sat[:,3]*(np.cos(satarr_2))

    Ne = (np.prod(np.shape(sat)) + 
          np.prod(np.shape(cen)) + 
          np.prod(np.shape(ns)))
    Ne += 2*N

    dt = (time.time()-t1)/60.
    if(params.rank==0 and params.verbose>0): 
        report('Time to populate halos:  '+str(dt)+'minutes',2)

    sat[:,3] = mhsat
    
    return sat

def make_muofn():
    
    def integrand(lmu):
        mu=np.exp(lmu)
        dns_dm = jiang_shmf(mu,1.0)
        return dns_dm

    mu = np.logspace(np.log10(1e-6),0,1000)
    n  = np.zeros(len(mu))
    for i in range(len(mu)):
        n[i], foo = scipy.integrate.quad(integrand,np.log(mu[i]),0.)

    muofn = interpolate.interp1d(n,mu)
    return muofn

