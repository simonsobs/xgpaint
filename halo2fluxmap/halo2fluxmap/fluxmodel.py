import flux2map
import globals

from utils   import *

def printbounds(a):
	print a.min(),a.max()

def sigma_cen(m):
    
    arg  = np.log10(m)
    arg -= np.log10(params.shang_Mpeak)
    arg  = arg**2
    arg /= 2*params.shang_sigmaM
    arg  = np.exp(-arg)
    arg *= m
    arg *= 1/np.sqrt(2*np.pi*params.shang_sigmaM)

    return arg

def sigma_sat(m):

    # Make table of <L_sat> vs. M_halo
    x = np.linspace(np.log(params.shang_Mmin),np.log(m.max()),100)
    L_mean     = np.zeros(len(x))
    for i in range(len(x)):
        L_mean[i],err = quad(integrand_L,np.log(params.shang_Mmin),
                     x[i],args=x[i])
    f_L = interpolate.interp1d(x,L_mean)

    return f_L(np.log(m))


def broken_power_law_sehgal(m_in):
    Lb = params.sehgal_L_b
    m = params.sehgal_m
    n = params.sehgal_n
    p0 = (1+m) * (1+n) / Lb / (m-n)
    
    def F_inv_lower(F_in):
        return (Lb**n * (1.0+n) * F_in / p0)**(1.0 / (1.0+n))
    def F_inv_upper(F_in):
        return (- Lb**m * (1+m) * (Lb * p0 * ( 1/(1.0+n) - 1/(1.0+m) ) - F_in) / p0 )**(1.0 / (1.0+m))
    Fb = Lb * p0 / (1+n)
    def F_inv(F_in):
        return np.piecewise(F_in, [F_in < Fb, F_in >= Fb], [F_inv_lower, F_inv_upper])
    
    unif_x = np.random.uniform(size=m_in.shape[0])
    return F_inv(unif_x)



def jiang_shmf(m,M_halo):	
    gamma_1    = 0.13
    alpha_1    = -0.83
    gamma_2    = 1.33
    alpha_2    = -0.02
    beta_2     = 5.67
    zeta       = 1.19 	

    dndm = (((gamma_1 * ((m/M_halo)**alpha_1)) +
            (gamma_2 * ((m/M_halo)**alpha_2))) *
            (np.exp(-(beta_2) * ((m / M_halo)**zeta))))

    return dndm

# <L_sat> interpolation values
def integrand_L(lm,lM_halo):
    m      = np.exp(lm)
    M_halo = np.exp(lM_halo)

    dns_dm = jiang_shmf(m,M_halo)
    dns_dm_sigma = sigma_cen(m) * dns_dm

    return dns_dm_sigma

def f2t(intensity):
    x = h*params.nu_obs/k/2.726
    T = ((np.exp(x)-1)**2/np.exp(x)*intensity).astype('float32')

    return T.astype('float32')

def nu2theta(nu):
    xnu     = globals.h*nu/globals.k/params.shang_Td
    Thetanu = (1/(nu*params.shang_I0) *
               xnu**(4.+params.shang_beta)/(np.exp(xnu)-1.))
    return Thetanu

def sehgal_B(n_obj):
    costh = np.random.uniform(size=n_obj)
    beta = np.sqrt(1 - 1./(params.sehgal_gamma**2))
    return ((1 - beta*costh)**(-2.0) + (1 + beta*costh)**(-2.0)) / 2.0

def sehgal_sed(n_obj):
    a0, (a1l, a1r), (a2l, a2r), (a3l ,a3r) = params.sehgal_a_coeff
    print(a0, (a1l, a1r), (a2l, a2r), (a3l ,a3r))
    a1 = np.random.uniform(low=a1l, high=a1r, size=n_obj)
    a2 = np.random.uniform(low=a2l, high=a2r, size=n_obj)
    a3 = np.random.uniform(low=a3l, high=a3r, size=n_obj)
    lnu = np.log10(params.nu_obs/1e9)
    sed = 10**( a0 + a1*lnu + a2*lnu**2 + a3*lnu**3 )
    return sed
    
def LF(M,x,y,z,gtype):

    r  = x**2 
    r += y**2
    r += z**2
    r  = r**0.5
    r *= globals.Rf

    z  = r2z(r)
    r  = (1+z)*params.nu_obs      #observed frequency in Hz

    if (params.LM=="Planck2013"):
        if(gtype == 'cen'): L = sigma_cen(M)
        if(gtype == 'sat'): L = sigma_sat(M)
    elif (params.LM=="Planck2015"):              #const * M_500/1e14M_sun
        L = M*np.sqrt(200./500) / 1.e14 #sqrt(200/500) to convert M_200 to M_500
    elif (params.LM=="sehgal_radio"):
        # just return, do not do Shang stuff
        if(gtype == 'cen'): return M * 0.0
        if(gtype == 'sat'):  
            L = (broken_power_law_sehgal(M) * 
                 (1 * sehgal_sed(M.shape[0]) +
                  params.sehgal_R_int * sehgal_B(M.shape[0]) * (params.nu_obs/145e9)**(-0.8))
                 / (1+params.sehgal_R_int))
            return L * (1+z) 

    L *= nu2theta(r)
    L *= (1+z)**params.shang_eta

    return L

def l2f(L,x,y,z):

    r  = x**2 
    r += y**2
    r += z**2
    r  = np.sqrt(np.abs(r))
    r *= globals.Rf
    return L / (r**2 / (1+r2z(r)) * globals.Rf**2)

def dimensionless(cen,sat,i):

    # cen[:,3]   *= globals.Mf**i

    cen[:,0:3] *= globals.Rf**i
    sat[:,0:3] *= globals.Rf**i

    return cen, sat
