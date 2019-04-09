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
	Thetanu = (1/(nu*params.shang_I0)*
		   xnu**(4.+params.shang_beta)/(np.exp(xnu)-1.))
	return Thetanu

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
	if (params.LM=="Planck2015"):              #const * M_500/1e14M_sun
		L = M*np.sqrt(200./500) / 1.e14 #sqrt(200/500) to convert M_200 to M_500

	L *= nu2theta(r)


	L *= (1+z)**params.shang_eta

	return L

def l2f(L,x,y,z):	

	r  = x**2 
	r += y**2
	r += z**2
	r  = r**0.5
	r *= globals.Rf	

	return L / r**2 / (1+r2z(r)) * globals.Rf**2

def dimensionless(cen,sat,i):

#	cen[:,3]   *= globals.Mf**i

	cen[:,0:3] *= globals.Rf**i
	sat[:,0:3] *= globals.Rf**i
	
	return cen, sat
