# List of parameters and their defaults
folder       = './'
inputfile    = 'input_halos.pksc'
outbase      = 'output_map'

verbose      = 2
scramble     = False
LightIO      = True
Nreads       = 1
format       = 'pksc'
octx         = 2
octy         = 2
octz         = 2

justify      = 25*' '
serial       = 0

flat         = 1            #0 for all sky, 1 for flatsky
freq_list    = ([545., 100., 143., 217., 353., 857.]) #all maps will be normalized w.r.t to the mean of the first map          
Inu_norm     = 0.0299988959602  # z<1.25 and z< 4.5 for M > 2.62e13Msun = 0.019840695547 0.0299988959602 
nside        = 1024         #nside for healpix, npix for flatsky
fov	     = 10
numdens      = 0            #0 for cib, 1 for optical 

hod          = 'shang'      #"shang" for cib, "manera" for optical 
LM           = 'Planck2013' #'Planck2015'

min_redshift = 0.0
max_redshift = 4.5
num_redshift = 1       # If memory limited, splits the map into redshift slices 
min_mass     = 2.62e13
box_size     = 40000

iwantcat     = 0     #save catalogue of satellites and centrals 
H0           = 70
omegab       = 0.043
omegac       = 0.207
omegal       = 0.75 
scalar_index = 0.96
sigma8       = 0.8   

shang_Td     = 24.4          #Planck 2013 values
shang_alpha  = 0.36
shang_beta   = 1.75
shang_eta    = 3.2           #same as delta_CIB in paper
shang_I0     = 46

#shang_Td     = 20           # old values
#shang_beta   = 1.4
#shang_eta    = 2.0
#shang_I0     = 46
shang_Mmin   = 1e10
shang_Msmin  = 1e11
shang_zplat  = 20.0
shang_Mpeak  = 10.**12.2
shang_sigmaM = 0.4
shang_nc2ns  = 10

#Satellite mass function parameters
M1           = 1e+14
M0           = 11939881044642.73
sigma_logM   = 0.569
logM_min     = 13.09
alpha        = 1.0127
mass_type    = 'nonlog'
