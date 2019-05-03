import numpy as np

from numpy.random import *

N      = 10000 # number of objects
chimin = 1     # minimum distance [Mpc]
chimax = 5e3   # maximum distance [Mpc]
latmin = -5.   # minimum latitute [deg]
latmax =  5.   # maximum latitute [deg]
lngmin = -5    # minimum longitute [deg]
lngmax =  5    # maximum longitute [deg]
Mmin   = 11.   # minimum log10 halo mass [Msun]
Mmax   = 15.   # maximum log10 halo mass [Msun]

fileout = open('./catalogs/example-catalog.pksc','w')

# mean matter density
omegam = 0.3
h      = 0.7
rho = 2.775e11 * omegam * h**2

# from latitude and longitude to spherical coordinates
thetamin = 90. - latmin
thetamax = 90. - latmax
phimin   = lngmin
phimax   = lngmax

# deg to rad
thetamin = 2*np.pi * thetamin / 360. 
thetamax = 2*np.pi * thetamax / 360. 
phimin   = 2*np.pi * phimin / 360.
phimax   = 2*np.pi * phimax / 360. 

mumin = np.cos(thetamin)
mumax = np.cos(thetamax)

# generate random positions and halo masses
chi = uniform(chimin, chimax, N)
mu  = uniform(mumin,   mumax, N)
phi = uniform(phimin, phimax, N) 
M   = uniform(Mmin,     Mmax, N)
z   = chi * mu
r   = chi * np.sqrt(1.-mu**2)
x   = r   * np.cos(phi)
y   = r   * np.sin(phi)

# zeros for velocities and lagrangian positions
ph = np.zeros(N)

# convert from log10 M to M
M = 10**M

# convert from M to Lagrangian radius RTH
RTH = (3*M/4./np.pi/rho)**(1./3.)

# write data to pksc file oriented such that longitude = latitude = 0 is along z-axis
np.asarray([N]).astype(np.int32).tofile(fileout)
np.asarray([RTH.max()]).astype(np.float32).tofile(fileout)
np.asarray([0]).astype(np.float32).tofile(fileout)
np.asarray([x,y,z,RTH]).transpose().astype(np.float32).tofile(fileout)

fileout.close()


