import flux2map
import params

import numpy      as np

from   fluxmodel  import *

def halos2sources(cen,nsat,nsatmean,sat):
	import time

	if params.flat==0: 
		fi = 12.*params.nside**2/4/np.pi               # dimensionless
	elif params.flat==1:
		fi = params.nside**2 / np.radians(params.fov)**2  # dimensionless

	t2=time.time()

	N = np.shape(sat)[0]

	# Convert to dimensionless units
	cen, sat = dimensionless(cen,sat,-1)

	mcen = cen[:,3] # Mass
	msat = flux2map.cen2sat_masses(np.column_stack((mcen,nsat)),nsat,nsatmean) # Parent Mass and Nsat

	#CALCULATE FLUX FROM CENTRALS AND SATELLITES
	if params.numdens==0:
		if (params.LM=="Planck2013"):
			lcen = LF(mcen,cen[:,0],cen[:,1],cen[:,2],'cen')
			lsat = LF(msat,sat[:,0],sat[:,1],sat[:,2],'cen') 
		if (params.LM=="Planck2015"):
		#divide total luminosity evenly between centrals and satellites
			lcen = LF(mcen,cen[:,0],cen[:,1],cen[:,2],'cen') / (nsat+1)
			lsat = LF(msat,sat[:,0],sat[:,1],sat[:,2],'sat') / (S[:,1]+1)
		
		fcen = l2f(lcen,cen[:,0],cen[:,1],cen[:,2])
		fsat = l2f(lsat,sat[:,0],sat[:,1],sat[:,2])

	elif params.numdens==1:
		fcen[:] = 1.0
		fsat[:] = 1.0

        pcen = cen[:,0:3]
        psat = sat[:,0:3]

        t1=time.time()
        dt = (t1-t2)/60.
        ng = np.shape(sat)[0]+np.shape(cen)[0]
        if(params.rank==0 and params.verbose>0 and
           params.nu_obs_GHz == params.norm_freq):
                report('Time to assign fluxes: '+str(dt)+' minutes',2)
        
        return (pcen,mcen,fcen,nsat,
                psat,msat,fsat,
                lcen,lsat)
