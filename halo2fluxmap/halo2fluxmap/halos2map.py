import flux2map
import params

import numpy      as np

from   fluxmodel  import *

def halos2map(cen,ns,nsmean,sat):
	import time

	if params.flat==0: 
		fi = 12.*params.nside**2/4/np.pi               # dimensionless
	elif params.flat==1:
		fi = params.nside**2 / np.radians(params.fov)**2  # dimensionless

	t2=time.time()

	N = np.shape(sat)[0]

	# Convert to dimensionless units
	cen, sat = dimensionless(cen,sat,-1)

	C = np.zeros((np.shape(cen)[0],1),dtype='float32')
	S = np.zeros((np.shape(sat)[0],2),dtype='float32')
	
	C[:,0] = cen[:,3] # Mass
	S[:,:] = flux2map.cen2sat_masses(np.column_stack((C[:,0],ns)),ns,nsmean) # Parent Mass and Nsat

	#CALCULATE FLUX FROM CENTRALS AND SATELLITES
	if params.numdens==0:
		if (params.LM=="Planck2013"):
#			C[:,0] = LF(C[:,0],cen[:,0],cen[:,1],cen[:,2],'cen')
#			S[:,0] = LF(S[:,0],sat[:,0],sat[:,1],sat[:,2],'sat') / S[:,1]
			C[:,0] = LF(C[:,0],cen[:,0],cen[:,1],cen[:,2],'cen')
			S[:,0] = LF(S[:,0],sat[:,0],sat[:,1],sat[:,2],'cen') 
		if (params.LM=="Planck2015"):
		#divide total luminosity evenly between centrals and satellites
			C[:,0] = LF(C[:,0],cen[:,0],cen[:,1],cen[:,2],'cen') / (ns+1)
			S[:,0] = LF(S[:,0],sat[:,0],sat[:,1],sat[:,2],'sat') / (S[:,1]+1)
		
		C[:,0] = l2f(C[:,0],cen[:,0],cen[:,1],cen[:,2])
		S[:,0] = l2f(S[:,0],sat[:,0],sat[:,1],sat[:,2])

	elif params.numdens==1:
		C[:,0] = 1.0
		S[:,0] = 1.0

	cen_fluxes = C[:,0]
	sat_fluxes = S[:,0]

	total_fluxl=np.zeros(1)
	total_flux=np.zeros(1)

	if params.flat==1:				
		thetaxc = np.abs(np.arctan(cen[:,1]/cen[:,0]))*2
		thetayc = np.abs(np.arctan(cen[:,2]/cen[:,0]))*2	
		dmc = [(thetaxc < np.radians(params.fov)) & (thetayc < np.radians(params.fov))
		       & (cen[:,0]>0)]

		thetaxs = np.abs(np.arctan(sat[:,1]/sat[:,0]))*2
		thetays = np.abs(np.arctan(sat[:,2]/sat[:,0]))*2	
		dms = [(thetaxs < np.radians(params.fov)) & (thetays < np.radians(params.fov))
		       & (sat[:,0]>0)]

		total_fluxl = cen_fluxes[dmc].sum() + sat_fluxes[dms].sum() # fiducial units
	else:
		total_fluxl = cen_fluxes.sum() + sat_fluxes.sum() # fiducial units
		
	total_fluxl = total_fluxl * globals.ff*globals.fs       # erg/sec/cm^2/Hz

	#MAKE MAP
	#To save memory no longer makes T_cen,T_sat,T_tot maps
	if params.flat==0:
		omega_map = 4.*np.pi		
		T_tot  = flux2map.makemap(
			cen[:,0],cen[:,1],cen[:,2],C[:,0],params.nside)
		T_tot += flux2map.makemap(
			sat[:,0],sat[:,1],sat[:,2],S[:,0],params.nside)
	elif params.flat==1:
		omega_map = np.radians(params.fov)**2
		T_tot  = flux2map.makemapflat(
			cen[:,0],cen[:,1],cen[:,2],C[:,0],params.nside,params.fov)
		T_tot += flux2map.makemapflat(
			sat[:,0],sat[:,1],sat[:,2],S[:,0],params.nside,params.fov)

	# Convert back to physical units
	cen, sat = dimensionless(cen,sat,1)
	# Convert from dimensionless temperature to muK	
	T_tot  = T_tot.astype("float32")

	if params.numdens==0:
		T_tot *= fi*globals.ff*globals.fs

	# Sum over all processes
	T = np.zeros(np.shape(T_tot)[0],dtype='float32')

        if params.parallel:
                from mpi4py import MPI 
                params.comm.Reduce([T_tot,MPI.FLOAT],[T,MPI.FLOAT],
                                   op = MPI.SUM,root=0)
                params.comm.Reduce(total_fluxl,total_flux,op=MPI.SUM,root=0)
        else:
                T          = T_tot
                total_flux = total_fluxl

	if(params.rank==0):
		# normalize final map values to specified frequency
		if params.nu_obs_GHz  == params.norm_freq:
			params.norm_value = np.mean(T)

		T = T/params.norm_value
		T = T*params.Inu_norm

                # get total flux of map by multiplying mean intensity by map solid angle
		map_flux = (T).mean()*omega_map

                # normalize sources to have the same total flux as map
                tot_source_flux = cen_fluxes.sum() + sat_fluxes.sum()
                cen_fluxes *= map_flux / tot_source_flux
                sat_fluxes *= map_flux / tot_source_flux
                tot_source_flux = cen_fluxes.sum() + sat_fluxes.sum()

                report('',2)
                report('Finished making map nu:     '+str(params.nu_obs_GHz)+' GHz',2)
                report('     mean source intensity: '+str((cen_fluxes.sum() + sat_fluxes.sum())/omega_map),2)
                report('        mean map intensity: '+str(T.mean()),2)
                report('    mean central intensity: '+str(cen_fluxes.sum() / omega_map),2)
                report('  mean satellite intensity: '+str(sat_fluxes.sum() / omega_map),2)
                report('           min, max of map: '+str(np.min(T))+', '+str(np.max(T)),2)

		nonzeropix = np.shape(T[T>0])[0]
		
                if params.iwantcat==1:   #save centrals and satellites
			np.savez(params.inputfile+"_censat_"+str(params.nu_obs_GHz),
                                 x_cen=cen[:,0],y_cen=cen[:,1],z_cen=cen[:,2],
				 M_cen=cen[:,3],ns_cen=ns,
                                 x_sat=sat[:,0],y_sat=sat[:,1],z_sat=sat[:,2],
				 M_sat=S[:,0],cen_fluxes=cen_fluxes,sat_fluxes=sat_fluxes)

		t1=time.time()
		dt = (t1-t2)/60.
		ng = np.shape(sat)[0]+np.shape(cen)[0]
		if(params.rank==0 and params.verbose>0 and 
		   params.nu_obs_GHz == params.norm_freq):
			report('Time to project galaxies: '+str(dt)+' minutes',2)
			report('Number of nonzero pixels:   '+str(nonzeropix),-1)
		return T, cen_fluxes, sat_fluxes
	else:
		return 0, 0, 0 
