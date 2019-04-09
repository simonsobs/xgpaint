import flux2map
import params

import numpy      as np

from   fluxmodel  import *

def sources2map(pcen,mcen,fcen,nsat,
                psat,msat,fsat,
                lcen,lsat):
	import time

	if params.flat==0:
		fi = 12.*params.nside**2/4/np.pi               # dimensionless
	elif params.flat==1:
		fi = params.nside**2 / np.radians(params.fov)**2  # dimensionless

	t2=time.time()

	#MAKE MAP
	#To save memory no longer makes T_cen,T_sat,T_tot maps
        xcen=pcen[:,0]
        ycen=pcen[:,1]
        zcen=pcen[:,2]
        xsat=psat[:,0]
        ysat=psat[:,1]
        zsat=psat[:,2]
	if params.flat==0:
		omega_map = 4.*np.pi		
		I_mapl  = flux2map.makemap(xcen,ycen,zcen,fcen,params.nside)
		I_mapl += flux2map.makemap(xsat,ysat,zsat,fsat,params.nside)
	elif params.flat==1:
		omega_map = np.radians(params.fov)**2
		I_mapl  = flux2map.makemapflat(xcen,ycen,zcen,fcen,params.nside,params.fov)
		I_mapl += flux2map.makemapflat(xsat,ysat,zsat,fsat,params.nside,params.fov)

	# Convert back to physical units
	pcen, psat = dimensionless(pcen,psat,1)
        I_mapl     = I_mapl.astype("float32")

	if params.numdens==0:
		I_mapl *= fi * globals.ff * globals.fs

        if params.flat==1:
                thetaxc = np.abs(np.arctan(pcen[:,0]/pcen[:,2]))*2
                thetayc = np.abs(np.arctan(pcen[:,1]/pcen[:,2]))*2
                dmc = [(thetaxc < np.radians(params.fov)) & (thetayc < np.radians(params.fov))
                       & (pcen[:,2]>0)]

                thetaxs = np.abs(np.arctan(psat[:,0]/psat[:,2]))*2
                thetays = np.abs(np.arctan(psat[:,1]/psat[:,2]))*2
                dms = [(thetaxs < np.radians(params.fov)) & (thetays < np.radians(params.fov))
                       & (psat[:,2]>0)]

                tot_src_fluxl = fcen[dmc].sum() + fsat[dms].sum() # fiducial units

        else:

                tot_src_fluxl = fcen.sum()      + fsat.sum()      # fiducial units

	# Sum over all processes
        tot_src_flux = 0.0
	I_map        = np.zeros(np.shape(I_mapl)[0],       dtype='float32')

        if params.parallel:

                from mpi4py import MPI 

                tot_src_flux_buff  = np.array(0,            dtype='float32')
                tot_src_flux_buffl = np.array(tot_src_fluxl,dtype='float32')

                params.comm.Reduce([I_mapl,             MPI.FLOAT],[I_map,             MPI.FLOAT],
                                   op = MPI.SUM, root = 0)
                params.comm.Reduce([tot_src_flux_buffl, MPI.FLOAT],[tot_src_flux_buff, MPI.FLOAT],
                                   op = MPI.SUM, root = 0)

                tot_src_flux       = tot_src_flux_buff

        else:

                I_map        = I_mapl
                tot_src_flux = tot_src_fluxl

	if(params.rank==0):
		# normalize final map values to specified frequency
		if params.nu_obs_GHz  == params.norm_freq:
			params.norm_value = np.mean(I_map)

                I_map = I_map / params.norm_value
		I_map = I_map * params.Inu_norm

                # get total flux of map by multiplying mean intensity by map solid angle
		map_flux = (I_map).mean()*omega_map

                # normalize sources to have the same total flux as map
                fcen *= map_flux / tot_src_flux
                fsat *= map_flux / tot_src_flux
                tot_src_flux = fcen.sum() + fsat.sum()

                report('',2)
                report('Finished making map nu:     '+str(params.nu_obs_GHz)+' GHz',2)
                report('     mean source intensity: '+str((tot_src_flux)/omega_map),2)
                report('        mean map intensity: '+str(I_map.mean()),2)
                report('    mean central intensity: '+str(fcen.sum() / omega_map),2)
                report('  mean satellite intensity: '+str(fsat.sum() / omega_map),2)
                report('           min, max of map: '+str(np.min(I_map))+', '+str(np.max(I_map)),2)

		nonzeropix = np.shape(I_map[I_map>0])[0]
		
                if params.iwantcat==1:   #save centrals and satellites
			np.savez(params.inputfile+"_censat_"+str(params.nu_obs_GHz),
                                 xcen=pcen[:,0],ycen=pcen[:,1],zcen=pcen[:,2],
				 mcen=mcen,fcen=fcen,nsat=nsat,
                                 xsat=psat[:,0],ysat=psat[:,1],zsat=psat[:,2],
				 msat=msat,fsat=fsat,lcen=lcen,lsat=lsat)

		t1=time.time()
		dt = (t1-t2)/60.
		ng = np.shape(psat)[0]+np.shape(pcen)[0]
		if(params.rank==0 and params.verbose>0 and 
		   params.nu_obs_GHz == params.norm_freq):
			report('Time to project galaxies: '+str(dt)+' minutes',2)
			report('Number of nonzero pixels:   '+str(nonzeropix),-1)
		return I_map, fcen, fsat
	else:
		return 0, 0, 0 
