#!/usr/bin/env python
import sys
import numpy             as np
import healpy            as hp
import matplotlib.pyplot as plt
import halo2fluxmap      as h2fm

import mpi4py.rc
import datetime
import os
import psutil

from   scipy.ndimage.filters import gaussian_filter

h2fm.init.initialize('param/param_II.txt')

if len(sys.argv) > 1:
    h2fm.params.octx = int(sys.argv[1])
    h2fm.params.octy = int(sys.argv[2])
    h2fm.params.octz = int(sys.argv[3])
    h2fm.params.outbase = h2fm.params.outbase+'-'+sys.argv[1]+sys.argv[2]+sys.argv[3]
    print h2fm.params.octx,h2fm.params.octy,h2fm.params.octz,h2fm.params.outbase

h2fm.utils.report('Reading, Shuffling, Culling Catalogue',2)

# Get catalog of centrals
cen = h2fm.io.read_catalog()

# Get number of satellites for each halo
ns, nsmean  = h2fm.hod.cen2ns(cen)

# Populate halos with satellites
sat = h2fm.hod.populate_centrals(cen,ns)

# Write time
h2fm.utils.write_time("HOD completed", h2fm.params.rank)

#loop over frequencies
for inu_map in range(len(h2fm.params.freq_list)):

    # Put halos in map
    h2fm.params.nu_obs     = h2fm.params.freq_list[inu_map] * 1.e9
    h2fm.params.nu_obs_GHz = h2fm.params.freq_list[inu_map]
        
     # Give flux to gals and put in map
#    intensity, flux_cen, flux_sat = h2fm.halos2map.halos2map(cen,ns,nsmean,sat) 
 
    # Give flux to gals and put in map
    #    intensity, flux_cen, flux_sat = h2fm.halos2map.halos2map(cen,ns,nsmean,sat)
    pcen,mcen,fcen,nsat,psat,msat,fsat,lcen,lsat = h2fm.halos2sources.halos2sources(cen,ns,nsmean,sat)
    intensity, flux_cen, flux_sat                = h2fm.sources2map.sources2map(pcen,mcen,fcen,nsat,
                                                                                psat,msat,fsat,lcen,lsat)       

    # Write map
    ns_str   = str(h2fm.params.nside)
    nu_str   = str(int(h2fm.params.freq_list[inu_map])).zfill(4)
    base = 'maps/'+h2fm.params.outbase 
    base += ('_ns'+ns_str+'_nu'+nu_str)

    h2fm.io.writemap(base,intensity)

# Write time
h2fm.utils.write_time("Halo projection completed",h2fm.params.rank)

