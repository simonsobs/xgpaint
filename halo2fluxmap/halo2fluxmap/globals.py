import numpy as np

                         # constants in cgs 
h     = 6.62606957e-27   # erg.s
c     = 3e+10            # cm/s
k     = 1.3806488e-16    # erg/K
Msun  = 2e33             # g
Mpc   = 3.086e24         # cm

                         # fiducial values 
Rf    = 1e3              # fiducial radius in Mpc
Rfg   = Rf*Mpc           # fiducial radius in cm
Mf    = 1                # fiducial mass value in Msun
Mfg   = Mf*Msun          # fiducial mass value in g

                         # conversion factors from 
                         # dimensionless to cgs units
ff = 1/Rfg**2/4/np.pi    # 1/cm^2
fs = Mfg                 # g



