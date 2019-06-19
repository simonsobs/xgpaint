import matplotlib.pyplot as plt
import numpy as np
from   scipy.ndimage.filters import gaussian_filter
import sys
mapfile = 'maps/example_ns4096_nu545.map'

npix = np.fromfile(mapfile,count=2,dtype=np.int32)
fov  = np.fromfile(mapfile,count=2,dtype=np.float32)
flux = np.fromfile(mapfile,count=npix[0]*npix[1],dtype=np.float32)

spix = 2

flux = np.reshape(flux,(npix[0],npix[1]))
flux = gaussian_filter(flux,spix)
fmin=1e-1
fmax=1e1
vmin=np.log10(fmin*flux.mean())
vmax=np.log10(fmax*flux.mean())
flux[flux==0]=fmin*flux.mean()

plt.imshow(np.log10(flux),vmin=vmin,cmap='jet')

plt.show()



