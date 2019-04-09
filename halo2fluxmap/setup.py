from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension("halo2fluxmap.flux2map"
              ,["halo2fluxmap/flux2map.pyx"],
              include_dirs=[np.get_include()]
    ),]

setup(name='halo2fluxmap',
      ext_modules=cythonize(ext_modules),
      include_dirs=[np.get_include()],
      version='0.1',
      description='for converting halo catalogs to maps for a given hod',
      url='http://github.com/marcelo-alvarez/halo2fluxmap',
      author='Marcelo Alvarez',
      author_email='marcelo.alvarez@berkeley.edu',
      license='GNU General Public License v3.0',
      packages=['halo2fluxmap'],
      zip_safe=False)

