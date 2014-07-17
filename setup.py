import numpy
from distutils.core import setup, Extension

qspline = Extension('qspline',
                    sources = ['qspline.c'])

setup(name = 'pyqspline',
      version = '1.0',
      description = 'Python version of the qspline C program',
      include_dirs = [numpy.get_include()],
      ext_modules = [qspline])
