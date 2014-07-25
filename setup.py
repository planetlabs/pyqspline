from distutils.core import setup, Extension

qspline = Extension('qspline',
                    sources = ['qspline.c'])

setup(name = 'pyqspline',
      version = '1.0',
      description = 'Python version of qspline, which produces a quaternion spline interpolation of sparse data.',
      ext_modules = [qspline])
