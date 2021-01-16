from distutils.core import setup, Extension

setup(name = 'pyqspline',
      version = '0.2',
      author = 'James McEnnan, Planet Labs',
      description = 'Python version of qspline, which produces a quaternion spline interpolation of sparse data.',
      py_modules = ['pyqspline'],
      ext_modules = [Extension('qspline', sources=['qspline.c'])])
