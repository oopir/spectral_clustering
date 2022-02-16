from setuptools import setup, Extension
import numpy

setup(name="mykmeanssp",
      version="1.0",
      description="our c sp`kmeans implementation",
      ext_modules=[Extension('mykmeanssp', sources=['spkmeans.c'])],
      include_dirs=[numpy.get_include()])