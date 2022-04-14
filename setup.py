from setuptools import setup, Extension, find_packages
import numpy

setup(name="mykmeanssp",
      version="1.0",
      description="our c spkmeans implementation",
      packages=find_packages(),
      ext_modules=[Extension('mykmeanssp', sources=['spkmeansmodule.c', 'spkmeans.c'])])