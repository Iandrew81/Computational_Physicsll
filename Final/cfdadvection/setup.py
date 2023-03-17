from setuptools import setup, find_packages

setup(name='cfdadvection',
      description='Module to solve 1D advection problems using finite volume approach and different slope limiters',
      url='https://github.com/Iandrew81/MACV/tree/main/Final/cfdadvection',
      author='Andres Villares, Mateo Carop',
      author_email='andres.villares@yachaytech.edu.ec, mateo.carpio@yachaytech.edu.ec',
      packages=find_packages(),
      install_requires=['numpy', 'matplotlib'])
