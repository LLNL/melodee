from setuptools import setup

setup(name='melodee',
      version='0.1',
      description='Differential Equation Compiler for HPC',
      url='http://github.com/llnl/melodee',
      author='Robert Clayton Blake III',
      author_email='blake14@llnl.gov',
      license='Apache-2.0',
      packages=['melodee'],
      install_requires=[
          'sympy',
          'ply'
      ],
      zip_safe=False)
