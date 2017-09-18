from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
      name='melodee',
      version='0.1',
      description='Differential Equation Compiler for HPC',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python 2',
          'Programming Language :: Python 3',
      ],
      url='https://github.com/llnl/melodee',
      author='Robert Clayton Blake III',
      author_email='blake14@llnl.gov',
      license='Apache-2.0',
      packages=['melodee'],
      install_requires=[
          'sympy',
          'ply',
          #'xml',
      ],
      zip_safe=False,
      entry_points = {
          'console_scripts' : ['cardioidGenerator=melodee.cardioidGenerator:main',
                               'continuityGenerator=melodee.continuityGenerator:main',
                               'cellmlConverter=melodee.cellmlConverter:main',
                               'matlabGenerator=melodee.matlabGenerator:main',
          ]
      },
)
