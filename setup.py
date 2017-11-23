from setuptools import setup

setup(name='phaseimaging',
      version='0.3',
      description='Library of functions for phase imaging',
      url='http://github.com/zac-k/phaseimaging',
      author='Zachary DC Kemp',
      author_email='zachary.kemp@monash.edu',
      license='GPLv3',
      packages=['phaseimaging'],
      install_requires=['numpy',
                        'matplotlib',
                        'warnings>=0.0.dev0',
                        'copy'])