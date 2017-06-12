#! /usr/bin/python3

from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='gfa',
      version='1.3.3',
      description='Python tools for the geospatial analysis of a velocity\
field, using Davis & Titus (2011) and Bevis & Brown (2014) papers.',
      long_description=readme(),
      url='https://github.com/VicenteYanez/GFA',
      keywords='gnss gps trajectory vorticity geospatial',
      author='Vicente Yanez',
      author_email='vicenteyanez@protonmail.com',
      license='GLP-3',
      packages=find_packages(),
      install_requires=[
          'numpy', 'scipy', 'matplotlib', 'cartopy', 'click'
      ],
      entry_points={
        'console_scripts': ['gfa=gfa.scripts.gfa:cli']
        },
      zip_safe=False)
