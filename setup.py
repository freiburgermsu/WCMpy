# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'wcmpy',      
  package_dir = {'wcm':'wcmpy'},
  packages = find_packages(),
#   package_data = {
# 	'test':['databases/*'],
#   },
  version = '0.0.1',
  license = 'MIT',
  description = "Simulates a range of biochemical processes of a cocci bacterium over time.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/chemw',   
  keywords = ['chemistry', 'molecular', 'biology', 'whole-cell', 'computational', 'modeling'],
  install_requires = ['scipy', 'codons', 'pandas']
)
