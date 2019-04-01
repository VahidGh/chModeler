#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='chModeler',
    version='0.9.0',
    python_requires='>=3.2',
    packages=find_packages(),
    description='Modeling ion channels from patch clamp data',
    author='Vahid Ghayoomie',
    author_email='vahidghayoomi@gmail.com',
    url='https://github.com/VahidGh/chModeler',
    install_requires=['argparse', 'matplotlib', 'numpy', 'scipy', 'requests'],
)