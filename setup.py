#!/usr/bin/env python
# -*- encoding: utf-8 -*-

from setuptools import find_packages
from setuptools import setup

setup(
    name='args_oap',
    license='MIT',
    description='ARGs-OAP: Online Analysis Pipeline for Anti-biotic Resistance Genes Detection from Meta-genomic Data Using an Integrated Structured ARG-database',
    author='Xiaole Yin',
    author_email='xiaole99@gmail.com',
    packages=find_packages(),
    classifiers=[
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3'],
    keywords=['antibiotics', 'antibiotic resistant genes'],
    package_data={
      'args_oap': ['*','*/*'],
    },
    entry_points={
        'console_scripts': [
            'args_oap = args_oap.args_oap:main',
        ]
    }
)