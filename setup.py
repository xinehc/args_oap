#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import os
import shutil
import subprocess
import re

from setuptools import find_packages
from setuptools import setup

setup(
    name='args_oap',
    license='MIT',
    description='test version',
    author='',
    author_email='',
    url='',
    packages=find_packages(),
    classifiers=[
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3'],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    package_data={
      'args_oap': ['*','*/*','bin/*/*', 'example/*/*'],
    },
    # install_requires=requirements,
    entry_points={
        'console_scripts': [
            'args_oap = args_oap.args_oap:main',
        ]
    }
)
