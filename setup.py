# -*- coding: utf-8 -*-
"""

Tis file is part of PySBS.
    PySBS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    PySBS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    
    If you use it in research please cite an appropriate paper from README file.
    
    If you have any questions email: marcinjmalinowski@gmail.com
    
    author: Marcin Malinowski
"""

#!/usr/bin/env python

from distutils.core import setup, find_packages


install_requires=[
          "fenics-ffc>=2018.1.0",
          "numpy>=1.13.3",
          "scipy>=0.19.1",
          "matplotlib>=2.1.1"
]

setup(name='PySBS',
      version='2018.0.0',
      description='Python Stimulated Brillouin Scattering calculations',
      author='Marcin Malinowski',
      author_email='marcinjmalinowski@gmail.com',
      install_requires=install_requires,
      #dependency_links=['add github url']
      #url=' add github url',
      package_dir = {'': 'pysbs'},
      packages=find_packages(),
      
     )