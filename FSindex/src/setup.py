"""Distutils based setup script for pfmfind.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install
"""

from distutils.core import setup, Extension
import sys
from os.path import join

PACKAGES = [
    'pfmfind',
    'pfmfind.GUI',
    'pfmfind.plugins',
    'pfmfind.search',
    'pfmfind.setup',
    ]

libs=["gcc"]

EXTENSIONS = [
    Extension('pfmfind.search.FS',
              ['swig/FS_wrap.c',
               'lib/bioseq.c',
               'lib/FSindex.c',
               'lib/fastadb.c',
               'lib/hit_list.c',
               'lib/misclib.c',
               'lib/smatrix.c',
               'lib/partition.c',
               ],
              include_dirs=['include'],
              define_macros=[('GCC_INLINE', None),
                             ('THREADS', 2)
                             ],
              libraries=libs,
              ),
    ]

SCRIPTS = [
    # 'pfmfind/scripts/FSsearchs.py',
    'pfmfind/scripts/FSsearchc.py',
    'pfmfind/scripts/FSsearchd.py',
    'pfmfind/scripts/PFMFind.pyw',
    'pfmfind/scripts/PFMFsetupdb.py',
    'pfmfind/scripts/PFMFsetupix.py',
    ]

PACKAGE_DATA = {'pfmfind.setup': ['sql-schema/*.sql',
                                  'sql-schema/README.txt',
                                  ]}

setup(
    name='PFMFind',
    version='0.5',
    author='Aleksandar Stojmirovic',
    author_email='aleksandar@stojmirovic.org',
    url='http://www.vuw.ac.nz/biodiscovery/publications/centre/pfmfind.aspx',
    package_dir={'pfmfind': 'pfmfind'},
    packages=PACKAGES,
    package_data=PACKAGE_DATA,
    scripts=SCRIPTS,
    ext_modules=EXTENSIONS,
    description='System for discovery of peptide homology and function',
    long_description='Set of routines/scripts supporting similarity search'\
    ' of datasets of short protein fragments of fixed length.',
    license='GNU GPL',
    platforms='POSIX',
    classifiers = [
      'Development Status :: 4 - Beta',
      'Environment :: X11 Applications',
      'Environment :: No Input/Output (Daemon)',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License (GPL)',
      'Operating System :: POSIX',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
    )
