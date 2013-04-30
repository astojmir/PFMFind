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
    'pfmfind.Expt',
    'pfmfind.Setup',
    'pfmfind.plugins',
    ]

libs=["gcc"]

EXTENSIONS = [
    Extension('pfmfind.Expt.FS',
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
              extra_compile_args=['-O3'],
              ),
    ]

SCRIPTS = [
    'pfmfind/scripts/FSsearchc.py',
    'pfmfind/scripts/PFMFind.pyw',
    'pfmfind/scripts/PFMFsetupdb.py',
    'pfmfind/scripts/PFMFsetupix.py',
    ]

from pfmfind import CONFIG_DATA_DIR, SQL_DATA_DIR

cfg_files = ['PFMFcf.dtd', 'PFMFdb.dtd', 'PFMFix.dtd',
             'dbsetup_sample01.xml', 'ixsetup_sample01.xml',
             ]
full_cfg_files = [join('pfmfind/data/setup_config', f) \
                  for f in cfg_files]

sql_files = ['biosqldb-pg-cnstr.sql',
             'biosqldb-pg-fk.sql',
             'biosqldb-pg-nocnstr.sql',
             'biosqldb-pg.sql',
             ]
full_sql_files = [join('pfmfind/data/sql-schema', f) \
                  for f in sql_files]

PACKAGE_DATA = {}


DATA_FILES = [
    (CONFIG_DATA_DIR, full_cfg_files),
    (SQL_DATA_DIR, full_sql_files),
    ]


SCRIPTS.append('pfmfind/scripts/FSsearchd.py')


setup(
    name='PFMFind',
    version='0.5',
    author='Aleksandar Stojmirovic',
    author_email='stojmira@ncbi.nlm.nih.gov',
    url='http://www.vuw.ac.nz/biodiscovery/publications/centre/pfmfind.aspx',
    package_dir={'pfmfind': '.'},
    packages=PACKAGES,
    package_data=PACKAGE_DATA,
    scripts=SCRIPTS,
    ext_modules=EXTENSIONS,
    data_files=DATA_FILES,
    description='System for discovery of peptide homology and function',
    long_description='Set of routines/scripts supporting similarity search'\
    ' of datasets of short protein fragments of fixed length.',
    license='GNU GPL',
    platforms='POSIX',
    classifiers = [
      'Development Status :: 3 - Alpha',
      'Environment :: X11 Applications',
      'Environment :: No Input/Output (Daemon)',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License (GPL)',
      'Operating System :: POSIX',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
    )
