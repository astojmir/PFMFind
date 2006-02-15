"""Distutils based setup script for ShortFrags.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install
"""

from distutils.core import setup, Extension
import sys
from os.path import join

PACKAGES = [
    'ShortFrags',
    'ShortFrags.GUI',
    'ShortFrags.Expt',
    'ShortFrags.Setup',
    'ShortFrags.plugins',
    ]

libs=["gcc"]
if sys.platform == 'win32':
    libs.append("gw32c")

EXTENSIONS = [
    Extension('ShortFrags.Expt.FS',
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
                             ('THREADS', 1)
                             ],
              libraries=libs,
              extra_compile_args=['-O3'],
              ),
    ]

SCRIPTS = [
    'ShortFrags/scripts/FSsearchd.py',
    'ShortFrags/scripts/PFMFind.py',
    'ShortFrags/scripts/PFMFsetupdb.py',
    'ShortFrags/scripts/PFMFsetupix.py',
    'ShortFrags/scripts/killslaves.py',
    'ShortFrags/scripts/searchclient.py',
    ]

from ShortFrags import CONFIG_DATA_DIR, SQL_DATA_DIR

cfg_files = ['PFMFcf.dtd', 'PFMFdb.dtd', 'PFMFix.dtd',
             'dbsetup_sample01.xml', 'ixsetup_sample01.xml',
             ]
full_cfg_files = [join('ShortFrags/data/setup_config', f) \
                  for f in cfg_files] 

sql_files = ['biosqldb-pg-cnstr.sql',
             'biosqldb-pg-fk.sql',
             'biosqldb-pg-ix.sql',
             'biosqldb-pg-nocnstr.sql',
             'biosqldb-pg-nofk.sql',
             'biosqldb-pg-plain.sql',
             'biosqldb-pg.sql',
             ]
full_sql_files = [join('ShortFrags/data/sql-schema', f) \
                  for f in sql_files] 

DATA_FILES = [
    (CONFIG_DATA_DIR, full_cfg_files),
    (SQL_DATA_DIR, full_sql_files),
    ]


setup(
    name='PFMFind',
    version='0.46.1',
    author='Aleksandar Stojmirovic',
    author_email='astojmir@uottawa.ca',
    url='http://aix1.uottawa.ca/~astojmir',
    packages=PACKAGES,
    scripts=SCRIPTS,
    ext_modules=EXTENSIONS,
    data_files=DATA_FILES,
    summary='System for discovery of peptide homology and function',
    description='Set of routines/scripts supporting similarity search'\
    ' of datasets of short protein fragments of fixed length.',
    license='GNU GPL',
    platform='POSIX, Windows',
    classifiers = [
      'Development Status :: 3 - Alpha',
      'Environment :: X11 Applications',
      'Environment :: No Input/Output (Daemon)',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License (GPL)',
      'Operating System :: POSIX',
      'Operating System :: Microsoft :: Windows',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
    )
