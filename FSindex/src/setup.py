"""Distutils based setup script for ShortFrags.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install
"""

from distutils.core import setup, Extension
import sys

PACKAGES = [
    'ShortFrags',
    'ShortFrags.GUI',
    'ShortFrags.Expt',
    ]

libs=["gcc"]
if sys.platform == 'win32':
    libs.append("gw32c")

EXTENSIONS = [
    Extension('ShortFrags._FS',
              ['swig/FS_wrap.c',
               #'ShortFrags/FS.i',
               'lib/bioseq.c',
               'lib/FSindex.c',
               'lib/avl.c',
               'lib/fastadb.c',
               'lib/hit_list.c',
               'lib/misclib.c',
               'lib/smatrix.c',
               'lib/partition.c',
               'sarray/lcp.c',
               'sarray/sarray.c',
               'sarray/scode.c',
               'sarray/ssarray.c',
               ],
              include_dirs=['include'],
              define_macros=[('GCC_INLINE', None),
                             ('THREADS', 4)
                             ],
              libraries=libs,
              extra_compile_args=['-O3'],
              ),
    ]

SCRIPTS = [
    'ShortFrags/frag_toolbox.pyw'
    ]

setup(
    name='ShortFrags',
    version='0.3',
    author='Aleksandar Stojmirovic',
    author_email='aleksand@mcs.vuw.ac.nz',
    url='http://www.mcs.vuw.ac.nz/~aleksand',
    packages=PACKAGES,
    scripts=SCRIPTS,
    ext_modules=EXTENSIONS,
#    data_files=DATA_FILES,
    )
