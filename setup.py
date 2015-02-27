#!/home/inaki/cvs/sage/local/bin/sage-python
from distutils.core import setup, Extension

scan_module = Extension('scan',
                        include_dirs = ['/home/inaki/cvs/sage/local/include',
                                        '/Developer/SDKs/MacOSX10.6.sdk/usr/include/'],
                        sources = ['scan.c'])

setup (name = 'Scan',
       version = '1.0',
       description = 'Scan over the list of 4d reflexive polytopes',
       ext_modules = [scan_module])
