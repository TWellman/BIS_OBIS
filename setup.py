# flake8: noqa

from setuptools import setup, find_packages
import pkg_resources
from codecs import open
from os import path
import json

# read package readme information
fpath = pkg_resources.resource_filename(__name__, 'README.rst')
with open(fpath, encoding='utf-8') as f:
    long_description = f.read()

# read package build information
fpath = pkg_resources.resource_filename(__name__, 'PKG_ID.json')
with open(fpath, encoding='utf-8') as f:
    _pkginfo = json.load(f)

# read Darwin Core standard vocabulary
resource_path = '/'.join(('data', 'DarwinCore_vocab_2018-03-13T17-03-40Z.json'))
fpath = pkg_resources.resource_filename(__name__, resource_path)
with open(fpath) as f:
    dc_vocab = json.load(f)


setup(
    name='obis',
    version=_pkginfo['version'],
    license='Alpha pre-release',
    description='A set of obis processing functions for Biogeographic Information System projects',
    long_description=long_description,
    url='https://www.sciencebase.gov/catalog/item/54540d80e4b0dc779374504a',
    author='T. Wellman and J. Long',
    maintainer='USGS BCB-Sciencebase Team',
    maintainer_email='bcb@usgs.gov',
    packages=find_packages(exclude=['docs', 'misc', 'tests', 'pkg_environments']),
    package_data={'data': ['data/*']},
    include_package_data=True,
    zip_safe=True,
    platforms='any',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Framework :: BIS Analytics',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Utilities',
        'Topic :: Software Development :: Libraries :: Python Modules'],
    keywords='BIS biogeography obis processor',
)
