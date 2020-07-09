import setuptools
from distutils.core import setup

try:
    import numpy as np
except:
    raise ImportError('NumPy must be installed to build GeoWombat.')


# Parse the version from the module.
# Source: https://github.com/mapbox/rasterio/blob/master/setup.py
with open('geosample/version.py') as f:

    for line in f:

        if line.find("__version__") >= 0:

            version = line.split("=")[1].strip()
            version = version.strip('"')
            version = version.strip("'")

            continue

pkg_name = 'geosample'
maintainer = 'Jordan Graesser'
maintainer_email = ''
description = 'Geo-sampling'
git_url = 'https://github.com/jgrss/geosample'
download_url = 'https://github.com/jgrss/geosample/archive/{VERSION}.tar.gz'.format(VERSION=version)
keywords = []

with open('README.md') as f:
    long_description = f.read()

with open('LICENSE.txt') as f:
    license_file = f.read()

with open('requirements.txt') as f:
    required_packages = f.readlines()


def get_packages():
    return setuptools.find_packages()


def get_package_data():
    return {'': ['*.md', '*.txt']}


def setup_package():

    include_dirs = [np.get_include()]

    metadata = dict(name=pkg_name,
                    maintainer=maintainer,
                    maintainer_email=maintainer_email,
                    description=description,
                    license=license_file,
                    version=version,
                    long_description=long_description,
                    packages=get_packages(),
                    package_data=get_package_data(),
                    zip_safe=False,
                    keywords=' '.join(keywords),
                    url=git_url,
                    download_url=download_url,
                    install_requires=required_packages,
                    include_dirs=include_dirs,
                    classifiers=['Intended Audience :: Science/Research',
                                 'License :: MIT',
                                 'Topic :: Scientific :: GIS',
                                 'Programming Language :: Python :: 3.6',
                                 'Programming Language :: Python :: 3.7',
                                 'Programming Language :: Python :: 3.8'])

    setup(**metadata)


if __name__ == '__main__':
    setup_package()
