try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Computes the log-likelihood of a single- or multi-epoch velocity distribution for an intrinsic Gaussian velocity distribution including the effect of binary orbital motions and Gaussian measurement uncertainties',
    'author': 'Michiel Cottaar',
    'url': 'https://github.com/MichielCottaar/velbin',
    'download_url': 'https://github.com/MichielCottaar/velbin',
    'author_email': 'MichielCottaar@gmail.com',
    'version': '0.1',
    'install_requires': ['scipy'],
    'packages': ['velbin'],
    'scripts': [],
    'name': 'velbin'
}

setup(**config)