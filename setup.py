try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Computes the log-likelihood of a single- or multi-epoch velocity distribution including the effect of binary orbital motions',
    'author': 'Michiel Cottaar',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'MichielCottaar@gmail.com',
    'version': '0.1',
    'install_requires': ['scipy'],
    'packages': ['NAME'],
    'scripts': [],
    'name': 'velbin'
}

setup(**config)