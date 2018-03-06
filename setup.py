try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'Andre_Macedo',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'andre.lopes.macedo@gmail.com',
    'version': '0.1pre',
    'install_requires': ['nose'],
    'packages': ['NAME'],
    'scripts': [],
    'name': 'projectname'
}

setup(**config)
