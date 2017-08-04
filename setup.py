import os

from setuptools import setup
from pip.req import parse_requirements

REQUIREMENTS_TXT = os.path.join(os.path.dirname(__file__), "requirements.txt")

setup(
    name='esse',
    version='0.1.0',
    description='Example and Schema Sources for Exabyte',
    url='https://github.com/Exabyte-io/exabyte-materials-json',
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    py_modules=["esse"],
    install_requires=[str(i.req) for i in parse_requirements(REQUIREMENTS_TXT, session=False)],
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development'
    ]
)
