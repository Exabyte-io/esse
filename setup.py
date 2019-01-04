from setuptools import setup

setup(
    name='esse',
    version='0.1.0',
    description='Example and Schema Sources for Exabyte',
    url='https://github.com/Exabyte-io/exabyte-esse',
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    py_modules=["esse"],
    install_requires=[
        "pyyaml==3.12",
        "jsonschema==2.6.0",
        "ordereddict==1.1",
    ],
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development'
    ]
)
