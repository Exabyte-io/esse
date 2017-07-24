from setuptools import setup, find_packages

setup(
    name='exabyte-materials-json',
    version='0.1.0',
    description='Contains schemas and examples for materials and simulations related data in JSON representation.',
    url='https://github.com/Exabyte-io/exabyte-materials-json',
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    packages=find_packages(exclude=['example', 'lib', 'tests']),
    install_requires=[
        "pyyaml==3.12",
        "jsonschema==2.6.0",
    ],
    dependency_links=[
        "git+git://github.com/Exabyte-io/json_include.git@master#egg=json_include"
    ],
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development'
    ]
)
