from setuptools import setup

setup(
    name='esse',
    version='0.1.0',
    description='Example and Schema Sources for Exabyte',
    url='https://github.com/Exabyte-io/exabyte-materials-json',
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    py_modules=["esse"],
    install_requires=[
        "pyyaml==3.12",
        "jsonschema==2.6.0",
        "json_include==0.2.9"
    ],
    dependency_links=[
        "git+https://git@github.com/Exabyte-io/json_include.git@3501e40dedbed4dc175411a05c8f1fa0d717e974#egg=json_include-0.2.9"
    ],
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development'
    ]
)
