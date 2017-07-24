from setuptools import setup, find_packages

setup(
    name='exabyte_materials_json',
    version='0.1.0',
    description='Contains schemas and examples for materials and simulations related data in JSON representation.',
    url='https://github.com/Exabyte-io/exabyte-materials-json',
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    packages=find_packages(exclude=['examples', 'docs', 'tests*']),
    install_requires=[
    ],
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development'
    ]
)
