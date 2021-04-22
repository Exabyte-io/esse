import os
from setuptools import find_packages, setup


def get_files_by_path(path):
    files = []
    for root, dirs_, files_ in os.walk(path, followlinks=True):
        for file_ in files_:
            files.append(os.path.join(root, file_))
    return files


with open('./README.md', 'r') as f:
    long_description = f.read()

DIR = os.path.abspath(os.path.dirname(__file__))
EXTRA_FILES = get_files_by_path(os.path.join(DIR, "example"))
EXTRA_FILES.extend(get_files_by_path(os.path.join(DIR, "schema")))
EXTRA_FILES.extend(get_files_by_path(os.path.join(DIR, "manifest")))

setup(
    name='esse',
    setup_requires=['setuptools_scm'],
    use_scm_version={
        'version_scheme': 'post-release',
    },
    description='Exabyte Source of Schemas and Examples',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Exabyte-io/exabyte-esse',
    author='Exabyte Inc.',
    author_email='info@exabyte.io',
    packages=find_packages(where='src/py', exclude=['src/py/tests*']),
    package_dir={'': 'src/py'},
    package_data={'esse': [f.replace(DIR, "data") for f in EXTRA_FILES]},
    install_requires=[
        "pyyaml>=4.2b1,<6",
        "jsonschema==2.6.0",
        "python-slugify==2.0.1",
        "exabyte_json_include>=2020.10.19"
    ],
    extras_require={
        "test": [
            "coverage[toml]>=5.3",
        ]
    },
    entry_points={
        'console_scripts': [
            'generate_dft_unit_functionals=esse.functionals:generate_dft_unit_functionals'
        ],
    },
    python_requires=">=3.6",
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development',
        'License :: OSI Approved :: Apache Software License'
    ]
)
