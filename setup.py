from setuptools import find_packages, setup

from DigCNVlib import digCNV_logger

setup(
    name='DigCNVlib',
    packages=find_packages(include=['DigCNVlib']),
    version='0.1.1',
    description='DigCNV library: Discriminates true from false CNVs called by at least two algorithm (QuantiSNP and PennCNV)',
    author='Thomas Renne',
    license='GPL >= 3.0',
    install_requires=['pandas', 'imblearn', 'sklearn', 'matplotlib', 'seaborn'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests'
)
