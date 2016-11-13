from setuptools import setup, find_packages
import os

os.system("pip install git+https://github.com/ComparativeGenomicsToolkit/sonLib@toil")

setup(
    name="progressiveCactus",
    version="1.0",
    author="Benedict Paten",
    package_dir = {'': 'src'},
    packages=find_packages(where='src'),

    install_requires=[
        'psutil',
        'networkx',
        'toil_lib==1.2.0a1.dev119'],
    
    entry_points={
        'console_scripts': [
            'progressiveCactus = cactus.progressive.cactus_progressive:main']},)

