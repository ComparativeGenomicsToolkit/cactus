from setuptools import setup, find_packages

setup(
    name="progressiveCactus",
    version="1.0",
    author="Benedict Paten",
    package_dir = {'': 'src'},
    packages=find_packages(where='src'),

    install_requires=[
        'psutil',
        'networkx'],
    
    entry_points={
        'console_scripts': [
            'progressiveCactus = cactus.progressive.cactus_progressive:main']},)

