from setuptools import setup, find_packages

setup(
    name="progressiveCactus",
    version="1.0",
    author="Benedict Paten",

    packages=find_packages(),

    install_requires=[
        'psutil',
        'toil==3.5.0a1.dev251',
        'toil-lib==1.1.0a1.dev100',
        'networkx'],
    
    entry_points={
        'console_scripts': [
        'progressiveCactus = cactus.progressive.cactus_progressive:main']},)

