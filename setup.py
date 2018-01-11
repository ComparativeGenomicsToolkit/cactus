from setuptools import setup, find_packages
import os
import subprocess

os.system("pip install git+https://github.com/ComparativeGenomicsToolkit/sonLib@toil")

versionFile = "src/cactus/shared/version.py"
if os.path.exists(versionFile):
    os.remove(versionFile)
git_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
with open(versionFile, 'w') as versionFH:
    versionFH.write("cactus_commit = '%s'" % git_commit)


setup(
    name="progressiveCactus",
    version="1.0",
    author="Benedict Paten",
    package_dir = {'': 'src'},
    packages=find_packages(where='src'),
    include_package_data=True,
    package_data={'cactus': ['*_config.xml']},
    # We use the __file__ attribute so this package isn't zip_safe.
    zip_safe=False,

    install_requires=[
        'decorator',
        'subprocess32',
        'psutil',
        'networkx==1.11'],
    
    entry_points={
        'console_scripts': ['cactus = cactus.progressive.cactus_progressive:main']},)
