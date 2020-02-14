from setuptools import setup, find_packages
from setuptools.command.install import install
import os
import subprocess

# FIXME this is duplicated in makefile, we need to sort this out
versionFile = "src/cactus/shared/version.py"
if os.path.exists(versionFile):
    os.remove(versionFile)
git_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], encoding="ascii").strip()
with open(versionFile, 'w') as versionFH:
    versionFH.write("cactus_commit = '%s'\n" % git_commit)

class PostInstallCommand(install):
    """Post-installation customization.  Ensure sonLib submodule in the tree is installed virtual env."""

    def run(self):
        subprocess.run(["pip", "install", "submodules/sonLib"], check=True)
        install.run(self)

setup(
    name = "Cactus",
    version = "1.0",
    author = "Benedict Paten",
    package_dir = {'': 'src'},
    packages = find_packages(where='src'),
    include_package_data = True,
    package_data = {
        'cactus': ['*_config.xml']
    },
    # We use the __file__ attribute so this package isn't zip_safe.
    zip_safe = False,

    python_requires = '>=3.6',

    install_requires = [
        'toil>=3.24',
        'decorator',
        'psutil',
        'networkx>=2,<3',
        'cython',
        'pytest'],

    cmdclass = {
        'install': PostInstallCommand,
    },
    entry_points= {
        'console_scripts': ['cactus = cactus.progressive.cactus_progressive:main',
                            'cactus_preprocess = cactus.preprocessor.cactus_preprocessor:main']},)
