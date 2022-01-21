from setuptools import setup, find_packages
from setuptools.command.install import install
import os, sys
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
        subprocess.run([sys.executable, "-m", "pip", "install", "submodules/sonLib"], check=True)
        install.run(self)

setup(
    name = "Cactus",
    version = "1.0",
    author = "Benedict Paten",
    package_dir = {'': 'src'},
    packages = find_packages(where='src'),
    include_package_data = True,
    package_data = {
        'cactus': ['*_config.xml', '*.knm']
    },
    # We use the __file__ attribute so this package isn't zip_safe.
    zip_safe = False,

    python_requires = '>=3.6',

    install_requires = [
        'decorator',
        'networkx>=2,<3',
        'pytest',
        'cigar',
        'biopython'], 

    cmdclass = {
        'install': PostInstallCommand,
    },
    entry_points= {
        'console_scripts': ['cactus = cactus.progressive.cactus_progressive:main',
                            'cactus-preprocess = cactus.preprocessor.cactus_preprocessor:main',
                            'cactus-prepare = cactus.progressive.cactus_prepare:main',
                            'cactus-prepare-toil = cactus.progressive.cactus_prepare:main_toil',
                            'cactus-blast = cactus.blast.cactus_blast:main',
                            'cactus-refmap = cactus.refmap.cactus_refmap:main',
                            'cactus-graphmap = cactus.refmap.cactus_graphmap:main',
                            'cactus-graphmap-split = cactus.refmap.cactus_graphmap_split:main',
                            'cactus-graphmap-join = cactus.refmap.cactus_graphmap_join:main',
                            'cactus-align = cactus.setup.cactus_align:main',
                            'cactus-align-batch = cactus.setup.cactus_align:main_batch']},)
