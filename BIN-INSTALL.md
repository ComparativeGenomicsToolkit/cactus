# Installation of the Cactus binary distribution 

This describes the steps require to install the Cactus
pre-compile binary, static linked distribution.

## Extracting
If you have not already extract the distribution and cd into the cactus directory:
```
tar -xzf cactus-bin-v2.7.2.tar.gz
cd cactus-bin-v2.7.2
```

## Setup

To build a python virtualenv and activate, do the following steps. This requires Python version >= 3.7 (so Ubuntu 18.04 users should use `-p python3.8` below):
```
virtualenv -p python3 venv-cactus-v2.7.2
printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\n" >> venv-cactus-v2.7.2/bin/activate
source venv-cactus-v2.7.2/bin/activate
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
```

Some tools required for `hal2assemblyHub.py`, `cactus-hal2chains` and `cactus-maf2bigmaf` are not included and must be downloaded separately.
They are `wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary`.  More information
can be found [here](https://hgdownload.cse.ucsc.edu/admin/exe/).  Note that some may require
a license for commercial use.  Static binaries are not available, but the following command
should set them up successfully on many 64 bit Linux systems:
```
cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod +x ${i}; done
```

## Testing

To test Cactus, the following will run a tiny sumulated alignment.
```
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal --realTimeLogging
``
