# Installation of the Cactus binary distribution 

This describes the steps require to install the Cactus
pre-compile binary, static linked distribution.

## Extracting
If you have not already extract the distribution and cd into the cactus directory:
```
tar -xzf cactus-bin-v2.9.8.tar.gz
cd cactus-bin-v2.9.8
```

## Setup

To build a python virtualenv and activate, do the following steps. This requires Python version >= 3.9 (so Ubuntu 18.04 users should use `-p python3.9` below):
```
virtualenv -p python3 venv-cactus-v2.9.8
printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\nexport LD_LIBRARY_PATH=$(pwd)/lib:\$LD_LIBRARY_PATH\n" >> venv-cactus-v2.9.8/bin/activate
source venv-cactus-v2.9.8/bin/activate
python3 -m pip install -U setuptools pip wheel
python3 -m pip install -U .
python3 -m pip install -U -r ./toil-requirement.txt
```

Some tools required for `hal2assemblyHub.py`, `cactus-hal2chains` and `cactus-maf2bigmaf` are not included and must be downloaded separately.
They are `wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary hgLoadChain`.  More information
can be found [here](https://hgdownload.cse.ucsc.edu/admin/exe/).  Note that some may require
a license for commercial use.  Static binaries are not available, but the following command
should set them up successfully on many 64 bit Linux systems:
```
cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary hgLoadChain; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod +x ${i}; done
```

`vcfwave` isn't included in the release binaries (but is in the docker image).  You can can try building it and adding it to `bin/` with the following command
```
build-tools/downloadVCFWave
```

## Testing

To test Cactus, the following will run a tiny sumulated alignment.
```
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal
``
