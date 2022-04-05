# Installation of the Cactus binary distribution 

This describes the steps require to install the Cactus
pre-compile binary, static linked distribution.

## Extracting
If you have not already extract the distribution:
```
tar -xzf "cactus-bin-${REL_TAG}.tar.gz"
cd "cactus-bin-${REL_TAG}"
```

## Setup

To build a python virtualenv and activate, do the following steps:
```
virtualenv -p python3.8 cactus_env
echo "export PATH=$(pwd)/bin:\$PATH" >> cactus_env/bin/activate
echo "export PYTHONPATH=$(pwd)/lib:\$PYTHONPATH" >> cactus_env/bin/activate
source cactus_env/bin/activate
python3 -m pip install -U setuptools pip==21.3.1
python3 -m pip install -U -r ./toil-requirement.txt
python3 -m pip install -U .
```

Some tools required for `hal2assemblyHub.py` are not included and must be downloaded separately.
They are `wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent`.  More information
can be found [here](https://hgdownload.cse.ucsc.edu/admin/exe/).  Note that some may require
a license for commercial use.  Static binaries are not available, but the following command
should set them up successfully on many 64 bit Linux systems:
```
cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod ugo+x ${i}; done
```

## Testing

To test Cactus, the following will run a moderately sized alignment.  It may
take several hours, depending on your system.
```
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal --realTimeLogging
``
