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
virtualenv -p python3.6 venv
source venv/bin/activate
pip install -U setuptools pip
pip install -U -r ./toil-requirement.txt -f ./wheels
pip install -U . -f ./wheels
export PATH=$(pwd)/bin:$PATH
```

## Testing

To test Cactus, the following will run a moderately sized alignment.  It may
take several hours, depending on your system.
```
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal --realTimeLogging
``
