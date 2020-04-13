
```
tar -xzf "cactus-bin-${REL_TAG}.tar.gz"
cd "cactus-bin-${REL_TAG}"
virtualenv -p python3.6 venv
source venv/bin/activate
pip install -U setuptools pip
pip install -U -r ./toil-requirement.txt
pip install -U . -f ./wheels
export PATH=$(pwd)/bin:$PATH
```
# test:
```
cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal --realTimeLogging
``
