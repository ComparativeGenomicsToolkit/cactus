# definitions and functions for release bash programs

PYTHON=python3.6
PIP="${PYTHON} -m pip"

# get the tag for the lastest release, in the form v1.2.3, from git

getLatestReleaseTag() {
    git describe --tags $(git rev-list --tags --max-count=10) | egrep -e '^v[0-9]+\.[0-9]+\.[0-9]+$' | head -1
}
