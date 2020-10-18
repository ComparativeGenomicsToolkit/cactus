FROM quay.io/glennhickey/cactus-ci-base:latest as builder

# apt dependencies for build
RUN apt-get update && apt-get install -y build-essential git python3 python3-dev python3-pip zlib1g-dev wget libbz2-dev pkg-config libhdf5-dev liblzo2-dev libtokyocabinet-dev wget

# build cactus binaries
RUN mkdir -p /home/cactus
COPY . /home/cactus

# compile with nehalem architecture target to improve portablity
ENV CFLAGS -march=nehalem
ENV CXXFLAGS -march=nehalem

# clean out stuff before build.
RUN find /home/cactus -name include.local.mk -exec rm -f {} \;
RUN cd /home/cactus && make clean -j $(nproc)
RUN cd /home/cactus && make -j $(nproc)

# download kent binaries used by hal for assembly hubs
# bedSort and hgGcPercent are part of a more restricted licence: https://hgdownload.cse.ucsc.edu/admin/exe/
# if you agree to it, add them by replacing the following line with this one:
# RUN cd /home/cactus/bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed bedSort hgGcPercent; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod ugo+x ${i}; done
RUN cd /home/cactus/bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod ugo+x ${i}; done

# download hal2vg
RUN cd /home/cactus/bin && ../build-tools/downloadHal2vg

# remove test executables
RUN cd /home/cactus && rm -f ${binPackageDir}/bin/*test ${binPackageDir}/bin/*tests ${binPackageDir}/bin/*Test ${binPackageDir}/bin/*Tests ${binPackageDir}/bin/cactus_runEndAlignment

# make the binaries smaller by removing debug symbols 
RUN cd /home/cactus && strip -d bin/* 2> /dev/null || true

# build cactus python3
RUN ln -fs /usr/bin/python3 /usr/bin/python
RUN mkdir -p /wheels && cd /wheels && python3 -m pip install -U pip && python3 -m pip wheel -r /home/cactus/toil-requirement.txt && python3 -m pip wheel /home/cactus

# Create a thinner final Docker image in which only the binaries and necessary data exist.
FROM ubuntu:bionic-20200112

# apt dependencies for runtime
RUN apt-get update && apt-get install -y --no-install-recommends git python3 python3-pip python3-distutils zlib1g libbz2-1.0 net-tools libhdf5-100 liblzo2-2 libtokyocabinet9 rsync libkrb5-3 libk5crypto3 time

# copy temporary files for installing cactus
COPY --from=builder /home/cactus /tmp/cactus
COPY --from=builder /wheels /wheels

# install the cactus binaries
RUN cd /tmp/cactus && cp -rP bin /usr/local/

# install the hal python modules
RUN rsync -avm --include='*.py' -f 'hide,! */' /tmp/cactus/submodules/hal /usr/local/lib
ENV PYTHONPATH /usr/local/lib:${PYTHONPATH}

# install the python3 binaries then clean up
RUN python3 -m pip install -U pip wheel setuptools && \
    python3 -m pip install -f /wheels -r /tmp/cactus/toil-requirement.txt && \
    python3 -m pip install -f /wheels /tmp/cactus && \
    rm -rf /wheels /root/.cache/pip/* /tmp/cactus && \
    apt-get remove -y git python3-pip rsync && \
    apt-get auto-remove -y

# check the linking on all our binaries (those kent tools above aren't static)
RUN for i in /usr/local/bin/* ; do if [ -f ${i} ] && [ $(ldd ${i} | grep "not found" | wc -l) -ge 1 ]; then exit 1; fi; done

# wrapper.sh is used when running using the docker image with --binariesMode docker
RUN mkdir /opt/cactus/
COPY runtime/wrapper.sh /opt/cactus/
RUN chmod 777 /opt/cactus/wrapper.sh

# log the memory usage (with --realTimeLogging) for local commands
ENV CACTUS_LOG_MEMORY 1

# remember where we came from
ARG CACTUS_COMMIT

# UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

