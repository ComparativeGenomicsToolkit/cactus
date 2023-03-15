FROM quay.io/comparative-genomics-toolkit/ubuntu:22.04 AS builder

# apt dependencies for build
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential git python3 python3-dev python3-pip zlib1g-dev wget libbz2-dev pkg-config libhdf5-dev liblzo2-dev libtokyocabinet-dev wget liblzma-dev libxml2-dev libssl-dev libpng-dev uuid-dev libcurl4-gnutls-dev libffi-dev python3-virtualenv rsync python-is-python3 libdeflate-dev

# build cactus binaries
RUN mkdir -p /home/cactus
COPY . /home/cactus

# Make sure abpoa doesn't build with -march=native, but something more portable
# Todo: It would be more portable to use "sse41", but that leads to segfaults in rare cases
# https://github.com/yangao07/abPOA/issues/26
ENV avx2 1

# install Phast and enable halPhyloP compilation
RUN cd /home/cactus && ./build-tools/downloadPhast
ENV ENABLE_PHYLOP 1

# Install UCSC browser libraries to compile UDC
# remote access.  The browser common.mk file checks for
RUN cd /home/cactus && ./build-tools/downloadUcscLib
ENV ENABLE_UDC 1
ENV KENTSRC /home/cactus/submodules/kent/src

# clean and build
RUN find /home/cactus -name include.local.mk -exec rm -f {} \; && \
	 cd /home/cactus && rm -rf bin/* && make clean -j $(nproc) && \
	 make -j $(nproc)

# download open-licenses kent binaries used by hal for assembly hubs and / or chains
RUN cd /home/cactus/bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod ugo+x ${i}; done

# download tools used for pangenome pipeline
RUN cd /home/cactus && ./build-tools/downloadPangenomeTools

# download tools used for working with MAF
RUN cd /home/cactus && ./build-tools/downloadMafTools

# remove test executables
RUN cd /home/cactus && rm -f ${binPackageDir}/bin/*test ${binPackageDir}/bin/*tests ${binPackageDir}/bin/*Test ${binPackageDir}/bin/*Tests

# make the binaries smaller by removing debug symbols (but leave them in cactus_consolidated)
RUN /bin/bash -O extglob -c "cd /home/cactus && strip -d bin/!(cactus_consolidated) 2> /dev/null || true"

# check the linking on all our binaries (those kent tools above aren't static)
RUN for i in /usr/local/bin/* ; do if [ -f ${i} ] && [ $(ldd ${i} | grep "not found" | wc -l) -ge 1 ]; then exit 1; fi; done

# build cactus python3
RUN cd /home/cactus && rm -rf cactus_env && \
	 python3 -m virtualenv -p python3 cactus_env  && \
	 . cactus_env/bin/activate && \
	 python3 -m pip install -U setuptools pip && \
	 python3 -m pip install -U -r ./toil-requirement.txt && \
	 python3 -m pip install -U .
	 
# prep the hal python install which is not part of the setup
RUN rm -rf /home/cactus/hal_lib && \
	 rsync -avm --include='*.py' -f 'hide,! */' /home/cactus/submodules/hal /home/cactus/hal_lib

# Create a thinner final Docker image in which only the binaries and necessary data exist.
FROM quay.io/comparative-genomics-toolkit/ubuntu:22.04

# apt dependencies for runtime
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends git python3 python3-pip python3-distutils zlib1g libbz2-1.0 net-tools libhdf5-103 liblzo2-2 libtokyocabinet9 libkrb5-3 libk5crypto3 time liblzma5 libcurl4 libcurl4-gnutls-dev libxml2 libgomp1 libffi7 parallel libdeflate0

# required for ubuntu22 but won't work anywhere else
RUN bash -c "if ! command -v catchsegv > /dev/null; then apt-get install glibc-tools; fi"

# copy cactus runtime essentials (note: important cactus_env keeps its path)
RUN mkdir /home/cactus
COPY --from=builder /home/cactus/cactus_env /home/cactus/cactus_env
COPY --from=builder /home/cactus/hal_lib/hal /home/cactus/hal
COPY --from=builder /home/cactus/bin /home/cactus/bin

# update the environment
ENV PATH="/home/cactus/cactus_env/bin:/home/cactus/bin:$PATH"
ENV PYTHONPATH="/home/cactus"

# sanity check to make sure cactus at least runs
RUN cactus --help

# wrapper.sh is used when running using the docker image with --binariesMode docker
RUN mkdir /opt/cactus/
COPY runtime/wrapper.sh /opt/cactus/
RUN chmod 777 /opt/cactus/wrapper.sh

# log the memory usage (with --realTimeLogging) for local commands
ENV CACTUS_LOG_MEMORY 1

# UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

