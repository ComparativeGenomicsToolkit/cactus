FROM ubuntu:16.04


RUN apt-get update && apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev

ARG CACTUS_COMMIT=fa2e170ceff7f102b8253d30b839af60605e40ef
ARG PINCHES_AND_CACTI_COMMIT=9c16563158d8d510767ca7ca2ed4f0b1e47e4277
ARG MATCHING_AND_AND_ORDERING_COMMIT=200a1b0d65ce816f7d83ca8d0a19bc807d4bde99
ARG CPECAN_COMMIT=86d5de5c6915ec0a9963d75583fe60062d734268
ARG SON_LIB_COMMIT=48e9613b0f884aaa9a2313556774a2934a5f109e

RUN git clone https://github.com/ComparativeGenomicsToolkit/sonLib.git /home/sonLib
RUN git clone https://github.com/ComparativeGenomicsToolkit/matchingAndOrdering.git /home/matchingAndOrdering
RUN git clone https://github.com/ComparativeGenomicsToolkit/pinchesAndCacti.git /home/pinchesAndCacti

RUN git clone https://github.com/benedictpaten/cPecan /home/cPecan
RUN git clone https://github.com/ComparativeGenomicsToolkit/cactus.git /home/cactus

RUN git clone https://github.com/UCSCReconGroup/hdf5 /home/hdf5
RUN git clone https://github.com/glennhickey/hal /home/hal
RUN git clone https://github.com/glennhickey/cactus2hal /home/cactus2hal

ENV sonLibRootPath /home/sonLib
RUN cd /home/sonLib && git pull && git checkout $SON_LIB_COMMIT && make
RUN cd /home/matchingAndOrdering && git pull && git checkout $MATCHING_AND_ORDERING_COMMIT && make
RUN cd /home/pinchesAndCacti && git pull && git checkout $PINCHES_AND_CACTI_COMMIT && make
RUN cd /home/cPecan && git pull && git checkout $CPECAN_COMMIT && make 2> /dev/null

RUN cd /home/cactus && git pull && git checkout $CACTUS_COMMIT
RUN cd /home/cactus && make

WORKDIR /home/hdf5
RUN ./configure --prefix=/home/hdf5 --enable-cxx 2> /dev/null
RUN CFLAGS=-std=c99 make -e 2> /dev/null
RUN make install 2> /dev/null
ENV PATH $PATH:/home/hdf5/bin

WORKDIR /home/hal
RUN make 2> /dev/null

WORKDIR /home/cactus2hal
RUN make 2> /dev/null


ENV PATH $PATH:/home/cactus/bin:/home/cactus2hal/bin

RUN mkdir /data
WORKDIR /data
