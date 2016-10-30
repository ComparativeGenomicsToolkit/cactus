FROM ubuntu:16.04


RUN apt-get update && apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev

ARG CACTUS_COMMIT=48398f54b864a3cac1526a28e7f5a2350a6aed63
ARG PINCHES_AND_CACTI_COMMIT=9c16563158d8d510767ca7ca2ed4f0b1e47e4277
ARG MATCHING_AND_AND_ORDERING_COMMIT=200a1b0d65ce816f7d83ca8d0a19bc807d4bde99
ARG CPECAN_COMMIT=86d5de5c6915ec0a9963d75583fe60062d734268
ARG SON_LIB_COMMIT=48e9613b0f884aaa9a2313556774a2934a5f109e

RUN git clone https://github.com/ComparativeGenomicsToolkit/sonLib.git /home/sonLib
RUN git clone https://github.com/ComparativeGenomicsToolkit/matchingAndOrdering.git /home/matchingAndOrdering
RUN git clone https://github.com/ComparativeGenomicsToolkit/pinchesAndCacti.git /home/pinchesAndCacti

RUN git clone https://github.com/benedictpaten/cPecan /home/cPecan
RUN git clone https://github.com/ComparativeGenomicsToolkit/cactus.git /home/cactus

ENV sonLibRootPath /home/sonLib
RUN cd /home/sonLib && git pull && git checkout $SON_LIB_COMMIT && make
RUN cd /home/matchingAndOrdering && git pull && git checkout $MATCHING_AND_ORDERING_COMMIT && make
RUN cd /home/pinchesAndCacti && git pull && git checkout $PINCHES_AND_CACTI_COMMIT && make
RUN cd /home/cPecan && git pull && git checkout $CPECAN_COMMIT && make

RUN cd /home/cactus && git pull && git checkout $CACTUS_COMMIT
RUN cd /home/cactus && make


#HACK - eventually all python code should be pip-installed and not run
#in containers, but for now there are various python executables
#that depend on sonLib
ENV PYTHONPATH $PYTHONPATH:/home/

COPY ./wrapper.sh /home/wrapper.sh
RUN mkdir /data
WORKDIR /data
ENTRYPOINT ["bash", "/home/wrapper.sh"]
