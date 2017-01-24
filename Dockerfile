FROM ubuntu:16.04


RUN apt-get update && apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev

RUN mkdir -p /home/cactus
COPY . /home/cactus

ENV sonLibRootPath /home/cactus/submodules/sonLib
RUN cd /home/cactus/submodules/sonLib && make
RUN cd /home/cactus/submodules/pinchesAndCacti && make
RUN cd /home/cactus/submodules/matchingAndOrdering && make
RUN cd /home/cactus/submodules/cPecan && make
RUN cd /home/cactus && make
RUN cd /home/cactus/submodules/hdf5 && ./configure --prefix=/home/cactus/submodules/hdf5 --enable-cxx 2> /dev/null
RUN cd /home/cactus/submodules/hdf5 && CFLAGS=-std=c99 make -e 2> /dev/null
RUN cd /home/cactus/submodules/hdf5 && make install 2> /dev/null
ENV PATH /home/cactus/submodules/hdf5/bin:$PATH
RUN cd /home/cactus/submodules/hal && make
RUN cd /home/cactus/submodules/cactus2hal && make


ENV PATH $PATH:/home/cactus/bin:/home/cactus2hal/bin

RUN mkdir /data
WORKDIR /data
