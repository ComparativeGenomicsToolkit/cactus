FROM ubuntu:12.04

MAINTAINER Alden Deran, adderan@ucsc.edu

RUN apt-get update && apt-get install -y build-essential zlib1g-dev git

ARG CPECAN_COMMIT=86d5de5c6915ec0a9963d75583fe60062d734268
ARG SONLIB_COMMIT=48e9613b0f884aaa9a2313556774a2934a5f109e

WORKDIR /home
RUN git clone https://github.com/benedictpaten/sonLib /home/sonLib
WORKDIR /home/sonLib
RUN git checkout $SON_LIB_COMMIT
RUN make
RUN git clone https://github.com/benedictpaten/cPecan /home/cPecan
RUN cd /home/cPecan && git checkout $CPECAN_COMMIT

WORKDIR /home/cPecan

RUN make

ENV PATH $PATH:/home/sonLib/bin

COPY ./wrapper.sh /home/wrapper.sh
RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/home/wrapper.sh"]