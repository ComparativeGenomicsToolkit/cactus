FROM ubuntu:12.04

MAINTAINER Alden Deran, adderan@ucsc.edu

RUN apt-get update && apt-get install -y build-essential zlib1g-dev git

WORKDIR /home
RUN git clone https://github.com/benedictpaten/sonLib /home/sonLib
WORKDIR /home/sonLib
RUN make
RUN git clone https://github.com/benedictpaten/cPecan /home/cPecan

WORKDIR /home/cPecan

RUN make
