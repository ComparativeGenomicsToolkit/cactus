FROM ubuntu:12.04

MAINTAINER Alden Deran, adderan@ucsc.edu

RUN mkdir /opt/cPecanRealign

COPY cPecanRealign /usr/local/bin/
COPY wrapper.sh /opt/cPecanRealign/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/cPecanRealign/wrapper.sh"]