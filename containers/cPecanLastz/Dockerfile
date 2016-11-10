FROM ubuntu:12.04

MAINTAINER Alden Deran, adderan@ucsc.edu

RUN mkdir /opt/cPecanLastz

COPY cPecanLastz /usr/local/bin/
COPY wrapper.sh /opt/cPecanLastz/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/cPecanLastz/wrapper.sh"]