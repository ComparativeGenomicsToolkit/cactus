FROM ubuntu:12.04

MAINTAINER Alden Deran, adderan@ucsc.edu

RUN mkdir /opt/lastz/

COPY lastz /usr/local/bin/
COPY wrapper.sh /opt/lastz/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/lastz/wrapper.sh"]
