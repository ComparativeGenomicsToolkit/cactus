FROM ubuntu:12.04

COPY ./tools/* /usr/local/bin/

RUN mkdir /opt/cpecan/
COPY ./wrapper.sh /opt/cpecan/wrapper.sh

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/opt/cpecan/wrapper.sh"]