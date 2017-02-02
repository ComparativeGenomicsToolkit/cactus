FROM ubuntu:16.04

RUN apt-get update && apt-get install -y libkyotocabinet-dev libkyototycoon-dev libtokyocabinet-dev python zlib1g-dev python-dev libbz2-dev build-essential python-pip git kyototycoon valgrind
RUN mkdir /opt/cactus/
COPY ./tools/* /usr/local/bin/

COPY wrapper.sh /opt/cactus/

ARG CACTUS_COMMIT

RUN pip install --pre toil
RUN pip install git+https://github.com/ComparativeGenomicsToolkit/sonLib@toil

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/opt/cactus/wrapper.sh"]
