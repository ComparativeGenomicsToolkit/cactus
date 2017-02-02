FROM ubuntu:16.04


RUN apt-get update && apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev

RUN mkdir -p /home/cactus

COPY . /home/cactus

RUN cd /home/cactus && make

RUN mkdir /data
WORKDIR /data
