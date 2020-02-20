FROM ubuntu:bionic-20200112 AS builder


RUN apt-get update
RUN apt-get install -y git gcc g++ build-essential python3 python3-dev zlib1g-dev wget valgrind libbz2-dev libhiredis-dev pkg-config libhdf5-dev liblzo2-dev libtokyocabinet-dev

RUN mkdir -p /home/cactus

COPY . /home/cactus

RUN cd /home/cactus && make -j $(nproc) clean
RUN cd /home/cactus && make -j $(nproc)

# Create a thinner final Docker image in which only the binaries and necessary data exist.
FROM ubuntu:bionic-20200112

RUN apt-get update

RUN apt-get install -y python3 zlib1g-dev python3-dev libbz2-dev build-essential python3-pip git net-tools redis-server libhiredis-dev libhdf5-100 liblzo2-2 libtokyocabinet-dev
COPY --from=builder /home/cactus/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/lib/* /usr/local/lib/
COPY --from=builder /home/cactus/submodules/sonLib/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/cactus2hal/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/sonLib /tmp/sonLib/

RUN mkdir /opt/cactus/
COPY runtime/wrapper.sh /opt/cactus/

ARG CACTUS_COMMIT

RUN pip3 install /tmp/sonLib
RUN rm -rf /tmp/sonLib

RUN mkdir /data
WORKDIR /data

ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

ENTRYPOINT ["bash", "/opt/cactus/wrapper.sh"]
