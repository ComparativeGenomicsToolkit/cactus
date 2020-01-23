FROM ubuntu:16.04 AS builder

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository "deb http://mirrors.kernel.org/ubuntu/ xenial universe"
RUN apt-get install -y git gcc g++ build-essential python3 python3-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev libhiredis-dev pkg-config libhdf5-cpp-11 libhdf5-dev

ENV kyotoTycoonIncl -I/usr/include -DHAVE_KYOTO_TYCOON=1
ENV kyotoTycoonLib -L/usr/lib -Wl,-rpath,/usr/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
RUN mkdir -p /home/cactus

COPY . /home/cactus

RUN cd /home/cactus && make -j 10 clean
RUN cd /home/cactus && make -j 10

# Create a thinner final Docker image in which only the binaries and necessary data exist.
FROM ubuntu:16.04
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository "deb http://mirrors.kernel.org/ubuntu/ xenial universe"
RUN apt-get install -y libkyotocabinet-dev libkyototycoon-dev libtokyocabinet-dev python3 zlib1g-dev python3-dev libbz2-dev build-essential python3-pip git kyototycoon net-tools redis-server libhiredis-dev libhdf5-cpp-11
COPY --from=builder /home/cactus/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/sonLib/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/cactus2hal/bin/* /usr/local/bin/

RUN mkdir /opt/cactus/
COPY runtime/wrapper.sh /opt/cactus/

ARG CACTUS_COMMIT

# FIXME: install from git until new release
RUN pip3 install --pre git+https://github.com/DataBiosphere/toil.git
RUN pip3 install /home/cactus/submodules/sonLib/src

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/opt/cactus/wrapper.sh"]
