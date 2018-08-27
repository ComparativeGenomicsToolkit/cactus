FROM ubuntu:16.04 AS builder

RUN apt-get update && apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev libhiredis-dev pkg-config

ENV kyotoTycoonIncl -I/usr/include -DHAVE_KYOTO_TYCOON=1
ENV kyotoTycoonLib -L/usr/lib -Wl,-rpath,/usr/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
RUN mkdir -p /home/cactus

COPY . /home/cactus

RUN cd /home/cactus && make clean && rm -f /home/cactus/submodules/hdf5/bin/h5c++
RUN cd /home/cactus && make deps
RUN cd /home/cactus && make

# Create a thinner final Docker image in which only the binaries and necessary data exist.
FROM ubuntu:16.04
RUN apt-get update && apt-get install -y libkyotocabinet-dev libkyototycoon-dev libtokyocabinet-dev python zlib1g-dev python-dev libbz2-dev build-essential python-pip git kyototycoon valgrind net-tools redis-server libhiredis-dev
COPY --from=builder /home/cactus/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/sonLib/bin/* /usr/local/bin/
COPY --from=builder /home/cactus/submodules/cactus2hal/bin/* /usr/local/bin/

RUN mkdir /opt/cactus/
COPY runtime/wrapper.sh /opt/cactus/

ARG CACTUS_COMMIT

RUN pip install --pre toil
RUN pip install git+https://github.com/ComparativeGenomicsToolkit/sonLib@toil

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["bash", "/opt/cactus/wrapper.sh"]
