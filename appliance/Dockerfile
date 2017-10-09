FROM quay.io/ucsc_cgl/toil:3.11.0
ENV CACTUS_BINARIES_MODE=local
ENV TOIL_AWS_AMI=ami-6357a31b
RUN apt-get update
RUN apt-get install -y git
RUN virtualenv --system-site-packages venv
COPY runtime/tools/* /bin/
RUN mkdir cactus
COPY *.py cactus/
COPY src cactus/src
RUN bash -c 'source /venv/bin/activate && cd cactus && pip install .'
CMD bash -c '. /venv/bin/activate; exec /bin/bash -i'