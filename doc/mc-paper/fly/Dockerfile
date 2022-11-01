FROM ubuntu:20.04
MAINTAINER vgteam

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

## python, snakemake and awscli
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    git screen wget curl gcc less nano \
    sudo \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    make \
    pigz \
    tabix \
    bzip2 \
    tzdata \
    gfortran-8 \
    g++ \
    pkg-config \
    build-essential cmake doxygen \
    apt-transport-https software-properties-common dirmngr gpg-agent \ 
    libncurses5-dev libncursesw5-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libxml2-dev libssl-dev libmariadbclient-dev libcurl4-openssl-dev \ 
    && rm -rf /var/lib/apt/lists/*

ENV TZ=America/Los_Angeles

RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir requests awscli snakemake==5.8.2 biopython pyfaidx pyvcf \
        pandas boto3 numpy matplotlib pyvcf seaborn cyvcf2

WORKDIR /build

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
        tar -xjf bcftools-1.10.2.tar.bz2 && \
        cd bcftools-1.10.2 && \
        ./configure && make && make install && \
        cd .. && rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

## vcftools
RUN wget --no-check-certificate https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
        tar -xzvf vcftools-0.1.16.tar.gz

ENV PERL5LIB=/build/vcftools-0.1.16/src/perl

WORKDIR /build/vcftools-0.1.16
RUN ./configure && make && make install && \
        cd .. && rm -rf vcftools-0.1.16 vcftools-0.1.16.tar.gz

## R and pandoc
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
        && apt-get update \
        && apt-get install -y r-base r-base-dev

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz')"

RUN R -e "install.packages(c('dplyr', 'knitr', 'ggplot2', 'rmarkdown', 'DT', 'ggrepel', 'tidyr', 'BiocManager', 'plotly', 'DiagrammeR', 'visNetwork', 'shiny'))"

RUN R -e "BiocManager::install(c('VariantAnnotation'))"

RUN apt-get update \
    && wget https://github.com/jgm/pandoc/releases/download/2.10.1/pandoc-2.10.1-1-amd64.deb \
    && dpkg -i pandoc-2.10.1-1-amd64.deb \
    && apt-get install -f \
    && rm pandoc-2.10.1-1-amd64.deb

## latex and packages
RUN R -e "install.packages('tinytex')"

RUN R -e "tinytex::install_tinytex()"

RUN R -e "tinytex::tlmgr_install(c('beamer', 'fp', 'ms', 'pgf', 'translator', 'caption'))"

RUN R -e "install.packages(c('R.utils'))"

## vg
WORKDIR /bin
RUN wget --no-check-certificate https://github.com/vgteam/vg/releases/download/v1.37.0/vg && \
        chmod +x vg

## jq
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    jq \
    && rm -rf /var/lib/apt/lists/*

## libbdsg
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pybind11

WORKDIR /build
RUN git clone --recursive https://github.com/vgteam/libbdsg.git && \
    cd libbdsg && \
    # git checkout cd51123736d13af17b5b4c43c762232a6c17c5c1 && \
    # git submodule update --init --recursive && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 4

RUN pip3 install --no-cache-dir bdsg

RUN R -e "BiocManager::install('remotes')"
RUN R -e "BiocManager::install('jmonlong/sveval@dev')"

## init home
WORKDIR /home
