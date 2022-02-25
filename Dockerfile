
FROM ubuntu:20.04

## Environment variables
ENV SAMTOOLS_VERSION=1.9

# set a directory for the app
WORKDIR /docker

# install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends && \
    apt-get install -y python3-pip && \
    apt-get install -y python3-venv && \
    apt-get install -y curl

## Create virtual environment
ENV VIRTUAL_ENV=$HOME/.virtualenvs/km
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

## Download and install jellyfish
RUN curl -L https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz --output jellyfish-2.2.6.tar.gz
RUN tar zxvf jellyfish-2.2.6.tar.gz && \
    cd jellyfish-2.2.6 && \
    ./configure --prefix=$VIRTUAL_ENV --enable-python-binding && \
    make -j 4 && \
    make install 

# Clean up
RUN rm -rf jellyfish-2.2.6 && \
    rm -rf jellyfish-2.2.6.tar.gz

# install km through a virtualenv - how the hell is this even possible?
RUN curl -L https://github.com/iric-soft/km/archive/refs/tags/2.0.2.tar.gz --output km-v2.0.2.tar.gz
RUN tar zxvf km-v2.0.2.tar.gz && \
    cd km-2.0.2 && \
    $VIRTUAL_ENV/bin/pip install pip setuptools wheel numpy --upgrade && \
    $VIRTUAL_ENV/bin/python setup.py install

# Clean up
RUN rm -rf km-2.0.2 && \
    rm -rf km-v2.0.2.tar.gz

# Install R by Add PPA for R 4.0
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && \
    apt-get -y upgrade && \ 
    apt -y install software-properties-common dirmngr && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -sc)-cran40/" && \
    apt-get -y install r-base && \
    Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")' && \
    Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' && \
    Rscript -e 'BiocManager::install(c("tidyverse", "Biostrings", "rtracklayer"))' && \
    apt-get clean && rm -r /var/cache/

# install samtools
RUN curl -SL -o samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR /docker/samtools-${SAMTOOLS_VERSION}
RUN ./configure --prefix=/usr/local
RUN make && \
    make install

# cleanup
WORKDIR ${HOME}/
RUN rm -rf /docker
ENV PATH="${PATH}:/usr/local/bin"