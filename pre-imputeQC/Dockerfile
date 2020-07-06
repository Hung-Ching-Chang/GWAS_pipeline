FROM ubuntu:18.04
MAINTAINER rick

ENV PRE_IMPUTEQC_VERSION=1.0.0
# set package dir.
ENV PATH="$PATH:/usr/local/bin" 
# without interactive
ENV DEBIAN_FRONTEND=noninteractive
# Install ubuntu packages
RUN apt update \
 && apt-get install -y git git-lfs perl gzip unzip zip tabix plink1.9 samtools bcftools \
# install python
 && apt-get install -y python3.7 python3-pip \
 && pip3 install fpdf \
# R & R packages 
 && apt-get install -y r-base libgdal-dev \
 && R -e 'install.packages("data.table"); install.packages("scales"); install.packages("ggplot2"); install.packages("plyr")' \
# git clone
 && cd /tmp/ \
 && git clone https://gitlab+deploy-token-4:WjMMZ33rcmz17RHfqszz@gitlab.corp.ailabs.tw/engine/3rd-party.git \ 
 && cd /tmp/3rd-party/ \
 && git lfs pull --include "pre-imputeQC-${PRE_IMPUTEQC_VERSION}" \
 && cd /tmp/3rd-party/pre-imputeQC-${PRE_IMPUTEQC_VERSION}/ \
 && cp -r . /opt/ \
# Clean up
 && rm -rf /tmp/*

# Load shell script
COPY run.sh /app/
WORKDIR /app

CMD ["bash", "run.sh"]
