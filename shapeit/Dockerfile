FROM ubuntu:18.04
MAINTAINER rick

ENV SHAPEIT_VERSION=2.17
# set package dir.
ENV PATH="$PATH:/usr/local/bin:/opt/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin" 
# without interactive
ENV DEBIAN_FRONTEND=noninteractive
# Install ubuntu packages
RUN apt update \
 && apt-get install -y git git-lfs perl gzip unzip zip tabix plink1.9 \
# git clone
 && cd /tmp/ \
 && git clone https://gitlab+deploy-token-4:WjMMZ33rcmz17RHfqszz@gitlab.corp.ailabs.tw/engine/3rd-party.git \ 
 && cd /tmp/3rd-party/ \
 && git lfs pull --include "shapeit-${SHAPEIT_VERSION}" \
 && cd /tmp/3rd-party/shapeit-${SHAPEIT_VERSION}/ \
# Install shapeit
 && tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz \
 && cp -r shapeit.v2.904.3.10.0-693.11.6.el7.x86_64 /opt/ \
# Clean up
 && rm -rf /tmp/*

# Load shell script
COPY run.sh /app/
WORKDIR /app

CMD ["bash", "run.sh"]
