FROM ubuntu:18.04
MAINTAINER rick

ENV MINIMAC_VERSION=4.0.0
ENV SHAPEIT_VERSION=2.17
# set package dir.
ENV PATH="$PATH:/usr/local/bin" 
ENV PATH="$PATH:/opt/Minimac4/release-build/" 
ENV PATH="$PATH:/opt/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin" 
# Install ubuntu packages
RUN apt update \
 && apt-get install -y git git-lfs unzip zip bash-completion libgomp1 plink1.9 bcftools \
# Download minimac from git lfs
 && cd /tmp/ \
 && git clone https://gitlab+deploy-token-4:WjMMZ33rcmz17RHfqszz@gitlab.corp.ailabs.tw/engine/3rd-party.git \ 
 && cd /tmp/3rd-party/ \
 && git lfs pull --include "minimac-${MINIMAC_VERSION}" \
 && git lfs pull --include "shapeit-${SHAPEIT_VERSION}" \
 && cd /tmp/3rd-party/minimac-${MINIMAC_VERSION}/ \
# minimac4 will be installed to /opt/Minimac4release-build/ 
 && unzip ./minimac4_install.zip -d /opt/ \
 && unzip ./plink2_linux_x86_64.zip -d /usr/local/bin \
 && cd /tmp/3rd-party/shapeit-${SHAPEIT_VERSION}/ \
 && tar -zxvf shapeit.v2.r904.glibcv2.17.linux.tar.gz \
# shapeit will be installed to /opt/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin
 && cp -r shapeit.v2.904.3.10.0-693.11.6.el7.x86_64 /opt/ \
# Clean up
 && rm -rf /tmp/*

# Load shell script
COPY run.sh /app/
WORKDIR /app

CMD ["bash", "run.sh"]
