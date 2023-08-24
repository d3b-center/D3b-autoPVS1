FROM ubuntu:20.04
LABEL maintainer="Miguel Brown <brownm28@chop.edu>"
LABEL description="Docker file with modified version of autopvs1"

ENV REPO_RELEASE=1.0.1

RUN apt update -y && apt install -y curl wget python3 python3-pip pigz libbz2-dev liblzma-dev zlib1g-dev
RUN pip3 install pysam==0.17.0 pyfaidx
RUN mkdir autopvs1
ADD ./ autopvs1/
RUN apt remove -y wget curl

COPY Dockerfile .