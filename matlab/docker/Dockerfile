# Build from the MATLAB base image
FROM mathworks/matlab:r2022b
ARG DEBIAN_FRONTEND=noninteractive
RUN sudo apt-get update -y && \
    sudo apt-get install -y \
    git wget

RUN git clone https://github.com/choilab-jefferson/LungCancerScreeningRadiomics.git
WORKDIR LungCancerScreeningRadiomics
RUN git submodule init && git submodule update
