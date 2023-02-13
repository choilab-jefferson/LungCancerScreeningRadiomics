# Build from the MATLAB base image
FROM mathworks/matlab-deep-learning:r2022b
ARG DEBIAN_FRONTEND=noninteractive
RUN sudo apt-get update -y && \
    sudo apt-get install -y --no-install-recommends apt-utils build-essential sudo git ca-certificates
RUN git clone https://github.com/mkazhdan/ConformalizedMCF.git && cd ConformalizedMCF && make || true
RUN sudo mv ConformalizedMCF/Bin/Linux/ConformalizedMCF /usr/bin
RUN git clone -b MSKLung https://github.com/choilab-jefferson/LungCancerScreeningRadiomics.git
WORKDIR LungCancerScreeningRadiomics
RUN git submodule init && git submodule update
