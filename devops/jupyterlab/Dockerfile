# Base Image
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

# Set Label
LABEL maintainer="Kailash Lohar <kailash@aganitha.ai>"

# Set Working Directory
WORKDIR /home

# Install Utilities & Miniconda Conda Environment
RUN apt-get update && \
    apt-get install -y make cmake wget curl git parallel gcc g++ ocl-icd-libopencl1 opencl-headers clinfo gnupg python3 python3-pip \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda && \
    rm -rf /var/lib/apt/lists/* && rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /miniconda/bin:$PATH
COPY environment.yaml .
RUN conda env create -f environment.yaml && conda clean --all --yes
SHELL ["conda", "run", "-n", "unidock", "/bin/bash", "-c"]

# Install Boost 1.77.0
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.gz \
    && tar -xzf boost_1_77_0.tar.gz && rm boost_1_77_0.tar.gz && cd /home/boost_1_77_0 && ./bootstrap.sh && ./b2 -d2 install

# Install NVIDIA Container Runtime
RUN curl -s -L https://nvidia.github.io/nvidia-container-runtime/gpgkey | apt-key add - \
    && distribution=$(. /etc/os-release; echo $ID$VERSION_ID) \
    && curl -s -L https://nvidia.github.io/nvidia-container-runtime/$distribution/nvidia-container-runtime.list | tee /etc/apt/sources.list.d/nvidia-container-runtime.list 
RUN mkdir -p /etc/OpenCL/vendors && echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# Set IOPub message rate limit
ENV JUPYTER_IOPUB_MSG_RATE_LIMIT=10000

# Copy and Compile Uni-Dock
COPY /unidock /home
RUN cmake -B build && cmake --build build -j`nprocs` && cmake --install build 

# Set permissions for files
RUN chmod -R 777 /home

WORKDIR /home

# Expose Port for JupyterLab
EXPOSE 4000

# Set Default Command for starting JupyterLab
CMD ["conda", "run", "-n", "unidock", "jupyter", "lab", "--ip=0.0.0.0", "--port=4000", "--no-browser", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]
