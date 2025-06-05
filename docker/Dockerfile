FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    git \
    cmake \
    build-essential \
    libeigen3-dev \
    libspdlog-dev \
    libsuitesparse-dev \
    qtdeclarative5-dev \
    qt5-qmake \
    libqglviewer-dev-qt5 \
    libmetis-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone and build g2o
WORKDIR /opt
RUN git clone https://github.com/RainerKuemmerle/g2o && \
    cd g2o && \
    git checkout master && \
    sed -i 's/^  void update(const double\* update);/  virtual void update(const double* update);/' g2o/core/sparse_optimizer.h && \
    mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DG2O_USE_OPENMP=ON -DCMAKE_CXX_FLAGS="-DEIGEN_USE_THREADS -DEIGEN_USE_OPENMP" .. && \
    make -j$(nproc) && \
    make install

# Copy startup script into container
COPY docker/startup.sh /usr/local/bin/startup.sh
RUN chmod +x /usr/local/bin/startup.sh

# Run this script when container starts
WORKDIR /workspace
ENTRYPOINT ["/usr/local/bin/startup.sh"]

