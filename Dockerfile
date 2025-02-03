FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    cmake \
    g++ \
    clang \
    libomp-dev \
    libeigen3-dev \
    libgsl-dev \
    wget \
    unzip \
    openjdk-17-jdk \
    maven \
    nano \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy project files
COPY src/ src/
COPY ApacheFlink_x_Java/ ApacheFlink_x_Java/
COPY CMakeLists.txt .
COPY main_benchmarking.cpp .

# Download and install Apache Flink
RUN FLINK_VERSION=1.20.0 && \
    FLINK_URL="https://archive.apache.org/dist/flink/flink-${FLINK_VERSION}/flink-${FLINK_VERSION}-bin-scala_2.12.tgz" && \
    wget $FLINK_URL && \
    tar -xzf flink-${FLINK_VERSION}-bin-scala_2.12.tgz && \
    mv flink-${FLINK_VERSION} /opt/flink && \
    rm flink-${FLINK_VERSION}-bin-scala_2.12.tgz

# Set FLINK_HOME and PATH environment variables
ENV FLINK_HOME=/opt/flink
ENV PATH=$FLINK_HOME/bin:$PATH


# Create build directory
RUN mkdir -p build
WORKDIR /app/build

# Configure and build project
RUN cmake .. && make -j$(nproc)

# Set default command
CMD ["/bin/bash"]