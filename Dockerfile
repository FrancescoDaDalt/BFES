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
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy project files
COPY src/ src/
COPY CMakeLists.txt .
COPY main_benchmarking.cpp .

# Create build directory
RUN mkdir -p build
WORKDIR /app/build

# Configure and build project
RUN cmake .. && make -j$(nproc)

# Set default command
CMD ["./BFES"]