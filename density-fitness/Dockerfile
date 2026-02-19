# --- STAGE 1: Build ---
FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install all build-time dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    catch2 \
    nlohmann-json3-dev \
    libeigen3-dev \
    libccp4-dev \
    libclipper-dev \
    libgsl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    libboost-regex-dev \
    libicu-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp/build

# Clone and build density-fitness
RUN git clone --depth 1 --branch v1.1.4 https://github.com/PDB-REDO/density-fitness.git && \
    cd density-fitness && \
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build -j$(nproc) && \
    cmake --install build && \
    strip /usr/local/bin/density-fitness

# --- STAGE 2: Runtime ---
FROM ubuntu:24.04

LABEL org.opencontainers.image.description="density-fitness: Calculate electron density statistics for X-ray structures"

ENV DEBIAN_FRONTEND=noninteractive

# Install only essential runtime libraries in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    libccp4c0 \
    libclipper2 \
    libgsl27 \
    libxml2 \
    libicu74 \
    libboost-regex1.83.0 \
    zlib1g \
    libbz2-1.0 \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy binary from builder
COPY --from=builder /usr/local/bin/density-fitness /usr/local/bin/density-fitness

# Copy only necessary dictionary files from builder
COPY --from=builder /usr/local/share/libcifpp /usr/local/share/libcifpp

# Download and setup data files
RUN mkdir -p /var/cache/libcifpp && \
    curl -fsSL -o /var/cache/libcifpp/components.cif \
        https://files.wwpdb.org/pub/pdb/data/monomers/components.cif && \
    curl -fsSL -o /var/cache/libcifpp/mmcif_pdbx.dic \
        https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic && \
    chmod 644 /var/cache/libcifpp/* && \
    # Remove curl after downloading to save space
    apt-get purge -y curl ca-certificates && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

ENV LIBCIFPP_DATA_DIR=/var/cache/libcifpp

WORKDIR /data
ENTRYPOINT ["density-fitness"]
