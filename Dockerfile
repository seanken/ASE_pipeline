FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    openjdk-11-jdk \
    r-base \
    git \
    wget \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip3 install --no-cache-dir \
    pandas>=1.2.3 \
    pysam>=0.16.0.1 \
    numpy>=1.19.5 \
    scipy>=1.6.3 \
    pillow>=8.2.0

# Install R packages
RUN R -e "install.packages(c('data.table', 'e1071', 'ggplot2'), repos='https://cloud.r-project.org/')"

# Clone and install Monopogen
WORKDIR /opt
RUN git clone https://github.com/KChen-lab/Monopogen.git
WORKDIR /opt/Monopogen
RUN pip3 install -e .

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
