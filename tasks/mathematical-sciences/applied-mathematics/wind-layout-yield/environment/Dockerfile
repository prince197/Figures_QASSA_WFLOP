FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*
RUN pip3 install --no-cache-dir --break-system-packages numpy==2.1.3 scipy==1.14.1 pandas==2.2.3
COPY data /app/data
RUN mkdir -p /app/output
WORKDIR /app
