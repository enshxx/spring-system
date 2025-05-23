FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

RUN chmod +x build.sh

RUN mkdir -p build && cd build

RUN ./build.sh

WORKDIR /app/build

CMD ["./gauss_newton"]