#!/usr/bin/env bash
sudo apt-get install --no-install-recommends    \
  build-essential                               \
  cmake                                         \
  git                                           \
  wget                                          \
  libatlas-base-dev   `#For Caffe`              \
  libboost-all-dev                              \
  libgflags-dev                                 \
  libgoogle-glog-dev                            \
  libhdf5-serial-dev                            \
  libleveldb-dev                                \
  liblmdb-dev                                   \
  libopencv-dev                                 \
  libprotobuf-dev                               \
  libsnappy-dev                                 \
  protobuf-compiler                             \
  python3-pip                                   \
  python3-dev                                   \
  libmatio-dev        `#MAT File I/O`           \
  libpugixml-dev      `#XML processing`         \
  libalglib-dev       `#Numerical analysis`     \
  libgsl-dev          `#GNU scientific library` \
  libopenblas-dev     `#Linear algebra library` \
  libeigen3-dev                                 \
  libpcap-dev         `#Packet capture library`

sudo pip3 install numpy matplotlib pandas scipy pybind11
