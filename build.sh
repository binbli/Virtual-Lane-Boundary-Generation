#!/usr/bin/env bash
cp 3rdParty/ICNet/PSPNet/Makefile.config.example 3rdParty/ICNet/PSPNet/Makefile.config
patch -p0 -i 3rdParty/patches/ICNet.patch

cd 3rdParty/ICNet/PSPNet
make
cd -

#patch -p0 -i 3rdParty/patches/PythonRobotics.patch

mkdir build && cd build
cmake ..
make
