## Camera-LIDAR Bimodal Lane Detection, Recognition and Tracking

## License

Our code is released under a [GPLv3 license](http://telerobot.cs.tamu.edu/lane/). In this code, we aim to recognize, generate and track virtual lane boundaries to help navigate autonomous vehicles in urban scenarios. It has been tested on Ubuntu 16.04LTS and should work on 16.04 or newer versions. We extend our previous work in

    @INPROCEEDINGS{li2018lane, 
        author={B. {Li} and D. {Song} and H. {Li} and A. {Pike} and P. {Carlson}}, 
        booktitle={2018 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)}, 
        title={Lane Marking Quality Assessment for Autonomous Driving}, 
        year={2018}, 
        pages={1-9}, 
        doi={10.1109/IROS.2018.8593855}, 
        ISSN={2153-0866}, 
        month={Oct}
    }
    
If you use our code in an academic work, please cite the aforementioned paper.

## Clone
```
git clone --recurse-submodules git@github.com:bli-tamu/LDRT.git
```

## Dependencies

You can use ``./install-dependencies.sh`` for convenience, or do it manually as follows:
    
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

    pip3 install numpy matplotlib pandas scipy pybind11

## Build Instruction

You can use ```./build.sh``` for default build, or do it manually as follows:

* ICNet

    Copy the example configuration file and apply patch.
    ```
    cp 3rdParty/ICNet/PSPNet/Makefile.config.example 3rdParty/ICNet/PSPNet/Makefile.config
    patch -p0 -i 3rdParty/patches/ICNet.patch
    ```
    Modified the configuration file before compile if needed.
    ```
    cd 3rdParty/ICNet/PSPNet
    make
    ```
* This project uses CMake (http://www.cmake.org), a cross-platform build system. 
    ```
    mkdir build && cd build
    cmake ..
    make
    ```

## External component
* Install [VINS-Fusion](https://github.com/HKUST-Aerial-Robotics/VINS-Fusion.git) to generate visual odometry ``vio.txt`` from the dataset, and put it under ``demo/data/datasetName`` folder.
* Download GPS way points from Google maps, utilize the script ``matlab/process_KITTIdata.m`` to generate the map, and put the data in the ``demo/data/datasetName/priorMap.txt``. 

## Prepare to run a demo

1. We use [KITTI dataset](http://www.cvlibs.net/datasets/kitti/raw_data.php?type=city) to demostrate our code, ``prepDemo.sh`` is a script that downloads and prepares the demo. For example, ``./prepDemo.sh 2011_09_26_drive_0056_sync`` will automatically prepare the PSPNet models, 2011_09_26_drive_0056_sync dataset and its calibration files.
2. Run the dataset with VINS-Fusion and place the ``vio.txt`` in ``demo/data/datasetName``.

## Run a demo
* Run with ``runDemo.sh ``, for example, ``./runDemo.sh 2011_09_26_drive_0056_sync`` will run the code and keep a log file in ``demo/data/2011_09_26_drive_0056_sync/``.


## Contact

1. Binbin Li <binbinli@tamu.edu>
2. Ankit Ramchandani <ankit61@tamu.edu>
3. Di Wang <ivanwang@tamu.edu>
4. Aaron Kingery <aaronmkingery@tamu.edu>
5. Aaron Angert <adangert@tamu.edu>
6. Dezhen Song <dzsong@cse.tamu.edu>
