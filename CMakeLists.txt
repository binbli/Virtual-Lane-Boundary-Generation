cmake_minimum_required(VERSION 3.5.1)
project(Evalutionproject)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/3rdParty/ICNet/PSPNet/cmake/Modules)

set(CMAKE_CXX_STANDARD 14)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-g0")


find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Glog REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(ALGLIB REQUIRED)
find_package(PugiXML REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(OpenCV 2.4.9.1 REQUIRED)
find_package(Boost COMPONENTS regex system filesystem REQUIRED)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(PYTHON_INCLUDE /usr/local/include/python3.5)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)
link_directories(${PROJECT_SOURCE_DIR}/3rdParty/ICNet/PSPNet/build/lib)

include_directories(
        ${GSL_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${GLOG_INCLUDE_DIR}
        ${EIGEN3_INCLUDE_DIR}
        ${ALGLIB_INCLUDES}
        ${PUGIXML_INCLUDE_DIR}
        ${OpenCV_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${PYTHON_INCLUDE}
        ${PYTHON_INCLUDE_DIRS}
        ${PROJECT_SOURCE_DIR}/3rdParty/ICNet/PSPNet/include
        ${PROJECT_SOURCE_DIR}/3rdParty/ICNet/PSPNet/build/src
        ${PROJECT_SOURCE_DIR}/3rdParty/LSD
        ${PROJECT_SOURCE_DIR}/3rdParty/TLinkage
        ${PROJECT_SOURCE_DIR}/3rdParty/DBScan
        ${PROJECT_SOURCE_DIR}/include
)

SET(GCC_COVERAGE_COMPILE_FLAGS "-Wno-sign-compare -Wno-uninitialized -MMD -MP -pthread -fPIC -g -O2")

add_definitions(
        -DCAFFE_VERSION=1.0.0-rc3
        -DUSE_OPENCV
        -DUSE_LEVELDB
        -DUSE_LMDB
        -DCPU_ONLY
)

add_executable(
        LDRT
        src/main.cpp
        src/LaneDetector.cpp
        src/RoadSegment.cpp
        src/LMsfromCam.cpp
        src/KPercentExtractor.cpp
        src/VeloPtsProjectCam.cpp
        src/LMsintersection.cpp
        src/SurfaceDataMaker.cpp
        src/Utilities.cpp
        src/PointsVisualizer.cpp
        src/LaneQualityChecker.cpp
        src/MapGenerator.cpp
        src/LaneMarkAssessor.cpp
        3rdParty/TLinkage/TLinkage.cpp
        3rdParty/TLinkage/BSplineTLinkage.cpp
        3rdParty/TLinkage/KLargestOR.cpp
        3rdParty/TLinkage/MaxDiffOR.cpp
        3rdParty/TLinkage/UniformSampler.cpp
        3rdParty/TLinkage/DistBasedSampler.cpp
        3rdParty/DBScan/DBScan.cpp
        3rdParty/LSD/LSDLines.cpp
)

target_link_libraries(LDRT
        ${OpenCV_LIBS}
        ${PUGIXML_LIBRARY}
        ${ALGLIB_LIBRARIES}
        ${GSL_LIBRARIES}
        ${Boost_LIBRARIES}
        ${GLOG_LIBRARY}
        ${PYTHON_LIBRARIES}
        caffe
)
