#!/usr/bin/env bash
#f=2011_09_26_drive_0035_sync
# prepKITTI.sh $f

if [ ! -f 3rdParty/ICNet/PSPNet/evaluation/model/*.caffemodel ]; then 
  cd 3rdParty/ICNet/PSPNet/evaluation/model/
  #wget https://github.tamu.edu/IvanWang/LaneRecognitionTracking/releases/download/Release/PSPNet-models.tar.gz
  #tar xzvf PSPNet-models.tar.gz
  cd -
fi


if [ ! -d demo/calib/${1:0:10}/ ]; then 
  cd demo/calib
  wget https://s3.eu-central-1.amazonaws.com/avg-kitti/raw_data/${1:0:10}_calib.zip
  unzip ${1:0:10}_calib.zip
  mv ${1:0:10} ${1}_calib
  cd -
fi


cd demo/data
mkdir $1 && cd $1
touch priorMap.txt
touch vio.txt

mkdir -p output/Images/segment_img
mkdir -p output/Images/overlayed_img
mkdir -p output/Images/lane_from_img
mkdir -p output/Images/roadBoundary_img
mkdir -p output/Images/intersect_img_lidar
mkdir -p output/Images/splineLane
mkdir -p output/Lidar/SplineModel
mkdir -p output/Lidar/SplineCluster
mkdir -p output/Lidar/curbsRT
mkdir -p output/Lanes

wget -c https://s3.eu-central-1.amazonaws.com/avg-kitti/raw_data/${1:0:(-5)}/$1.zip
unzip $1.zip
ls ${1:0:10}/$1/image_02/data|head -n -20 > seq.txt
sed -e "s/\$f0/${1:0:10}/g" -e "s/\$f1/${1}/g" ../../Parameter-template.xml > $1.xml

