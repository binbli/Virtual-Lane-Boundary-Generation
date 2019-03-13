#ifndef LANE_QUALITY_CHECKER
#define LANE_QUALITY_CHECKER

#include"VeloPtsProjectCam.h"
#include<fstream>

namespace LD {

    class LaneQualityChecker : public VeloPtsProjectCam {
    public:

        LaneQualityChecker(string _xmlFile);

        virtual void
        operator()(const Eigen::ArrayXXf &_intersectedPts, const Eigen::ArrayXf &_clusters, const cv::Mat &_veloPoints,
                   const Mat &_reflectivity, const Eigen::MatrixXf &_veloImg, const Mat &_inputImg,
                   const Mat &_segmentedImg, const Mat &_refinedImg, float &_reflectivityRatio,
                   float &_brightnessRatio);

    protected:

        virtual void ParseXML() override;

        virtual void
        ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints, Mat &_reflectivity,
                                 Mat &_inputImg) override {}

        string m_clustersFile;
        string m_ptsFile;
        int m_noiseID;
        std::ifstream finClusters;
        std::ifstream finPts;

    };

}

#endif
