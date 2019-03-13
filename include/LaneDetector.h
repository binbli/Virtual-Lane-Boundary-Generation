#ifndef LANE_DETECTOR_H_
#define LANE_DETECTOR_H_

#include"Solver.h"
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/ml/ml.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<eigen3/Eigen/Dense>
#include <chrono>
#include"RoadSegment.h"
#include"KPercentExtractor.h"
#include"LMsintersection.h"
#include"../3rdParty/TLinkage/BSplineTLinkage.h"
#include"LaneQualityChecker.h"
#include"PointsVisualizer.h"
#include "MapGenerator.h"

namespace LD {

    struct LaneAssessmentParams {
        float brightnessRatio = 0;
        float reflectivityRatio = 0;
        float areaDiff = 0;
        float avgDist = 0;
    };


    class LaneDetector : public Solver {
    public:
        void
        operator()(const cv::Mat &_inputImg, cv::Mat &segImg, const cv::Mat &_veloPoints, const string _imgBaseName,
                   vector<Eigen::ArrayXXf> &_ctrlPts,
                   vector<bool> &_isVirtual, LaneAssessmentParams &_assessParams, Eigen::Matrix3f &RotationG,
                   Eigen::RowVector3f &translationG, Eigen::MatrixXf &RoadEdgePts,
                   vector<Eigen::ArrayXXf> &previousCtrlPts, bool &isCenterVirtual, ArrayXXf &centerCtrlPts,
                   Eigen::ArrayXXf &GpsPts, Eigen::ArrayXXf &mapgpsPtsL, cv::Mat &roadMaskImg);

        virtual void Run() override;

        LaneDetector(string _xmlFile);

    protected:
        RoadSegment m_RoadSegment;
        KPercentExtractor m_LMsfromCam;
        LMsintersection m_LMsintersection;
        BSplineTLinkage m_bSplineTLinkage;
        LaneQualityChecker m_laneQualityChecker;
        PointsVisualizer m_visualizer;
        MapGenerator m_mapGenerator;
        KPercentExtractor m_KPercentExtractor;

        string m_imgFile, m_imgRoot, m_veloRoot;
        string m_ratiosFile;
        string m_imgBaseName;
        string m_isCenterVirtual;
        bool m_saveVizImg;
        string m_vizImgPrefix;

        virtual void ParseXML() override;
    };
}

#endif