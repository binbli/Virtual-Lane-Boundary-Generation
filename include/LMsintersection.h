#ifndef LMs_INTERSECTION_H_
#define LMs_INTERSECTION_H_

#include"../include/VeloPtsProjectCam.h"
#include"../include/Utilities.h"
#include "../include/KPercentExtractor.h"

namespace LD {
    class LMsintersection : public VeloPtsProjectCam{

    public:
        using BaseLD::m_debug;

        static Eigen::MatrixXf H_p;

        LMsintersection(string _xmlFile) : VeloPtsProjectCam(_xmlFile) {
            ParseXML();
        }

        void
        operator()(const Mat &_veloPoints, const Mat &_segImg, const Mat &_refinedImg, Eigen::ArrayXXf &_intersectedPts,
                   Mat &_reflectivity, Eigen::MatrixXf &_veloImgPoints, const Mat &_inputImg,
                   const string m_imgBaseName, Eigen::Matrix3f &RotationL2G, Eigen::RowVector3f &translationL2G,
                   Eigen::MatrixXf &refineCurbSet, cv::Mat &roadMaskImg);

        bool isMode2D() { return m_printOnly2D; }

        void RoadEdgePtsDetector(const Eigen::MatrixXf _veloImg, const Mat _veloPoints, const Mat _segImg,
                                 Eigen::MatrixXf &refineCurbSet, Eigen::Matrix3f &RotationL2G,
                                 Eigen::RowVector3f &translationL2G, cv::Mat &roadMaskImg);

    protected:

        string m_segRoot;
        string m_refinedRoot;
        string m_refImgPrefix;
        string m_segImgPrefix;
        string m_vizImgPrefix;
        string m_outputFile;
        string m_CurbRTOutFile;
        string m_CurbRTRoot;
        std::ofstream m_fout;
        bool m_saveVizImg;
        bool m_printOnly2D;
        int m_maxWidth, m_maxLength;
        float m_dThrehold;
        int m_nIter;
        int m_minSamples;
        float m_percentInliers;
        float m_radius;
        float m_height;
        float m_curbH;
        float m_surfCurvePara;

        virtual void
        ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints, Mat &_reflectivity,
                                 Mat &_inputImg) override;

        void IntersectIn3D(const Eigen::MatrixXf _veloImg, const Mat &_veloPoints, const Mat &_reflectivity,
                           const Mat &_refinedImg, const double &_thresh, Eigen::ArrayXXf &_intersectedImg,
                           Mat &_vizImg);

        void
        PrintEdgePtsRT(Eigen::MatrixXf refineCurbSet, Eigen::Matrix3f RotationL2G, Eigen::RowVector3f translationL2G);

        void PrintToFile(Eigen::ArrayXXf &_intersectedPts);

        double OtsuThresholdRoad(const Eigen::MatrixXf _veloImg, const Mat &_segImg, const Mat &_reflectivity);

        virtual void ParseXML() override;

        void  getRoadEdgePoints(Eigen::MatrixXf roadlidarpts, Eigen::MatrixXf _model, Eigen::MatrixXd &_curbPointSet);

        void edgePtsRefiner(const Eigen::MatrixXd _curbPointSet, Eigen::MatrixXd &_refineCurbSet, float radius);

        void computeRoadEdgePts(Eigen::MatrixXf roadlidarpts, Eigen::MatrixXf _model, float radius,
                                Eigen::MatrixXd &refineCurbSet);

        void shortestDistance(Eigen::Vector3d point, Eigen::MatrixXd _set, float radius, Eigen::MatrixXd &_nearestPoints, std::vector<double> &distance);

    };
}

#endif
