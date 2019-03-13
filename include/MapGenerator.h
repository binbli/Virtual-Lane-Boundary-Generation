#ifndef EVALUTIONPROJECT_MAPGENERATOR_H
#define EVALUTIONPROJECT_MAPGENERATOR_H

#include "../include/Solver.h"
#include "../3rdParty/TLinkage/BSplineTLinkage.h"
#include "../include/Utilities.h"
#include "../include/LaneMarkAssessor.h"
#include "../include/PointsVisualizer.h"

namespace LD {
    class MapGenerator : public Solver {

    public:

        enum LaneType {
            LEFT, CENTER, RIGHT
        };

        MapGenerator(string _file);

        virtual void Run() override;

        void Project(const Eigen::ArrayXXf &_veloPoints, const Eigen::MatrixXf &_rotation,
                     const Eigen::VectorXf &_translation, Eigen::ArrayXXf &_newVeloPoints);

        void Project(const Eigen::ArrayXf &_veloPoints, const Eigen::MatrixXf &_rotation,
                     const Eigen::VectorXf &_translation, Eigen::ArrayXf &_newVeloPoints);

        void GetWorldCtrlPts(const string &_imgName, const vector<ArrayXf> &_models, vector<ArrayXXf> &_worldCtrlPts);

        void
        TransformCtrlPtsToNext(const Eigen::ArrayXf &_model, const string &_curImgFileName, Eigen::ArrayXXf &_ctrlPts);

        void TransformCtrlPtsToNext(const Eigen::ArrayXXf &_curCtrlPts, const string &_curImgFileName,
                                    Eigen::ArrayXXf &_nextCtrlPts);

        void GenerateCenterCtrlPts(ArrayXf _centerLaneModel, const string &_imgFileName, ArrayXXf &_ctrlPts,
                                   ArrayXXf &_prectrilPts, bool &_isVirtual, Matrix3f Rotation2G,
                                   RowVector3f translation2G, MatrixXf RoadEdgePts, ArrayXXf &GpsPts,
                                   ArrayXXf &mapgpsPtsL);

        void PrintPassedWorldCtrlPts(const ArrayXXf &_ctrlPts, const string &_imgBaseName, LaneType _type);

        void vizLCCurve(ArrayXXf _centerCtrlPts, ArrayXXf GpsPts, ArrayXXf mapgpsPtsL, bool _isCenterVirtual,
                        const cv::Mat _inputImg, string _imgBaseName);

        void loadPriorGPS(std::string _path, ArrayXXf &_rawMapGPSMat);

        void priorMapGPsInL(ArrayXXf &_rawMapGPSMat, Eigen::MatrixXd gTruthGPS, Eigen::ArrayXXf &priorGPSInL);

        class noLCCgenerated {
        };

    protected:

        virtual void ParseXML();

        void ParseVioFile(vector<Eigen::ArrayXf> &_transformationInfo);

        void
        GetRotationTranslation(const string &_imgFileName, Eigen::MatrixXf &_rotation, Eigen::VectorXf &_translation);

        void GetRotationTranslation(const Eigen::ArrayXf &_trasformationVec, Eigen::MatrixXf &_rotation,
                                    Eigen::VectorXf &_translation);

        void GetWorldCtrlPts(const ArrayXf &_model, const MatrixXf &_rotation, const VectorXf &_translation,
                             ArrayXXf &_worldCtrlPts);

        void PrintWorldCtrlPts(const ArrayXXf &_ctrlPts, LaneType _type);

        void printCostoFile(string imgBaseName, int costMetricMat, RowVectorXf standardizedCostrow);

        double degreeToRadian(const double degree) { return (degree * PI / 180); };

        double radianToDegree(const double radian) { return (radian * 180 / PI); };

        bool loadGPSFile(std::string _path, std::string _varname, Eigen::MatrixXd &_outputMat, int _nums);

        double CoordinatesToMeters(double latitude1, double longitude1, double latitude2, double longitude2);

        double CoordinatesToAngle(double latitude1, const double longitude1, double latitude2, const double longitude2);

        std::pair<double, double>
        CoordinateToCoordinate(double latitude, double longitude, double angle, double meters);

        bool GPSCoordidateTran(std::string _path, std::string _fimeName, int _nums, double &_cVhcesd,
                               Eigen::MatrixXd &_rawGPSfMat,
                               Eigen::MatrixXd &_GPSCoordinatInL);

        void
        GPScuve(const std::string _GPSDirRoot, const std::string _fName, const int _nCtrlPts, double &_cVhcesd,
                Eigen::MatrixXd &_rawGPSfMat,
                Eigen::MatrixXd &_GPSCoordinateInL, Eigen::MatrixXf &_GPSpline, double &_goalPosAngle);

        bool CandidateCurvSelector(Eigen::MatrixXf _GPSpline, double _goalPosAngle, double _NumfLCCCurves,
                                   double _minLCCCurveLen, Eigen::MatrixXf _RoadEdgePts,
                                   Eigen::ArrayXXf &_LCCostSum, Eigen::MatrixXf &_OptLCCSpline, double _cVhcesd,
                                   ArrayXXf infeasibleLCC, RowVectorXf &optimalLCCost, ArrayXXf &_prectrlPts);

        void GetPassedWorldCtrlPts(const ArrayXXf &_curCtrlPts, const string &_curImgFileName,
                                   ArrayXXf &_passedWorldCtrlPts);


        string m_vioFile;
        string m_splineModelsFile;
        string m_gpsimudir;
        string m_CurbRTRoot;
        string m_CurbRTinFile;
        string m_leftCtrlPtsFile, m_centerCtrlPtsFile, m_rightCtrlPtsFile;
        string m_dataRoot, m_dataFile, m_mapgpsdir;
        string m_vizImgPrefix;
        string m_costFile;
        double m_carWidth;
        double m_minX;
        double m_maxX;
        double m_LCC2orign;
        bool m_saveVizImg;
        int m_NumfCurves;
        int m_minSamples;
        int m_IntersecGPSample;
        float m_startAngleThreshold;
        float m_UpperElevation;
        float m_steplen;
        float m_laneCenterWidth;
        float m_InterSectionTrigger, m_IntersectionLen;
        float m_maxAngle, m_minAngle, m_TargetAngleDiff;
        float m_weightfs, m_weightfo, m_weightfg, m_weightfc, m_weightfa;
        const double m_earthDiameterMeters = 6365.998 * 2 * 1000; //TODO: table for earthdiameter

        vector<ArrayXf> m_transformationInfo;
        ArrayXXf m_rawMapGPSMat;
        BSplineTLinkage m_bSplineTLinkage;
        LaneMarkAssessor m_LaneMarkAssessor;
        PointsVisualizer m_visualizer;
        std::ofstream m_foutCenterCtrlPts, m_foutLeftCtrlPts, m_foutRightCtrlPts, m_fcostCenterMetric;

    };

}

#endif //EVALUTIONPROJECT_MAPGENERATOR_H
