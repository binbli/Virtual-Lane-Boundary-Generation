#include "../include/MapGenerator.h"
#include "../include/Utilities.h"
#include <algorithm>

namespace LD {

    using namespace Eigen;

    MapGenerator::MapGenerator(string _file) : Solver(_file), m_bSplineTLinkage(_file), m_LaneMarkAssessor(_file),
                                               m_visualizer(_file) {
        ParseXML();
        ParseVioFile(m_transformationInfo);
        loadPriorGPS(m_mapgpsdir, m_rawMapGPSMat);
    }

    void MapGenerator::ParseXML() {

        if (m_debug)
            cout << "Entering MapGenerator::ParseXML()" << endl;

        m_xml = m_xml.child("MapGenerator");
        pugi::xml_node solverInstance = m_xml.child("SolverInstance");

        m_dataRoot = solverInstance.attribute("dataRoot").as_string();
        m_dataFile = solverInstance.attribute("dataFile").as_string();
        m_saveVizImg = solverInstance.attribute("saveVizImg").as_bool();
        m_vioFile = solverInstance.attribute("vioFile").as_string();
        m_splineModelsFile = solverInstance.attribute("splineModelsFile").as_string();
        m_mapgpsdir = solverInstance.attribute("mapgpsdir").as_string();
        m_gpsimudir = solverInstance.attribute("gpsimudir").as_string();
        m_vizImgPrefix = solverInstance.attribute("vizImgPrefix").as_string();

        m_CurbRTRoot = solverInstance.attribute("CurbRTRoot").as_string();
        m_CurbRTinFile = solverInstance.attribute("CurbRTinFile").as_string();

        m_leftCtrlPtsFile = solverInstance.attribute("leftCtrlPtsFile").as_string();
        m_rightCtrlPtsFile = solverInstance.attribute("rightCtrlPtsFile").as_string();
        m_centerCtrlPtsFile = solverInstance.attribute("centerCtrlPtsFile").as_string();
        m_costFile = solverInstance.attribute("costFile").as_string();

        m_minSamples = solverInstance.attribute("minSamples").as_int();
        m_minX = solverInstance.attribute("minX").as_double();
        m_maxX = solverInstance.attribute("maxX").as_double();
        m_steplen = solverInstance.attribute("steplen").as_double();
        m_LCC2orign = solverInstance.attribute("LCC2orign").as_double();

        m_NumfCurves = solverInstance.attribute("NumfCurves").as_int();
        m_carWidth = solverInstance.attribute("carWidth").as_double();
        m_laneCenterWidth = solverInstance.attribute("laneCenterWidth").as_float();
        m_startAngleThreshold = solverInstance.attribute("startAngleThreshold").as_float();
        m_UpperElevation = solverInstance.attribute("UpperElevation").as_float();
        m_InterSectionTrigger = solverInstance.attribute("InterSectionTrigger").as_float();
        m_IntersectionLen = solverInstance.attribute("IntersectionLen").as_float();
        m_IntersecGPSample = solverInstance.attribute("IntersecGPSample").as_int();

        m_maxAngle = solverInstance.attribute("maxAngle").as_float();
        m_minAngle = solverInstance.attribute("minAngle").as_float();
        m_TargetAngleDiff = solverInstance.attribute("TargetAngleDiff").as_float();

        m_weightfs = solverInstance.attribute("weightfs").as_float();
        m_weightfo = solverInstance.attribute("weightfo").as_float();
        m_weightfg = solverInstance.attribute("weightfg").as_float();
        m_weightfc = solverInstance.attribute("weightfc").as_float();
        m_weightfa = solverInstance.attribute("weightfa").as_float();

        if (m_dataRoot.empty() || m_dataFile.empty() || m_vioFile.empty() || m_splineModelsFile.empty() ||
            m_mapgpsdir.empty() || m_vizImgPrefix.empty() || m_leftCtrlPtsFile.empty() ||
            m_CurbRTRoot.empty() || m_CurbRTinFile.empty() || m_rightCtrlPtsFile.empty() || m_costFile.empty() ||
            m_gpsimudir.empty() ||
            m_centerCtrlPtsFile.empty() || !m_minSamples || !m_startAngleThreshold || !m_maxX || !m_NumfCurves ||
            !m_carWidth || !m_laneCenterWidth || !m_IntersectionLen || !m_IntersecGPSample)

            throw runtime_error(
                    "missing values in one of the following attributes: vioFile, splineModelFile, gpsimudir, minSamples, rightCtrlPtsFile, leftCtrlPtsFile, centerCtrlPtsFile, startAngleThreshold, maxCurvatureThreshold, minX, maxX, minSamples, numfCurves, carWidth, curbRTRoot, curbRTOutFile");

        if (m_debug)
            cout << "Exiting MapGenerator::ParseXML()" << endl;

    }

    void
    MapGenerator::Project(const ArrayXXf &_veloPoints, const Eigen::MatrixXf &_rotation, const VectorXf &_translation,
                          ArrayXXf &_newVeloPoints) {
        if (m_debug)
            cout << "Entering MapGenerator::Project()" << endl;

        //_veloPoints expected to have shape of 3 x n, where  n are the number of points
        _newVeloPoints = (((_rotation) * _veloPoints.matrix()).colwise() + _translation).array();

        if (m_debug)
            cout << "Exiting MapGenerator::Project()" << endl;
    }

    void MapGenerator::Project(const ArrayXf &_veloPoint, const MatrixXf &_rotation, const VectorXf &_translation,
                               ArrayXf &_newVeloPoint) {
        if (m_debug)
            cout << "Entering MapGenerator::Project()" << endl;

        _newVeloPoint = ((_rotation * _veloPoint.matrix()) + _translation).array();

        if (m_debug)
            cout << "Exiting MapGenerator::Project()" << endl;
    }

    void MapGenerator::Run() {

        if (m_debug)
            cout << "Entering MapGenerator::Run()" << endl;

        Matrix3f RotationG;
        RowVector3f translationG;
        MatrixXf roadEdgePtsSet;
        ArrayXXf centerCtrlPts, prectrlPts;
        ArrayXXf GpsPts, mapgpsPtsL;
        ArrayXf DummyMat;
        vector<ArrayXf> models;
        vector<ArrayXXf> ctrlPts;
        bool isCenterVirtual;
        vector<bool> isVirtual;
        string _imgBaseName;
        std::ifstream inlist(m_dataFile.c_str());
        if (!inlist)
            throw runtime_error("Cannot open " + m_dataFile);
        string line;
        DummyMat.setZero();

        while ((std::getline(inlist, _imgBaseName))) {

            cv::Mat inputImg = cv::imread(m_dataRoot + "/" + _imgBaseName);
            if (inputImg.empty())
                std::cerr << "could not open or find " << m_dataRoot + "/" + _imgBaseName << endl;
            if (m_debug)
                cout << "Successfully opened " << _imgBaseName << endl;

            string SplineModalFilePah = m_splineModelsFile + _imgBaseName.substr(0, line.size() - 4) + ".txt";
            std::ifstream finModels(SplineModalFilePah.c_str());

            if (finModels.is_open()) {
                string temp;
                std::getline(finModels, temp);
                ReadEigenArrayVecFromFile(finModels, models);
            }

            if (!models.empty())
                GenerateCenterCtrlPts(models.back(), _imgBaseName, centerCtrlPts, prectrlPts, isCenterVirtual,
                                      RotationG,
                                      translationG, roadEdgePtsSet, GpsPts, mapgpsPtsL);
            else
                GenerateCenterCtrlPts(DummyMat, _imgBaseName, centerCtrlPts, prectrlPts, isCenterVirtual, RotationG,
                                      translationG, roadEdgePtsSet, GpsPts, mapgpsPtsL);

            vizLCCurve(centerCtrlPts, GpsPts, mapgpsPtsL, isCenterVirtual, inputImg, _imgBaseName);

            PrintPassedWorldCtrlPts(centerCtrlPts, _imgBaseName, LaneType::CENTER);

        }

        if (m_debug)
            cout << "Exiting MapGenerator::Run()" << endl;
    }

    void MapGenerator::PrintPassedWorldCtrlPts(const ArrayXXf &_ctrlPts, const string &_imgBaseName, LaneType _type) {
        ArrayXXf worldCtrlPts;
        GetPassedWorldCtrlPts(_ctrlPts, _imgBaseName, worldCtrlPts);
        PrintWorldCtrlPts(worldCtrlPts.transpose(), _type);
    }

    void
    MapGenerator::GenerateCenterCtrlPts(ArrayXf _centerLaneModel, const string &_imgFileName, ArrayXXf &_ctrlPts,
                                        ArrayXXf &_prectrilPts,
                                        bool &_isVirtual, Matrix3f Rotation2G,
                                        RowVector3f translation2G, MatrixXf RoadEdgePts, ArrayXXf &gpsPtsL,
                                        ArrayXXf &priorGPSfcarInL) {
        if (m_debug)
            cout << "Entering MapGenerator::GenerateCenterCtrlPts()" << endl;

        MatrixXf GPSpline, _OptfLCCurveInG;
        MatrixXd rawGPsMat, GPScoordinateInL;
        ArrayXf CenterModel;
        RowVectorXf optimalLCCost;
        int CenterVirtul = 1;
        bool isLCCgenerated;
        double _cVhcesd = 0.0, _StartPosAngle, _goalPosAngle, _LCCurveLen = 0.0, _minLCCCurveLen, _LCCCurvatureMax;
        ArrayXXf localCenterCurveL, LCCost, GpsPts, priorGPSInL, mapgpsPtsL, CenterCurve, infeasibleLCC, costMetricMat;

        std::string curbRTPath =
                m_CurbRTRoot + m_CurbRTinFile + _imgFileName.substr(0, _imgFileName.size() - 4) + ".txt";
        std::string _imgName = _imgFileName.substr(0, _imgFileName.size() - 3);


        // if has center LCC
        if (!_centerLaneModel.isZero()) {
            m_LaneMarkAssessor.operator()(_centerLaneModel, _StartPosAngle, _LCCurveLen,
                                          _LCCCurvatureMax,
                                          CenterCurve); //computer LCC _StartPosAngle, _LCCurveLen, _LCCCurvatureMax
        }

        loadGPSFile(m_gpsimudir, _imgFileName.substr(0, _imgFileName.size() - 3), rawGPsMat, 5);

        //load map gps in L
        priorMapGPsInL(m_rawMapGPSMat, rawGPsMat, priorGPSInL);


        if (rawGPsMat.rows() != 0)
            _cVhcesd = rawGPsMat.col(4).mean();

        if (_cVhcesd > 8)
            _minLCCCurveLen = std::pow(_cVhcesd, 2) / (1.5 * 9.8094); // safe stopping distance
        else
            _minLCCCurveLen = m_minX;

        if (m_debug) {
            cout << "The LCC has _StartPosAngle: " << fabs(_StartPosAngle) << ", the maximum allowed angle: "
                 << m_startAngleThreshold << "\nLCC curve lengh: " << _LCCurveLen
                 << ", and the minimum allowed length: " << _minLCCCurveLen << "\nMax LCC curvature: "
                 << _LCCCurvatureMax << ", the maximum allowed:  "
                 << (((0.2 + m_UpperElevation) * 9.8) / pow(_cVhcesd, 2)) << endl;

            if (!_centerLaneModel.isZero()) {
                cout << "The Center LCC has distance to the orgin : " << fabs(_centerLaneModel[3]) << endl;
            }
        }

        if (!_centerLaneModel.isZero() && (fabs(_StartPosAngle) < m_startAngleThreshold) &&
            _LCCurveLen > _minLCCCurveLen &&
            (_LCCCurvatureMax >= (((0.2 + m_UpperElevation) * 9.8) / pow(_cVhcesd, 2))) &&
            (fabs(_centerLaneModel[3]) >= m_LCC2orign)) {
            infeasibleLCC = CenterCurve;
        } else infeasibleLCC.setZero();

        if (Rotation2G.isZero() || translation2G.isZero())
            readCurbRT(curbRTPath, Rotation2G, translation2G, RoadEdgePts);

        // inverse the Rotation2G in case of nan value occurs
        MatrixXf invRotation2G;
        cv::Mat cvRotation, pinvRotation;
        cv::eigen2cv(Rotation2G, cvRotation);
        cv::invert(cvRotation, pinvRotation, DECOMP_SVD);
        cv::cv2eigen(pinvRotation, invRotation2G);
        if (invRotation2G(0, 0) < 0) invRotation2G(0, 0) = std::fabs(invRotation2G(0, 0));
        if (invRotation2G(1, 1) < 0) invRotation2G(1, 1) = std::fabs(invRotation2G(1, 1));

        GPScuve(m_gpsimudir, _imgName, m_minSamples, _cVhcesd, rawGPsMat, GPScoordinateInL, GPSpline,
                _goalPosAngle);

        if (fabs(_goalPosAngle) > m_InterSectionTrigger)
            GPScuve(m_gpsimudir, _imgName, m_IntersecGPSample, _cVhcesd, rawGPsMat, GPScoordinateInL, GPSpline,
                    _goalPosAngle);

        GpsPts = GPScoordinateInL.cast<float>();

        gpsPtsL = ((invRotation2G * (GpsPts.matrix()).transpose()).colwise() +
                   translation2G.transpose()).array(); //3 x n array


        mapgpsPtsL = ((invRotation2G * (priorGPSInL.matrix()).transpose()).colwise() +
                      translation2G.transpose()).array();

        PtsRemoveBackfCar(mapgpsPtsL, m_minX, m_maxX, priorGPSfcarInL);

        // ALB
        if (!_centerLaneModel.isZero() && fabs(_StartPosAngle) < m_startAngleThreshold &&
            _LCCCurvatureMax < (((0.2 + m_UpperElevation) * 9.8) / pow(_cVhcesd, 2)) && _LCCurveLen > _minLCCCurveLen &&
            fabs(_centerLaneModel[3]) < m_LCC2orign) {

            if (m_debug)
                cout << "Enter ALB Generation..." << endl;

            m_bSplineTLinkage.GetControlPts(_centerLaneModel, _ctrlPts);
            _isVirtual = false;
            CenterVirtul = 0;

        } else {

            if (m_debug)
                cout << "Enter VLB Generation..." << endl;

            CenterVirtul = 0;

            isLCCgenerated = CandidateCurvSelector(GPSpline, _goalPosAngle, m_NumfCurves, _minLCCCurveLen, RoadEdgePts,
                                                   LCCost,
                                                   _OptfLCCurveInG, _cVhcesd, infeasibleLCC, optimalLCCost,
                                                   _prectrilPts);
            if (isLCCgenerated == false)
                throw noLCCgenerated();

            localCenterCurveL = ((invRotation2G * (_OptfLCCurveInG.transpose())).colwise() +
                                 translation2G.transpose()).array();//3 x n array

            m_bSplineTLinkage.FitModel(localCenterCurveL.transpose(), CenterModel);

            m_bSplineTLinkage.GetControlPts(CenterModel, _ctrlPts);

            _isVirtual = true;

        }

        if (m_debug)
            cout << "Exiting MapGenerator::GenerateCenterCtrlPts()" << endl;
    }


    bool
    MapGenerator::CandidateCurvSelector(Eigen::MatrixXf _GPSpline, double _goalPosAngle, double _NumfLCCCurves,
                                        double _minLCCCurveLen,
                                        Eigen::MatrixXf _RoadEdgePts, Eigen::ArrayXXf &_LCCostSum,
                                        Eigen::MatrixXf &_OptLCCSpline, double _cVhcesd, ArrayXXf infeasibleLCC,
                                        RowVectorXf &optimalLCCost, ArrayXXf &_prectrlPts) {

        if (m_debug)
            cout << "Entering MapGenerator::CandidateCurvSelector" << endl;

        int _CurveIdx, totalfCurves;
        float ptsDis, lcc2GPSplineCost, shortTempDis, lcc2albtemp;
        float shortDistoCurbPts = std::numeric_limits<float>::max(), roadEdgeCost;
        float lcc2pretemp, lcc2precost;
        double y, dy, d2y, smoothCost;
        double dycost, d2ycost;
        double lcc2albcost;
        std::vector<int> zeroIdx;
        RowVectorXf LatticeParameter(11);
        MatrixXf pointKXY, OptimalXY, GpsfrontSpline, RoadfEdgePoints, res, preCurve;
        ArrayXf curLCCmodel;
        ArrayXXf costMetricMat;
        vector<RowVector3f> GPSfcurve, roadfedgepts;
        vector<Eigen::RowVectorXf> tempvec;
        alglib::spline1dinterpolant yOfX, zOfX;

        //remove any road edge points larger than _minLCCCurveLen
        for (int i = 0; i < _RoadEdgePts.rows(); ++i) {
            if (_RoadEdgePts(i, 0) <= _minLCCCurveLen) {
                roadfedgepts.push_back(_RoadEdgePts.row(i));
            }
        }

        RoadfEdgePoints = Eigen::MatrixXf(roadfedgepts.size(), 3);
        for (int j = 0; j < roadfedgepts.size(); ++j) {
            RoadfEdgePoints.row(j) = roadfedgepts[j];
        }

        if (m_debug)
            cout << "Start to generate candidate LCC, the _goalPosAngle is : " << _goalPosAngle << endl;

        if (_cVhcesd == 0) _cVhcesd = 0.01;

        try {
            pybind11::module latticexy = pybind11::module::import("Call_LatticePlanner");
            LatticeParameter << 0.0, _NumfLCCCurves, 1, _minLCCCurveLen, m_minAngle, m_maxAngle, _goalPosAngle -
                                                                                                 m_TargetAngleDiff,
                    _goalPosAngle + m_TargetAngleDiff, 100, _goalPosAngle, _cVhcesd;

            pybind11::object result = latticexy.attr("callLatticePlanner")(LatticeParameter);
            pybind11::object numfCurves = latticexy.attr("call_numberofline")();
            res = result.cast<Eigen::MatrixXf>();
            totalfCurves = numfCurves.cast<int>();
        }
        catch (pybind11::error_already_set const &pythonErr) {
            std::cout << "PythonRobotics Error founded : " << pythonErr.what() << endl;
        }

        if (m_debug)
            cout << "Finish generating candidate LCC, total " << totalfCurves << " candidates." << endl;

        if (totalfCurves == 0) {
            return false;
        }

        Eigen::RowVectorXf xPostion = res.row(0);
        for (int i = 0; i < res.cols(); ++i) {
            if (xPostion(i) == 0)
                zeroIdx.push_back(i);
        }

        if (!preCurve.isZero())
            m_bSplineTLinkage.CubicSplineModelToCurve(_prectrlPts, preCurve, m_steplen);

        if (fabs(_goalPosAngle) < m_InterSectionTrigger) {
            if (infeasibleLCC.isZero()) {
                _LCCostSum = Eigen::ArrayXXf(totalfCurves, 4);
            } else {
                _LCCostSum = Eigen::ArrayXXf(totalfCurves, 5);
            }

            for (int k = 0; k < totalfCurves; ++k) {
                lcc2GPSplineCost = 0;
                smoothCost = 0;
                dycost = 0, d2ycost = 0;
                lcc2albcost = 0;
                lcc2precost = 0;
                if (k == 0) {
                    pointKXY = res.block(0, zeroIdx[k], 2, zeroIdx[k + 1]);
                };
                if (k > 0 && k < totalfCurves - 1) {
                    pointKXY = res.block(0, zeroIdx[k], 2, zeroIdx[k + 1] - zeroIdx[k] + 1);
                }
                if (k == totalfCurves - 1) {
                    pointKXY = res.block(0, zeroIdx[k], 2, res.cols() - zeroIdx[k]);
                }

                // cost to GPS spline
                for (int l = 0; l < pointKXY.cols(); ++l) {
                    shortestDistance(pointKXY.col(l).transpose(), GpsfrontSpline, ptsDis);
                    lcc2GPSplineCost += ptsDis;
                }

                //smooth cost
                m_bSplineTLinkage.Fit2DModel(pointKXY.transpose(), curLCCmodel, yOfX);
                for (double x = curLCCmodel(1); x < curLCCmodel(curLCCmodel(0) - 5); x += m_steplen) {
                    alglib::spline1ddiff(yOfX, x, y, dy, d2y);
                    smoothCost += sqrt(pow(dy, 2) + pow(d2y, 2)) * m_steplen;
                }

                //avoid road edges
                if (!RoadfEdgePoints.isZero()) {
                    for (int l = 0; l < pointKXY.cols(); ++l) {
                        shortestDistance(pointKXY.col(l), RoadfEdgePoints, shortTempDis);
                        shortDistoCurbPts = std::min(shortDistoCurbPts, shortTempDis);
                    }
                }
                if (shortDistoCurbPts <= m_carWidth)
                    roadEdgeCost = std::numeric_limits<float>::max();
                else if (shortDistoCurbPts > m_carWidth && shortDistoCurbPts < 5.12)
                    roadEdgeCost = -(pow(10, 4) / 4.12) * shortDistoCurbPts + pow(10, 4) * 5.12 / 4.12;
                else roadEdgeCost = 0;

                //compute cost for infeasible ALB;
                if (!infeasibleLCC.isZero()) {
                    for (int l = 0; l < pointKXY.cols(); ++l) {
                        MatrixXf VLBmat = infeasibleLCC.matrix();
                        shortestDistance(pointKXY.col(l), VLBmat, lcc2albtemp);
                        lcc2albcost += lcc2albtemp;
                    }
                }

                //compute cost for previous lcc
                if (!preCurve.isZero()) {
                    for (int l = 0; l < pointKXY.cols(); ++l) {
                        shortestDistance(pointKXY.col(l), preCurve, lcc2pretemp);
                        lcc2precost += lcc2pretemp;
                    }
                }

                if (infeasibleLCC.isZero()) {
                    _LCCostSum.row(k) << smoothCost, roadEdgeCost, lcc2GPSplineCost, lcc2precost;
                } else {
                    _LCCostSum.row(k) << smoothCost, roadEdgeCost, lcc2GPSplineCost, lcc2precost, lcc2albcost;
                }
            }
        } else { // Road intersection
            if (m_debug)
                cout << "Road Intersection: Generate candidate LCC..." << endl;
            _LCCostSum = Eigen::ArrayXXf(totalfCurves, 3);
            for (int k = 0; k < totalfCurves; ++k) {
                lcc2GPSplineCost = 0;
                smoothCost = 0;
                dycost = 0, d2ycost = 0;
                lcc2albcost = 0;

                if (k == 0) {
                    pointKXY = res.block(0, zeroIdx[k], 2, zeroIdx[k + 1]);
                };
                if (k > 0 && k < totalfCurves - 1) {
                    pointKXY = res.block(0, zeroIdx[k], 2, zeroIdx[k + 1] - zeroIdx[k] + 1);
                }
                if (k == totalfCurves - 1) {
                    pointKXY = res.block(0, zeroIdx[k], 2, res.cols() - zeroIdx[k]);
                }

                //smooth cost
                m_bSplineTLinkage.Fit2DModel(pointKXY.transpose(), curLCCmodel, yOfX);
                for (double x = curLCCmodel(1); x < curLCCmodel(curLCCmodel(0) - 5); x += m_steplen) {
                    alglib::spline1ddiff(yOfX, x, y, dy, d2y);
                    smoothCost += sqrt(pow(dy, 2)) * m_steplen;
                }

                //avoid road edges
                if (!RoadfEdgePoints.isZero()) {
                    for (int l = 0; l < pointKXY.cols(); ++l) {
                        shortestDistance(pointKXY.col(l), RoadfEdgePoints, shortTempDis);
                        shortDistoCurbPts = std::min(shortDistoCurbPts, shortTempDis);
                    }
                }
                if (shortDistoCurbPts <= m_carWidth)
                    roadEdgeCost = std::numeric_limits<float>::max();
                else if (shortDistoCurbPts > m_carWidth && shortDistoCurbPts < 5.12)
                    roadEdgeCost = -(pow(10, 4) / 4.12) * shortDistoCurbPts + pow(10, 4) * 5.12 / 4.12;
                else roadEdgeCost = 0;

                // cost to GPS spline
                for (int l = 0; l < pointKXY.cols(); ++l) {
                    shortestDistance(pointKXY.col(l).transpose(), GpsfrontSpline, ptsDis);
                    lcc2GPSplineCost += ptsDis;
                }

                _LCCostSum.row(k) << smoothCost, roadEdgeCost, lcc2GPSplineCost;

            }
        }

        //standardize cost
        for (int i = 0; i < _LCCostSum.cols(); ++i) {
            if (_LCCostSum.col(i).maxCoeff() != 0) {
                _LCCostSum.col(i) = _LCCostSum.col(i) / (_LCCostSum.col(i).maxCoeff());
            }
        }

        if (fabs(_goalPosAngle) < m_InterSectionTrigger) {
            _LCCostSum.col(0) = m_weightfs * _LCCostSum.col(0);
            _LCCostSum.col(1) = m_weightfo * _LCCostSum.col(1);
            _LCCostSum.col(2) = m_weightfg * _LCCostSum.col(2);
            _LCCostSum.col(3) = m_weightfc * _LCCostSum.col(3);
            if (!infeasibleLCC.isZero()) {
                _LCCostSum.col(4) = m_weightfa * _LCCostSum.col(4);
            }
        } else {
            _LCCostSum.col(0) = m_weightfs * _LCCostSum.col(0);
            _LCCostSum.col(1) = m_weightfo * _LCCostSum.col(1);
            _LCCostSum.col(2) = m_weightfg * _LCCostSum.col(2);
        }

        ArrayXf TotalCost = _LCCostSum.rowwise().sum();
        float minC = TotalCost.minCoeff();

        for (int i = 0; i < TotalCost.rows(); ++i) {
            if (TotalCost(i) == minC)
                _CurveIdx = i;
        }

        optimalLCCost = (_LCCostSum.row(_CurveIdx)).matrix();

        if (_CurveIdx == 0) {
            OptimalXY = res.block(0, zeroIdx[_CurveIdx], 2, zeroIdx[_CurveIdx + 1]);
        };
        if (_CurveIdx > 0 && _CurveIdx < totalfCurves - 1) {
            OptimalXY = res.block(0, zeroIdx[_CurveIdx], 2, zeroIdx[_CurveIdx + 1] - zeroIdx[_CurveIdx] + 1);
        }
        if (_CurveIdx == totalfCurves - 1) {
            OptimalXY = res.block(0, zeroIdx[_CurveIdx], 2, res.cols() - zeroIdx[_CurveIdx]);
        }

        if (OptimalXY.rows() == 2) {
            OptimalXY.conservativeResize(OptimalXY.rows() + 1, OptimalXY.cols());
            Eigen::RowVectorXf vec = Eigen::MatrixXf::Zero(1, OptimalXY.cols());
            OptimalXY.row(OptimalXY.rows() - 1) = vec;
        }

        _OptLCCSpline = OptimalXY.transpose();

        if (m_debug)
            cout << "Existing MapGenerator::CandidateCurvSelector()" << endl;

        return true;
    }

    void MapGenerator::printCostoFile(string imgBaseName, int costMetricMat, RowVectorXf standardizedCostrow) {
        if (m_debug)
            cout << "Entering MapGenerator::printCostoFile()" << endl;

        if (!m_fcostCenterMetric.is_open())
            m_fcostCenterMetric.open(m_costFile.c_str());

        m_fcostCenterMetric << imgBaseName << "\t" << costMetricMat << "\t" << standardizedCostrow << endl;

        if (m_debug)
            cout << "Existing MapGenerator::printCostoFile()" << endl;

    }


    void MapGenerator::GetWorldCtrlPts(const ArrayXf &_model, const MatrixXf &_rotation, const VectorXf &_translation,
                                       ArrayXXf &_worldCtrlPts) {
        if (m_debug)
            cout << "Entering MapGenerator::GetWorldCtrlPts()" << endl;

        ArrayXXf localCtrlPts;
        m_bSplineTLinkage.GetControlPts(_model, localCtrlPts);
        Project(localCtrlPts, _rotation, _translation, _worldCtrlPts);

        if (m_debug)
            cout << "Exiting MapGenerator::GetWorldCtrlPts()" << endl;
    }

    void MapGenerator::GetWorldCtrlPts(const string &_imgName, const vector<ArrayXf> &_models,
                                       vector<ArrayXXf> &_worldCtrlPts) {
        MatrixXf rotation;
        VectorXf translation;

        GetRotationTranslation(_imgName, rotation, translation);
        _worldCtrlPts = vector<ArrayXXf>(_models.size());
        for (int i = 0; i < _models.size(); i++)
            GetWorldCtrlPts(_models[i], rotation, translation, _worldCtrlPts[i]);

    }


    void MapGenerator::GetPassedWorldCtrlPts(const ArrayXXf &_curCtrlPts, const string &_curImgFileName,
                                             ArrayXXf &_passedWorldCtrlPts) {
        if (m_debug)
            cout << "Entering MapGenerator::GetPassedCtrlPts()" << endl;

        MatrixXf rotation;
        VectorXf translation;
        GetRotationTranslation(_curImgFileName, rotation, translation);

        if (ImgFile2Int(_curImgFileName) == m_transformationInfo.size() - 1) {
            _passedWorldCtrlPts = _curCtrlPts;
            Project(_passedWorldCtrlPts, rotation, translation, _passedWorldCtrlPts);
            return;
        }

        ArrayXXf nextPts;
        TransformCtrlPtsToNext(_curCtrlPts, _curImgFileName, nextPts);
        _passedWorldCtrlPts.resizeLike(_curCtrlPts);
        int passedPts = 0;

        for (int c = 0; c < nextPts.cols(); c++) {
            if (nextPts(0, c) < m_minX)
                _passedWorldCtrlPts.col(passedPts++) = _curCtrlPts.col(c);
        }

        _passedWorldCtrlPts.conservativeResize(NoChange, passedPts);

        Project(_passedWorldCtrlPts, rotation, translation, _passedWorldCtrlPts);

        if (m_debug)
            cout << "Exiting MapGenerator::GetPassedCtrlPts()" << endl;

    }

    void
    MapGenerator::TransformCtrlPtsToNext(const ArrayXf &_model, const string &_curImgFileName, ArrayXXf &_ctrlPts) {
        m_bSplineTLinkage.GetControlPts(_model, _ctrlPts);
        TransformCtrlPtsToNext(_ctrlPts, _curImgFileName, _ctrlPts);
    }

    void MapGenerator::TransformCtrlPtsToNext(const Eigen::ArrayXXf &_curCtrlPts, const string &_curImgFileName,
                                              Eigen::ArrayXXf &_nextCtrlPts) {
        int curImg = ImgFile2Int(_curImgFileName);
        if (curImg == m_transformationInfo.size() - 1)
            throw runtime_error("no next image.  Current image must not be the very last one");

        MatrixXf curRotation, nextRotation, relRotation;

        VectorXf curTranslation, nextTranslation, relTranslation;

        GetRotationTranslation(m_transformationInfo[curImg], curRotation, curTranslation);

        GetRotationTranslation(m_transformationInfo[curImg + 1], nextRotation, nextTranslation);

        relRotation = nextRotation.inverse() * curRotation;

        relTranslation = nextTranslation - relRotation * curTranslation;

        Project(_curCtrlPts, relRotation, relTranslation, _nextCtrlPts);

    }


    void
    MapGenerator::vizLCCurve(ArrayXXf _centerCtrlPts, ArrayXXf GpsPts, ArrayXXf mapgpsPtsL, bool _isCenterVirtual,
                             const cv::Mat _inputImg,
                             string _imgBaseName) {
        if (m_debug)
            cout << "Entering MapGenerator::vizLCCurve()" << endl;

        ArrayXXf CtrlPts = _centerCtrlPts.transpose();
        ArrayXXf GPSLpts = GpsPts.transpose();
        ArrayXXf mapGPSLpts = mapgpsPtsL.transpose();

        Mat vizImg = _inputImg.clone();

        if (m_saveVizImg) {

            m_visualizer(vizImg, CtrlPts, (_isCenterVirtual ? Scalar(255, 0, 0) : Scalar(0, 255, 0)));

            m_visualizer(vizImg, GPSLpts, Scalar(0, 0, 255));

            m_visualizer(vizImg, mapGPSLpts, Scalar(255, 0, 255));

            imwrite(m_vizImgPrefix + _imgBaseName, vizImg);
            if (m_debug)
                cout << "The LCC is in " << m_vizImgPrefix + _imgBaseName << endl;
        }

        if (m_debug)
            cout << "Existing MapGenerator::vizLCCurve()" << endl;
    }


    void MapGenerator::ParseVioFile(vector<ArrayXf> &_transformationInfo) {
        if (m_debug)
            cout << "Entering MapGenerator::ParseVioFile()" << endl;

        std::ifstream fin(m_vioFile);
        string line;

        if (!fin)
            throw runtime_error("couldn't open " + m_vioFile);

        while (std::getline(fin, line)) {
            std::istringstream iss(line);
            ArrayXf cur(12); //num points per line in vio file = 12
            for (int i = 0; i < cur.size(); i++)
                iss >> cur(i);
            _transformationInfo.push_back(cur);
        }

        fin.close();

        if (m_debug)
            cout << "Exiting MapGenerator::ParseVioFile()" << endl;

    }

    void MapGenerator::PrintWorldCtrlPts(const ArrayXXf &_ctrlPts, LaneType _type) {

        if (m_debug)
            cout << "Entering MapGenerator::PrintCtrlPts()" << endl;

        switch (_type) {
            case CENTER:
                if (!m_foutCenterCtrlPts.is_open())
                    m_foutCenterCtrlPts.open(m_centerCtrlPtsFile.c_str());

                m_foutCenterCtrlPts << _ctrlPts << endl;
                break;
            case LEFT:
                if (!m_foutLeftCtrlPts.is_open())
                    m_foutLeftCtrlPts.open(m_leftCtrlPtsFile.c_str());

                m_foutLeftCtrlPts << _ctrlPts << endl;
                break;
            case RIGHT:
                if (!m_foutRightCtrlPts.is_open())
                    m_foutRightCtrlPts.open(m_rightCtrlPtsFile.c_str());

                m_foutRightCtrlPts << _ctrlPts << endl;
                break;
        }

        if (m_debug)
            cout << "Exiting MapGenerator::PrintCtrlPts()" << endl;
    }

    void MapGenerator::GetRotationTranslation(const string &_imgName, MatrixXf &_rotation, VectorXf &_translation) {

        if (m_debug)
            cout << "Entering MapGenerator::GetRotationTranslation()" << endl;

        //This function is highly questionable

        int imgNum = ImgFile2Int(_imgName);
        GetRotationTranslation(m_transformationInfo[imgNum], _rotation, _translation);

        if (m_debug)
            cout << "Exiting MapGenerator::GetRotationTranslation()" << endl;
    }

    void MapGenerator::GetRotationTranslation(const ArrayXf &_transformationVec, MatrixXf &_rotation,
                                              VectorXf &_translation) {

        if (m_debug)
            cout << "Entering MapGenerator::GetRotationTranslation()" << endl;

        assert(_transformationVec.size() == 12);
        _rotation = ArrayXXf(3, 3);
        _translation = ArrayXf(3);
        _rotation << _transformationVec(0), _transformationVec(1), _transformationVec(2),
                _transformationVec(4), _transformationVec(5), _transformationVec(6),
                _transformationVec(8), _transformationVec(9), _transformationVec(10);

        _translation << _transformationVec(11), _transformationVec(3), _transformationVec(7);

        if (m_debug)
            cout << "Exiting MapGenerator::GetRotationTranslation()" << endl;
    }


    void
    MapGenerator::GPScuve(const std::string _GPSDirRoot, const std::string _fName, const int _nCtrlPts,
                          double &_cVhcesd,
                          MatrixXd &_rawGPSfMat, Eigen::MatrixXd &_GPSCoordinateInL, Eigen::MatrixXf &_GPSpline,
                          double &_goalPosAngle) {

        if (m_debug)
            cout << "Entering MapGenerator::GPScuve" << endl;

        Eigen::ArrayXf model;
        std::vector<Eigen::RowVector3f> gpscuves;
        double sx, dsx, d2sx;
        alglib::spline1dinterpolant yOfX, zOfX;

        // get gps coordinates in font of the vehicle
        GPSCoordidateTran(_GPSDirRoot, _fName, _nCtrlPts, _cVhcesd, _rawGPSfMat, _GPSCoordinateInL);

        Eigen::ArrayXXf sample = (_GPSCoordinateInL.array()).cast<float>(); // cast to float
        double minx = _GPSCoordinateInL.col(0).minCoeff(), maxx = _GPSCoordinateInL.col(0).maxCoeff();

        m_bSplineTLinkage.FitModel(sample, model);
        m_bSplineTLinkage.CreateSplineInterpolants(model, yOfX, zOfX);

        // gps curves
        for (float x = minx; x <= maxx; x += m_steplen) {
            Eigen::RowVector3f ptOnSpline;
            ptOnSpline << x, spline1dcalc(yOfX, x), spline1dcalc(zOfX, x);
            gpscuves.push_back(ptOnSpline);
        }

        _GPSpline = Eigen::MatrixXf(gpscuves.size(), 3);
        for (int i = 0; i < gpscuves.size(); i++)
            _GPSpline.row(i) = gpscuves[i];

        // goal angle to help select LCC according to GPS spline curve
        spline1ddiff(yOfX, maxx, sx, dsx, d2sx);
        _goalPosAngle = atan(dsx) * 180 / PI;

        if (m_debug)
            cout << "Existing MapGenerator::GPScuve" << endl;

    }

    void MapGenerator::loadPriorGPS(std::string _path, ArrayXXf &_rawMapGPSMat) {

        if (m_debug)
            cout << "Entering MapGenerator::loadPriorGPS" << endl;

        vector<Eigen::RowVector2f> priorMap;
        string gpsline, gps, tempstr;
        vector<vector<float>> gpsMat;
        vector<float> gpsCline;

        std::ifstream fin(_path);
        if (!fin.is_open())
            fin.open(_path);

        while (std::getline(fin, gpsline)) {
            std::istringstream instr(gpsline);
            while (instr >> gps) {
                gpsCline.emplace_back(stof(gps));
            }
            gpsMat.emplace_back(gpsCline);
            gpsCline.clear();
        }

        _rawMapGPSMat = ArrayXXf(gpsMat.size(), 2);
        for (int i = 0; i < gpsMat.size(); ++i) {
            for (int j = 0; j < gpsMat[i].size(); ++j) {
                _rawMapGPSMat(i, j) = gpsMat[i][j];
            }
        }

        if (m_debug)
            cout << "Entering MapGenerator::loadPriorGPS" << endl;
    }


    void
    MapGenerator::priorMapGPsInL(ArrayXXf &_rawMapGPSMat, Eigen::MatrixXd gTruthGPS, Eigen::ArrayXXf &priorGPSInL) {
        //m_rawMapGPSMat and priorGPSInL are n x 2 array,

        if (m_debug)
            cout << "Entering MapGenerator::priorMapGPsInL" << endl;

        Matrix2d rotMat;
        Vector2d original, transformed;
        priorGPSInL = ArrayXXf(_rawMapGPSMat.rows(), _rawMapGPSMat.cols());
        Eigen::RowVectorXf curCarGPs = gTruthGPS.row(1).cast<float>();

        rotMat << cos(curCarGPs(0, 3)), sin(curCarGPs(0, 3)),
                -sin(curCarGPs(0, 3)), cos(curCarGPs(0, 3));

        for (int i = 1; i < _rawMapGPSMat.rows(); ++i) {
            double angle = CoordinatesToAngle(curCarGPs(0, 0), curCarGPs(0, 1), _rawMapGPSMat(i, 0),
                                              _rawMapGPSMat(i, 1));
            double meter = CoordinatesToMeters(curCarGPs(0, 0), curCarGPs(0, 1), _rawMapGPSMat(i, 0),
                                               _rawMapGPSMat(i, 1));

            original << meter * sin(degreeToRadian(angle)),
                    meter * cos(degreeToRadian(angle));
            transformed = rotMat * original;
            priorGPSInL(i, 0) = transformed(0, 0);
            priorGPSInL(i, 1) = transformed(1, 0);
        }

        priorGPSInL.conservativeResize(priorGPSInL.rows(), priorGPSInL.cols() + 1);
        priorGPSInL.col(priorGPSInL.cols() - 1) = Eigen::VectorXf::Zero(priorGPSInL.rows(), 1);

        if (m_debug)
            cout << "Entering MapGenerator::priorMapGPsInL" << endl;
    }

    bool MapGenerator::loadGPSFile(std::string _path, std::string _varname, Eigen::MatrixXd &_outputMat, int _nums) {

        if (m_debug)
            cout << "Entering MapGenerator::loadGPSFile" << endl;

        std::string value, temp;
        int s = 0;
        std::vector<double> IMUMat;
        _outputMat = Eigen::MatrixXd(_nums, 5);
        string gpsFileName;

        for (int i = 0; i < _nums; ++i) {
            // process file name and path

            std::string _name = _varname.substr(0, _varname.length() - 1);
            gpsFileName = _name.substr(0, _name.length() - std::to_string(stoi(_name)).length()) +
                          std::to_string(stoi(_name) + i);

            if (gpsFileName.size() > 10)
                gpsFileName = gpsFileName.substr(gpsFileName.size() - 10, gpsFileName.size());

            std::string _filename = _path + gpsFileName + ".txt";

            std::ifstream in(_filename);
            if (!in.is_open())
                in.open(_filename);
            std::getline(in, value);
            std::istringstream ss(value);

            while (getline(ss, temp, ' ')) {
                s++;
                std::cout.precision(15);
                if (s == 1) _outputMat(i, 0) = stod(temp);
                else if (s == 2) _outputMat(i, 1) = stod(temp);
                else if (s == 3) _outputMat(i, 2) = stod(temp);
                else if (s == 4) continue;
                else if (s == 5) continue;
                else if (s == 6) _outputMat(i, 3) = stod(temp);
                else if (s == 7) continue;
                else if (s == 8) continue;
                else if (s == 9) _outputMat(i, 4) = stod(temp);
                else break;

            }
            s = 0;
            in.close();

        }

        if (m_debug)
            cout << "Existing MapGenerator::loadGPSFile" << endl;

        if (_outputMat.rows() != 0 && _outputMat.cols() != 0)
            return true;
        else
            return false;
    }

    // get a set of control points from the GPS coordinates
    bool MapGenerator::GPSCoordidateTran(std::string _path, std::string _fimeName, int _nums, double &_cVhcesd,
                                         Eigen::MatrixXd &_rawGPSfMat,
                                         Eigen::MatrixXd &_GPSCoordinatInL) {


        if (m_debug)
            cout << "Entering MapGenerator::GPSCoordidateTran" << endl;

        Eigen::Matrix2d rotMat;
        Eigen::Vector2d original, transformed;
        _GPSCoordinatInL = Eigen::MatrixXd(_nums, 3);
        _GPSCoordinatInL(0, 0) = 0;
        _GPSCoordinatInL(0, 1) = 0;
        _GPSCoordinatInL(0, 2) = 0;
        bool status = loadGPSFile(_path, _fimeName, _rawGPSfMat, _nums);

        if (_rawGPSfMat.rows() != 0)
            _cVhcesd = _rawGPSfMat.col(4).mean();
        else return false;

        if (status == false) return false;
        for (int i = 1; i < _rawGPSfMat.rows(); ++i) {
            double angle = CoordinatesToAngle(_rawGPSfMat(0, 0), _rawGPSfMat(0, 1), _rawGPSfMat(i, 0),
                                              _rawGPSfMat(i, 1));
            double meter = CoordinatesToMeters(_rawGPSfMat(0, 0), _rawGPSfMat(0, 1), _rawGPSfMat(i, 0),
                                               _rawGPSfMat(i, 1));
            rotMat << cos(_rawGPSfMat(i, 3)), sin(_rawGPSfMat(i, 3)),
                    -sin(_rawGPSfMat(i, 3)), cos(_rawGPSfMat(i, 3));
            original << meter * sin(degreeToRadian(angle)),
                    meter * cos(degreeToRadian(angle));
            transformed = rotMat * original;
            _GPSCoordinatInL(i, 0) = transformed(0, 0);
            _GPSCoordinatInL(i, 1) = -transformed(1, 0);
            _GPSCoordinatInL(i, 2) = _rawGPSfMat(i, 2) - _rawGPSfMat(0, 2);
        }

        if (m_debug)
            cout << "Existing MapGenerator::GPSCoordidateTran" << endl;

        if (_GPSCoordinatInL.isZero(0)) return false;
        else return true;
    }

    double MapGenerator::CoordinatesToAngle(double latitude1,
                                            const double longitude1,
                                            double latitude2,
                                            const double longitude2) {

        const auto longitudeDifference = degreeToRadian(longitude2 - longitude1);
        latitude1 = degreeToRadian(latitude1);
        latitude2 = degreeToRadian(latitude2);

        using namespace std;
        const auto x = (cos(latitude1) * sin(latitude2)) -
                       (sin(latitude1) * cos(latitude2) * cos(longitudeDifference));
        const auto y = sin(longitudeDifference) * cos(latitude2);
        auto degree = radianToDegree(atan2(y, x));

        //degree = degree - 90;
        return (degree >= 0) ? degree : (degree + 360);
    }

    double MapGenerator::CoordinatesToMeters(double latitude1, double longitude1, double latitude2, double longitude2) {

        latitude1 = degreeToRadian(latitude1);
        longitude1 = degreeToRadian(longitude1);
        latitude2 = degreeToRadian(latitude2);
        longitude2 = degreeToRadian(longitude2);

        auto x = sin((latitude2 - latitude1) / 2), y = sin((longitude2 - longitude1) / 2);

        return m_earthDiameterMeters * asin(sqrt((x * x) + (cos(latitude1) * cos(latitude2) * y * y)));
    }

    std::pair<double, double>
    MapGenerator::CoordinateToCoordinate(double latitude, double longitude, double angle, double meters) {
        latitude = degreeToRadian(latitude);
        longitude = degreeToRadian(longitude);
        angle = degreeToRadian(angle);
        meters *= 2 / m_earthDiameterMeters;
        using namespace std;
        pair<double, double> coordinate;
        coordinate.first = asin((sin(latitude) * cos(meters))
                                + (cos(latitude) * sin(meters) * cos(angle)));
        coordinate.second = longitude + atan2((sin(angle) * sin(meters) * cos(latitude)),
                                              cos(meters) - (sin(latitude) * sin(coordinate.first)));
        coordinate.first = radianToDegree(coordinate.first);
        coordinate.second = radianToDegree(coordinate.second);

        return coordinate;
    }

}