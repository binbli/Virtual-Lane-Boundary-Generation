#include"../include/LaneDetector.h"
#include "../3rdParty/TLinkage/Sampler.h"
#include<fstream>
#include<libgen.h>
#include <chrono>
#include <ratio>
#include <thread>
#include"../include/Utilities.h"

namespace LD {

    LaneDetector::LaneDetector(string _xmlFile) : Solver(_xmlFile), m_RoadSegment(_xmlFile),
                                                  m_KPercentExtractor(_xmlFile), m_LMsfromCam(_xmlFile),
                                                  m_LMsintersection(_xmlFile), m_bSplineTLinkage(_xmlFile),
                                                  m_laneQualityChecker(_xmlFile), m_visualizer(_xmlFile),
                                                  m_mapGenerator(_xmlFile) {
        ParseXML();
        assert(!m_LMsintersection.isMode2D());
    }

    void LaneDetector::ParseXML() {
        if (m_debug)
            cout << "Entering LaneDetector::ParseXML()" << endl;

        m_xml = m_xml.child("LaneDetector");

        m_imgRoot = m_xml.attribute("dataRoot").as_string();
        m_imgFile = m_xml.attribute("dataFile").as_string();
        m_veloRoot = m_xml.attribute("veloRoot").as_string();
        m_ratiosFile = m_xml.attribute("ratiosFile").as_string();
        m_saveVizImg = m_xml.attribute("saveVizImg").as_bool(true);
        m_vizImgPrefix = m_xml.attribute("vizImgPrefix").as_string();
        m_isCenterVirtual = m_xml.attribute("isCenterVirtual").as_string();
        m_laneAssess = m_xml.attribute("lassAsses").as_bool();

        if (m_imgRoot.empty() || m_imgFile.empty() || m_veloRoot.empty() || m_ratiosFile.empty() ||
            (m_saveVizImg && m_vizImgPrefix.empty()))
            throw runtime_error(
                    "at least one of the following attributes is missing: imgRoot, imgFile, veloRoot, ratiosFile, vizImgPrefix, saveVizImg");

        if (m_debug)
            cout << "Existing LaneDetector::ParseXML()" << endl;
    }

    void LaneDetector::operator()(const cv::Mat &_inputImg, cv::Mat &segImg, const cv::Mat &_veloPoints,
                                  const string _imgBaseName,
                                  vector<Eigen::ArrayXXf> &_ctrlPts, vector<bool> &_isVirtual,
                                  LaneAssessmentParams &_assessParams,
                                  Eigen::Matrix3f &RotationG, Eigen::RowVector3f &translationG,
                                  Eigen::MatrixXf &RoadEdgePts, vector<Eigen::ArrayXXf> &_previousCtrlPts,
                                  bool &isCenterVirtual, ArrayXXf &centerCtrlPts, Eigen::ArrayXXf &GpsPts,
                                  Eigen::ArrayXXf &mapgpsPtsL, cv::Mat &roadMaskImg) {
        if (m_debug)
            cout << "Entering LaneDetector::()" << endl;

        cv::Mat refinedImg, reflectivity, savedImg, newsegImg;
        Eigen::ArrayXXf intersectedPts;
        Eigen::ArrayXf clusters;
        Eigen::MatrixXf veloImg;
        Eigen::ArrayXf centerLaneModel;
        vector<ArrayXf> models;
        ArrayXXf _updatedNewCtrolPts;
        ArrayXXf _nextLeftCtrlPts, _nextRightCtrlPts, _nextCenterCtrlPts;
        ArrayXf DummyMat;
        DummyMat.setZero();

        m_RoadSegment(_inputImg, segImg);
        m_LMsfromCam(_inputImg, segImg, refinedImg);
        m_LMsintersection(_veloPoints, segImg, refinedImg, intersectedPts, reflectivity, veloImg,
                          _inputImg, _imgBaseName, RotationG, translationG, RoadEdgePts, roadMaskImg);
        m_bSplineTLinkage(intersectedPts, clusters, models, _isVirtual);


        if (!models.empty()) {

            m_mapGenerator.GenerateCenterCtrlPts(models.back(), _imgBaseName, centerCtrlPts, _previousCtrlPts[2], isCenterVirtual, RotationG,
                                                 translationG, RoadEdgePts, GpsPts, mapgpsPtsL);
        } else {
            m_mapGenerator.GenerateCenterCtrlPts(DummyMat, _imgBaseName, centerCtrlPts, _previousCtrlPts[2], isCenterVirtual, RotationG,
                                                 translationG, RoadEdgePts, GpsPts, mapgpsPtsL);
        }

        if (!isCenterVirtual) {
            _ctrlPts.resize(models.size());
            for (int i = 0; i < _ctrlPts.size(); i++)
                m_bSplineTLinkage.GetControlPts(models[i], _ctrlPts[i]);
        } else {
            ShiftCenterCtrlPts(centerCtrlPts, m_bSplineTLinkage.GetAvgCenterLaneWidth(),
                               m_bSplineTLinkage.GetAvgCenterLaneHeight(), m_bSplineTLinkage.IsRightLaneHigher(),
                               _ctrlPts);
            _ctrlPts.push_back(centerCtrlPts);
            _isVirtual = vector<bool>(_ctrlPts.size(), true);
        }


        if (m_saveVizImg) {
            Mat vizImg = _inputImg.clone();

            if (!GpsPts.isZero()) {
                ArrayXXf GPSLpts = GpsPts.transpose();
                m_visualizer(vizImg, GPSLpts, Scalar(0, 0, 255));
            }

            if (!mapgpsPtsL.isZero()) {
                ArrayXXf mapGPSLpts = mapgpsPtsL.transpose();
                m_visualizer(vizImg, mapGPSLpts, Scalar(255, 0, 255));
            }


            for (int i = 0; i < _ctrlPts.size(); i++) {
                ArrayXXf curCtrlPts = _ctrlPts[i].transpose();
                m_visualizer(vizImg, curCtrlPts, ((i == 2) ? Scalar(255, 0, 0) : Scalar(0, 255, 0)));
            }

            m_KPercentExtractor.noiseRemover(roadMaskImg);
            m_RoadSegment.OverlayMask(roadMaskImg, newsegImg);
            cv::addWeighted(vizImg, 0.8, newsegImg, 0.2, 0.0, savedImg);
            imwrite(m_vizImgPrefix + m_imgBaseName, savedImg);


            if (m_debug)
                cout << "The LBs and LCC are displayed in " << m_vizImgPrefix + _imgBaseName << endl;
        }

        m_mapGenerator.PrintPassedWorldCtrlPts(_ctrlPts[0], m_imgBaseName,
                                               _ctrlPts[0](1, 0) < 0 ? MapGenerator::LaneType::RIGHT
                                                                     : MapGenerator::LaneType::LEFT);
        m_mapGenerator.PrintPassedWorldCtrlPts(_ctrlPts[1], m_imgBaseName,
                                               _ctrlPts[0](1, 0) < 0 ? MapGenerator::LaneType::LEFT
                                                                     : MapGenerator::LaneType::RIGHT);
        m_mapGenerator.PrintPassedWorldCtrlPts(_ctrlPts[2], m_imgBaseName,
                                               MapGenerator::LaneType::CENTER); //center is always last

        for (int i = 0; i < _ctrlPts.size(); ++i) {
            if (i == 0 || i == 1) {
                m_mapGenerator.TransformCtrlPtsToNext(_ctrlPts[i], m_imgBaseName,
                                                      _ctrlPts[i](1, 0) < 0 ? _nextRightCtrlPts : _nextLeftCtrlPts);
            } else m_mapGenerator.TransformCtrlPtsToNext(_ctrlPts[i], m_imgBaseName, _nextCenterCtrlPts);
        }

        _previousCtrlPts.clear();
        _previousCtrlPts.push_back(_nextLeftCtrlPts);
        _previousCtrlPts.push_back(_nextRightCtrlPts);
        _previousCtrlPts.push_back(_nextCenterCtrlPts);

        if (m_laneAssess) {
            m_laneQualityChecker(intersectedPts, clusters, _veloPoints, reflectivity, veloImg, _inputImg, segImg,
                                 refinedImg, _assessParams.brightnessRatio, _assessParams.reflectivityRatio);
            int npt = 1;
            double sum = 0;
            for (int i = 0; i < models.size(); i++) {
                for (int j = 0; j < intersectedPts.cols(); ++j) {
                    double PtsToSpline = m_bSplineTLinkage.Distance(intersectedPts.col(j), models[i]);
                    if (PtsToSpline <= 0.3) {
                        sum += PtsToSpline;
                        npt++;
                    }
                }

            }
            _assessParams.avgDist = (float) sum / npt; // the shape metric
        }


        if (m_debug)
            cout << "Exiting LaneDetector::()" << endl;
    }

    void LaneDetector::Run() {

        if (m_debug)
            cout << "Entering LaneDetector::Run()" << endl;

        std::ifstream fin(m_imgFile.c_str());
        std::ofstream fout(m_ratiosFile.c_str());
        std::ofstream fvirout(m_isCenterVirtual);
        string line;
        float AreaDeriva;
        ulli index = 0;

        Eigen::Matrix3f RotationG;
        Eigen::RowVector3f translationG;
        Eigen::MatrixXf RoadEdgePts;

        vector<ArrayXf> models;
        vector<Eigen::ArrayXXf> previousCtrlPts;
        ArrayXXf centerCtrlPts, updatedNewCtrolPts, GpsPts, mapgpsPtsL;
        ArrayXXf nextLeftCtrlPts, nextRightCtrlPts, nextCenterCtrlPts;
        ArrayXf DummyMat;
        cv::Mat roadMaskImg;

        bool isCenterVirtual;
        previousCtrlPts.clear();
        DummyMat.setZero();

        while (std::getline(fin, line)) {
            LaneAssessmentParams assessParams;
            m_imgBaseName = basename(const_cast<char *>(line.c_str()));
            cv::Mat inputImg, veloPoints, savedImg, segImg, newsegImg;
            vector<Eigen::ArrayXXf> ctrlPts;
            vector<bool> isVirtual;
            inputImg = cv::imread(m_imgRoot + "/" + line);
            if (!inputImg.data) {
                throw runtime_error("invalid path: " + m_imgRoot + "/" + line);
            }
            m_LMsintersection.ReadVeloData(m_veloRoot + "/" + line.substr(0, line.size() - 3) + "bin", veloPoints);

            try {
                auto start = std::chrono::steady_clock::now();
                this->operator()(inputImg, segImg, veloPoints, m_imgBaseName, ctrlPts, isVirtual, assessParams,
                                 RotationG,
                                 translationG, RoadEdgePts, previousCtrlPts, isCenterVirtual, centerCtrlPts, GpsPts,
                                 mapgpsPtsL, roadMaskImg);
                auto finish = std::chrono::steady_clock::now();

                if (m_debug)
                    cout << "It takes " << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                         << " milliseconds to process the Image " << m_imgBaseName << endl;

                if (isCenterVirtual) {
                    fvirout << m_imgBaseName << "\t 1" << endl;
                } else {
                    fvirout << m_imgBaseName << "\t 0" << endl;
                }

                fout << m_imgBaseName << endl;

                if (m_laneAssess) {
                    if (index == 0)
                        AreaDeriva = assessParams.areaDiff;

                    index += 1;
                    fout << max(assessParams.brightnessRatio, assessParams.reflectivityRatio) << endl;
                    fout << assessParams.avgDist << endl;
                    fout << (float) assessParams.areaDiff / AreaDeriva << endl;
                }
            }
            catch (Sampler::MinSamplesNotFound) {

                if (m_debug)
                    cout << "// BSplineTLinkage lacks enough samples for spline" << endl;

                m_mapGenerator.GenerateCenterCtrlPts(DummyMat, m_imgBaseName, centerCtrlPts, previousCtrlPts[2], isCenterVirtual, RotationG,
                                                     translationG, RoadEdgePts, GpsPts, mapgpsPtsL);

                if (!isCenterVirtual) {
                    ctrlPts.resize(models.size());
                    for (int i = 0; i < ctrlPts.size(); i++)
                        m_bSplineTLinkage.GetControlPts(models[i], ctrlPts[i]);
                } else {
                    ShiftCenterCtrlPts(centerCtrlPts, m_bSplineTLinkage.GetAvgCenterLaneWidth(),
                                       m_bSplineTLinkage.GetAvgCenterLaneHeight(),
                                       m_bSplineTLinkage.IsRightLaneHigher(), ctrlPts);
                    ctrlPts.push_back(centerCtrlPts);
                    isVirtual = vector<bool>(ctrlPts.size(), true);
                }

                if (m_saveVizImg) {
                    Mat vizImg = inputImg.clone();

                    if (!GpsPts.isZero()) {
                        ArrayXXf GPSLpts = GpsPts.transpose();
                        m_visualizer(vizImg, GPSLpts, Scalar(0, 0, 255));
                    }

                    if (!mapgpsPtsL.isZero()) {
                        ArrayXXf mapGPSLpts = mapgpsPtsL.transpose();
                        m_visualizer(vizImg, mapGPSLpts, Scalar(255, 0, 255));
                    }

                    for (int i = 0; i < ctrlPts.size(); i++) {
                        ArrayXXf curCtrlPts = ctrlPts[i].transpose();
                        m_visualizer(vizImg, curCtrlPts, ((i == 2) ? Scalar(255, 0, 0) : Scalar(0, 255, 0)));
                    }

                    m_KPercentExtractor.noiseRemover(roadMaskImg);
                    m_RoadSegment.OverlayMask(roadMaskImg, newsegImg);
                    cv::addWeighted(vizImg, 0.8, newsegImg, 0.2, 0.0, savedImg);

                    imwrite(m_vizImgPrefix + m_imgBaseName, savedImg);

                    if (m_debug)
                        cout << "The LBs and LCC are displayed in " << m_vizImgPrefix + m_imgBaseName << endl;
                }

                m_mapGenerator.PrintPassedWorldCtrlPts(ctrlPts[0], m_imgBaseName,
                                                       ctrlPts[0](1, 0) < 0 ? MapGenerator::LaneType::RIGHT
                                                                            : MapGenerator::LaneType::LEFT);
                m_mapGenerator.PrintPassedWorldCtrlPts(ctrlPts[1], m_imgBaseName,
                                                       ctrlPts[0](1, 0) < 0 ? MapGenerator::LaneType::LEFT
                                                                            : MapGenerator::LaneType::RIGHT);
                m_mapGenerator.PrintPassedWorldCtrlPts(ctrlPts[2], m_imgBaseName,
                                                       MapGenerator::LaneType::CENTER); //center is always last


                for (int i = 0; i < ctrlPts.size(); ++i) {
                    if (i == 0 || i == 1) {
                        m_mapGenerator.TransformCtrlPtsToNext(ctrlPts[i], m_imgBaseName,
                                                              ctrlPts[i](1, 0) < 0 ? nextRightCtrlPts
                                                                                   : nextLeftCtrlPts);
                    } else m_mapGenerator.TransformCtrlPtsToNext(ctrlPts[i], m_imgBaseName, nextCenterCtrlPts);
                }

                previousCtrlPts.clear();
                previousCtrlPts.push_back(nextLeftCtrlPts);
                previousCtrlPts.push_back(nextRightCtrlPts);
                previousCtrlPts.push_back(nextCenterCtrlPts);

                if (isCenterVirtual) {
                    fvirout << m_imgBaseName << "\t 1" << endl;
                } else {
                    fvirout << m_imgBaseName << "\t 0" << endl;
                }


                fout << m_imgBaseName << endl;
                fout << "---------------------WARNING: Could not find enough samples-----------" << endl;
            }
            catch(MapGenerator::noLCCgenerated){
                fout << m_imgBaseName << endl;
                fout << "---------------------WARNING: Could not generate enough LCC -----------" << endl;
            }
            catch (alglib::ap_error &e) {
                fout << m_imgBaseName << endl;
                fout << "---------------------WARNING: ALGLIB Error: " << e.msg << "-----------" << endl;
            }
            catch (std::runtime_error &e) {
                fout << m_imgBaseName << endl;
                fout << "---------------------WARNING: runtime_error: " << e.what() << "-----------" << endl;
            }
            catch (...) {
                fout << m_imgBaseName << endl;
                fout << "----------------------------Warning: Unknown Error------------------------------" << endl;
            }
        }

        if (m_debug)
            cout << "Exiting LaneDetector::Run()" << endl;
    }
}