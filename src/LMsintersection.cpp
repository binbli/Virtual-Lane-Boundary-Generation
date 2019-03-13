#include"../include/LMsintersection.h"
#include"../include/Utilities.h"


namespace LD {


    Eigen::MatrixXf LMsintersection::H_p = Eigen::MatrixXf::Zero(4, 1);

    void LMsintersection::ParseXML() {

        m_xml = m_xml.child("LMintersection");
        m_segRoot = m_xml.attribute("segRoot").as_string();
        m_refinedRoot = m_xml.attribute("refinedRoot").as_string();
        m_segImgPrefix = m_xml.attribute("segImgPrefix").as_string();
        m_refImgPrefix = m_xml.attribute("refImgPrefix").as_string();
        m_outputFile = m_xml.attribute("outputFile").as_string();
        m_saveVizImg = m_xml.attribute("saveVizImg").as_bool();
        m_vizImgPrefix = m_xml.attribute("vizImgPrefix").as_string();
        m_maxWidth = m_xml.attribute("maxWidth").as_int();
        m_maxLength = m_xml.attribute("maxLength").as_int();
        m_printOnly2D = m_xml.attribute("printOnly2D").as_bool(false);
        m_dThrehold = m_xml.attribute("dThrehold").as_float();
        m_nIter = m_xml.attribute("nIter").as_int();
        m_minSamples = m_xml.attribute("minSamples").as_int();
        m_percentInliers = m_xml.attribute("percentInliers").as_float();
        m_radius = m_xml.attribute("radius").as_float();
        m_CurbRTRoot = m_xml.attribute("CurbRTRoot").as_string();
        m_CurbRTOutFile = m_xml.attribute("CurbRTOutFile").as_string();
        m_height = m_xml.attribute("height").as_float();
        m_curbH = m_xml.attribute("curbH").as_double();
        m_surfCurvePara = m_xml.attribute("surfCurvePara").as_double();

        if (m_segRoot.empty() || m_refinedRoot.empty() || m_segImgPrefix.empty() || m_refImgPrefix.empty() ||
            m_outputFile.empty() || (m_saveVizImg && m_vizImgPrefix.empty()) || !m_maxWidth || !m_maxLength ||
            !m_nIter || !m_minSamples || !m_dThrehold || !m_radius || !m_percentInliers || !m_height || !m_curbH ||
            !m_surfCurvePara || m_CurbRTOutFile.empty() || m_CurbRTRoot.empty())
            throw runtime_error(
                    "at least one of the following attributes are missing in LMsintersection node: segRoot, refinedRoot, "
                    "outputFile, segImgPrefix, refImgPrefix, vizImgPrefix, maxWidth, maxHeight, m_nIter, m_minSamples, m_dThrehold, m_nIter, m_minSamples, m_radius, m_percentInliers, m_height, m_CurbRTinFile, m_surfCurvePara, m_surfCurvePara, m_CurbRTRoot");
    }

    void LMsintersection::ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints,
                                                   Mat &_reflectivity, Mat &_inputImg) {

        if (m_debug)
            cout << "Enter LMintersection::ProcessProjectedLidarPts()" << endl;

        Eigen::ArrayXXf intersectedPts, roadboundaryPts;
        Eigen::MatrixXf _curbptsImg;
        Mat _roadBvizImg, _vizcurbImg;
        Eigen::MatrixXf refineCurbSet;
        Eigen::Matrix3f RotationL2G;
        Eigen::RowVector3f translationL2G;

        string segImgName = m_segRoot + "/" + m_segImgPrefix + m_imgBaseName;
        string refImgName = m_refinedRoot + "/" + m_refImgPrefix + m_imgBaseName;

        Mat segImg = imread(segImgName, IMREAD_GRAYSCALE);
        Mat refinedImg = imread(refImgName, IMREAD_GRAYSCALE);
        Mat roadMaskImg;
        _roadBvizImg = _inputImg.clone();
        _vizcurbImg = _inputImg.clone();

        if (segImg.empty())
            throw std::runtime_error("Can't open " + segImgName);
        else if (m_debug)
            cout << "Successfully read segmented image: " << segImgName << endl;

        if (refinedImg.empty())
            throw std::runtime_error("Can't open " + refImgName);
        else if (m_debug)
            cout << "Successfully read refined image: " << refImgName << endl;

        double thresh = OtsuThresholdRoad(_veloImg, segImg, _reflectivity);
        if (m_debug)
            cout << "Threshold set to " << thresh << endl;

        this->operator()(_veloPoints, segImg, refinedImg, intersectedPts, _reflectivity, _veloImg,
                         _inputImg, m_imgBaseName, RotationL2G, translationL2G, refineCurbSet, roadMaskImg);

        if (m_debug)
            cout << "Exiting LMintersection::ProcessProjectedLidarPts()" << endl;

    }

    void LMsintersection::operator()(const Mat &_veloPoints, const Mat &_segImg, const Mat &_refinedImg,
                                     Eigen::ArrayXXf &_intersectedPts, Mat &_reflectivity,
                                     Eigen::MatrixXf &_veloImgPts, const Mat &_inputImg,
                                     const string m_imgBaseName, Eigen::Matrix3f &RotationL2G,
                                     Eigen::RowVector3f &translationL2G, Eigen::MatrixXf &refineCurbSet,
                                     cv::Mat &roadMaskImg) {
        if (m_debug)
            cout << "Entering LMintersection::()" << endl;

        cv::Mat _origImg = _inputImg.clone(), _edgeRoadImg, _roadBvizImg;
        Eigen::ArrayXXf _roadboundaryPts;
        _roadBvizImg = _inputImg.clone();
        Project(_veloPoints, _veloImgPts, _reflectivity);

        double _thresh = OtsuThresholdRoad(_veloImgPts, _segImg, _reflectivity);

        if (m_debug)
            cout << "Threshold set to " << _thresh << endl;

        // get the 3D lane marking points
        IntersectIn3D(_veloImgPts, _veloPoints, _reflectivity, _refinedImg, _thresh, _intersectedPts, _origImg);

        if (m_debug)
            PrintToFile(_intersectedPts);

        RoadEdgePtsDetector(_veloImgPts, _veloPoints, _segImg, refineCurbSet, RotationL2G, translationL2G, roadMaskImg);

        if (m_debug)
            PrintEdgePtsRT(refineCurbSet, RotationL2G, translationL2G);

        if (m_saveVizImg && m_debug)
            imwrite(m_outputRoot + "/" + m_vizImgPrefix + m_imgBaseName, _origImg);

        if (m_debug)
            cout << "LMsintersection::() Successfully saved " << m_outputRoot + "/" + m_vizImgPrefix + m_imgBaseName
                 << endl;


        if (m_debug)
            cout << "Entering LMintersection::()" << endl;
    }

    void
    LMsintersection::IntersectIn3D(const Eigen::MatrixXf _veloImg, const Mat &_veloPoints, const Mat &_reflectivity,
                                   const Mat &_refinedImg, const double &_thresh, Eigen::ArrayXXf &_intersectedPts,
                                   Mat &_vizImg) {
        if (m_debug)
            cout << "Entering LMintersection::IntersectIn3D()" << endl;

        vector<vector<float> > intersectedPtsVec;

        for (ulli i = 0; i < _veloImg.rows(); i++) {
            int xImg = _veloImg(i, 0), yImg = _veloImg(i, 1);
            float xLidar = _veloPoints.at<float>(i, 0), yLidar = _veloPoints.at<float>(i, 1);
            float reflect = _reflectivity.at<float>(i, 0);
            if (isValid(yImg, xImg, _refinedImg.rows, _refinedImg.cols) && std::abs(yLidar) <= m_maxWidth &&
                std::abs(xLidar) <= m_maxLength &&
                _refinedImg.at<unsigned char>(yImg, xImg) && reflect >= _thresh) {
                if (m_printOnly2D)
                    intersectedPtsVec.push_back({xLidar, yLidar});
                else {
                    intersectedPtsVec.push_back({xLidar, yLidar, _veloPoints.at<float>(i, 2)});
                }
                if (m_debug) {
                    if (m_saveVizImg && !_vizImg.empty())
                        circle(_vizImg, Point(xImg, yImg), 3, Scalar(0, 255, 0), CV_FILLED);
                }
            }
        }

        _intersectedPts.resize(intersectedPtsVec.size() ? intersectedPtsVec[0].size() : 0, intersectedPtsVec.size());

        for (ulli r = 0; r < _intersectedPts.rows(); r++)
            for (ulli c = 0; c < _intersectedPts.cols(); c++)
                _intersectedPts(r, c) = intersectedPtsVec[c][r];

        if (m_debug)
            cout << "Exiting LMintersection::IntersectIn3D()" << endl;
    }

    void LMsintersection::PrintToFile(Eigen::ArrayXXf &_intersectedPts) {
        if (m_debug)
            cout << "Entering LMintersection::PrintToFile()" << endl;
        string intersectedPtsPath =
                m_outputRoot + "/" + m_outputFile + m_imgBaseName.substr(0, m_imgBaseName.size() - 4) + ".txt";
        m_fout.open(intersectedPtsPath);
        if (!m_fout.is_open())
            m_fout.open(intersectedPtsPath);

        m_fout << m_imgBaseName << endl;
        m_fout << _intersectedPts.cols() << "\t" << _intersectedPts.rows() << endl;

        for (ulli c = 0; c < _intersectedPts.cols(); c++) {
            for (ulli r = 0; r < _intersectedPts.rows(); r++)
                m_fout << _intersectedPts(r, c) << "\t";

            m_fout << endl;
        }
        m_fout.close();
        if (m_debug)
            cout << "Exiting LMintersection::PrintToFile()" << endl;
    }

    void LMsintersection::PrintEdgePtsRT(Eigen::MatrixXf refineCurbSet, Eigen::Matrix3f RotationL2G,
                                         Eigen::RowVector3f translationL2G) {
        if (m_debug)
            cout << "Entering LMintersection::PrintToFile()" << endl;

        string filepath = m_CurbRTRoot + m_CurbRTOutFile + m_imgBaseName.substr(0, m_imgBaseName.size() - 4) + ".txt";
        std::ofstream file(filepath);
        if (file.is_open()) {
            //print R
            for (int i = 0; i < RotationL2G.rows(); ++i) {
                for (int j = 0; j < RotationL2G.cols(); ++j)
                    file << RotationL2G(i, j) << "\t";
            }

            file << "\n";

            // print T
            for (int k = 0; k < translationL2G.cols(); ++k) {
                file << translationL2G(0, k) << "\t";
            }

            file << "\n";

            for (int i = 0; i < refineCurbSet.rows(); ++i) {
                for (int j = 0; j < refineCurbSet.cols(); ++j) {
                    file << refineCurbSet(i, j) << "\t";
                }
            }

            if (m_debug)
                cout << " Successfully printed to file LMsintersection::PrintEdgePtsRT" << endl;
        } else {
            if (m_debug)
                cout << " Cannot print to file LMsintersection::PrintEdgePtsRT" << endl;
        }

        if (file.is_open()) file.close();

        if (m_debug)
            cout << "Existing LMintersection::PrintToFile()" << endl;

    }


    double
    LMsintersection::OtsuThresholdRoad(const Eigen::MatrixXf _veloImg, const Mat &_segImg, const Mat &_reflectivity) {
        if (m_debug)
            cout << "Entering LMintersection::OtsuThresholdRoad()" << endl;

        //Find Otsu thresholding for points that are on road and have positive reflectivity
        Mat scaledReflectivity = 255 * _reflectivity;
        scaledReflectivity.convertTo(scaledReflectivity, CV_8UC1);
        vector<unsigned char> onRoadRef;
        onRoadRef.reserve(scaledReflectivity.rows);

        for (ulli i = 0; i < scaledReflectivity.rows; i++) {
            int x = _veloImg(i, 0), y = _veloImg(i, 1);
            int reflect = scaledReflectivity.at<unsigned char>(i, 0);
            if (isValid(y, x, _segImg.rows, _segImg.cols) && _segImg.at<unsigned char>(y, x))
                onRoadRef.push_back(scaledReflectivity.at<unsigned char>(i, 0));
        }

        double thresh = threshold(onRoadRef, onRoadRef, 1, 255, THRESH_TOZERO | THRESH_OTSU) / 255;

        if (m_debug)
            cout << "Exiting LMintersection::OtsuThresholdRoad()" << endl;

        return thresh;
    }

    void LMsintersection::RoadEdgePtsDetector(const Eigen::MatrixXf _veloImg, const Mat _veloPoints,
                                              const Mat _segImg, Eigen::MatrixXf &refineCurbSet,
                                              Eigen::Matrix3f &RotationL2G, Eigen::RowVector3f &translationL2G,
                                              cv::Mat &roadMaskImg) {

        if (m_debug)
            cout << "Entering LMintersection:: RoadEdgePtsDetector" << endl;

        Mat roadsurfLidarpts;
        int upIter = 500;
        int nsumforplane = 3;
        Eigen::MatrixXf roadlidarpts, roadwayPts, roadwayPtsImg;
        Eigen::MatrixXd CurbsPointSet;
        Eigen::ArrayXXf roadwayImgPts;
        vector<vector<int> > roadwayPtsVec;
        roadMaskImg = cv::Mat::zeros(_segImg.rows, _segImg.cols, CV_8UC1);

        for (int i = 0; i < _veloImg.rows(); i++) {
            int x = _veloImg(i, 0), y = _veloImg(i, 1);
            if (isValid(y, x, _segImg.rows, _segImg.cols) && abs(_veloPoints.at<float>(i, 0)) <= m_maxLength &&
                abs(_veloPoints.at<float>(i, 1)) <= m_maxWidth) {
                cv::Mat m = _veloPoints.row(i).colRange(0, 3);
                roadsurfLidarpts.push_back(m);
            }
        }

        cv2eigen(roadsurfLidarpts, roadlidarpts);
        int ninlierThresholdforPlane = std::floor(m_percentInliers * roadlidarpts.rows());
        m_nIter = std::min(m_nIter, upIter);
        retresult coeffplane = runRansacplanefit(roadlidarpts, nsumforplane, m_nIter, ninlierThresholdforPlane,
                                                 m_dThrehold);
        H_p = coeffplane.modelcoff;
        retresult area_plane = mcdPlaneComput(coeffplane.modelcoff, roadlidarpts, m_dThrehold, false);
        coordinatetransL2G(area_plane.inliers, RotationL2G, translationL2G);

        //get road edges and obstacle points
        computeRoadEdgePts(roadlidarpts, H_p, m_radius, CurbsPointSet);
        refineCurbSet = CurbsPointSet.cast<float>();

        //obtain roadway points
        cv2eigen(_veloPoints(Range(0, _veloPoints.rows - 1), Range(0, 3)), roadwayPts);
        retresult area_road = mcdPlaneComput(coeffplane.modelcoff, roadwayPts, m_height, true);
        Eigen::ArrayXXf roadway3DPts = area_road.inliers.array();
        Project(roadway3DPts, roadwayPtsImg);
        for (int i = 0; i < roadwayPtsImg.rows(); ++i) {
            int xImg = roadwayPtsImg(i, 1), yImg = roadwayPtsImg(i, 0);
            if (isValid(xImg, yImg, _segImg.rows, _segImg.cols)) {
                circle(roadMaskImg, Point(yImg, xImg), 5, Scalar(255, 255, 0), CV_FILLED);
            }
        }

        morphologyEx(roadMaskImg, roadMaskImg, MORPH_DILATE,
                     getStructuringElement(MORPH_CROSS, Size(10, 10)));

        for (int i = 0; i < _segImg.rows; i++) {
            for (int j = 0; j < _segImg.cols; j++) {
                if (_segImg.at<unsigned char>(i, j) > 0 && roadMaskImg.at<unsigned char>(i, j) > 0) {
                    roadMaskImg.at<unsigned char>(i, j) = 255;
                }
            }
        }

        if (m_debug && refineCurbSet.rows() < m_minSamples)
            cout << "Warning: Not enough refined curb points (LMintersection:: RoadEdgePtsDetector)" << endl;

        if (m_debug) {
            cout << "Exiting LMintersection::RoadEdgePtsDetector()" << endl;
        }
    }

    void
    LMsintersection::computeRoadEdgePts(Eigen::MatrixXf roadlidarpts, Eigen::MatrixXf _model, float radius,
                                        Eigen::MatrixXd &refineCurbSet) {
        if (m_debug)
            cout << "Entering LMintersection:: computeRoadEdgePts" << endl;

        Eigen::MatrixXd curbPointSet;
        getRoadEdgePoints(roadlidarpts, _model, curbPointSet);
        if (curbPointSet.rows() > 0) {
            edgePtsRefiner(curbPointSet, refineCurbSet, radius);
        }

        if (m_debug)
            cout << "Existing LMintersection:: computeRoadEdgePts" << endl;
    }


    void
    LMsintersection::getRoadEdgePoints(Eigen::MatrixXf roadlidarpts, Eigen::MatrixXf _model,
                                       Eigen::MatrixXd &_curbPointSet) {
        if (m_debug)
            cout << "Entering LMintersection:: getRoadEdgePoints" << endl;

        cv::Mat _curbpts;
        if (_model.rows() == 4) {
            roadlidarpts.conservativeResize(roadlidarpts.rows(), roadlidarpts.cols() + 1);
            roadlidarpts.col(roadlidarpts.cols() - 1) = Eigen::VectorXf::Ones(roadlidarpts.rows(), 1);
            float coeffsum = sqrt(std::pow(_model(0, 0), 2) + std::pow(_model(1, 0), 2) + std::pow(_model(2, 0), 2));
            if (fabs(coeffsum) != 0 && roadlidarpts.cols() != 0) {
                Eigen::VectorXf computeZ = (roadlidarpts * _model).cwiseAbs() / coeffsum;
                Eigen::ArrayXf distance =
                        computeZ.array() + 1.2 * pow(10, -3) * (((roadlidarpts.col(0)).array()).abs2() +
                                                                ((roadlidarpts.col(1)).array()).abs2() +
                                                                ((roadlidarpts.col(
                                                                        2)).array()).abs2()).sqrt()
                        - 0.1 * pow(10, -3);

                for (int i = 0; i < distance.rows(); ++i) {
                    if (distance(i, 0) > m_dThrehold &&
                        distance(i, 0) < m_curbH) { //distance assumption for curb points
                        cv::Mat inlierpts = (cv::Mat_<float>(1, 3) << roadlidarpts(i, 0), roadlidarpts(i,
                                                                                                       1), roadlidarpts(
                                i, 2));
                        _curbpts.push_back(inlierpts);
                    }
                }

            }
        } else {
            int lrows = roadlidarpts.rows(), lcols = _model.rows();
            Eigen::MatrixXf inlierroadpts, computeZ = Eigen::MatrixXf::Ones(lrows, lcols);
            if (lrows != 0 && lcols != 0) {
                computeZ.col(1) = roadlidarpts.col(0);
                computeZ.col(2) = roadlidarpts.col(1);
                computeZ.col(3) = (roadlidarpts.col(0)).array().pow(2);
                computeZ.col(4) = (roadlidarpts.col(0)).array() * (roadlidarpts.col(1)).array();
                computeZ.col(5) = ((roadlidarpts.col(1)).array()).pow(2);
                computeZ.col(6) = ((roadlidarpts.col(0)).array()).pow(3);
                computeZ.col(7) = ((roadlidarpts.col(0)).array()).pow(2) * (roadlidarpts.col(1)).array();
                computeZ.col(8) = (roadlidarpts.col(0)).array() * (roadlidarpts.col(1)).array().pow(2);
                computeZ.col(9) = ((roadlidarpts.col(1)).array()).pow(3);
                Eigen::MatrixXf distance = computeZ * _model; //TODO: use vertical distance

                for (int i = 0; i < distance.rows(); ++i) {
                    if (fabs(distance(i)) > m_dThrehold && fabs(distance(i)) < m_curbH) {
                        cv::Mat inlierpts = (cv::Mat_<float>(1, 3) << roadlidarpts(i, 0), roadlidarpts(i,
                                                                                                       1), roadlidarpts(
                                i, 2));
                        _curbpts.push_back(inlierpts);
                    }
                }
            }
        }

        cv::cv2eigen(_curbpts, _curbPointSet);

        if (m_debug)
            cout << "Existing LMintersection:: getRoadEdgePoints" << endl;
    }


    void LMsintersection::shortestDistance(Eigen::Vector3d point, Eigen::MatrixXd _set, float radius,
                                           Eigen::MatrixXd &_nearestPoints, std::vector<double> &distance) {

        // _set is nX2 or nX3 matrix
        Eigen::MatrixXd _sets = Eigen::MatrixXd(_set.rows(), _set.cols());
        std::vector<Eigen::RowVector3d> radiusMat;
        _sets = _set;

        if (_set.cols() < 3) {
            Eigen::VectorXd pointZ;
            pointZ = Eigen::MatrixXd::Zero(_sets.rows(), 1);
            _sets.conservativeResize(_sets.rows(), _sets.cols() + 1);
            _sets.col(_sets.cols() - 1) = pointZ;
        }

        for (int i = 0; i < _sets.rows(); ++i) {
            float dis = sqrt(
                    pow(point(0) - _sets(i, 0), 2) + pow(point(1) - _sets(i, 1), 2) + pow(point(2) - _sets(i, 2), 2));
            if (fabs(dis) <= fabs(radius)) {
                radiusMat.push_back(_sets.row(i));
                distance.push_back(fabs(dis) * fabs(dis));
            }
        }

        _nearestPoints = Eigen::MatrixXd(radiusMat.size(), 3);

        for (int i = 0; i < radiusMat.size(); ++i) {
            _nearestPoints.row(i) = radiusMat[i];
        }

    }

    void
    LMsintersection::edgePtsRefiner(const Eigen::MatrixXd _curbPointSet, Eigen::MatrixXd &_refineCurbSet,
                                    float radius) {

        if (m_debug)
            cout << "Entering LMintersection:: edgePtsRefiner" << endl;

        double angle, PointX = 0, PointY = 0, PointZ = 0, distance, anglesum, distancesum;
        double minVal, maxVal, mineigenvalue, maxeigenvalue;
        cv::Mat _refinecurbset, _newcurbset, _weightMat;
        Eigen::MatrixXd newcurbset;
        Eigen::Vector3d CenterPts, NeighborPts;
        Eigen::Array33d weightMat;
        vector<unsigned char> tempCluster;
        vector<float> curvature;
        double meanx = 0, meany = 0, meanz = 0, avgdis = 0;

        for (int j = 0; j < _curbPointSet.rows(); ++j) {

            double depth = sqrt(
                    pow(_curbPointSet(j, 0), 2) + pow(_curbPointSet(j, 1), 2) + pow(_curbPointSet(j, 2), 2));
            Eigen::MatrixXd nearestPoints;
            std::vector<double> distance;
            shortestDistance(_curbPointSet.row(j).transpose(), _curbPointSet, radius, nearestPoints, distance);

            if (nearestPoints.rows() > 0) {
                Eigen::VectorXf weight(nearestPoints.rows());
                //get the center point
                for (size_t i = 0; i < nearestPoints.rows(); ++i) {
                    meanx += nearestPoints(i, 0);
                    meany += nearestPoints(i, 1);
                    meanz += nearestPoints(i, 2);
                    avgdis += distance[i];
                }
                if (nearestPoints.rows() != 0) {
                    avgdis = distancesum / nearestPoints.rows();
                    meanx = PointX / nearestPoints.rows();
                    meany = PointY / nearestPoints.rows();
                    meanz = PointZ / nearestPoints.rows();
                }

                CenterPts << meanx,
                        meany,
                        meanz;

                // get the metric
                for (size_t i = 0; i < nearestPoints.rows(); ++i) {
                    if (distance[i] >= avgdis)
                        weight(i) = exp(-distance[i] / pow(avgdis, 2));
                    else weight(i) = 1.0;
                    NeighborPts = nearestPoints.row(i).transpose();
                    weightMat +=
                            weight(i) * ((NeighborPts - CenterPts) * ((NeighborPts - CenterPts).transpose())).array();
                }

                Eigen::MatrixXf _tempMat = (weightMat.matrix()).cast<float>();


                cv::eigen2cv(_tempMat, _weightMat);
                cv::PCA _PCAweight(_weightMat, cv::Mat(), CV_PCA_DATA_AS_ROW, 0);
                cv::Mat _weightEigenValue = _PCAweight.eigenvalues;
                curvature.push_back(_weightEigenValue.at<float>(2, 0) /
                                    (_weightEigenValue.at<float>(0, 0) + _weightEigenValue.at<float>(1, 0) +
                                     _weightEigenValue.at<float>(2, 0)));
            }
        }

        // uniform [0, 1]
        auto _minc = *std::min_element(curvature.begin(), curvature.end());
        auto _maxc = *std::max_element(curvature.begin(), curvature.end());
        for (auto &s : curvature)
            s = (s - _minc) / (_maxc - _minc);

        for (int j = 0; j < _curbPointSet.rows(); ++j) {
            if (curvature[j] > m_surfCurvePara) {
                cv::Mat temp = (cv::Mat_<double>(1, 3) << _curbPointSet(j, 0), _curbPointSet(j, 1), _curbPointSet(j,
                                                                                                                  2));
                _newcurbset.push_back(temp);
            }
        }
        cv2eigen(_newcurbset, newcurbset);

        if (newcurbset.rows() != 0) {
            Eigen::MatrixXd resMat(newcurbset.rows(), 2);
            for (int j = 0; j < newcurbset.rows(); ++j) {
                anglesum = 0, distancesum = 0; // reset the sum
                double depth = sqrt(pow(newcurbset(j, 0), 2) + pow(newcurbset(j, 1), 2) + pow(newcurbset(j, 2), 2));

                Eigen::MatrixXd newestPts;
                std::vector<double> pointRadiusSquaredDistance;

                shortestDistance(newcurbset.row(j).transpose(), newcurbset, radius, newestPts,
                                 pointRadiusSquaredDistance);

                if (newestPts.rows() > 0) {
                    for (size_t i = 0; i < newestPts.rows(); ++i) {
                        PointX = newestPts(i, 0);
                        PointY = newestPts(i, 1);
                        PointZ = newestPts(i, 2);
                        angle = atan2(abs(newcurbset(j, 1) - PointY), abs(newcurbset(j, 0) - PointX)) * 180 / PI +
                                atan2(abs(newcurbset(j, 2) - PointZ), abs(newcurbset(j, 0) - PointX)) * 180 / PI;;
                        distance = pointRadiusSquaredDistance[i];
                        anglesum += angle;
                        distancesum += distance;
                    }
                }

                // get the average distance and angle for one search point
                if (newcurbset.rows() != 0) {
                    resMat(j, 0) = distancesum / (newcurbset.rows() * depth);
                    resMat(j, 1) = anglesum / newcurbset.rows();
                } else {
                    resMat(j, 0) = 0;
                    resMat(j, 1) = 0;
                }
                if (std::isnan(resMat(j, 0))) resMat(j, 0) = 0;
                if (std::isnan(resMat(j, 1))) resMat(j, 1) = 0;
            }

            // get the threshold value for all the curb points
            double mind = resMat.col(0).minCoeff(), maxd = resMat.col(0).maxCoeff();
            double mina = resMat.col(1).minCoeff(), maxa = resMat.col(1).maxCoeff();

            if (maxd - mind != 0)
                resMat.col(0) = ((resMat.col(0).array() - mind) / (maxd - mind)).matrix();

            if (maxa - mina != 0)
                resMat.col(1) = ((resMat.col(1).array() - mina) / (maxa - mina)).matrix();

            Eigen::VectorXd disMat = resMat.col(0) + resMat.col(1);
            double minm = disMat.minCoeff(), maxm = disMat.maxCoeff();
            disMat = ((disMat.array() - minm) * 255 / (maxm - minm)).matrix();

            for (int j = 0; j < disMat.rows(); ++j)
                tempCluster.push_back(static_cast<unsigned char>(disMat(j, 0)));

            minMaxLoc(tempCluster, &minVal, &maxVal);
            double thres = threshold(tempCluster, tempCluster, minVal, maxVal, CV_THRESH_TOZERO | CV_THRESH_OTSU);


            for (int j = 0; j < newcurbset.rows(); ++j) {
                if (disMat(j, 0) >= m_surfCurvePara * 255.0) {
                    cv::Mat temp = (cv::Mat_<double>(1, 3) << newcurbset(j, 0), newcurbset(j, 1), newcurbset(j, 2));
                    _refinecurbset.push_back(temp);
                }
            }
        } else cout << "Not curb points in LMsintersection::edgePtsRefiner!" << endl;

        if (newcurbset.rows() != 0)
            cv::cv2eigen(_refinecurbset, _refineCurbSet);

        if (m_debug)
            cout << "Existing LMintersection:: edgePtsRefiner" << endl;

    }
}