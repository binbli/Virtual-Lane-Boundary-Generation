#include"../include/LaneQualityChecker.h"
#include<unordered_set>
#include"../include/Utilities.h"
#include<sstream>
#include<iomanip>

namespace LD {

    LaneQualityChecker::LaneQualityChecker(string _xmlFile) : VeloPtsProjectCam(_xmlFile),
                                                              finClusters(m_clustersFile.c_str()),
                                                              finPts(m_ptsFile.c_str()) {
		ParseXML();
    }
	
	void LaneQualityChecker::ParseXML() {
        if (m_debug)
			cout << "Entering LaneQualityChecker::ParseXML()" << endl;

        m_xml = m_xml.child("LaneQualityChecker");
		m_noiseID = m_xml.attribute("noiseID").as_int(-1);

        if (m_debug)
			cout << "Exiting LaneQualityChecker::ParseXML()" << endl;
	}

    string toString(const float &_f1, const float &_f2, const float &_f3) {
		std::ostringstream os;
		os << std::fixed << std::setprecision(2) << _f1 << "-" << _f2 << "-" << _f3;
		return os.str();
	}

    void LaneQualityChecker::operator()(const Eigen::ArrayXXf &_intersectedPts, const Eigen::ArrayXf &_clusters,
                                        const cv::Mat &_veloPoints, const Mat &_reflectivity,
                                        const Eigen::MatrixXf &_veloImg, const Mat &_inputImg, const Mat &_segmentedImg,
                                        const Mat &_refinedImg, float &_brightnessRatio, float &_reflectivityRatio) {

        if (m_debug)
			cout << "Entering LaneQualityChecker::()" << endl;

        assert(_intersectedPts.rows() == 3);
		assert(_veloImg.cols() == 2);
		assert(_veloPoints.cols == 4);

		std::unordered_set<string> lanePtsSet;
		Mat grayImg;
		cvtColor(_inputImg, grayImg, CV_BGR2GRAY);

        for (ulli i = 0; i < _intersectedPts.cols(); i++)
            if (_clusters(i) != m_noiseID)
				lanePtsSet.insert(toString(_intersectedPts(0, i), _intersectedPts(1, i), _intersectedPts(2, i)));

		float laneReflectivityMean = 0, roadReflectivityMean = 0;
		float laneBrightnessMean = 0, roadBrightnessMean = 0;
		ulli roadPts = 0, lanePts = 0;
        for (ulli i = 0; i < _veloPoints.rows; i++) {
			int c = _veloImg(i, 0), r = _veloImg(i, 1);
            if (isValid(r, c, grayImg.rows, grayImg.cols)) {
                if (lanePtsSet.find(toString(_veloPoints.at<float>(i, 0), _veloPoints.at<float>(i, 1),
                                             _veloPoints.at<float>(i, 2))) != lanePtsSet.end()) {
					laneReflectivityMean += _reflectivity.at<float>(i, 0);
					lanePts++;
                } else if (_segmentedImg.at<unsigned char>(r, c)) {
					roadReflectivityMean += _reflectivity.at<float>(i, 0);
					roadPts++;
				}
			}
		}

        Mat lanes;
		bitwise_and(grayImg, _refinedImg, lanes);
		bitwise_and(grayImg, _segmentedImg, grayImg);
		int lanePtsImg = countNonZero(lanes);
		laneBrightnessMean = (float) sum(lanes)[0];
		roadBrightnessMean = (sum(grayImg)[0] - laneBrightnessMean) / (countNonZero(grayImg) - lanePtsImg);
		laneBrightnessMean /= lanePtsImg;
		laneReflectivityMean /= lanePts;
        roadReflectivityMean /= roadPts;
		_reflectivityRatio = laneReflectivityMean / roadReflectivityMean;
        _brightnessRatio = laneBrightnessMean / roadBrightnessMean;

        if (m_debug)
			cout << "Exiting LaneQualityChecker::()" << endl;
	}

}
