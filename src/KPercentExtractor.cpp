#include"../include/KPercentExtractor.h"
#include "../include/Utilities.h"

namespace LD {

    void KPercentExtractor::ParseXML() {
        m_xml = m_xml.child("KPercentExtractor");
        m_kpercent = m_xml.attribute("kpercent").as_float();
        m_dilateSize = m_xml.attribute("dilateSize").as_int();
        if (!m_kpercent && !m_dilateSize)
            throw runtime_error("at least one of the following attribute is missing: k or m_dilateSize");
    }

    void KPercentExtractor::runBFS(cv::Mat &_noisefreeImg, int r, int c, int &_nClusters, vector<cv::Point> &_pts) {
        std::queue<Point> q;
        q.push(Point(c, r));
        _nClusters = 0;
        Rect imgRect(0, 0, _noisefreeImg.cols, _noisefreeImg.rows);
        while (q.size()) {
            Point cur = q.front();
            q.pop();
            if (imgRect.contains(cur) && _noisefreeImg.at<uchar>(cur) > 0) {
                _noisefreeImg.at<uchar>(cur) = 0;
                _pts.push_back(cur);
                _nClusters++;
                q.push(Point(cur.x - 1, cur.y));
                q.push(Point(cur.x + 1, cur.y));
                q.push(Point(cur.x, cur.y - 1));
                q.push(Point(cur.x, cur.y + 1));
            }
        }
    }

    int KPercentExtractor::noiseRemover(cv::Mat &_segImg) {

        if (m_debug)
            cout << "Entering Refiner::noiseRemover()" << endl;

        if (_segImg.rows == 0 || _segImg.cols == 0)
            return 0;

        int res = 0, nCluster = 0;
        cv::Mat _noisefreeImg = _segImg.clone(), tempImg = _segImg.clone();
        vector<cv::Point> pts;
        vector<int> cluster;
        vector<unsigned char> tempCluster;
        vector<vector<cv::Point>> PtsMat;
        int width = _segImg.cols;
        int height = _segImg.rows;

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (_segImg.at<uchar>(i, j)) {
                    ++res;
                    runBFS(_noisefreeImg, i, j, nCluster, pts);
                    if (nCluster > 0) {
                        cluster.push_back(nCluster);
                        PtsMat.push_back(pts);
                    }
                    nCluster = 0;
                    pts.clear();
                }
            }
        }

        vector<int> sortedCluster(cluster.begin(), cluster.end());
        std::sort(sortedCluster.begin(), sortedCluster.end(), std::greater<int>());

        if (!sortedCluster.empty()) {
            for (int k = 0; k < cluster.size(); ++k) {
                if (cluster[k] < sortedCluster[0] && cluster[k] > 0 && !PtsMat[k].empty()) {
                    for (auto s : PtsMat[k]) {
                        tempImg.at<uchar>(s) = 0;
                    }
                }
            }
        }
        _segImg = tempImg.clone();

        if (m_debug)
            cout << "Exiting Refiner::noiseRemover()" << endl;

        return res;
    }


    void KPercentExtractor::Preprocess(const Mat &_original, const Mat &_segImg, Mat &_preprocessed) {
        if (m_debug)
            cout << "Entering KPercentExtractor::Preprocess()" << endl;
        if (_original.channels() == 3)
            cvtColor(_original, _preprocessed, COLOR_RGB2GRAY);
        else _preprocessed = _original;
        bitwise_and(_preprocessed, _segImg, _preprocessed);

        if (m_debug)
            cout << "Exiting KPercentExtractor::Preprocess()" << endl;
    }

    void KPercentExtractor::LMsfromImg(const Mat &_original, const Mat &_segImg, Mat &_refinedImg) {

        if (m_debug)
            cout << "Entering KPercentExtractor::Refine()" << endl;

        Mat extractedImg;
        Preprocess(_original, _segImg, extractedImg);
        //get top k% pixels
        Mat flattened = extractedImg.reshape(1, 1).clone();
        if (flattened.isContinuous()) {
            std::sort(flattened.data, flattened.data + flattened.total());
            int numZeros = std::distance(flattened.datastart,
                                         std::upper_bound(flattened.datastart, flattened.dataend, 0));
            int threshIndex = ((flattened.total() - numZeros) * (100 - m_kpercent)) / 100;
            int thresh = flattened.data[numZeros + threshIndex] - 1;

            if (m_debug)
                cout << "Threshold set to " << thresh << endl;

            threshold(extractedImg, _refinedImg, thresh, 255, THRESH_BINARY);

            morphologyEx(_refinedImg, _refinedImg, MORPH_OPEN,
                         getStructuringElement(MORPH_CROSS, Size(m_dilateSize, m_dilateSize))); // remove noise
        } else
            throw runtime_error("Matrix is not continuous");//ideally it should always be continuous


        if (m_debug)
            cout << "Exiting KPercentExtractor::Refine()" << endl;
    }

}