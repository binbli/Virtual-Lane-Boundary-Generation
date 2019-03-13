#include "LMsfromCam.h"
#include<fstream>
#include<stdlib.h>

namespace LD {

    using namespace cv;

    LMsfromCam::LMsfromCam(string _xmlFile) : Solver(_xmlFile) {

        if (m_debug)
            cout << "Entering LMsfromCam::LMsfromCam()" << endl;

        ParseXML();

        if (m_debug)
            cout << "Exiting LMsfromCam::LMsfromCam()" << endl;
    }

    void LMsfromCam::ParseXML() {
        m_xml = m_xml.child("LMsfromCam");
        pugi::xml_node solverInstance = m_xml.child("SolverInstance");

        m_dataRoot = solverInstance.attribute("dataRoot").as_string();
        m_dataFile = solverInstance.attribute("dataFile").as_string();
        m_segRoot = solverInstance.attribute("segRoot").as_string();
        m_refinedRoot = solverInstance.attribute("refinedRoot").as_string();
        m_vizImgPrefix = solverInstance.attribute("vizImgPrefix").as_string();
        m_saveVizImg = solverInstance.attribute("saveVizImg").as_bool(true);
        m_refinedImgPrefix = solverInstance.attribute("refinedImgPrefix").as_string();
        m_roadbRoot = solverInstance.attribute("roadbRoot").as_string();

        if (m_dataRoot.empty() || m_dataFile.empty() || m_segRoot.empty() || m_refinedRoot.empty() ||
            (m_saveVizImg && m_vizImgPrefix.empty()) || m_refinedImgPrefix.empty() || m_roadbRoot.empty())
            throw runtime_error(
                    "One of the following attributes are missing in SolverInstance node of Segmenter: dataRoot, dataFile, segRoot, refinedRoot, vizImgPrefix, saveVizImg, refinedImgPrefix, m_roadbRoot");
    }

    void LMsfromCam::operator()(const Mat &_original, const Mat &_segImg, Mat &_refinedImg) {
        if (m_debug)
            cout << "Entering LMsfromCam::()" << endl;

        Mat preprocessed;
        LMsfromImg(_original, _segImg, _refinedImg);

        if (m_debug)
            cout << "Exiting LMsfromCam::()" << endl;
    }

    void LMsfromCam::Run() {

        if (m_debug)
            cout << "Entering LMsfromCam::Run()" << endl;

        std::ifstream fin(m_dataFile.c_str());
        string line;
        cv::Mat edgeImg, refinedImg;
        while (std::getline(fin, line)) {

            m_imgBaseName = string(basename(const_cast<char *>(line.c_str())));
            string fullImgName = m_dataRoot + "/" + line, fullSegName = m_segRoot + "/segmented_" + line;

            Mat original = imread(fullImgName);
            Mat segImg = imread(fullSegName, IMREAD_GRAYSCALE);

            if (original.empty()) {
                cout << "could not open or find " << fullImgName << endl;
                exit(1);
            }
            if (segImg.empty()) {
                cout << "could not open or find " << fullSegName << endl;
                exit(1);
            }

            if (m_debug)
                cout << "Successfully opened " << fullImgName << " and " << fullSegName << endl;

            this->operator()(original, segImg, refinedImg);

            if (m_debug) {
                imwrite(m_refinedRoot + "/" + m_refinedImgPrefix + m_imgBaseName, refinedImg);
                cout << "Successfully saved LM's image mask at "
                     << m_refinedRoot + "/" + m_refinedImgPrefix + m_imgBaseName << endl;

            }

            //save overlayed image
            if (m_saveVizImg) {
                Mat regions, overlayed;
                string overlayedName = m_refinedRoot + "/" + m_vizImgPrefix + m_imgBaseName;
                vector<Mat> channels;
                channels.push_back(Mat::zeros(refinedImg.size(), CV_8UC1));
                channels.push_back(refinedImg);
                channels.push_back(Mat::zeros(refinedImg.size(), CV_8UC1));
                cv::merge(channels, regions);
                cv::addWeighted(regions, 0.5, original, 0.5, 0, overlayed);
                imwrite(overlayedName, overlayed);
            }

        }

        if (m_debug)
            cout << "Exiting LMsfromCam::Run()" << endl;
    }


}
