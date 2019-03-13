#include"../include/SurfaceDataMaker.h"

namespace LD {
    void SurfaceDataMaker::ParseXML() {
        m_xml = m_xml.child("SurfaceDataMaker");

        m_RoadSegRoot = m_xml.attribute("segRoot").as_string();
        m_segImgPrefix = m_xml.attribute("segImgPrefix").as_string();
        m_outputFilePrefix = m_xml.attribute("outputFilePrefix").as_string();
        m_saveVizImg = m_xml.attribute("saveVizImg").as_bool(true);
        m_vizImgPrefix = m_xml.attribute("vizImgPrefix").as_string();
        m_minPoints = m_xml.attribute("minPoints").as_int();
        m_printProjectedPts = m_xml.attribute("printProjectedPts").as_bool();

        if (m_RoadSegRoot.empty() || m_segImgPrefix.empty() || m_outputFilePrefix.empty() ||
            (m_saveVizImg && m_vizImgPrefix.empty()) || !m_minPoints)
            throw runtime_error(
                    "at least one of the following attributes are missing in SurfaceDataMaker node: segRoot, saveVizImg, outputFile, segImgPrefix, vizImgPrefix, minPoints");
    }

    void SurfaceDataMaker::ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints,
                                                    Mat &_reflectivity, Mat &_inputImg) {

        if (m_debug)
            cout << "Entering ProcessProjectedLidarPts()" << endl;

        string segImgName = m_RoadSegRoot + "/" + m_segImgPrefix + m_imgBaseName;
        Mat segImg = imread(segImgName, IMREAD_GRAYSCALE);

        if (segImg.empty())
            throw std::runtime_error("Can't open " + segImgName);
        else if (m_debug)
            cout << "Successfully read segmented image: " << segImgName << endl;

        string outFile =
                m_outputRoot + "/" + m_outputFilePrefix + m_imgBaseName.substr(0, m_imgBaseName.size() - 3) + "txt";
        std::ofstream fout(outFile);

        for (int i = 0; i < _veloImg.rows(); i++) {
            int x = _veloImg(i, 0), y = _veloImg(i, 1);
            float reflect = _reflectivity.at<float>(i, 0);
            if (isValid(y, x, segImg.rows, segImg.cols) && segImg.at<unsigned char>(y, x)) {
                if (m_printProjectedPts)
                    fout << _veloImg(i, 0) << "\t" << _veloImg(i, 1) << endl;
                else
                    fout << _veloPoints.at<float>(i, 0) << "\t" << _veloPoints.at<float>(i, 1) << "\t"
                         << _veloPoints.at<float>(i, 2) << endl;

                if (m_saveVizImg)
                    circle(_inputImg, Point(x, y), 5, Scalar(255 * reflect, 0, 0));
            }
        }
        fout.flush();

        if (m_saveVizImg)
            imwrite(m_outputRoot + "/" + m_vizImgPrefix + m_imgBaseName, _inputImg);

        if (m_debug)
            cout << "Exiting ProcessProjectedLidarPts()" << endl;
    }
}
