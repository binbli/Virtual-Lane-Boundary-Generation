#include "VeloPtsProjectCam.h"

namespace LD {

    VeloPtsProjectCam::VeloPtsProjectCam(string _xmlFile) : Solver(_xmlFile), m_calibDataLoader(_xmlFile) {

        if (m_debug)
            cout << "Entering VeloPtsProjectCam::VeloPtsProjectCam() " << endl;

        ParseXML();
        string camCalibFile = m_calibRoot + "/calib_cam_to_cam.txt";
        string veloCalibFile = m_calibRoot + "/calib_velo_to_cam.txt";

        bool isSuccess = true;

        Eigen::MatrixXf R, T;

        isSuccess &=
                m_calibDataLoader.ReadVariable(camCalibFile, "P_rect_0" + std::to_string(m_camNum), 3, 4, m_PRect) &&
                m_calibDataLoader.ReadVariable(camCalibFile, "R_rect_00", 3, 3, m_RRect) &&
                m_calibDataLoader.ReadVariable(veloCalibFile, "R", 3, 3, R) &&
                m_calibDataLoader.ReadVariable(veloCalibFile, "T", 3, 1, T);

        if (!isSuccess)
            throw std::runtime_error("Incorrect format of calibration files");

        m_Tr = Eigen::MatrixXf(R.rows() + 1, R.cols() + T.cols());

        m_Tr << R, T,
                0, 0, 0, 1;

        ComputeProjMat();

        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::VeloPtsProjectCam() " << endl;
    }

    void VeloPtsProjectCam::ParseXML() {
        if (m_debug)
            cout << "Entering VeloPtsProjectCam::ParseXML()" << endl;

        m_xml = m_xml.child("VeloPtsProjectCam");

        pugi::xml_node solverInstance;
        solverInstance = m_xml.child("SolverInstance");
        m_dataRoot = solverInstance.attribute("dataRoot").as_string();
        m_dataFile = solverInstance.attribute("dataFile").as_string();
        m_veloRoot = solverInstance.attribute("veloRoot").as_string();
        m_calibRoot = solverInstance.attribute("calibRoot").as_string();
        m_outputRoot = solverInstance.attribute("outputRoot").as_string();
        m_retentionFrequency = solverInstance.attribute("retentionFrequency").as_int();
        m_camNum = solverInstance.attribute("camNum").as_int(-1);
        m_minX = solverInstance.attribute("minX").as_int();

        if (m_dataRoot.empty() || m_dataFile.empty() || m_veloRoot.empty() || m_calibRoot.empty() ||
            m_outputRoot.empty() || !m_retentionFrequency || (m_camNum == -1) || !m_minX)
            throw runtime_error(
                    "at least one of the following attributes are missing in SolverInstance node of VeloProject: dataRoot, dataFile, segRoot, refinedRoot, veloRoot, calibRoot, outputRoot, retentionFrequency, camNum, minX, outputFile, segImgPrefix, refImgPrefix");

        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::ParseXML()" << endl;
    }


    void VeloPtsProjectCam::Run() {
        if (m_debug)
            cout << "Entering VeloPtsProjectCam::Run() " << endl;

        std::ifstream fin(m_dataFile.c_str());
        string line;
        Eigen::MatrixXf projectionMat, veloImgPts;
        Mat veloPoints, reflectivity;
        while (std::getline(fin, line)) {

            m_imgBaseName = string(basename(const_cast<char *>(line.c_str())));
            string inputImgName = m_dataRoot + "/" + line;

            Mat inputImg = imread(inputImgName);

            if (inputImg.empty())
                throw std::runtime_error("Can't open " + inputImgName);
            else if (m_debug)
                cout << "Successfully read input image: " << inputImgName << endl;

            ReadVeloData(m_veloRoot + "/" + line.substr(0, line.size() - 3) + "bin", veloPoints);
            Project(veloPoints, veloImgPts, reflectivity);
            ProcessProjectedLidarPts(veloImgPts, veloPoints, reflectivity, inputImg);
        }

        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::Run() " << endl;
    }


    void VeloPtsProjectCam::ReadVeloData(string _binFile, Mat &_veloPoints) {
        //taken from KITTI website
        if (m_debug)
            cout << "Entering VeloPtsProjectCam::ReadVeloData() " << endl;

        _veloPoints.release();
        // allocate 4 MB buffer (only ~130*4*4 KB are needed)
        int32_t num = 1000000;
        float *data = (float *) malloc(num * sizeof(float));

        // pointers
        float *px = data + 0;

        // load point cloud
        FILE *stream;
        stream = fopen(_binFile.c_str(), "rb");
        num = fread(data, sizeof(float), num, stream) / 4;
        for (int32_t i = 0; i < num; i++) {
            if (*px >= m_minX && (i % m_retentionFrequency) == 0) {
                float x = *px;
                float y = *(px +1);
                float z = *(px +2);
                float intensity = *(px +3 );
                Mat m = (Mat_<float>(1, 4)<< x, y, z, intensity);
                _veloPoints.push_back(m);
            }
            px += 4;
        }
        fclose(stream);
        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::ReadVeloData() " << endl;
    }

    void VeloPtsProjectCam::ComputeProjMat() {
        if (m_debug)
            cout << "Entering VeloPtsProjectCam::ComputeProjMat() " << endl;

        Eigen::MatrixXf RCamToRect = Eigen::MatrixXf::Identity(4, 4);

        RCamToRect.topLeftCorner<3, 3>() = m_RRect;

        m_projectionMat = m_PRect * RCamToRect * m_Tr;
        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::ComputeProjMat() " << endl;
    }

    void VeloPtsProjectCam::Project(const Mat &_veloPoints, Eigen::MatrixXf &_veloImg, Mat &_reflectivity) {

        if (m_debug)
            cout << "Entering VeloPtsProjectCam::Project() " << endl;

        int dimNorm = m_projectionMat.rows();
        int dimProj = m_projectionMat.cols();

        if (_veloPoints.cols == dimProj - 1) {
            Mat col(_veloPoints.rows, 1, _veloPoints.type(), 1);
            hconcat(_veloPoints, col, _veloPoints);
        } else if (_veloPoints.cols < dimProj)
            throw std::runtime_error("incorrect dimensions to multiply");

        if (!_veloPoints.isContinuous())
            throw std::runtime_error("matrix is not continuous");

        Eigen::MatrixXf veloPtsEg;
        cv::cv2eigen(_veloPoints, veloPtsEg);
        _reflectivity = _veloPoints.col(dimProj - 1).clone();
        veloPtsEg.col(dimProj - 1) = Eigen::MatrixXf::Ones(veloPtsEg.rows(), 1);

        Eigen::MatrixXf newPoints = (m_projectionMat * veloPtsEg.transpose()).transpose();

        for (int i = 0; i < dimNorm - 1; i++)
            newPoints.col(i).array() = newPoints.col(i).array() / newPoints.col(dimNorm - 1).array();

        _veloImg = newPoints.topLeftCorner(newPoints.rows(), newPoints.cols() - 1);

        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::Project() " << endl;
    }

    void VeloPtsProjectCam::Project(Eigen::ArrayXXf &_veloPoints, Eigen::MatrixXf &_veloImg) {
        if (m_debug)
            cout << "Entering VeloPtsProjectCam::Project() " << endl;

        int dimNorm = m_projectionMat.rows();
        int dimProj = m_projectionMat.cols();

        if (_veloPoints.cols() == dimProj - 1) {
            _veloPoints.conservativeResize(Eigen::NoChange, _veloPoints.cols() + 1);
            _veloPoints.col(_veloPoints.cols() - 1) = 1;

        } else if (_veloPoints.cols() < dimProj)
            throw std::runtime_error("incorrect dimensions to multiply");

        Eigen::MatrixXf newPoints = (m_projectionMat * _veloPoints.matrix().transpose()).transpose();

        for (int i = 0; i < dimNorm - 1; i++)
            newPoints.col(i).array() = newPoints.col(i).array() / newPoints.col(dimNorm - 1).array();

        _veloImg = newPoints.topLeftCorner(newPoints.rows(), newPoints.cols() - 1);

        if (m_debug)
            cout << "Exiting VeloPtsProjectCam::Project() " << endl;

    }

}
