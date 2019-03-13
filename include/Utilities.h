#ifndef UTILITIES_H_
#define UTILITIES_H_

#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/SVD>
#include <eigen3/Eigen/Geometry>
#include"BaseLD.h"
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/core/eigen.hpp>
#include <opencv2/core/eigen.hpp>
#include"interpolation.h"
#include<omp.h>

namespace LD {

    struct retresult {
        Eigen::MatrixXf modelcoff;
        Eigen::MatrixXf inliers;
        long int ninliers;
    };

    bool isValid(long long int r, long long int c, long long int rows, long long int cols);

    void ReadEigenMatFromFile(std::ifstream &_fin, Eigen::ArrayXXf &_data, bool _shouldTranspose = false);

    void ReadEigenMatFromFile(std::ifstream &_fin, Eigen::ArrayXf &_data, bool _shouldTranspose);

    void ReadEigenMatFromFile(const string &_fileName, Eigen::ArrayXXf &_data, bool _shouldTranspose = false);

    void ReadEigenArrayVecFromFile(std::ifstream &_fin, vector<Eigen::ArrayXf> &_data);

    void readCurbRT(const string filepath, Eigen::Matrix3f &Rotation, Eigen::RowVector3f &translation, Eigen::MatrixXf &curbs);

    void CreateAlglibArray(const vector<Eigen::ArrayXf> &_samples, vector<alglib::real_1d_array> &_coordinates);

    void CreateAlglibArray(const Eigen::ArrayXXf _samples, vector<alglib::real_1d_array> &_coordinates);

    void ShiftCenterCtrlPts(const Eigen::ArrayXXf& _centerCtrlPts, const float& _horizontalDist, const float& _upwardDist, const bool& _isRightHigher, vector<Eigen::ArrayXXf>& _laneCtrlPts);

    Eigen::MatrixXf roadptsview(cv::Mat roadpts, Eigen::MatrixXf projection_mat);

    cv::Mat getroadPts(retresult model, cv::Mat velo_points_, float dThreshold, float dx, float dy);

    void newlaneptsget(Eigen::MatrixXf road_img_newpts, cv::Mat roadpts, cv::Mat &_refined_img, double _thresh,
                       cv::Mat &_img);

    Eigen::MatrixXf computeModel(Eigen::MatrixXf &fitdata);

    retresult mcdcomputeres(Eigen::MatrixXf modelcoff, Eigen::MatrixXf lidarpts, float thresdistance, bool status);

    retresult runRansacRoad(Eigen::MatrixXf &lidarpts, int nSamNumForModel, int maxtrials, int ninlierThreshold,
                            float dThreshold);

    Eigen::MatrixXf computePlaneModel(Eigen::MatrixXf &fitdata);

    retresult mcdPlaneComput(Eigen::MatrixXf modelcoff, Eigen::MatrixXf lidarpts, float thresdistance, bool status);

    retresult runRansacplanefit(Eigen::MatrixXf &lidarpts, int nSamNumForModel, int maxtrials, int ninlierThreshold,
                                float dThreshold);

    Eigen::MatrixXf ptSelect(cv::Mat roadsurfLidarpts, float dx = 30.0, float dy = 3.0);

    void coordinatetransL2G(Eigen::MatrixXf roadlidarpts, Eigen::Matrix3f &RotationL2G, Eigen::RowVector3f &translationL2G);

    void PtsRemoveBackfCar(Eigen::ArrayXXf gpsPts, float xmin, float m_maxX, Eigen::ArrayXXf &gpsPtsInL);

    class AverageMeter {
    private:
        float m_sum;
        vector<float> m_vals;
        size_t m_index;

    public:

        AverageMeter(int _numToTrack);

        void Update(const float& _val);

        float Average();
    };

    void CurrentOptCtrlPts(Eigen::ArrayXXf _PrectrlPts, Eigen::ArrayXXf _CurctrlPts, Eigen::ArrayXXf &_updatedNewCtrolPts);

    void shortestDistance(Eigen::VectorXf point, Eigen::MatrixXf _set, float &_dis);

    void shortestDistance(Eigen::Vector3f point, Eigen::MatrixXf _set, float radius,
                          Eigen::MatrixXf &_nearestPoints, std::vector<float> &distance);

    int ImgFile2Int(const string& _imgFileName);

}

#endif
