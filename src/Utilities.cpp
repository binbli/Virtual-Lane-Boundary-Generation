#include"../include/Utilities.h"

namespace LD {

    bool isValid(long long int r, long long int c, long long int rows, long long int cols) {
        return r >= 0 && c >= 0 && r < rows && c < cols;
    }

    void ReadEigenMatFromFile(std::ifstream &_fin, Eigen::ArrayXXf &_data, bool _shouldTranspose) {

        unsigned long long int rows, cols;
        _fin >> rows >> cols;

        _data.resize(rows, cols);
        for (unsigned long long int r = 0; r < rows; r++)
            for (unsigned long long int c = 0; c < cols; c++)
                _fin >> _data(r, c);

        if (_shouldTranspose)
            _data.transposeInPlace();

    }

    void ReadEigenMatFromFile(std::ifstream &_fin, Eigen::ArrayXf &_data, bool _shouldTranspose) {

        unsigned long long int rows, cols;
        _fin >> rows;

        _data = Eigen::ArrayXf(rows, 1);
        for (unsigned long long int r = 0; r < rows; r++)
            _fin >> _data(r, 0);

        if (_shouldTranspose)
            _data.transposeInPlace();

    }

    void ReadEigenArrayVecFromFile(std::ifstream &_fin, vector<Eigen::ArrayXf> &_data) {
        unsigned long long int rows, cols;
        _fin >> rows >> cols;
        _data = vector<Eigen::ArrayXf>(rows, Eigen::ArrayXf(cols));
        for (unsigned long long int r = 0; r < rows; r++)
            for (unsigned long long int c = 0; c < cols; c++)
                _fin >> _data[r](c);
    }

    void ReadEigenMatFromFile(const string &_fileName, Eigen::ArrayXXf &_data, bool _shouldTranspose) {
        typedef unsigned long long int ulli;
        std::ifstream fin(_fileName);
        if (!fin)
            throw runtime_error("Can't open " + _fileName);
        string line;
        vector<vector<float> > els;
        ulli rows = 0, cols = 0;
        while (std::getline(fin, line)) {
            std::istringstream is(line);
            float el;
            els.push_back(vector<float>());
            while (is >> el)
                els[rows].push_back(el);
            if (!cols)
                cols = els[rows].size();
            else if (cols != els[rows].size())
                throw runtime_error("Number of columns are inconsistent");
            rows++;
        }

        _data.resize(rows, cols);
        for (ulli r = 0; r < rows; r++)
            for (ulli c = 0; c < cols; c++)
                _data(r, c) = els[r][c];

        if (_shouldTranspose)
            _data.transposeInPlace();
    }

    void readCurbRT(const string filepath, Eigen::Matrix3f &Rotation, Eigen::RowVector3f &translation,
                    Eigen::MatrixXf &curbs) {

        std::ifstream fin(filepath);
        std::string tempstr, value;
        vector<float> valuMat;
        int index = 0;
        if (!fin) {
            cout << "Cannot open file in Utilizies::readCurbRT!" << endl;
            return;
        }
        while (std::getline(fin, tempstr)) {
            ++index;
            std::istringstream instr(tempstr);
            if (index == 1) {
                while (instr >> value) {
                    valuMat.push_back(stof(value));
                }
                if (valuMat.size() == 9) {
                    Rotation << valuMat[0], valuMat[1], valuMat[2],
                            valuMat[3], valuMat[4], valuMat[5],
                            valuMat[6], valuMat[7], valuMat[8];
                    valuMat.clear();
                } else {
                    return;
                }
            } else if (index == 2) {
                while (instr >> value) {
                    valuMat.push_back(stof(value));
                }
                if (valuMat.size() == 3) {
                    translation << valuMat[0], valuMat[1], valuMat[2];
                    valuMat.clear();
                } else {
                    return;
                }
            } else {
                while (instr >> value) {
                    valuMat.push_back(stof(value));
                }
                if (valuMat.size() != 0) {
                    curbs = Eigen::MatrixXf(valuMat.size() / 3, 3);
                    for (int i = 0; i < valuMat.size(); ++i) {
                        int row = i / 3;
                        int col = i % 3;
                        curbs(row, col) = valuMat[i];
                    }
                }
            }
        }

        if (Rotation.rows() != 3 || Rotation.cols() != 3 || translation.cols() != 3) {
            throw runtime_error("Cannot get Rotation2G and translation2G");
        }
        if (curbs.rows() == 0) {
            cout << "No curb Points TLinkage::Run()...." << endl;
        }

    }

    void CreateAlglibArray(const vector<Eigen::ArrayXf> &_samples, vector<alglib::real_1d_array> &_coordinates) {
        if (!_samples.size())
            return;

        _coordinates.resize(_samples[0].size());
        for (int i = 0; i < _coordinates.size(); i++)
            _coordinates[i].setlength(_samples.size());

        for (int i = 0; i < _coordinates.size(); i++)
            for (int j = 0; j < _samples.size(); j++)
                _coordinates[i](j) = _samples[j](i);
    }


    void CreateAlglibArray(const Eigen::ArrayXXf _samples, vector<alglib::real_1d_array> &_coordinates) {
        if (!_samples.cols())
            return;

        _coordinates.resize(_samples.cols());
        for (int i = 0; i < _coordinates.size(); i++)
            _coordinates[i].setlength(_samples.rows());

        for (int i = 0; i < _coordinates.size(); i++)
            for (int j = 0; j < _samples.rows(); j++)
                _coordinates[i](j) = _samples(j, i);
    }


// get the model coeff using SVD deccomposition
    Eigen::MatrixXf computeModel(Eigen::MatrixXf &fitdata) {
        // fitdata: The minimum number of data points to fit for the model
        cv::Mat A(fitdata.rows(), fitdata.rows(), CV_32FC1); // A*x = b, using SVD to get x
        cv::Mat b;
        Eigen::MatrixXf modelcoff;
        for (int i = 0; i < fitdata.rows(); ++i) {
            for (int j = 0; j < fitdata.rows(); ++j) {
                switch (j) {
                    case 0:
                        A.at<float>(i, j) = 1;
                        break;
                    case 1:
                        A.at<float>(i, j) = fitdata(i, 0);
                        break;
                    case 2:
                        A.at<float>(i, j) = fitdata(i, 1);
                        break;
                    case 3:
                        A.at<float>(i, j) = fitdata(i, 0) * fitdata(i, 0);
                        break;
                    case 4:
                        A.at<float>(i, j) = fitdata(i, 0) * fitdata(i, 1);
                        break;
                    case 5:
                        A.at<float>(i, j) = fitdata(i, 1) * fitdata(i, 1);
                        break;
                    case 6:
                        A.at<float>(i, j) = pow(fitdata(i, 0), 3);
                        break;
                    case 7:
                        A.at<float>(i, j) = pow(fitdata(i, 0), 2) * fitdata(i, 1);
                        break;
                    case 8:
                        A.at<float>(i, j) = fitdata(i, 0) * pow(fitdata(i, 1), 2);
                        break;
                    case 9:
                        A.at<float>(i, j) = pow(fitdata(i, 1), 3);
                        break;
                    default:
                        std::cerr << "wrong seg coeff result" << std::endl;
                }
            }

            b.push_back(fitdata(i, 2));
        }
        Eigen::MatrixXf _A;
        cv2eigen(A, _A);
        Eigen::MatrixXf _b;
        cv2eigen(b, _b);
        Eigen::JacobiSVD<Eigen::MatrixXf> svd(_A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        modelcoff = svd.solve(_b);
        return modelcoff;
    }


    retresult mcdcomputeres(Eigen::MatrixXf modelcoff, Eigen::MatrixXf lidarpts, float thresdistance, bool status) {
        // modelcoeff : parameters of the model
        // lidarpts : all the lidar 3D points
        retresult res;
        cv::Mat inliers;
        int lrows = lidarpts.rows(), lcols = modelcoff.rows();
        Eigen::MatrixXf inlierroadpts, computeZ = Eigen::MatrixXf::Ones(lrows, lcols);
        if (lrows != 0 && lcols != 0) {
            computeZ.col(1) = lidarpts.col(0);
            computeZ.col(2) = lidarpts.col(1);
            computeZ.col(3) = (lidarpts.col(0)).array().pow(2);
            computeZ.col(4) = (lidarpts.col(0)).array() * (lidarpts.col(1)).array();
            computeZ.col(5) = ((lidarpts.col(1)).array()).pow(2);
            computeZ.col(6) = ((lidarpts.col(0)).array()).pow(3);
            computeZ.col(7) = ((lidarpts.col(0)).array()).pow(2) * (lidarpts.col(1)).array();
            computeZ.col(8) = (lidarpts.col(0)).array() * (lidarpts.col(1)).array().pow(2);
            computeZ.col(9) = ((lidarpts.col(1)).array()).pow(3);
            Eigen::MatrixXf distance = computeZ * modelcoff; //TODO: use vertical distance
            if (status) {
                for (int i = 0; i < distance.rows(); ++i) {
                    if (fabs(distance(i)) <= thresdistance) {
                        cv::Mat inlierpts = (cv::Mat_<float>(1, 4) << lidarpts(i, 0), lidarpts(i, 1), lidarpts(i,
                                                                                                               2), fabs(
                                distance(i)));
                        inliers.push_back(inlierpts);
                    }
                }
            } else {
                for (int i = 0; i < distance.rows(); ++i) {
                    cv::Mat inlierpts = (cv::Mat_<float>(1, 4) << lidarpts(i, 0), lidarpts(i, 1), lidarpts(i,
                                                                                                           2), fabs(
                            distance(i)));
                    inliers.push_back(inlierpts);

                }
            }
        }
        cv::cv2eigen(inliers, inlierroadpts);
        res.modelcoff = modelcoff;
        res.inliers = inlierroadpts;
        res.ninliers = inlierroadpts.rows();
        return res;
    }


    retresult
    runRansacRoad(Eigen::MatrixXf &lidarpts, int nSamNumForModel, int maxtrials, int ninlierThreshold,
                  float dThreshold) {
        if (lidarpts.cols() != 3)
            std::cerr << "The input matrix should be 3d..." << std::endl;

        int histcountinliers = 0;
        long long int nIter = maxtrials;
        Eigen::MatrixXf retinliers;
        retresult finalres, retdata;

        for (int count = 0; count < nIter; ++count) {
            if (count > nIter)
                break;

            // 1. sampling
            Eigen::MatrixXi sampleMask = Eigen::ArrayXXi::Zero(1, nSamNumForModel);
            // Takes nSamNumForModel different samples
            if (sampleMask.sum() != nSamNumForModel) {
                Eigen::ArrayXXf ind =
                        lidarpts.rows() / 2 * (Eigen::MatrixXf::Random(1, abs(nSamNumForModel - sampleMask.sum())) +
                                               Eigen::MatrixXf::Constant(1, abs(nSamNumForModel - sampleMask.sum()),
                                                                         1.));
                for (int nidx = 0; nidx < ind.cols(); ++nidx) {
                    sampleMask(0, nidx) = abs(floor(ind(0, nidx)));
                }
            }
            cv::Mat fitdata;
            for (int fitindex = 0; fitindex < nSamNumForModel; ++fitindex) {
                int row_pts = abs(sampleMask(0, fitindex));
                float pts[1][3] = {lidarpts(row_pts, 0), lidarpts(row_pts, 1), lidarpts(row_pts, 2)};
                cv::Mat tempts(1, 3, CV_32FC1, &pts, 2);
                fitdata.push_back(tempts);
            }
            // 2. create the model
            Eigen::MatrixXf fitlidardata;
            cv2eigen(fitdata, fitlidardata);
            Eigen::MatrixXf curModel = computeModel(fitlidardata);
            // 3. model inlier estimation
            retdata = mcdcomputeres(curModel, lidarpts, dThreshold, true);
            //4. Check the size of inliers
            if (retdata.ninliers >= ninlierThreshold) {
                finalres.modelcoff = retdata.modelcoff;
                finalres.inliers = retdata.inliers;
                finalres.ninliers = retdata.ninliers;
                return finalres;
            }
            if (retdata.ninliers > histcountinliers) {
                histcountinliers = retdata.ninliers;
                finalres.modelcoff = retdata.modelcoff;
                finalres.inliers = retdata.inliers;
                finalres.ninliers = retdata.ninliers;
                float p = 0.99;
                float e = 1 - finalres.ninliers / lidarpts.rows();
                nIter = abs(ceil(log(1 - p) / log(1 - pow(1 - e, nSamNumForModel))));
                if (fabs(nIter) > maxtrials || fabs(nIter) < -maxtrials) {
                    nIter = maxtrials;
                }
            }
            if (count >= maxtrials) {
                std::cout << "The maximum trails has been reached" << std::endl;
                break;
            }
        }
        return finalres;
    }

    Eigen::MatrixXf computePlaneModel(Eigen::ArrayXXf fitdata) {
        Eigen::ArrayXf _model;
        Eigen::ArrayXf means = fitdata.colwise().mean();
        Eigen::MatrixXf zeroCentered = (fitdata.rowwise() - means.transpose()).matrix();
        Eigen::EigenSolver<Eigen::MatrixXf> solver(zeroCentered.transpose() * zeroCentered);
        unsigned long long int minEigenValueIndex;
        solver.eigenvalues().real().array().minCoeff(
                &minEigenValueIndex); //eigen values will be real as matrix is symmetric
        _model.resize(4);
        _model.head(3) = solver.eigenvectors().col(minEigenValueIndex).real();
        _model.head(3).matrix().normalize();
        _model(3) = -_model.head(3).matrix().dot(means.matrix());
        return _model.matrix();
    }

    retresult mcdPlaneComput(Eigen::MatrixXf modelcoff, Eigen::MatrixXf lidarpts, float thresdistance, bool status) {
        // modelcoeff : parameters of the model
        // lidarpts : all the lidar 3D points

        cv::Mat inliers;
        Eigen::MatrixXf inlierroadpts;
        lidarpts.conservativeResize(lidarpts.rows(), lidarpts.cols() + 1);
        lidarpts.col(lidarpts.cols() - 1) = Eigen::VectorXf::Ones(lidarpts.rows(), 1);
        float coeffsum = sqrt(
                std::pow(modelcoff(0, 0), 2) + std::pow(modelcoff(1, 0), 2) + std::pow(modelcoff(2, 0), 2));
        if (fabs(coeffsum) != 0 && lidarpts.rows() != 0) {
            Eigen::VectorXf computeZ = (lidarpts * modelcoff).cwiseAbs() / coeffsum;
            Eigen::ArrayXf distance = computeZ.array() + 1.2 * pow(10, -3) * (((lidarpts.col(0)).array()).abs2() +
                                                                              ((lidarpts.col(1)).array()).abs2() +
                                                                              ((lidarpts.col(2)).array()).abs2()).sqrt()
                                      - 0.1 * pow(10, -3);
            if (status) {
                for (int i = 0; i < distance.rows(); ++i) {
                    if (distance(i, 0) <= thresdistance) {
                        cv::Mat inlierpts = (cv::Mat_<float>(1, 4) << lidarpts(i, 0), lidarpts(i, 1), lidarpts(i,
                                                                                                               2), distance(
                                i, 0));
                        inliers.push_back(inlierpts);
                    }
                }
            } else {
                for (int i = 0; i < distance.rows(); ++i) {
                    cv::Mat inlierpts = (cv::Mat_<float>(1, 4) << lidarpts(i, 0), lidarpts(i, 1), lidarpts(i,
                                                                                                           2), distance(
                            i, 0));
                    inliers.push_back(inlierpts);

                }
            }
        }
        cv::cv2eigen(inliers, inlierroadpts);
        return {modelcoff, inlierroadpts, inlierroadpts.rows()};
    }


    // fit the plane using ransac
    retresult
    runRansacplanefit(Eigen::MatrixXf &lidarpts, int nSamNumForModel, int maxtrials, int ninlierThreshold,
                      float dThreshold) {

        if (lidarpts.cols() != 3)
            std::cerr << "The input matrix should be 3d..." << std::endl;

        int histcountinliers = 0;
        long long int nIter = maxtrials;
        Eigen::MatrixXf retinliers;
        retresult finalres, retdata;

        for (int count = 0; count < nIter; ++count) {
            // 1. sampling
            Eigen::MatrixXi sampleMask = Eigen::ArrayXXi::Zero(1, nSamNumForModel);
            // Takes nSamNumForModel different samples
            if (sampleMask.sum() != nSamNumForModel) {
                Eigen::ArrayXXf ind =
                        lidarpts.rows() / 2 * (Eigen::MatrixXf::Random(1, abs(nSamNumForModel - sampleMask.sum())) +
                                               Eigen::MatrixXf::Constant(1, abs(nSamNumForModel - sampleMask.sum()),
                                                                         1.));
                for (int nidx = 0; nidx < ind.cols(); ++nidx) {
                    sampleMask(0, nidx) = abs(floor(ind(0, nidx)));
                }
            }
            cv::Mat fitdata;
            for (int fitindex = 0; fitindex < nSamNumForModel; ++fitindex) {
                int row_pts = abs(sampleMask(0, fitindex));
                float pts[1][3] = {lidarpts(row_pts, 0), lidarpts(row_pts, 1), lidarpts(row_pts, 2)};
                cv::Mat tempts(1, 3, CV_32FC1, &pts, 2);
                fitdata.push_back(tempts);
            }
            Eigen::MatrixXf fitlidardata;
            cv2eigen(fitdata, fitlidardata);
            Eigen::MatrixXf curModel = computePlaneModel(fitlidardata.array()); // 2. create the model
            retdata = mcdPlaneComput(curModel, lidarpts, dThreshold, true); // 3. model inlier estimation
            //4. Check the size of inliers
            if (retdata.ninliers >= ninlierThreshold) {
                finalres.modelcoff = retdata.modelcoff;
                finalres.inliers = retdata.inliers;
                finalres.ninliers = retdata.ninliers;
                return finalres;
            }

            if (retdata.ninliers >= histcountinliers) {
                histcountinliers = retdata.ninliers;
                finalres.modelcoff = retdata.modelcoff;
                finalres.inliers = retdata.inliers;
                finalres.ninliers = retdata.ninliers;
                float p = 0.99;
                float e = 1 - finalres.ninliers / lidarpts.rows();
                nIter = abs(ceil(log(1 - p) / log(1 - pow(1 - e, nSamNumForModel))));
                if (fabs(nIter) > maxtrials || fabs(nIter) < -maxtrials) {
                    nIter = maxtrials;
                }
            }

            if (count >= maxtrials || count >= nIter) {
                std::cout << "The maximum trails has been reached" << std::endl;
                break;
            }
        }
        return finalres;
    }

    void
    coordinatetransL2G(Eigen::MatrixXf roadlidarpts, Eigen::Matrix3f &RotationL2G, Eigen::RowVector3f &translationL2G) {

        if (roadlidarpts.rows() == 0)
            std::cerr << " Not enough road lidar point on the current lane!" << endl;

        Eigen::MatrixXf RoadPts;
        if (roadlidarpts.cols() == 4)
            RoadPts = roadlidarpts.block(0, 0, roadlidarpts.rows(), 3);
        else RoadPts = roadlidarpts;

        Eigen::RowVector3f meanP = RoadPts.colwise().mean();
        Eigen::MatrixXf reducedroadlidarpts = (RoadPts.array() -
                                               (Eigen::MatrixXf::Ones(RoadPts.rows(), 1) * meanP).array()).matrix();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(
                reducedroadlidarpts.transpose() * reducedroadlidarpts);
        Eigen::Matrix3f egienvectors = eigensolver.eigenvectors();
        Eigen::RowVector3f n_z = egienvectors.col(0);
        Eigen::RowVector3f n_x;
        n_x << 1, 0, -n_z(0) / n_z(2);
        n_x = n_x / n_x.norm();
        Eigen::RowVector3f n_y = n_x.cross(n_z);
        n_y = n_y / n_y.norm();

        RotationL2G << n_x,
                n_y,
                n_z;

        translationL2G << 5, 0, ((meanP(0) - 5) * n_z(0) + n_z(1) * meanP(1)) / n_z(2) + meanP(2);
    }


    void PtsRemoveBackfCar(Eigen::ArrayXXf gpsPts, float xmin, float m_maxX, Eigen::ArrayXXf &gpsPtsInL){
        vector<Eigen::Vector3f> gpsfcar;
        for(int i = 0; i < gpsPts.cols(); ++i){
            if(gpsPts(0, i) > xmin && gpsPts(0, i) < m_maxX){
                gpsfcar.emplace_back(gpsPts.col(i));
            }
        }
        gpsPtsInL = Eigen::ArrayXXf(3, gpsfcar.size());
        for(int j = 0; j < gpsfcar.size(); ++j){
            gpsPtsInL.col(j) = gpsfcar[j];
        }
    }

    AverageMeter::AverageMeter(int _numToTrack) : m_vals(_numToTrack, 0), m_index(0), m_sum(0) {}

    void AverageMeter::Update(const float &_val) {
        if (m_index >= m_vals.size())
            m_sum -= m_vals[m_index % m_vals.size()];
        m_sum += (m_vals[m_index % m_vals.size()] = _val);
        m_index++;
    }

    float AverageMeter::Average() {
        return m_sum / std::max(size_t(1), std::min(m_index, m_vals.size()));
    }

    void
    CurrentOptCtrlPts(Eigen::ArrayXXf _PrectrlPts, Eigen::ArrayXXf _CurctrlPts, Eigen::ArrayXXf &_updatedNewCtrolPts) {

        std::vector<Eigen::RowVector3f> temp;
        Eigen::ArrayXXf _invPreCtrlPts = _PrectrlPts.transpose();

        if (!_invPreCtrlPts.isZero()) {
            for (int k = 0; k < _invPreCtrlPts.rows(); ++k)
                temp.emplace_back(_invPreCtrlPts.row(k));

/*            for (int i = 0; i < _CurctrlPts.rows(); ++i) {
                if (_CurctrlPts(i, 0) > _invPreCtrlPts(_invPreCtrlPts.rows() - 1, 0)) {
                    temp.emplace_back(_CurctrlPts.row(i));
                }
            }*/

            _updatedNewCtrolPts = Eigen::ArrayXXf(temp.size(), 3);
            for (int j = 0; j < temp.size(); ++j) {
                _updatedNewCtrolPts.row(j) = temp[j];
            }
        }
        else _updatedNewCtrolPts = _CurctrlPts;

    }

    void shortestDistance(Eigen::VectorXf point, Eigen::MatrixXf _set, float &_dis) {

        // _set is nX2 or nX3 matrix
        Eigen::MatrixXf _sets = Eigen::MatrixXf(_set.rows(), _set.cols());
        _sets = _set;
        float shortest = std::numeric_limits<float>::max();

        if (point.rows() == 3) {
            if (_set.cols() < 3) {
                Eigen::VectorXf pointZ;
                pointZ = Eigen::MatrixXf::Zero(_sets.rows(), 1);
                _sets.conservativeResize(_sets.rows(), _sets.cols() + 1);
                _sets.col(_sets.cols() - 1) = pointZ;
            }

            for (int i = 0; i < _sets.rows(); ++i) {
                float dis = sqrt(
                        pow(point(0) - _sets(i, 0), 2) + pow(point(1) - _sets(i, 1), 2) +
                        pow(point(2) - _sets(i, 2), 2));
                _dis = std::min(shortest, dis);
            }
        } else {
            for (int i = 0; i < _sets.rows(); ++i) {
                float dis = sqrt(pow(point(0) - _sets(i, 0), 2) + pow(point(1) - _sets(i, 1), 2));
                _dis = std::min(shortest, dis);
            }
        }
    }


    void shortestDistance(Eigen::Vector3f point, Eigen::MatrixXf _set, float radius,
                          Eigen::MatrixXf &_nearestPoints, std::vector<float> &distance) {
        // _set is nX2 or nX3 matrix
        Eigen::MatrixXf _sets = Eigen::MatrixXf(_set.rows(), _set.cols());
        std::vector<Eigen::RowVector3f> radiusMat;
        _sets = _set;
        float shortest = std::numeric_limits<float>::max();
        if (_set.cols() < 3) {
            Eigen::VectorXf pointZ;
            pointZ = Eigen::MatrixXf::Zero(_sets.rows(), 1);
            _sets.conservativeResize(_sets.rows(), _sets.cols() + 1);
            _sets.col(_sets.cols() - 1) = pointZ;
        }

        for (int i = 0; i < _sets.rows(); ++i) {
            float dis = sqrt(
                    pow(point(0) - _sets(i, 0), 2) + pow(point(1) - _sets(i, 1), 2) + pow(point(2) - _sets(i, 2), 2));
            if (dis <= radius) {
                radiusMat.push_back(_sets.row(i));
                distance.push_back(dis);
            }
        }
    }

    int ImgFile2Int(const string &_imgFileName) {
        return std::stoi(_imgFileName.substr(0, _imgFileName.size() - 3));
    }


    void
    ShiftCenterCtrlPts(const Eigen::ArrayXXf &_centerCtrlPts, const float &_horizontalDist, const float &_upwardDist,
                       const bool &_isRightHigher, vector<Eigen::ArrayXXf> &_laneCtrlPts) {
        _laneCtrlPts = vector<Eigen::ArrayXXf>(2, _centerCtrlPts);

        _laneCtrlPts[0].row(1) -= _horizontalDist;
        _laneCtrlPts[0].row(2) += _isRightHigher ? _upwardDist : -_upwardDist;

        _laneCtrlPts[1].row(1) += _horizontalDist;
        _laneCtrlPts[1].row(2) += _isRightHigher ? -_upwardDist : _upwardDist;
    }

}