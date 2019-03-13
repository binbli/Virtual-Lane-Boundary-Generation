#include"Line3DTLinkage.h"
#include<Eigen/Geometry>

namespace LD {

    void Line3DTLinkage::ParseXML() {
        if (m_debug)
            cout << "Entering Line3DTLinkage::ParseXML()" << endl;

        m_xml = m_xml.child("Models").child("Line3D");
        m_resolution = m_xml.attribute("resolution").as_float(0);
        m_maxShift = m_xml.attribute("maxShift").as_float(0);
        m_minShift = m_xml.attribute("minShift").as_float(0);
        m_errorThreshold = m_xml.attribute("errorThreshold").as_float(0);

        if (m_resolution <= 0 || m_minShift <= 0 || m_maxShift <= 0 || m_errorThreshold <= 0)
            throw runtime_error(
                    "at least one of the following attributes in Line3D node is missing/invalid: resolution, minShift, maxShift, m_errorThreshold");

        if (m_debug)
            cout << "Entering Line3DTLinkage::ParseXML()" << endl;
    };

    ArrayXf Line3DTLinkage::GenerateHypothesis(const vector<ArrayXf> &_samples) {

        ArrayXf hypothesis(m_modelParams);

        hypothesis.head<3>() = _samples[0];
        hypothesis.tail<3>() = (_samples[1] - _samples[0]).matrix().normalized();

        return hypothesis;
    }

    double Line3DTLinkage::Distance(ArrayXf _dataPoint, ArrayXf _model) {
        Array3f ptOnLine = _model.head<3>() + _model.tail<3>();
        Array3f diff1 = _dataPoint - ptOnLine;
        Array3f diff2 = _dataPoint - _model.head<3>();
        return diff1.matrix().cross(diff2.matrix()).norm() / _model.tail<3>().matrix().norm();
    }

    void Line3DTLinkage::FitModel(const ArrayXXf &_cluster, ArrayXf &_model) {
        if (m_debug)
            cout << "Entering Line3DTLinkage::FitModel()" << endl;

        //Find Least Squares Lines
        ArrayXf means = _cluster.colwise().mean();
        MatrixXf zeroCentered = (_cluster.rowwise() - means.transpose()).matrix();
        EigenSolver<MatrixXf> solver(zeroCentered.transpose() * zeroCentered);
        ulli maxEigenValueIndex;
        solver.eigenvalues().real().array().maxCoeff(
                &maxEigenValueIndex); //eigen values will be real as matrix is symmetric

        _model.resize(m_modelParams);
        _model.head<3>() = means;
        _model.tail<3>() = solver.eigenvectors().col(maxEigenValueIndex).real().normalized();

        if (m_debug)
            cout << "Exiting Line3DTLinkage::FitModel()" << endl;
    }

    virtual void ShiftModelBy(ArrayXf &_model, const float &_shiftBy) {
        _model(1) += _shiftBy;
    }

    bool Line3DTLinkage::IsModelOnRight(const ArrayXf &_model) {
        if (m_debug)
            cout << "Exiting Line3DTLinkage::IsModelOnRight()" << endl;

        float t = -_model(0) /
                  _model(3); //for line <x, y, z> + t * <a, b, c>, x' = x + t * b.  So t = (x' - x) / b.  (here, x' = 0)

        float yIntercept = _model(1) + t * _model(4);

        return yIntercept > 0;

        if (m_debug)
            cout << "Exiting Line3DTLinkage::IsModelOnRight()" << endl;

    }

}
