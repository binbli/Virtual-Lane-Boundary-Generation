#include"BSplineTLinkage.h"
#include"Utilities.h"
#include<stdlib.h>
#include<limits>

namespace LD {

    BSplineTLinkage::BSplineTLinkage(string _xmlFile) : TLinkage(0, 0, _xmlFile), m_params1dSpline(6) {
        ParseXML();
        m_modelParams = 2 * m_params1dSpline * (m_minSamples - 1) + 1;
    }

    void BSplineTLinkage::ParseXML() {
        if (m_debug)
            cout << "Entering BSplineTLinkage::ParseXML()" << endl;

        m_xml = m_xml.child("Models").child("BSpline");
        m_minSamples = m_xml.attribute("minSamples").as_int();
        m_regularizationConst = m_xml.attribute("regularizationConst").as_float(-15);
        m_resolution = m_xml.attribute("resolution").as_float();
        m_minX = m_xml.attribute("minX").as_float();
        m_maxX = m_xml.attribute("maxX").as_float();

        if (m_minSamples <= 5 || m_resolution <= 0 || m_minX <= 0 || m_maxX <= 0)
            throw runtime_error(
                    "at least one attribute is missing/invalid: minSamples, regularizationConst, resolution, minX, maxX");

        if (m_debug)
            cout << "Exiting BSplineTLinkage::ParseXML()" << endl;
    }

    ArrayXf BSplineTLinkage::GenerateHypothesis(const vector<ArrayXf> &_samples) {

        alglib::spline1dinterpolant yOfX, zOfX;
        alglib::real_2d_array yCoeffs, zCoeffs;
        vector<alglib::real_1d_array> coordinates;
        CreateAlglibArray(_samples, coordinates);

        alglib::ae_int_t ctrlPts = coordinates[0].length();

        spline1dbuildcubic(coordinates[0], coordinates[1], ctrlPts, 0, 0.0, 0, 0.0, yOfX);
        spline1dbuildcubic(coordinates[0], coordinates[2], ctrlPts, 0, 0.0, 0, 0.0, zOfX);

        spline1dunpack(yOfX, ctrlPts, yCoeffs);
        spline1dunpack(zOfX, ctrlPts, zCoeffs);

        ArrayXf model = ConvertCoeffsTable2Model(yCoeffs, zCoeffs);

        return model;

    }

    double BSplineTLinkage::Distance(ArrayXf _dataPoint, ArrayXf _model) {
        alglib::spline1dinterpolant yOfX, zOfX;

        //build from _model
        CreateSplineInterpolants(_model, yOfX, zOfX);

        //brute force to find minimum distance (shouldn't be very expensive for realistic m_minX and m_maxX)
        float minDist = std::numeric_limits<float>::max();

        for (float x = m_minX; x <= m_maxX; x += m_resolution) {
            Array3f ptOnSpline;
            ptOnSpline << x, spline1dcalc(yOfX, x), spline1dcalc(zOfX, x);
            minDist = std::min(minDist, (ptOnSpline - _dataPoint).matrix().norm());
        }

        return minDist;

    }

    float BSplineTLinkage::MeanModelHeight(const ArrayXf &_model, const vector<ulli> &_clusterIndices, const ArrayXXf &_data) {

        alglib::spline1dinterpolant yOfX, zOfX;

        //build from _model
        CreateSplineInterpolants(_model, yOfX, zOfX);

        float height = 0;
        for (int i = 0; i < _clusterIndices.size(); i++)
            height += spline1dcalc(zOfX, _data(0, _clusterIndices[i]));

        return height / _clusterIndices.size();

    }

    void BSplineTLinkage::FitModel(const ArrayXXf &_clusters, ArrayXf &_model) {
        if (m_debug)
            cout << "Entering BSplineTLinkage::FitModel()" << endl;

        alglib::spline1dinterpolant yOfX, zOfX;
        alglib::real_2d_array yCoeffs, zCoeffs;
        alglib::spline1dfitreport rep;
        alglib::ae_int_t info;
        vector<alglib::real_1d_array> coordinates;
        CreateAlglibArray(_clusters, coordinates);  //FIXME: Deep copy can be avoided, by doing a shallow copy
        spline1dfitpenalized(coordinates[0], coordinates[1], m_minSamples - 2, m_regularizationConst, info, yOfX, rep);

        if (info < 0)
            throw runtime_error("Can't fit least square spline");

        spline1dfitpenalized(coordinates[0], coordinates[2], m_minSamples - 2, m_regularizationConst, info, zOfX, rep);

        if (info < 0)
            throw runtime_error("Can't fit least square spline");

        alglib::ae_int_t ctrlPts = m_minSamples;

        spline1dunpack(yOfX, ctrlPts, yCoeffs);
        spline1dunpack(zOfX, ctrlPts, zCoeffs);
        _model = ConvertCoeffsTable2Model(yCoeffs, zCoeffs);

        if (m_debug)
            cout << "Exiting BSplineTLinkage::FitModel()" << endl;
    }

    void BSplineTLinkage::Fit2DModel(const ArrayXXf &_clusters, ArrayXf &_model, alglib::spline1dinterpolant &_yOfX) {
        if (m_debug)
            cout << "Entering BSplineTLinkage::Fit2DModel()" << endl;

        alglib::spline1dinterpolant yOfX, zOfX;
        alglib::real_2d_array yCoeffs, zCoeffs;
        alglib::spline1dfitreport rep;
        alglib::ae_int_t info;
        vector<alglib::real_1d_array> coordinates;
        alglib::real_1d_array x, y;

        x.setlength(m_minSamples);
        y.setlength(m_minSamples);
        _model = ArrayXf(m_modelParams);

        CreateAlglibArray(_clusters, coordinates);
        spline1dfitpenalized(coordinates[0], coordinates[1], m_minSamples - 2, m_regularizationConst, info, yOfX, rep);

        if (info < 0)
            throw runtime_error("Can't fit least square spline");

        alglib::ae_int_t ctrlPts = m_minSamples;
        spline1dunpack(yOfX, ctrlPts, yCoeffs);

        _model(0) = 1 + yCoeffs.rows() * yCoeffs.cols(); // first index where zCoeffs start

        for (int r = 0; r < yCoeffs.rows(); r++)
            for (int c = 0; c < yCoeffs.cols(); c++)
                _model(1 + c + r * yCoeffs.cols()) = yCoeffs(r, c);

        for (int i = 0; i < m_minSamples - 1; i++) {
            x(i) = _model(1 + m_params1dSpline * i);
            y(i) = _model(3 + m_params1dSpline * i);
        }

        int i = m_minSamples - 2;
        int ySegmentStart = 1 + m_params1dSpline * i;
        x(i + 1) = _model(ySegmentStart + 1);
        float t = x(i) - _model(ySegmentStart);
        y(i + 1) = _model(2 + ySegmentStart) +
                   t * (_model(3 + ySegmentStart) + t * (_model(4 + ySegmentStart) + t * _model(5 + ySegmentStart)));

        spline1dbuildcubic(x, y, ctrlPts, 0, 0.0, 0, 0.0, _yOfX);

        if (m_debug)
            cout << "Exiting BSplineTLinkage::Fit2DModel()" << endl;

    }


    void BSplineTLinkage::CubicSplineModelToCurve(Eigen::ArrayXf model, Eigen::MatrixXf &curve, float steplen){
        if (m_debug)
            cout << "Entering BSplineTLinkage::CubicSplineModelToCurve()" << endl;

        std::vector<Eigen::RowVector3f> Curves;
        alglib::spline1dinterpolant yOfX, zOfX;
        ArrayXXf _ctrlPts;

        CreateSplineInterpolants(model, yOfX, zOfX);
        GetControlPts(model, _ctrlPts);
        float minx = _ctrlPts.col(0).minCoeff();
        float maxx = _ctrlPts.col(0).maxCoeff();
        for (float x = minx; x <= maxx; x += steplen) {
            Eigen::RowVector3f ptOnSpline;
            ptOnSpline << x, spline1dcalc(yOfX, x), spline1dcalc(zOfX, x);
            Curves.push_back(ptOnSpline);
        }

        curve = Eigen::MatrixXf(Curves.size(), 3);
        for (int i = 0; i < Curves.size(); i++)
            curve.row(i) = Curves[i];

        if (m_debug)
            cout << "Exiting BSplineTLinkage::CubicSplineModelToCurve()" << endl;
    }



    void BSplineTLinkage::CreateSplineInterpolants(ArrayXf _model, alglib::spline1dinterpolant &_yOfX,
                                                   alglib::spline1dinterpolant &_zOfX) {

        alglib::real_1d_array x, y, z;
        x.setlength(m_minSamples);
        y.setlength(m_minSamples);
        z.setlength(m_minSamples);

        for (int i = 0; i < m_minSamples - 1; i++) {
            x(i) = _model(1 + m_params1dSpline * i);
            y(i) = _model(3 + m_params1dSpline * i);
            z(i) = _model(_model(0) + 2 + m_params1dSpline * i);
        }

        int i = m_minSamples - 2;
        int ySegmentStart = 1 + m_params1dSpline * i;
        int zSegmentStart = _model(0) + m_params1dSpline * i;
        x(i + 1) = _model(ySegmentStart + 1);
        float t = x(i) - _model(ySegmentStart);
        y(i + 1) = _model(2 + ySegmentStart) +
                   t * (_model(3 + ySegmentStart) + t * (_model(4 + ySegmentStart) + t * _model(5 + ySegmentStart)));
        z(i + 1) = _model(2 + zSegmentStart) +
                   t * (_model(3 + zSegmentStart) + t * (_model(4 + zSegmentStart) + t * _model(5 + zSegmentStart)));

        alglib::ae_int_t ctrlPts = m_minSamples;

        spline1dbuildcubic(x, y, ctrlPts, 0, 0.0, 0, 0.0, _yOfX);
        spline1dbuildcubic(x, z, ctrlPts, 0, 0.0, 0, 0.0, _zOfX);

    }

    ArrayXf
    BSplineTLinkage::ConvertCoeffsTable2Model(alglib::real_2d_array &_yCoeffs, alglib::real_2d_array &_zCoeffs) {

        ArrayXf hypothesis(m_modelParams);
        hypothesis(0) = 1 + _yCoeffs.rows() * _yCoeffs.cols(); // first index where zCoeffs start

        for (int r = 0; r < _yCoeffs.rows(); r++)
            for (int c = 0; c < _yCoeffs.cols(); c++)
                hypothesis(1 + c + r * _yCoeffs.cols()) = _yCoeffs(r, c);

        for (int r = 0; r < _zCoeffs.rows(); r++)
            for (int c = 0; c < _zCoeffs.cols(); c++)
                hypothesis(hypothesis(0) + c + r * _zCoeffs.cols()) = _zCoeffs(r, c);

        return hypothesis;

    }

    void BSplineTLinkage::PrintModelsToFile(vector<ArrayXf> _models, const string &_imgName) {
        if (m_debug)
            cout << "Entering BSplineTLinkage::PrintModelsToFile()" << endl;

        string ModalFilePath = m_LBmodelFile + _imgName.substr(0, _imgName.size() - 4) + ".txt";

        if(!m_foutModels.is_open())
            m_foutModels.open(ModalFilePath);

        m_foutModels << _imgName << endl;
        m_foutModels << _models.size() << "\t" << m_modelParams << endl;
        for (int i = 0; i < _models.size(); i++) {
            m_foutModels << _models[i](0) << endl;
            for (int r = 0; r < (_models[i].size() - 1) / m_params1dSpline; r++) {
                for (int c = 0; c < m_params1dSpline; c++) {
                    m_foutModels << _models[i](1 + c + r * m_params1dSpline) << "\t";
                }
                m_foutModels << endl;
            }
            m_foutModels << endl;
        }
        m_foutModels.close();

        if (m_debug)
            cout << "Exiting BSplineTLinkage::PrintModelsToFile()" << endl;
    }

    void BSplineTLinkage::VisualizeModel(ArrayXf &_model, ArrayXXf &_coordinates) {
        if (m_debug)
            cout << "Entering BSplineTLinkage::VizualizeModels()" << endl;

        int zCoeffsStart = _model(0);
        int numXTerms = (_model(zCoeffsStart - 5) - _model(1)) / m_resolution;
        _coordinates.resize(numXTerms, 3);
        _coordinates = 0;
        _coordinates.col(0).setLinSpaced(numXTerms, _model(1), _model(zCoeffsStart - 5));
        int j = 0;
        for (int i = 0; i < _coordinates.rows(); i++) {
            float x = _coordinates(i, 0);
            int yStart = m_params1dSpline * j + 1;
            int zStart = m_params1dSpline * j + zCoeffsStart;
            if((yStart + 6) <= _model.rows() && (zStart + 6) <= _model.rows()) {
                j += (x >= _model(yStart) && x <= _model(yStart + 1)) ? 0 : 1;
                float t = x - _model(yStart);
                _coordinates(i, 1) +=
                        _model(yStart + 2) +
                        t * (_model(yStart + 3) + t * (_model(yStart + 4) + t * _model(yStart + 5)));
                _coordinates(i, 2) +=
                        _model(zStart + 2) +
                        t * (_model(zStart + 3) + t * (_model(zStart + 4) + t * _model(zStart + 5)));
            }
        }


        if (m_debug)
            cout << "Exiting BSplineTLinkage::VizualizeModels()" << endl;
    }

    bool BSplineTLinkage::IsModelOnRight(const ArrayXf &_model) {
        return _model[3] < 0; //_models[3] = y value at minimum x
    }

    void BSplineTLinkage::ShiftModelHorizontallyBy(ArrayXf &_model, const float &_shiftBy) {
        for (int i = 2; i < _model(0); i += m_params1dSpline) {
            _model(1 + i) += _shiftBy;
        }
    }

    void BSplineTLinkage::ShiftModelUpwardBy(ArrayXf &_model, const float &_shiftBy) {
        for (int i = _model(0); i < _model.size(); i += m_params1dSpline) {
            _model(i + 2) += _shiftBy;
        }
    }

    float BSplineTLinkage::ModelHeightDiff(const ArrayXf& _model1, const ArrayXf& _model2) {
        float diff = 0;
        int num = 0;
        for (int i = _model1(0); i < _model1.size(); i += m_params1dSpline) {
            diff += _model2(i + 2) - _model1(i + 2);
            num++;
        }
        return diff / num;
    }

    float BSplineTLinkage::ModelWidthDiff(const ArrayXf& _model1, const ArrayXf& _model2) {
        float diff = 0;
        int num = 0;
        for (int i = 2; i < _model1(0); i += m_params1dSpline) {
            diff += _model2(i + 1) - _model1(i + 1);
            num++;
        }
        return diff / num;
    }
    void BSplineTLinkage::GetControlPts(const ArrayXf& _model, ArrayXXf& _ctrlPts) {
        _ctrlPts = ArrayXXf(3, m_minSamples);
        //get ys
        int ctrlIndex = 0;
        for(int i = 1; i < _model(0); i += m_params1dSpline) {
            _ctrlPts(0, ctrlIndex) = _model(i);                 // = x
            _ctrlPts(1, ctrlIndex) = _model(i + 2);             // = y
            _ctrlPts(2, ctrlIndex) = _model(_model(0) + i + 1); // = z
            ctrlIndex++;
        }
        int i = _model(0) - m_params1dSpline;
        _ctrlPts(0, ctrlIndex) = _model(i + 1);
        float t = _model(i + 1) - _model(i);
        _ctrlPts(1, ctrlIndex) = _model(i + 2) + t * (_model(i + 3) + t * (_model(i + 4) + t * (_model(i + 5))));
        _ctrlPts(2, ctrlIndex) = _model(_model(0) + i + 1) + t * (_model(_model(0) + i + 2) + t * (_model(_model(0) + i + 3) + t * (_model(_model(0) + i + 4))));

    }

}