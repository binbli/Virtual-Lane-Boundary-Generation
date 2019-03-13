#ifndef B_SPLINE_T_LINKAGE_
#define B_SPLINE_T_LINKAGE_

#include"TLinkage.h"
#include<unordered_map>
#include"interpolation.h"

namespace LD {

    class BSplineTLinkage : public TLinkage {
    public:
        virtual ArrayXf GenerateHypothesis(const vector<ArrayXf> &_samples);

        virtual double Distance(ArrayXf _dataPoint, ArrayXf _model) override;

        virtual void FitModel(const ArrayXXf &_clusters, ArrayXf &_model) override;

        void Fit2DModel(const ArrayXXf &_clusters, ArrayXf &_model, alglib::spline1dinterpolant &_yOfX);

        BSplineTLinkage(string _xmlFile);

        ArrayXf ConvertCoeffsTable2Model(alglib::real_2d_array &_yCoeffs, alglib::real_2d_array &_zCoeffs);

        void CreateSplineInterpolants(ArrayXf _model, alglib::spline1dinterpolant &_yOfX, alglib::spline1dinterpolant &_zOfX);

        virtual void PrintModelsToFile(vector<ArrayXf> _models, const string &_imgName) override;

        virtual void VisualizeModel(ArrayXf &_model, ArrayXXf &_coordinates) override;

        virtual bool IsModelOnRight(const ArrayXf &_model) override;

        virtual void ShiftModelHorizontallyBy(ArrayXf &_model, const float &_shiftBy) override;

        virtual void ShiftModelUpwardBy(ArrayXf &_model, const float &_shiftBy) override;

        virtual float ModelHeightDiff(const ArrayXf& _model1, const ArrayXf& _model2) override;

        virtual float ModelWidthDiff(const ArrayXf& _model1, const ArrayXf& _model2) override;

        virtual float MeanModelHeight(const ArrayXf &_model, const vector<ulli> &_clusterIndices, const ArrayXXf &_data);

        virtual void GetControlPts(const ArrayXf& _model, ArrayXXf& _cntrlPts);

        void CubicSplineModelToCurve(Eigen::ArrayXf model, Eigen::MatrixXf &curve, float steplen);

    protected:
        void ParseXML() override;

        float m_regularizationConst;
        float m_resolution;
        float m_minX, m_maxX;
        int m_params1dSpline;

    };

}
#endif