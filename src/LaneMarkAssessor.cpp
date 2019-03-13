#include"../include/LaneMarkAssessor.h"
#include"interpolation.h"

namespace LD {

    LaneMarkAssessor::LaneMarkAssessor(string _xmlFile) : BaseLD(_xmlFile), m_bSplineTLinkage(_xmlFile) {
        ParseXML();
    }

    void LaneMarkAssessor::ParseXML() {
        if (m_debug)
            cout << "Entering LaneMarkAssessor::ParseXML() " << endl;

        m_xml = m_xml.child("LaneMarkAssessor");
        m_precision = m_xml.attribute("precision").as_double();

        if(!m_precision)
            throw runtime_error("At least one of the following attributes is undefined in LaneMarkAssessor node: precision");

        if (m_debug)
            cout << "Existing LaneMarkAssessor::ParseXML() " << endl;

    }

    void LaneMarkAssessor::operator()(const Eigen::ArrayXf& _splineModel, double& _headingAngle, double& _length, double& _maxCurvature, Eigen::ArrayXXf &splineCurve) {

        if (m_debug)
            cout << "Entering LaneMarkAssessor::operator() " << endl;

        _headingAngle = GetHeadingAngle(_splineModel);
        _length = 0;
        _maxCurvature = 0;
        double y, z;
        double dy, dz;
        double d2y, d2z;
        double curvature;
        std::vector<Eigen::RowVector3f> Curves;

        alglib::spline1dinterpolant yOfX, zOfX;
        m_bSplineTLinkage.CreateSplineInterpolants(_splineModel, yOfX, zOfX);

        for(double x = _splineModel(1); x < _splineModel(_splineModel(0) - 5); x += m_precision) {
            alglib::spline1ddiff(yOfX, x, y, dy, d2y);
            alglib::spline1ddiff(zOfX, x, z, dz, d2z);
            Eigen::RowVector3f ptOnSpline;
            ptOnSpline << x, spline1dcalc(yOfX, x), spline1dcalc(zOfX, x);
            Curves.push_back(ptOnSpline);
            _length += sqrt(1 + pow(dy, 2) + pow(dz, 2)) * m_precision;
            Eigen::Vector3f firstDeriv, secondDeriv;
            firstDeriv << 1, dy, dz;
            secondDeriv << 0, d2y, d2z;
            if(x < 5){ // we only consider near distance curvature
                curvature = firstDeriv.cross(secondDeriv).norm() / pow(firstDeriv.norm(), 3);
                _maxCurvature = std::max(curvature, _maxCurvature);
            }
        }

        splineCurve = Eigen::ArrayXXf(Curves.size(), 3);
        for (int i = 0; i < Curves.size(); i++)
            splineCurve.row(i) = Curves[i];

        if (m_debug)
            cout << "Existing LaneMarkAssessor::operator() " << endl;

    }

    double LaneMarkAssessor::GetHeadingAngle(const Eigen::ArrayXf& _splineModel) {
        return atan(_splineModel(3) / _splineModel(1)) * 180 / PI;
    }

}
