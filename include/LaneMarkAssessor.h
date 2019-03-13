#ifndef EVALUTIONPROJECT_LANEMARKASSESSOR_H
#define EVALUTIONPROJECT_LANEMARKASSESSOR_H

#include"BaseLD.h"
#include"../3rdParty/TLinkage/BSplineTLinkage.h"

namespace LD {
    class LaneMarkAssessor : public BaseLD {
    protected:
        void ParseXML();
        double m_precision;
        BSplineTLinkage m_bSplineTLinkage;

    public:

        LaneMarkAssessor(string _xmlFile);

        void operator()(const Eigen::ArrayXf& _splineModel, double& _headingAngle, double& _length, double& _maxCurvature, Eigen::ArrayXXf &splineCurve);

        double GetHeadingAngle(const Eigen::ArrayXf& _splineModel);
    };
};

#endif //EVALUTIONPROJECT_LANEMARKASSESSOR_H
