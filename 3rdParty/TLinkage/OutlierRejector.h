#ifndef OUTLIER_REJECTOR_H_
#define OUTLIER_REJECTOR_H_

#include<eigen3/Eigen/Dense>
#include "BaseLD.h"
#include<unordered_map>

namespace LD {
    using namespace Eigen;

    class OutlierRejector : public BaseLD {
    public:
        virtual void
        operator()(const ArrayXf &_clusters, const std::unordered_map<int, vector<ulli> > &_clusterID2PtIndices,
                   ArrayXf &_out, const int &_noiseIndex = -1) = 0;

        virtual void ParseXML() override {
            m_xml = m_xml.child("Solvers").child("TLinkage").child("OutlierRejectors");
        }

        OutlierRejector(string _xmlFile) : BaseLD(_xmlFile) { ParseXML(); }
    };
}
#endif
