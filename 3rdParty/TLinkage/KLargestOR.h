#ifndef K_LARGEST_OR_H_
#define K_LARGEST_OR_H_

#include"OutlierRejector.h"

//TODO: Make the function that calculates clusters also calculate their respective sizes

namespace LD {
    using namespace Eigen;

    class KLargestOR : public OutlierRejector {
    public:
        virtual void
        operator()(const ArrayXf &_clusters, const std::unordered_map<int, vector<ulli> > &_clusterID2PtIndices,
                   ArrayXf &_out, const int &_noiseIndex = -1);

        virtual void ParseXML() override;

        KLargestOR(string _xmlFile) : OutlierRejector(_xmlFile) { ParseXML(); }

    protected:
        int m_k;
    };
}
#endif
