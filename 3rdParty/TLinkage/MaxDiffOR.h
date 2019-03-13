#ifndef MAX_DIFF_OR_H_
#define MAX_DIFF_OR_H_

#include "OutlierRejector.h"
#include<unordered_map>
#include<unordered_set>

namespace LD {

    class MaxDiffOR : public OutlierRejector {
    public:
        virtual void
        operator()(const ArrayXf &_clusters, const std::unordered_map<int, vector<ulli> > &_clusterID2PtIndices,
                   ArrayXf &_out, const int &_noiseIndex = -1) override;

        virtual void ParseXML() override {
            m_xml = m_xml.child("Max_size_change");
        }

        MaxDiffOR(int _minSamples, string _xmlFile) : OutlierRejector(_xmlFile),
                                                      m_minSamples(_minSamples) { ParseXML(); }

    protected:
        int m_minSamples;
    };
}
#endif
