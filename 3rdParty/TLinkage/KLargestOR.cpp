#include"KLargestOR.h"
#include<unordered_set>

namespace LD {

    void
    KLargestOR::operator()(const ArrayXf &_clusters, const std::unordered_map<int, vector<ulli> > &_clusterID2PtIndices,
                           ArrayXf &_out, const int &_noiseIndex) {
        if (m_debug)
            cout << "Entering KLargestOR::operator()" << endl;

        vector<std::pair<ulli, ulli> > count;

        for (auto it = _clusterID2PtIndices.begin(); it != _clusterID2PtIndices.end(); it++)
            if (it->first != _noiseIndex)
                count.push_back(std::make_pair(it->first, it->second.size()));

        int k = std::min(m_k, int(count.size()));

        std::partial_sort(count.begin(), count.begin() + k, count.end(),
                          [](std::pair<ulli, ulli> &a, std::pair<ulli, ulli> &b) {
                              return a.second > b.second;
                          });

        std::unordered_set<ulli> clustersToKeep;

        for (ulli i = 0; i < k; i++)
            clustersToKeep.insert(count[i].first);

        _out = _clusters;
        for (ulli i = 0; i < _clusters.size(); i++)
            if (clustersToKeep.find(_clusters(i)) == clustersToKeep.end())
                _out(i) = _noiseIndex;

        if (m_debug)
            cout << "Exiting KLargestOR::operator()" << endl;
    }

    void KLargestOR::ParseXML() {
        m_xml = m_xml.child("KLargest");
        m_k = m_xml.attribute("k").as_int();
        if (!m_k)
            throw runtime_error("at least one of the following attributes is missing in KLargest node: k");
    }
}
