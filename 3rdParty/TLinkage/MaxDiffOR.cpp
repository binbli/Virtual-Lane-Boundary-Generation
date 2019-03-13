#include"MaxDiffOR.h"
#include<limits>

namespace LD {
    void
    MaxDiffOR::operator()(const ArrayXf &_clusters, const std::unordered_map<int, vector<ulli> > &_clusterID2PtIndices,
                          ArrayXf &_out, const int &_noiseIndex) {
        vector<std::pair<int, ulli> > count;

        for (auto it = _clusterID2PtIndices.begin(); it != _clusterID2PtIndices.end(); it++)
            count.push_back(std::make_pair(it->first, it->second.size()));

        count.push_back(std::make_pair(std::numeric_limits<int>::min(), m_minSamples)); //add dummy

        sort(count.begin(), count.end(), [](std::pair<int, ulli> &a, std::pair<int, ulli> &b) {
            return a.second < b.second;
        });

        ulli maxChange = 0, maxChangeIndex = 0;

        for (ulli i = 0;
             i < count.size() - 2; i++) {  //intentionally skipping last value; FIXME: may be parallelized;
            ulli diff = count[i + 1].second - count[i].second;
            if (diff > maxChange) {
                maxChange = diff;
                maxChangeIndex = i + 1;
            }
        }

        std::unordered_set<ulli> clustersToDiscard;
        for (ulli i = 0; i < maxChangeIndex; i++)
            clustersToDiscard.insert(count[i].first);

        _out = _clusters;
        for (ulli i = 0; i < _clusters.size(); i++)
            if (clustersToDiscard.find(_out(i)) != clustersToDiscard.end())
                _out(i) = _noiseIndex;
    }
}
