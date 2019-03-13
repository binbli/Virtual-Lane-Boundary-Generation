#ifndef DIST_BASED_SAMPLER_H_
#define DIST_BASED_SAMPLER_H_

#include "Sampler.h"
#include<unordered_set>
#include<algorithm>
#include <boost/algorithm/string/predicate.hpp>

namespace LD {

    class DistBasedSampler : public Sampler {
    public:
        virtual void operator()(const ArrayXXf &_data, ArrayXXf &_sampleIndices,
                                const ulli &_numSamples, const ulli &_minSamples) override;

        virtual void ParseXML() override;

        DistBasedSampler(string _xmlFile) : Sampler(_xmlFile) { ParseXML(); }

        enum MeasurementWays {
            EUCLIDEAN,
            ABS_VERTICAL_DEGREE
        };

    protected:
        virtual void FindCDF(const ArrayXXf &_data, ArrayXXf &_cdfs);

        double Distance(const ArrayXf &pt1, const ArrayXf &pt2);

        double m_sigma, m_maxDiff;
        bool m_shouldXUnique;
        MeasurementWays m_measureBy;

    };
}
#endif
