#ifndef UNIFORM_SAMPLER_H_
#define UNIFORM_SAMPLER_H_

#include "Sampler.h"
#include<unordered_set>

namespace LD {

    class UniformSampler : public Sampler {
    public:
        virtual void operator()(const ArrayXXf &_data, ArrayXXf &_sampleIndices,
                                const ulli &_numSamples, const ulli &_minSamples) override;

        virtual void ParseXML() override {
            m_xml = m_xml.child("Uniform");
        }

        UniformSampler(string _xmlFile) : Sampler(_xmlFile) { ParseXML(); }
    };

}

#endif
