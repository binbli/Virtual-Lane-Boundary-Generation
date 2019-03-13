#ifndef SAMPLER_H_
#define SAMPLER_H_

#include<eigen3/Eigen/Dense>
#include "BaseLD.h"

namespace LD {
    using namespace Eigen;

    class Sampler : public BaseLD {
    public:
        virtual void operator()(const ArrayXXf &_data, ArrayXXf &_sampleIndices,
                                const ulli &_numSamples, const ulli &_minSamples) = 0;

        virtual void ParseXML() override {
            m_xml = m_xml.child("Solvers").child("TLinkage").child("Samplers");
            m_maxIterationsFactor = m_xml.attribute("maxIterationsFactor").as_int();
            m_maxTries = m_xml.attribute("maxTries").as_int();

            if (!m_maxIterationsFactor || !m_maxTries)
                throw runtime_error(
                        "At least one attribute is invalid/missing in Samplers node: maxIterationsFactor, maxTries");

        }

        Sampler(string _xmlFile) : BaseLD(_xmlFile) { ParseXML(); }

        class MinSamplesNotFound {
        };

        int m_maxIterationsFactor, m_maxTries;
    };
}
#endif