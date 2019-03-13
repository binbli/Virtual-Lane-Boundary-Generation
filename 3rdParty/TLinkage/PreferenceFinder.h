#ifndef PREFERENCE_FINDER_H_
#define PREFERENCE_FINDER_H_

#include<eigen3/Eigen/Dense>
#include "BaseLD.h"

namespace LD {

    using namespace Eigen;

    class PreferenceFinder : public BaseLD {
    public:
        virtual void operator()(const ArrayXXf &_residuals, ArrayXXf &_preferences) = 0;

        virtual void ParseXML() override {

            m_xml = m_xml.child("Solvers").child("TLinkage").child("PreferenceFinders");
            m_inlierThreshold = m_xml.attribute("inlierThreshold").as_float();

            if (!m_inlierThreshold)
                throw runtime_error("inlierThreshold attribute not specified in PreferenceFinder Node");
        }

        PreferenceFinder(string _xmlFile) : BaseLD(_xmlFile) { ParseXML(); }

    protected:
        float m_inlierThreshold;
    };

}

#endif
