#ifndef GAUSS_PREFERENCE_FINDER_H_
#define GAUSS_PREFERENCE_FINDER_H_

#include"PreferenceFinder.h"

namespace LD {

    class GaussPreferenceFinder : public PreferenceFinder {
    public:
        virtual void operator()(const ArrayXXf &_residuals, ArrayXXf &_preferences) override {

            _preferences = ArrayXXf::Zero(_residuals.rows(), _residuals.cols());
            double sigma = m_inlierThreshold / m_xml.attribute("divFactor").as_int();
            for (ulli r = 0; r < _residuals.rows(); r++)
                for (ulli c = 0; c < _residuals.cols(); c++)
                    if (_residuals(r, c) < m_inlierThreshold)
                        _preferences(r, c) = std::exp(-std::pow(_residuals(r, c), 2) / std::pow(sigma, 2));

        }

        virtual void ParseXML() override {
            m_xml = m_xml.child("Gauss");
            m_divFactor = m_xml.attribute("divFactor").as_int();
            if (!m_divFactor)
                throw runtime_error("divFactor attribute not specified in Gauss preference finder node");
        }

        GaussPreferenceFinder(string _xmlFile) : PreferenceFinder(_xmlFile) { ParseXML(); }

    protected:
        int m_divFactor;

    };

}

#endif
