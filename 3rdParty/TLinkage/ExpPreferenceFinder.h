#ifndef EXP_PREFERENCE_FINDER_H_
#define EXP_PREFERENCE_FINDER_H_

#include"PreferenceFinder.h"

namespace LD {

    class ExpPreferenceFinder : public PreferenceFinder {
    public:
        virtual void operator()(const ArrayXXf &_residuals, ArrayXXf &_preferences) override {

            if (m_debug)
                cout << "Entering ExpPreferenceFinder::operator()" << endl;

            _preferences = ArrayXXf::Zero(_residuals.rows(), _residuals.cols());
            double tau = m_inlierThreshold / m_divFactor;
            for (ulli r = 0; r < _residuals.rows(); r++)
                for (ulli c = 0; c < _residuals.cols(); c++)
                    if (_residuals(r, c) < m_inlierThreshold)
                        _preferences(r, c) = std::exp(-_residuals(r, c) / tau);

            if (m_debug)
                cout << "Exiting ExpPreferenceFinder::operator()" << endl;

        }

        virtual void ParseXML() override {
            m_xml = m_xml.child("Exp");
            m_divFactor = m_xml.attribute("divFactor").as_int();
            if (!m_divFactor)
                throw runtime_error("divFactor attribute not specified in Exp preference finder node");
        }

        ExpPreferenceFinder(string _xmlFile) : PreferenceFinder(_xmlFile) { ParseXML(); }

    protected:
        int m_divFactor;
    };

}
#endif
