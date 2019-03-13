#ifndef TUKEY_PREFERENCE_FINDER_H_
#define TUKEY_PREFERENCE_FINDER_H_

#include"PreferenceFinder.h"

namespace LD {

    class TukeyPreferenceFinder : public PreferenceFinder {
    public:
        virtual void operator()(const ArrayXXf &_residuals, ArrayXXf &_preferences) override {

            _preferences = ArrayXXf::Zero(_residuals.rows(), _residuals.cols());

            for (ulli r = 0; r < _residuals.rows(); r++)
                for (ulli c = 0; c < _residuals.cols(); c++)
                    if (_residuals(r, c) < m_inlierThreshold)
                        _preferences(r, c) = std::pow(1 - std::pow(_residuals(r, c) / m_inlierThreshold, 2), 3);

        }

        virtual void ParseXML() override {
            m_xml = m_xml.child("Tukey");
        }

        TukeyPreferenceFinder(string _xmlFile) : PreferenceFinder(_xmlFile) { ParseXML(); }
    };

}

#endif
