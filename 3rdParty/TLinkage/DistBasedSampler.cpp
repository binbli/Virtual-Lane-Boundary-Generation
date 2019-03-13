#include"DistBasedSampler.h"
#include<iomanip>

namespace LD {

    void DistBasedSampler::ParseXML() {
        m_xml = m_xml.child("DistBased");
        m_sigma = m_xml.attribute("sigma").as_float();
        m_maxDiff = m_xml.attribute("maxDiff").as_float();
        m_shouldXUnique = m_xml.attribute("shouldXUnique").as_bool(true);

        string measurementWay = m_xml.attribute("measurementWay").as_string();

        if (!m_sigma || !m_maxDiff || measurementWay.empty())
            throw runtime_error(
                    "at least one attibute missing in DistBased Sampler node:sigma, maxDiff, measurementWay");
        if (boost::iequals(measurementWay, "Euclidean"))
            m_measureBy = EUCLIDEAN;
        else if (boost::iequals(measurementWay, "AbsVerticalDegree"))
            m_measureBy = ABS_VERTICAL_DEGREE;
        else
            throw runtime_error(measurementWay + " is not a valid way to measure");

    }

    void DistBasedSampler::operator()(const ArrayXXf &_data, ArrayXXf &_sampleIndices,
                                      const ulli &_numSamples, const ulli &_minSamples) {

        if (m_debug)
            cout << "Entering DistBasedSampler::operator() " << endl;

        if (_data.cols() < _minSamples)
            throw MinSamplesNotFound();

        ArrayXXf cdfs;
        _sampleIndices.resize(_minSamples, _numSamples);

        FindCDF(_data, cdfs);
        for (ulli i = 0; i < _numSamples; i++) {

            _sampleIndices(0, i) = i % _data.cols();
            std::unordered_set<ulli> uniques;
            std::unordered_set<float> xs;
            uniques.insert(_sampleIndices(0, i));

            ulli iterations = 0;
            int tries = 0;

            while (uniques.size() != _minSamples) { //TODO: should be able to ensure 0 repetitions
                auto it = std::upper_bound(cdfs.col(_sampleIndices(0, i)).data(),
                                           cdfs.col(_sampleIndices(0, i)).data() + cdfs.rows(),
                                           (double) rand() / RAND_MAX);
                _sampleIndices(uniques.size(), i) = std::distance(cdfs.col(_sampleIndices(0, i)).data(), it);
                float x = _data(0, _sampleIndices(uniques.size(), i));
                if (_sampleIndices(uniques.size(), i) != _data.cols() &&
                    (m_shouldXUnique && xs.find(x) == xs.end())) {
                    uniques.insert(_sampleIndices(uniques.size(), i));
                    xs.insert(x);
                }
                if (++iterations > m_maxIterationsFactor * _minSamples) {
                    (tries < m_maxTries) ? tries++ : throw MinSamplesNotFound();
                    _sampleIndices(0, i) = rand() % _data.cols();
                    uniques.clear();
                    uniques.insert(_sampleIndices(0, i));
                    iterations = 0;
                }
            }

        }

        if (m_debug)
            cout << "Exiting DistBasedSampler::operator() " << endl;
    }

    void DistBasedSampler::FindCDF(const ArrayXXf &_data, ArrayXXf &_cdfs) {
        if (m_debug)
            cout << "Entering DistBasedSampler::FindCDF() " << endl;

        //create pdf
        double dist;
        _cdfs.resize(_data.cols(), _data.cols());
        for (ulli c = 0; c < _cdfs.cols(); c++) {
            for (ulli r = 0; r < _cdfs.rows(); r++) {
                dist = DistBasedSampler::Distance(_data.col(r), _data.col(c));
                if (dist < m_maxDiff && dist > 0)
                    _cdfs(r, c) = std::exp(-dist * m_sigma / m_maxDiff);
                else
                    _cdfs(r, c) = 0;
            }
        }

        //create cdf
        double colSum;
        for (ulli c = 0; c < _cdfs.cols(); c++) {
            colSum = _cdfs.col(c).sum();
            if (colSum > 0.001) {
                _cdfs(0, c) /= colSum;
                for (ulli r = 1; r < _cdfs.rows(); r++)
                    _cdfs(r, c) = _cdfs(r, c) / colSum + _cdfs(r - 1, c);
            } else {
                //create uniform distribution
                for (ulli r = 0; r < _cdfs.rows(); r++)
                    _cdfs(r, c) = (double) (r + 1) / _cdfs.rows();
            }
        }

        if (m_debug)
            cout << "Exiting DistBasedSampler::FindCDF() " << endl;
    }

    double DistBasedSampler::Distance(const ArrayXf &pt1, const ArrayXf &pt2) {
        switch (m_measureBy) {
            case EUCLIDEAN:
                return (pt1 - pt2).matrix().norm();
            case ABS_VERTICAL_DEGREE:
                return (std::abs(pt2[0] - pt1[0]) <= 0.01) ? 0 : (180 / M_PI) * (std::atan(
                        std::abs((double) (pt2[1] - pt1[1]) / (pt2[0] - pt1[0]))));
        }
    }
	
}
