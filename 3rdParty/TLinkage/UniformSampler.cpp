#include"UniformSampler.h"

namespace LD {

    void UniformSampler::operator()(const ArrayXXf &_data, ArrayXXf &_sampleIndices,
                                    const ulli &_numSamples, const ulli &_minSamples) {

        if (m_debug)
            cout << "Entering UniformSampler::operator() " << endl;


        if (_data.cols() < _minSamples)
            throw MinSamplesNotFound();

        _sampleIndices.resize(_minSamples, _numSamples);


        for (ulli i = 0; i < _numSamples; i++) {
            _sampleIndices(0, i) = i % _data.cols();
            std::unordered_set<ulli> uniques;
            uniques.insert(_sampleIndices(0, i));
            ulli iterations = 0;
            int tries = 0;
            while (uniques.size() != _minSamples) {
                _sampleIndices(uniques.size(), i) = rand() % _data.cols();
                uniques.insert(_sampleIndices(uniques.size(), i));

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
            cout << "Exiting UniformSampler::operator() " << endl;
    }

}
