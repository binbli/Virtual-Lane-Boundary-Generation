#ifndef LINE_3D_T_LINKAGE_H_
#define LINE_3D_T_LINKAGE_H_

#include"TLinkage.h"

namespace LD {

    class Line3DTLinkage : public TLinkage {

    protected:

        virtual ArrayXf GenerateHypothesis(const vector<ArrayXf> &_samples);

        virtual double Distance(ArrayXf _dataPoint, ArrayXf _model) override;

        virtual void FitModel(const ArrayXXf &_clusters, ArrayXf &_model) override;

        virtual void RefineModels(const vector<ArrayXf> &_models, const ArrayXXf &_data, const ArrayXf &_clusters,
                                  const std::unordered_map<int, ulli> &_clusterID2Index,
                                  const std::unordered_map<int, vector<ulli> > &_clusterID2PtIndices,
                                  vector<ArrayXf> &_refinedModels, const int &_noiseIndex = -1);

        virtual void ParseXML() override;

        virtual void ShiftModelBy(ArrayXf &_model, const float &_shiftBy);

        virtual bool IsModelOnRight(const ArrayXf &_model) override;


    public:

        Line3DTLinkage(string _xmlFile) : TLinkage(2, 6, _xmlFile) { ParseXML(); }
    };

}
#endif
