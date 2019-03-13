#ifndef K_PERCENT_EXTRACTOR_H_
#define K_PERCENT_EXTRACTOR_H_

#include"LMsfromCam.h"

namespace LD {

    class KPercentExtractor : public LMsfromCam {

    private:

        std::vector<int> temp;

    protected:

        int m_kpercent;
        int m_dilateSize;

        virtual void ParseXML() override;

    public:

        void runBFS(cv::Mat &_noisefreeImg, int x, int y, int &_nclusters, vector<cv::Point> &_pts);

        int noiseRemover(cv::Mat &_segImg);

        virtual void Preprocess(const Mat &_original, const Mat &_segImg, Mat &_preprocessed);

        virtual void LMsfromImg(const Mat &_original, const Mat &_segImg, Mat &_refinedImg) override;

        KPercentExtractor(string _xmlFile) : LMsfromCam(_xmlFile) {ParseXML(); };

    };

}

#endif
