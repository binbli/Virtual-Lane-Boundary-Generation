#ifndef ROADSEGMENT_H_
#define ROADSEGMENT_H_

#include<string>
//<libgen.h> is included in Linux; it's only used to extract base name from full file path
//user can easily code this functionality on other platforms manually
#include<libgen.h>
#include<vector>
#include<iostream>
#include<unordered_set>
#include<caffe/caffe.hpp>
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include <omp.h>
#include "Solver.h"
#include "../3rdParty/DBScan/DBScan.h"
#include <opencv2/core/eigen.hpp>
#include "eigen3/Eigen/Core"

namespace LD {

    using namespace caffe;
    using std::string;
    using std::vector;
    using std::cout;
    using std::endl;

    class RoadSegment : public LD::Solver {
    protected:

        /**< corresponds to the width of image that model takes as input */
        int m_cropSizeW;

        /**< corresponds to the height of image that model takes as input */
        int m_cropSizeH;

        double m_meanR;    /**< mean in red channel */
        double m_meanG;        /**< mean in green channel */
        double m_meanB;    /**< mean in blue channel */

        /**< stores 19 categories of cityscape on which the model is trained */
        vector<string> m_labels;

        /**< corresponds to the width of original image given by user */
        int m_originalW;

        /**< corresponds to the height of original image given by user */
        int m_originalH;

        /**< name of current image file being processed	*/
        string m_imgBaseName;

        /**< names of .caffemodel and .prototxt files respectively */
        string m_weightsFile, m_deployFile;

        /**< name of file which stores relative paths of all images to process */
        string m_dataFile;

        /**< directory in which images to segment are stored */
        string m_dataRoot;

        /**< directory in which results of road segmentation should be stored */
        string m_segRoot;

        /**< directory in which overlayed images should be stored */
        string m_overlayedRoot;

        /**< the neural network used to process */
        std::shared_ptr<Net < float> >
        m_net;

        /**< the colormap which would make output more readable */
        cv::Mat m_colormap;

        bool m_saveVizImg;
        string m_vizImgPrefix;
        string m_segImgPrefix;

        vector<cv::Mat> m_inputChannels;

        /**
         * \brief preprocesses the input image (resizing, zero centering) so
         * it can be fed into the model
         * \pre assumes Segmenter::WrapInputLayer() has been called before
         * \param _inputImg input image as given by user
         */
        void Preprocess(const cv::Mat &_inputImg);

        /**
         * \brief outputs raw output of ICNet, where every pixel value stores 0-18
         *  which is the index (in m_labels) of the class that that pixel belongs to
         * \pre assumes Preprocess() has been called before
         * \param _rawOp stores raw output after completion
         */
        void Segment(cv::Mat &_rawOp);

        /**
         * \brief refines _segImg into a viewable form
         * \pre assumes Segment() is called and its output is the first parameter
         * \param _rawOp same as outputted by Segmenter::Segment()
         * \param _segImg outputs b\w image where all white pixels correspond to road
         */
        void PostProcess(const cv::Mat &_rawOp, cv::Mat &_segImg);

        /**
         * \brief saves segmented and overlayed images at desired locations
         * \pre assumes PostProcess() is called and its output is the second parameter
         * \param _inputImg input image as given by user
         * \param _segImg as outputted by Segmenter::PostProcess()
         */
        void Save(cv::Mat &_inputImg, cv::Mat &_segmentedImg);

        /**
         * \brief makes elements of m_inputChannels point to input layer of network
         */
        void WrapInputLayer();

        virtual void ParseXML() override;

    public:
        /**
         * \brief initializes relevant member variables
         * \param _xmlFile xml file with all parameters
         */
        RoadSegment(string _xmlFile);

        /**
         * \brief coordinates calls to other member functions and saves the final input
         */
        virtual void Run() override;

        /**
         * \brief responsible for the core task of Segmenter: takes an input image
         *        and returns a b/w image where road pixels are white
          */
        virtual void operator()(const cv::Mat &_inputImg, cv::Mat &_segImg);


        void OverlayMask(cv::Mat &_segImg,  cv::Mat &_newSegImg);
    };

}

#endif
