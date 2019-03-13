#include "RoadSegment.h"

namespace LD {

    RoadSegment::RoadSegment(string _xmlFile) : Solver(_xmlFile),m_cropSizeH(1025), m_cropSizeW(2049),
                                            m_meanR(123.68),
                                            m_meanG(116.779), m_meanB(103.939),
                                            m_deployFile("3rdParty/ICNet/evaluation/prototxt/icnet_cityscapes.prototxt"),
                                            m_weightsFile("3rdParty/ICNet/evaluation/model/icnet_cityscapes_train_30k.caffemodel") {
        //initializes relevant member variables
        ParseXML();
/*        m_labels = {"road", "sidewalk", "building", "wall", "fence", "pole", "traffic light", "traffic sign",
                    "vegetation", "terrain", "sky", "person", "rider", "car", "truck", "bus", "train", "motorcycle",
                    "bicycle"};*/

        m_labels = {"road"};

        Caffe::set_mode(Caffe::CPU);
        m_net.reset(new Net<float>(m_deployFile, TEST));
        m_net->CopyTrainedLayersFrom(m_weightsFile);

        if (m_debug)
            cout << "net initialized with .prototxt and .caffemodel files" << endl;

        Blob<float> *inputLayer = m_net->input_blobs()[0];
        inputLayer->Reshape(1, 3, m_cropSizeH, m_cropSizeW); // numChannels=3
        m_net->Reshape();

        m_colormap = cv::imread("./src/colormap.png", cv::IMREAD_GRAYSCALE);

        WrapInputLayer();

        if (m_debug)
            cout << "Exiting RoadSegment::RoadSegment()" << endl;
    }

    void RoadSegment::ParseXML() {
        m_xml = m_xml.child("RoadSegment");
        pugi::xml_node solverInstance = m_xml.child("SolverInstance");

        m_dataRoot = solverInstance.attribute("dataRoot").as_string();
        m_dataFile = solverInstance.attribute("dataFile").as_string();
        m_segRoot = solverInstance.attribute("segRoot").as_string();
        m_overlayedRoot = solverInstance.attribute("overlayedRoot").as_string();
        m_saveVizImg = solverInstance.attribute("saveVizImg").as_bool(true);
        m_vizImgPrefix = solverInstance.attribute("vizImgPrefix").as_string();
        m_segImgPrefix = solverInstance.attribute("segImgPrefix").as_string();

        if (m_dataRoot.empty() || m_dataFile.empty() || m_segRoot.empty() || m_overlayedRoot.empty() ||
            m_segImgPrefix.empty() || (m_vizImgPrefix.empty() && m_saveVizImg))
            throw runtime_error(
                    "One of the following attributes are missing in SolverInstance node of RoadSegment: dataRoot, dataFile, segRoot, overlayedRoot, saveVizImg, vizImgPrefix, segImgPrefix");
    }

    void RoadSegment::operator()(const cv::Mat &_inputImg, cv::Mat &_segImg) {
        if (m_debug)
            cout << "Entering RoadSegment::()" << endl;

        if(_inputImg.empty())
            cout << "RoadSegment::operator() inputImg is empty!" << endl;

        m_originalW = _inputImg.cols, m_originalH = _inputImg.rows;
        Preprocess(_inputImg);
        Segment(_segImg);
        PostProcess(_segImg, _segImg);

        if (m_debug)
            cout << "Exiting Segmenter::()" << endl;
    }

    void RoadSegment::Run() {
        //coordinates calls to other member functions and saves the final input

        if (m_debug)
            cout << "Entering RoadSegment::Run()" << endl;

        std::ifstream fin(m_dataFile.c_str());
        string line;

        while (std::getline(fin, line)) {
            m_imgBaseName = basename(const_cast<char *>(line.c_str()));
            cv::Mat inputImg, segImg;

            inputImg = cv::imread(m_dataRoot + "/" + line);

            if(inputImg.empty())
                std::cerr << "could not open or find " << m_dataRoot + "/" + line << endl;

            if (m_debug)
                cout << "Successfully opened " << line << endl;

            this->operator()(inputImg, segImg);
            if (m_debug)
                Save(inputImg, segImg);
        }

        if (m_debug)
            cout << "Exiting RoadSegment::Run()" << endl;

    }

    void RoadSegment::WrapInputLayer() {
        //makes elements of m_inputChannels point to input layer of network

        if (m_debug)
            cout << "Entering RoadSegment::WrapInputLayer()" << endl;

        m_inputChannels = vector<cv::Mat>(3, cv::Mat(cv::Size(m_cropSizeW, m_cropSizeH), CV_32FC3));
        Blob<float> *inputLayer = m_net->input_blobs()[0];
        int width = inputLayer->width();
        int height = inputLayer->height();
        float *inputData = inputLayer->mutable_cpu_data();

        for (int i = 0; i < inputLayer->channels(); ++i) {
            cv::Mat channel(height, width, CV_32FC1, inputData);
            m_inputChannels[i] = channel;
            inputData += width * height;
        }

        if (m_debug)
            cout << "Exiting RoadSegment::WrapInputLayer()" << endl;

    }

    void RoadSegment::Preprocess(const cv::Mat &_inputImg) {
        //preprocesses the input image (resizing, zero centering) so it can be fed into the model

        if (m_debug)
            cout << "Entering RoadSegment::Preprocess()" << endl;

        if(_inputImg.empty())
            std::cerr << "RoadSegment::Preprocess() _inputImg is empty!" << endl;

        cv::Mat equalizedImg;
        cv::cvtColor(_inputImg, equalizedImg, cv::COLOR_BGR2YUV);
        vector<cv::Mat> channels;
        cv::split(equalizedImg, channels);
        cv::equalizeHist(channels[0], channels[0]);
        cv::merge(channels, equalizedImg);
        cv::cvtColor(equalizedImg, equalizedImg, CV_YUV2BGR);
        cv::Mat img;
        equalizedImg.convertTo(img, CV_32FC3);
        cv::resize(img, img, cv::Size(m_cropSizeW, m_cropSizeH));
        cv::subtract(img, cv::Scalar(m_meanB, m_meanG, m_meanR), img);
        cv::split(img, m_inputChannels); //puting image in input layer

        if (m_debug) {
            cout << "Image set to be input of network" << endl;
            cout << "Exiting RoadSegment::Preprocess()" << endl;
        }
    }

    void RoadSegment::Segment(cv::Mat &_rawOp) {
        //outputs raw output of ICNet, where every pixel value stores 0-18 which is the
        //index (in m_labels) of the class that that pixel belongs to

        if (m_debug)
            cout << "Entering RoadSegment::Segment()" << endl;


         //const std::vector<caffe::Blob < float> *> &results = m_net->ForwardPrefilled(&loss);
        const std::vector<caffe::Blob < float> *> &results = m_net->ForwardPrefilled();

        if (m_debug)
            cout << "Forward pass successfully finsihed " << endl;

        Blob<float> *outputLayer = m_net->output_blobs()[0];
        cv::Mat score(outputLayer->channels(), outputLayer->width() * outputLayer->height(), CV_32FC1,
                      outputLayer->mutable_cpu_data());
        cv::transpose(score, score);
        _rawOp = cv::Mat(outputLayer->height(), outputLayer->width(), CV_8UC1);

        double maxVal;
        cv::Point maxIndex;
        std::unordered_set<int> thingsIdentified;

        for (int i = 0; i < score.rows; i++) {
            minMaxLoc(score.row(i), 0, &maxVal, 0, &maxIndex);
            _rawOp.at<uchar>(i) = maxIndex.x;
/*            if (m_debug) {
                if (thingsIdentified.find(maxIndex.x) == thingsIdentified.end()) {
                    //cout << m_labels[maxIndex.x] << " found!" << endl;
                    thingsIdentified.insert(maxIndex.x);
                }
            }*/
        }

        cv::resize(_rawOp, _rawOp, cv::Size(m_originalW, m_originalH));

        if (m_debug) {
            cout << "Found segmented image" << endl;
            cout << "Exiting RoadSegment::Segment()" << endl;
        }
    }

    void RoadSegment::PostProcess(const cv::Mat &_rawOp, cv::Mat &_segImg) {
        //saves segmented and overlayed images at desired locations
        if (m_debug) {
            cout << "Entering RoadSegment::PostProcess()" << endl;
            cout << "The number of channels and depths of the segmented image and the colormap must be the same"
                 << endl;

            cout << "segmented image(rows, columns, channels, depth): " <<
                 _segImg.rows << " " << _segImg.cols << " " <<
                 _segImg.channels() << " " << _segImg.depth() << endl;

            cout << "Colormap(rows, columns, channels, depth): " <<
                 m_colormap.rows << " " << m_colormap.cols << " " <<
                 m_colormap.channels() << " " << m_colormap.depth() << endl;
        }

        cv::LUT(_rawOp, m_colormap, _segImg);

        if (m_debug)
            cout << "Exiting RoadSegment::PostProcess()" << endl;
    }

    void RoadSegment::Save(cv::Mat &_inputImg, cv::Mat &_segImg) {
        // saves segmented and overlayed images at desired locations
        cv::Mat origimg = _inputImg.clone();
        //create Segmentation
        if (m_debug)
            cout << "Entering RoadSegment::Save()" << endl;

        string segImgName = m_segImgPrefix + m_imgBaseName;
        vector<cv::Mat> channels;

        if (m_debug)
            cout << "The number of rows, columns, and depths must be the same for different matrices" << endl;

        imwrite(m_segRoot + "/" + segImgName, _segImg);

        if (m_debug)
            cout << "Road Image Mask " << segImgName << " is saved in " << m_segRoot << endl;

        if (m_saveVizImg) {

            //create overlay
            cv::Mat overlayed;
            string overlayedImgName = m_vizImgPrefix + m_imgBaseName;

            vector<cv::Mat> channels;
            cv::Mat black = cv::Mat::zeros(_segImg.rows, _segImg.cols, CV_8UC1);
            _segImg.convertTo(_segImg, CV_8UC1);
            channels.push_back(black);
            channels.push_back(_segImg);
            channels.push_back(black);
            cv::merge(channels, _segImg);

            cv::addWeighted(_segImg, 0.5, origimg, 0.5, 0.0, overlayed);

            imwrite(m_overlayedRoot + "/" + overlayedImgName, overlayed);
            if (m_debug) {
                cout << "Overlayed Road Image " << overlayedImgName << " is saved in " << m_overlayedRoot << endl;
                cout << "Exiting RoadSegment::Save()" << endl;
            }
        }

        if (m_debug)
            cout << "Exiting RoadSegment::Save()" << endl;

    }

    void RoadSegment::OverlayMask(cv::Mat &_segImg, cv::Mat &_newSegImg){
        _newSegImg = _segImg.clone();
        vector<cv::Mat> channels;
        cv::Mat black = cv::Mat::zeros(_newSegImg.rows, _newSegImg.cols, CV_8UC1);
        _newSegImg.convertTo(_newSegImg, CV_8UC1);
        channels.push_back(black);
        channels.push_back(_newSegImg);
        channels.push_back(black);
        cv::merge(channels, _newSegImg);
    }

}