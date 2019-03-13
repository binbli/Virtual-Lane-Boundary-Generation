//TODO: Check pre conditions of functions. Some only except CV_8U type images

#ifndef LMsfromCam_H_
#define LMsfromCam_H_

//<libgen.h> is included in Linux; it's only used to extract base name from full file path
//user can easily code this functionality on other platforms manually
#include<libgen.h>
#include<queue>
#include<fstream>
#include<unordered_set>
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/ml/ml.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<limits>
#include"Solver.h"

namespace LD {

	using namespace cv;

	class LMsfromCam : public Solver {
	protected:

		/**< directory which stores original images (not segmented)  */
		string m_dataRoot;

		/**< name of file which stores relative paths of all original images to process */
		string m_dataFile;

		/**< directory in which segmented images to refine are stored;
             segmentation results of original images must have "segmented_" prepended to their names
             Example: if uu_000000.png is name of original image, then segmented_uu_000000.png must be
             the name of its segmented version stored in segmented_root_    */
		string m_segRoot;

		/**< directory in which refined image results should be stored */
		string m_refinedRoot;

		/**< name of current image file being processed	*/
		string m_imgBaseName;

		string m_vizImgPrefix;
		bool m_saveVizImg;
		string m_refinedImgPrefix;
		string m_roadbRoot;
		/**
         * \brief extracts lanes from segmented image
         * \param _extractedImg image which is black except where there is road
         * \param _refinedImg stores output
         */
		virtual void LMsfromImg(const Mat &_original, const Mat &_segImg, Mat &_refinedImg) = 0;


		virtual void ParseXML() override;

	public:
		/**
         * \brief initializes relevant member variables
         * \param _dataRoot root of directory where all original images are stored
         * \param _dataFile name of file which lists relative paths (to _dataRoot) of all images original images
         * \param _segRoot directory where segmented images are stored
         * \param _refinedRoot directory where refined images should be stored
         */
		LMsfromCam(string _xmlFile);

		/**
         * \brief coordinates calls to other member functions and saves the final input
         */
		virtual void Run() override;

		virtual void operator()(const Mat &_original, const Mat &_segImg, Mat &_refinedImg);

	};

}

#endif
