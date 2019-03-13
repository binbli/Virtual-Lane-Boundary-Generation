#ifndef SURFACE_DATA_MAKER_H_
#define SURFACE_DATA_MAKER_H_

#include"VeloPtsProjectCam.h"
#include"Utilities.h"

namespace LD {
	class SurfaceDataMaker : public VeloPtsProjectCam {
	protected:
		string m_RoadSegRoot;
		string m_segImgPrefix;
		string m_outputFilePrefix;
		string m_vizImgPrefix;
		bool m_saveVizImg;
		bool m_printProjectedPts;
		int m_minPoints;

		virtual void ParseXML() override;

	public:
		virtual void
		ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints, Mat &_reflectivity,
								 Mat &_inputImg) override;

		SurfaceDataMaker(string _xmlFile) : VeloPtsProjectCam(_xmlFile) { ParseXML(); }
	};
}

#endif
