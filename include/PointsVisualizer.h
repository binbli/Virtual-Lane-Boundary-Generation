#ifndef MODEL_VISUALIZER_H_
#define MODEL_VISUALIZER_H_

#include"VeloPtsProjectCam.h"

namespace LD {

	class PointsVisualizer : public VeloPtsProjectCam {

	public:

		PointsVisualizer(const string &_xmlFile) : VeloPtsProjectCam(_xmlFile) {}

		virtual void operator()(Mat &_inputImg, Eigen::ArrayXXf &_coordinates, const Scalar& _color = Scalar(0, 255, 0));

		virtual void
		ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints, Mat &_reflectivity,
								 Mat &_inputImg) {}
	};

}

#endif