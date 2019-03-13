#include"../include/PointsVisualizer.h"
#include"../include/Utilities.h"

namespace LD {

	void PointsVisualizer::operator()(Mat &_inputImg, Eigen::ArrayXXf &_coordinates, const Scalar& _color) {
		Eigen::MatrixXf imgPoints;
		Project(_coordinates, imgPoints);
		for (ulli i = 0; i < imgPoints.rows(); i++) {
			int x = imgPoints(i, 0), y = imgPoints(i, 1);
			if (isValid(y, x, _inputImg.rows, _inputImg.cols))
				circle(_inputImg, Point(x, y), 5, _color, -1);
		}

	}

}