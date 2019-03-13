#ifndef VELO_PROJECTOR_H_
#define VELO_PROJECTOR_H_

#include"loadProjectionMat.h"
#include"Solver.h"

namespace LD {

	using namespace cv;

	class VeloPtsProjectCam : public Solver {
	protected:

		string m_dataRoot;
		string m_calibRoot;
		string m_veloRoot;
		string m_outputRoot;
		string m_dataFile;
		string m_imgBaseName;
		double m_minX;
		int m_retentionFrequency;
		int m_camNum;
		Eigen::MatrixXf m_PRect, m_Tr, m_RRect, m_projectionMat;
		loadProjectionMat m_calibDataLoader;

		virtual void Project(const Mat &_veloPoints, Eigen::MatrixXf &_veloImg, Mat &_reflectivity) final;

		virtual void Project(Eigen::ArrayXXf &_veloPoints, Eigen::MatrixXf &_veloImg) final;

		void ComputeProjMat();

		virtual void
		ProcessProjectedLidarPts(Eigen::MatrixXf &_veloImg, const Mat &_veloPoints, Mat &_reflectivity,
								 Mat &_inputImg) = 0;

	public:

		VeloPtsProjectCam(string _xmlFile);

		virtual void ParseXML() override;

		virtual void Run() override;

		virtual void ReadVeloData(string _binFile, Mat &_veloPoints);

		friend class roadEdgeDetector;
			
	};

}

#endif
