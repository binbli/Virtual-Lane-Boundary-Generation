#ifndef DBSCAN_H_
#define DBSCAN_H_

#include<eigen3/Eigen/Dense>
#include"Solver.h"
#include<unordered_set>

namespace LD {

    class DBScan : public Solver {
    protected:

        void ParseXML() override;

        virtual double Distance(const Eigen::ArrayXf &_pt1, const Eigen::ArrayXf &_pt2);

        virtual ulli GetNeighbors(const Eigen::ArrayXXf &_data, const ulli &_ptIndex, vector<ulli> &_neighborIndices,
                                  const std::unordered_set<ulli> &_uniques = std::unordered_set<ulli>()) final;

        void PrintOutputFile(const string &_imageFile, const Eigen::ArrayXXf &_data, const vector<int> &_labels);

        float m_eps;
        int m_minMinPts, m_minX, m_maxMinPts;
        string m_dataFile, m_dataRoot;
        string m_outputRoot, m_outputFilePrefix;
        string m_inputFilePrefix;
        bool m_shouldTranspose;
        float m_declineSlope;

    public:

        enum {
            UNDEFINED = -2,
            NOISE
        };

        //DBScan();
        DBScan(const string &_xmlFile) : Solver(_xmlFile) { ParseXML(); }

        virtual void Cluster(const Eigen::ArrayXXf &_data, vector<int> &_labels);

        virtual void Run();

        virtual void operator()(const Eigen::ArrayXXf &_data, vector<int> &_labels);
    };

}

#endif
