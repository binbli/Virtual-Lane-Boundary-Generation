#ifndef CALIB_DATA_LOADER_H_
#define CALIB_DATA_LOADER_H_

#include<fstream>
#include<string>
#include<iostream>
#include<vector>
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<stdio.h>
#include<stdlib.h>
#include<eigen3/Eigen/Dense>
#include"BaseLD.h"

namespace LD {

	class loadProjectionMat : BaseLD {
    protected:

        virtual void ParseXML() {}

    public:

        loadProjectionMat(string _xmlFile) : BaseLD(_xmlFile) {}

        bool ReadVariable(string _file, string _varName, int _rows, int _cols, Eigen::MatrixXf &_output) {

            if (m_debug)
                cout << "Entering loadProjectionMat::ReadVariable()" << endl;

            FILE *stream = fopen(_file.c_str(), "r");
            char word[80] = "";
            while (!feof(stream) && !ferror(stream) && strcmp(word, (_varName + ":").c_str()) != 0)
                int s = fscanf(stream, "%s", word);
            _output = Eigen::MatrixXf(_rows, _cols);

            if (feof(stream))
                return false;

            for (int r = 0; r < _rows; r++) {
                for (int c = 0; c < _cols; c++) {
                    if (ferror(stream))
                        return false;
                    int t = fscanf(stream, "%f", &_output(r, c));
                }
            }
            fclose(stream);

            if (m_debug)
                cout << "Exiting loadProjectionMat::ReadVariable()" << endl;

            return true;

        };
	};

}

#endif
