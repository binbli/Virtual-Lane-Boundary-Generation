#ifndef BaseLD_H_
#define BaseLD_H_

#include<iostream>
#include<vector>
#include<string>
#include<algorithm>
#include<math.h>
#include<pugixml.hpp>
#include<exception>
#include<libgen.h>
#include<unordered_set>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/SVD>
#include<eigen3/Eigen/Core>
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/core/eigen.hpp>
#include<stdio.h>
#include<exception>
#include<stdlib.h>
#include<omp.h>
#include<fstream>
#include<sstream>
#include<memory>
#include<iomanip>
#include <limits>
#include<pybind11/pybind11.h>
#include<pybind11/embed.h>
#include<pybind11/eigen.h>

namespace LD {

    using std::cout;
    using std::endl;
    using std::string;
    using std::vector;
    using std::runtime_error;

    class BaseLD {

    protected:
        string m_xmlFileName;
        pugi::xml_node m_xml;
        pugi::xml_document m_xmlDoc;

        virtual void ParseXML() = 0;

    public:

#define PRINT(name) print(#name, (name))

        template<class T>
        void print(std::string name, T value) {
            std::cout << "----------------- Variable : " << name << "   -----------------" << std::endl;
            std::cout << value << std::endl;
            std::cout << "------------------------------------------------------" << std::endl;
            getchar();
            std::cout << std::endl;
        }

#ifndef PI
#define PI 3.14159265358979323846
#endif

        bool m_debug; //to be changed

        typedef unsigned long long int ulli;

        virtual void ParseXML(string _file) final {
            pugi::xml_parse_result docStatus = m_xmlDoc.load_file(_file.c_str());

            if (!docStatus)
                throw runtime_error("Error with " + _file + ": " + docStatus.description());

            m_xml = m_xmlDoc.document_element();

            m_debug = m_xml.child("Main").attribute("debug").as_bool();
        }

        BaseLD(string _file) : m_xmlFileName(_file) {
            ParseXML(_file); //will call above function
        }

    };

}

#endif
