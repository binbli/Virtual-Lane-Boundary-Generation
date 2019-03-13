#include "../include/Solver.h"
#include"../include/LaneDetector.h"
#include"RoadSegment.h"
#include"LMsfromCam.h"
#include"../include/KPercentExtractor.h"
#include"VeloPtsProjectCam.h"
#include"LMsintersection.h"
#include"../3rdParty/TLinkage/TLinkage.h"
#include"../3rdParty/TLinkage/BSplineTLinkage.h"
#include"../3rdParty/DBScan/DBScan.h"
#include"MapGenerator.h"
#include<pugixml.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include<memory>

using namespace LD;

int main(int argc, char *argv[]) {
    if (argc != 2)
        throw runtime_error("xml file name expected as first argument");

    pybind11::initialize_interpreter();

    pugi::xml_document xml;
    pugi::xml_parse_result status = xml.load_file(argv[1]);
    if (!status)
        throw runtime_error(string("xml error") + status.description());

    string solver = xml.document_element().child("Main").attribute("solver").as_string();

    std::unique_ptr<LD::Solver> solverPtr = nullptr;

    if (boost::iequals(solver, "RoadSegment"))
        solverPtr = std::make_unique<RoadSegment>(argv[1]);
    else if (boost::iequals(solver, "KPercentExtractor"))
        solverPtr = std::make_unique<KPercentExtractor>(argv[1]);
    else if (boost::iequals(solver, "LMsintersection"))
        solverPtr = std::make_unique<LMsintersection>(argv[1]);
    else if (boost::iequals(solver, "BSplineTLinkage"))
        solverPtr = std::make_unique<BSplineTLinkage>(argv[1]);
    else if (boost::iequals(solver, "MapGenerator"))
        solverPtr = std::make_unique<MapGenerator>(argv[1]);
    else if (boost::iequals(solver, "LaneDetector"))
        solverPtr = std::make_unique<LaneDetector>(argv[1]);
    else
        throw runtime_error("No such solver implemented: " + solver);

    solverPtr->Run();

    pybind11::finalize_interpreter();

    return 0;
}
