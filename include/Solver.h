#ifndef SOLVER_H_
#define SOLVER_H_
#include"BaseLD.h"

namespace LD {

	class Solver : public BaseLD {
	protected:
		virtual void ParseXML() override { m_xml = m_xml.child("Solvers"); }

	public:

		virtual void Run() = 0;

		bool m_laneAssess;

		Solver(string _xmlFile) : BaseLD(_xmlFile) { ParseXML(); }
	};

}

#endif
