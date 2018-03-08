#include "stdafx.h"

#include "fedesc.h"

#include <iostream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////



FEDescH1::FEDescH1(int order, int dim)
{
	// generate name
	std::ostringstream oss;
	oss << "H1" << ";D=" << dim << ";P=" << order;
	name = oss.str();

	// 
	elements.fill(NULL);
	dofs.fill(0);


}

FEDescH1::~FEDescH1()
{
}

const char * FEDescH1::getName() const
{
	return name.c_str();
}

const FiniteElement * FEDescH1::getFiniteElementForGeom(int geom) const
{
	return nullptr;
}

int FEDescH1::numDofForGeom(int geom) const
{
	return 0;
}

////////////////////////////////////////////////////////////////////////////////



const char * LinearFEDescription::getName() const
{
	static const std::string name("Linear");
	return name.c_str();
}

const FiniteElement * LinearFEDescription::getFiniteElementForGeom(int geom) const
{
	if (geom == BasicGeomType::Line) {
		return &fe_segment;
	}
	else if (geom == BasicGeomType::Triangle) {
		return &fe_triangle;
	}
	else if (geom == BasicGeomType::Quadrilateral) {
		return &fe_quadrilateral;
	}
	else {
		std::cerr << __FUNCTION__ << ": invalid geometry=" << geom << std::endl;
		exit(1);
		return nullptr;
	}
}

int LinearFEDescription::numDofForGeom(int geom) const
{
	return 0;
}


////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


