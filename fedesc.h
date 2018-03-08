#pragma once

#include "sfem.h"
#include "geometry.h"
#include "fe.h"


////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


class FiniteElementDescription
{
protected:



	////////////////////////////////////////////////////////////
public:
	FiniteElementDescription() {}
	virtual ~FiniteElementDescription() {}

	////////////////////////////////////////////////////////////
public:

	virtual const char* getName() const = 0;

	virtual const FiniteElement* getFiniteElementForGeom(int geom) const = 0;

	virtual int numDofForGeom(int geom) const = 0;


};


////////////////////////////////////////////////////////////////////////////////

class FEDescH1 : public FiniteElementDescription
{
protected:

	static const int NumGeom = BasicGeomType::NumGeom;

	// should look like "H1;Dim=2;Order=1"
	std::string name;

	std::array<FiniteElement*, NumGeom> elements;
	std::array<int, NumGeom> dofs;

	////////////////////////////////////////////////////////////
public:
	FEDescH1(int order, int dim);

	virtual ~FEDescH1();

	////////////////////////////////////////////////////////////
public:

	virtual const char * getName() const override;

	virtual const FiniteElement * getFiniteElementForGeom(int geom) const override;

	virtual int numDofForGeom(int geom) const override;

};

////////////////////////////////////////////////////////////////////////////////

class LinearFEDescription : public FiniteElementDescription
{
protected:

	Linear1DFiniteElement fe_segment;
	Linear2DFiniteElement fe_triangle;
	Bilinear2DFiniteElement fe_quadrilateral;

public:
	LinearFEDescription() {}
	virtual ~LinearFEDescription() {}

	virtual const char * getName() const override;

	virtual const FiniteElement * getFiniteElementForGeom(int geom) const override;

	virtual int numDofForGeom(int geom) const override;

};




////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


