#include "stdafx.h"

#include "fe.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


FiniteElement::FiniteElement(int ndim_, int geom_, int ndof_, int order_, int fun_)
	: ndim(ndim_), geom(geom_), ndof(ndof_), order(order_), func(fun_)
	, nodes(ndof_)

{
	// allocate nodes
	// they must be set by subclass
	for (int i = 0; i < ndof; i++) {
		nodes[i].resize(ndim);
	}

	// TODO
}

////////////////////////////////////////////////////////////////////////////////


Linear1DFiniteElement::Linear1DFiniteElement()
	: SimpleFiniteElement(1, BasicGeomType::Line, 2, 1, FunctionType::Pk)
{
	// set nodes
	nodes[0] << 0;
	nodes[1] << 1;
}

void Linear1DFiniteElement::calcShape(const VectorXd & p, VectorXd & shape) const
{
	double x = p(0);

	shape.resize(2); // ndof=2
	shape(0) = 1.0 - x;
	shape(1) = x;
}

void Linear1DFiniteElement::calcDeriv(const VectorXd & p, MatrixXd & dshape) const
{
	double x = p(0);

	dshape.resize(2, 1); // ndof=2, ndim=1
	dshape << -1, +1;
}


////////////////////////////////////////////////////////////////////////////////

Linear2DFiniteElement::Linear2DFiniteElement()
	: SimpleFiniteElement(2, BasicGeomType::Triangle, 3, 1, FunctionType::Pk)
{
	// set nodes
	nodes[0] << 0, 0;
	nodes[1] << 1, 0;
	nodes[2] << 0, 1;
}

void Linear2DFiniteElement::calcShape(const VectorXd & p, VectorXd & shape) const {
	double x = p(0);
	double y = p(1);

	shape.resize(3); // ndof=3
	shape(0) = 1.0 - x - y;
	shape(1) = x;
	shape(2) = y;
}

void Linear2DFiniteElement::calcDeriv(const VectorXd & p, MatrixXd & dshape) const {
	dshape.resize(3, 2); // ndof=3, ndim=2
	dshape << 
		-1, -1, 
		+1,  0, 
		 0, +1;
}


////////////////////////////////////////////////////////////////////////////////


Bilinear2DFiniteElement::Bilinear2DFiniteElement()
	: SimpleFiniteElement(2, BasicGeomType::Quadrilateral, 4, 1, FunctionType::Qk)
{
	// set nodes
	nodes[0] << 0, 0;
	nodes[1] << 1, 0;
	nodes[2] << 1, 1;
	nodes[3] << 0, 1;
}

void Bilinear2DFiniteElement::calcShape(const VectorXd & p, VectorXd & shape) const
{
	double x = p(0);
	double y = p(1);

	shape.resize(4); // ndof=4
	shape << (1.0 - x)*(1.0 - y), x*(1.0 - y), x*y, (1.0 - x)*y;
}

void Bilinear2DFiniteElement::calcDeriv(const VectorXd & p, MatrixXd & dshape) const
{
	double x = p(0);
	double y = p(1);

	dshape.resize(4, 2); // ndof=4, ndim=2
	dshape <<
		-1.0 + y, -1.0 + x,
		1.0 - y,  -x,
		y,        x,
		-y,       1.0 - x;
}

////////////////////////////////////////////////////////////////////////////////


Quadratic2DFiniteElement::Quadratic2DFiniteElement()
	: SimpleFiniteElement(2, BasicGeomType::Triangle, 6, 2, FunctionType::Pk)
{
	// set nodes
	nodes[0] << 0, 0;
	nodes[1] << 1, 0;
	nodes[2] << 0, 1;
	nodes[3] << 0.5, 0;
	nodes[4] << 0.5, 0.5;
	nodes[5] << 0, 0.5;
}

void Quadratic2DFiniteElement::calcShape(const VectorXd & p, VectorXd & shape) const
{
	const double x = p(0);
	const double y = p(1);
	double a1 = 1.0 - x - y;
	double a2 = x;
	double a3 = y;

	shape.resize(6); // ndof=6
	shape <<
		a1 * (2 * a1 - 1),
		a2 * (2 * a2 - 1),
		a3 * (2 * a3 - 1),
		4 * a1 * a2,
		4 * a2 * a3,
		4 * a3 * a1;
}

void Quadratic2DFiniteElement::calcDeriv(const VectorXd & p, MatrixXd & dshape) const
{
	const double x = p(0);
	const double y = p(1);

	dshape.resize(6, 2); // ndof=6, ndim=2
	dshape <<
		4 * (x + y) - 3, 4 * (x + y) - 3,
		4 * x - 1, 0,
		0, 4 * y - 1,
		-4 * (2 * x + y - 1), -4 * x,
		4 * y, 4 * x,
		-4 * y, -4 * (x + 2 * y - 1);
}


////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

