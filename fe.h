#pragma once

#include "sfem.h"
#include "geometry.h"
#include "linalg.h"
#include "quadrature.h"

BEGIN_SFEM_NS;

enum FunctionType
{
	Pk,
	Qk,
};


class FiniteElement 
{

	////////////////////////////////////////////////////////////
protected:

	int ndim; // reference space dimension
	int geom; // basic geometry type
	int ndof; // number of DOFs
	int order; // order of function space
	int func; // function type, P/Q

	// 
	//QuadratureRule nodes;
	Array<VectorXd> nodes;

	////////////////////////////////////////////////////////////
public:

	FiniteElement(int ndim, int geom, int ndof, int order, int fun);

	virtual ~FiniteElement() {}

	////////////////////////////////////////////////////////////
public:

	// number of DOFs
	int numDof() const { return ndof; }

	// number of reference space dimension
	int numDim() const { return ndim; }

	// basic geometry type
	int getGeom() const { return geom; }

	// element order
	int getOrder() const { return order; }

	// function type
	int getFuncType() const { return func; }

	// get node position
	const VectorXd& getNodePos(int i) const { return nodes[i]; }

	// get all nodes positions in matrix (NDIM*NDOF)
	void getNodePositions(MatrixXd &xs) const {
		xs.resize(ndim, ndof);
		for (int i = 0; i < ndof; i++) {
			xs.col(i) = nodes[i];
		}
	}

	////////////////////////////////////////////////////////////
	// functions to be override
public:
	// Calculate shape function
	// at local position
	virtual void calcShape(const VectorXd &p, VectorXd &shape) const = 0;

	// Calcualte shape derivative, size should be (Ndof * Ndim)
	virtual void calcDeriv(const VectorXd &p, MatrixXd &deriv) const = 0;

	////////////////////////////////////////////////////////////
public:

	// Shape function at quadrature point
	virtual void calcShape(const QuadraturePoint &qp, VectorXd &shape) const {
		calcShape(qp.p, shape);
	}

	// Derivative function at quadrature point, size should be (Ndof * Ndim)
	virtual void calcDeriv(const QuadraturePoint &qp, MatrixXd &dshape) const {
		calcDeriv(qp.p, dshape);
	}


};


////////////////////////////////////////////////////////////////////////////////

class SimpleFiniteElement : public FiniteElement
{
public:
	SimpleFiniteElement(int ndim, int geom, int ndof, int order, int fun)
		: FiniteElement(ndim, geom, ndof, order, fun)
	{
		// nothing
	}
};

//
// linear element on 1D line
//
class Linear1DFiniteElement : public SimpleFiniteElement {
public:
	Linear1DFiniteElement();
public:
	virtual void calcShape(const VectorXd &p, VectorXd &shape) const override;
	virtual void calcDeriv(const VectorXd &p, MatrixXd &dshape) const override;
};


// 
// linear element on 2D triangle
// 
class Linear2DFiniteElement : public SimpleFiniteElement {
public:
	Linear2DFiniteElement();
public:
	virtual void calcShape(const VectorXd &p, VectorXd &shape) const override;
	virtual void calcDeriv(const VectorXd &p, MatrixXd &dshape) const override;
};

//
// bi-linear element on 2D quad.
//
class Bilinear2DFiniteElement : public SimpleFiniteElement 
{
public:
	Bilinear2DFiniteElement();
public:
	virtual void calcShape(const VectorXd &p, VectorXd &shape) const override;
	virtual void calcDeriv(const VectorXd &p, MatrixXd &dshape) const override;
};

//
// quadratic element on 1D line
//
class Quadratic1DFiniteElement : public SimpleFiniteElement {
public:
	Quadratic1DFiniteElement();
public:
	virtual void calcShape(const VectorXd &p, VectorXd &shape) const override;
	virtual void calcDeriv(const VectorXd &p, MatrixXd &dshape) const override;
};

// 
// quadratic element on 2D triangle
//
class Quadratic2DFiniteElement : public SimpleFiniteElement {
public:
	Quadratic2DFiniteElement();
public:
	virtual void calcShape(const VectorXd &p, VectorXd &shape) const override;
	virtual void calcDeriv(const VectorXd &p, MatrixXd &dshape) const override;
};


////////////////////////////////////////////////////////////////////////////////







END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////



