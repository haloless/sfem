#pragma once

#include "sfem.h"
#include "linalg.h"
#include "quadrature.h"
#include "fe.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

//
// Abstract element transformation
// Reference space -> Physical space
// 
class Transformation
{
protected:

	// geometry on which the transformation is defined
	// this determines the reference dimension
	int geom;

	// physical space dimensions
	int phys_dim;

	// target point
	const QuadraturePoint *qpoint;

	// physical point
	VectorXd xpos;

	// for Jacobian, (Nphys*Nref)
	// (x,y,z) <-> (r,s,t)
	// [dx/dr, dx/ds, dx/dt]
	// [dy/dr, dy/ds, dy/dt]
	// [dz/dr, dz/ds, dz/dt]
	MatrixXd Jmat; 

	// weight of this transformation
	// w = det(J) * 
	double wgt;
	// det(J)
	double detJ;

	// inverse of Jacobian
	// must be explicitly calculated
	MatrixXd invJ;


	////////////////////////////////////////////////////////////
public:

	Transformation();

	virtual ~Transformation() { /* nothing */ }

public:

	// reference space
	int numDim() const { return GeometryDimension(geom); }
	
	// physical space
	int numPhysDim() const { return phys_dim; }

	// whether reference space and physical space have the same dimensions
	bool isEqualDim() const { return numDim() == numPhysDim(); }

	// fundamental geometry
	int getGeom() const { return geom; }

	// set current point
	void setPoint(const QuadraturePoint &qp);
	
	// get current point, must be set previously
	const QuadraturePoint& getPoint() const { return *qpoint; }

	// whether the point has been set
	bool hasPoint() const { return qpoint != nullptr; }

	////////////////////////////////////////////////////////////
public:

	// the following 

	const VectorXd& evalTransformedPoint();

	double evalThickness();

	virtual const MatrixXd& evalJacobian() = 0;
	
	double evalWeight();

	const MatrixXd& evalInverseJacobian();

	// transform reference point to 
	virtual void transform(const QuadraturePoint &qp, VectorXd &x) = 0;
	//virtual void transform(const QuadratureRule &q, MatrixXd &m) = 0;

public:

	// TODO

	// 
	// Note the definition of Jacobian
	// whose dimension is [N_physical * N_reference]
	// If (x,y,z) <-> (r,s,t)
	// then 
	//		dx/dr dx/ds dx/dt
	// J =	dy/dr dy/ds dy/dt
	//		dz/dr dz/ds dz/dt
	//
	const MatrixXd& Jacobian() {
		return evalJacobian();
	}

	// W = det(J) * thickness
	double weight() {
		return evalWeight();
	}

	const MatrixXd& inverseJacobian() {
		return evalInverseJacobian();
	}

	double thickness() {
		return evalThickness();
	}
};

////////////////////////////////////////////////////////////////////////////////

class IsoparametricTransformation : public Transformation
{
protected:
	
	// shape function, (ndof)
	VectorXd shape;
	// derivative shape function, (ndof*ndim)
	MatrixXd dshape;

	// the corresponding finite element
	const FiniteElement *fe;

	// copy of node positions, (ndim*ndof)
	// [x1, x2, x3, ...]
	// [y1, y2, y3, ...]
	// [z1, z2, z3, ...]
	MatrixXd nodesMat;

	////////////////////////////////////////////////////////////
public:

	IsoparametricTransformation();

	virtual ~IsoparametricTransformation() { /* nothing */ }

public:

	// get node positions buffer, (ndim*ndof)
	MatrixXd& getNodesMat() { return nodesMat; }
	// get node positions buffer, (ndim*ndof)
	const MatrixXd& getNodesMat() const { return nodesMat; }

	// 
	const FiniteElement* getFE() const { return fe; }

	// set to finite element
	void setFE(const FiniteElement *pfe);

	// Set to corresponding finite element
	// This uses the singleton FEs pre-defined 
	void setSameAs(int geomType);

protected:

	// reset according to present finite element
	void resetFE();

	////////////////////////////////////////////////////////////
public:

	virtual const MatrixXd& evalJacobian() override;

	virtual void transform(const QuadraturePoint &qp, VectorXd &x) override;

	////////////////////////////////////////////////////////////
public:

	// Get default mapping from Geometry to Finite Element
	static FiniteElement* GetDefaultIsoparamFE(int geomType);
};











////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////



