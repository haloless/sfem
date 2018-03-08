#include "stdafx.h"

#include "transform.h"
#include "element.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;

////////////////////////////////////////////////////////////////////////////////

Transformation::Transformation()
	: geom(-1), phys_dim(-1)
	, qpoint(nullptr)
{
	// TODO
}



// set current point
void Transformation::setPoint(const QuadraturePoint & qp) { 
	qpoint = &qp; 
}

const VectorXd & Transformation::evalTransformedPoint()
{
	assert(qpoint);

	// perform transformation for current point
	// physical position is stored
	transform(*qpoint, xpos);

	return xpos;
}

double Transformation::evalThickness()
{
	if (sfem::IsRZCoord()) {
		// in R-Z coordinate, we use x as R
		evalTransformedPoint();
		return xpos(0);
	}
	else {
		// 
		return 1.0;
	}
}

// 
double Transformation::evalWeight()
{
	// force evaluation of J
	evalJacobian();

	// compute |J|
	detJ = Jmat.determinant();

	// weight
	wgt = detJ;
	// 
	wgt *= thickness();

	return wgt;
}

const MatrixXd & Transformation::evalInverseJacobian()
{
	// force evaluation of J
	evalJacobian();

	invJ = Jmat.inverse();

	return invJ;
}


////////////////////////////////////////////////////////////////////////////////

IsoparametricTransformation::IsoparametricTransformation()
	: Transformation()
	, fe(nullptr)
{
	// TODO
}


// set to finite element
void IsoparametricTransformation::setFE(const FiniteElement * pfe) {
	// bind present FE
	fe = pfe;

	// set transformation
	resetFE();
}

void IsoparametricTransformation::setSameAs(int geomType)
{
	// choose proper finite element
	fe = GetDefaultIsoparamFE(geomType);

	// set this transformation
	resetFE();
}

void IsoparametricTransformation::resetFE()
{
	assert(fe);

	geom = fe->getGeom();
	phys_dim = fe->numDim();
	fe->getNodePositions(nodesMat);
}




const MatrixXd & IsoparametricTransformation::evalJacobian()
{
	assert(fe);
	assert(qpoint);

	const int ndim = fe->numDim();
	const int ndof = fe->numDof();
	assert(nodesMat.rows() == ndim && nodesMat.cols() == ndof);

	// derivative shape function
	dshape.resize(ndof, ndim);
	fe->calcDeriv(*qpoint, dshape);

	// J = 
	Jmat = nodesMat * dshape;

	return Jmat;
}

void IsoparametricTransformation::transform(const QuadraturePoint & qp, VectorXd & x)
{
	assert(fe);

	// shape func.
	shape.resize(fe->numDof());
	fe->calcShape(qp, shape);

	// 
	x.resize(nodesMat.rows());
	x = nodesMat * shape;

}


FiniteElement * IsoparametricTransformation::GetDefaultIsoparamFE(int geomType)
{
	// choose proper finite element
	FiniteElement *fe = nullptr;

	switch (geomType) {
	case BasicGeomType::Line:
		fe = &TheSegmentFE;
		break;
	case BasicGeomType::Triangle:
		fe = &TheTriangleFE;
		break;
	case BasicGeomType::Quadrilateral:
		fe = &TheQuadrilateralFE;
		break;
	default:
		std::cerr << "Invalid geometry=" << geomType << std::endl;
		exit(1);
		break;
	}

	assert(fe);
	return fe;
}



////////////////////////////////////////////////////////////////////////////////



END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

