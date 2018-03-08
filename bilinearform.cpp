#include "stdafx.h"
#include "bilinearform.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////


void BilinearFormIntegrator::assembleElementMatrix(
	const FiniteElement & elem, Transformation & trans, MatrixXd & elmat)
{
	std::cerr << __FUNCTION__ << ": not implemented" << std::endl;
	exit(1);
}

void BilinearFormIntegrator::assembleElementVector(
	const FiniteElement & elem, Transformation & trans, const VectorXd & v, VectorXd & elvec)
{
	std::cerr << __FUNCTION__ << ": not implemented" << std::endl;
	exit(1);
}

const QuadratureRule * BilinearFormIntegrator::ensureQuadratureRule() const
{
	// the rule set previously
	const QuadratureRule *qr = getQuadratureRule();

	if (!qr) {
		// quadrature rule not set
		// TODO load proper rules here
		std::cerr << __FUNCTION__ << ": quadrature rule not set" << std::endl;
		exit(1);
	}

	return qr;
}

////////////////////////////////////////////////////////////////////////////////


void LaplacianIntegrator::assembleElementMatrix(
	const FiniteElement & elem, Transformation & trans, 
	MatrixXd & elmat)
{
	const int ndof = elem.numDof();
	const int ndim = elem.numDim();

	//const int spacedim = trans.numPhysDim();
	if (/*ndim != spacedim*/!trans.isEqualDim()) {
		std::cerr << __FUNCTION__ << ": dimension do not match" << std::endl;
		exit(1);
	}

	// local buffers
	MatrixXd dshape(ndof, ndim);
	MatrixXd dshapedxt(ndof, ndim);

	// result matrix
	elmat.resize(ndof, ndof);
	elmat.setZero();

	// get quadrature rule
	const QuadratureRule *qr = ensureQuadratureRule();
	assert(qr);

	// loop quadrature points
	for (int i = 0; i < qr->numPoints(); i++) {
		const QuadraturePoint &qp = qr->getPoint(i);

		// derivative shape function
		elem.calcDeriv(qp, dshape);

		// current transformation
		trans.setPoint(qp);

		double wgt = trans.weight();

		MatrixXd Jmat = trans.Jacobian();
		MatrixXd invJ = Jmat.inverse();

		dshapedxt = dshape * invJ;

		elmat += dshapedxt * dshapedxt.transpose() * wgt * qp.w;
	}
}

void LaplacianIntegrator::assembleElementVector(
	const FiniteElement & elem, Transformation & trans, 
	const VectorXd & v, VectorXd & elvec)
{
	const int ndof = elem.numDof();
	const int ndim = elem.numDim();

	// local buffers
	MatrixXd dshape(ndof, ndim);

	std::cerr << __FUNCTION__ << ": not implemented" << std::endl;
	exit(1);
}


////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

