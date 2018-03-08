#include "stdafx.h"
#include "linearform.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

void DomainLinearFormIntegrator::assembleElementVector(
	const FiniteElement & el, Transformation & trans, VectorXd & vec)
{
	const int ndof = el.numDof();

	vec.resize(ndof);
	vec.setZero();

	// get quadrature rule
	const QuadratureRule *qr = qrule;
	if (qr == nullptr) {
		//int a = 1;
		//int b = 1;
		//qr = QuadratureRule::GetRule(el.getGeom(), el.getOrder() * a + b);
		//assert(qr);
		std::cerr << __FUNCTION__ << ": quadrature rule not set" << std::endl;
		exit(1);
	}

	// shape function
	VectorXd shape(ndof);

	// loop quadrature points
	for (int i = 0; i < qr->numPoints(); i++) {
		const QuadraturePoint &qp = qr->getPoint(i);

		// shape function at this point
		el.calcShape(qp, shape);

		trans.setPoint(qp);

		double wgt = trans.weight(); // detJ

		double val = qfun.eval(trans, qp); // f(x)

		vec += qp.w * val * wgt * shape;
	}
}

////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

