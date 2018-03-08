#pragma once

#include "sfem.h"
#include "linalg.h"
#include "fe.h"
#include "quadrature.h"
#include "transform.h"
#include "coefficient.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////


class LinearFormIntegrator
{
protected:
	const QuadratureRule *qrule;

	////////////////////////////////////////////////////////////
public:

	LinearFormIntegrator(const QuadratureRule *qr = nullptr)
		: qrule(qr)
	{}

	virtual ~LinearFormIntegrator() {}

	void setQuadratureRule(const QuadratureRule *qr) { qrule = qr; }
	const QuadratureRule* getQuadratureRule() const { return qrule; }

	////////////////////////////////////////////////////////////
public:
	// for Element and Transformation
	// assemble element vector
	virtual void assembleElementVector(
		const FiniteElement &el, 
		Transformation &tr, 
		VectorXd &vec) = 0;

	//
};


////////////////////////////////////////////////////////////////////////////////

class DomainLinearFormIntegrator : public LinearFormIntegrator
{
protected:
	typedef LinearFormIntegrator super_type;

	//VectorXd shape;
	ScalarCoefficient &qfun;

public:
	DomainLinearFormIntegrator(ScalarCoefficient &coef, const QuadratureRule *qr)
		: super_type(qr)
		, qfun(coef)
	{}

public:

	// LinearFormIntegrator ÇâÓÇµÇƒåpè≥Ç≥ÇÍÇ‹ÇµÇΩ
	virtual void assembleElementVector(
		const FiniteElement & el, Transformation & tr, VectorXd & vec) override;

};






////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


