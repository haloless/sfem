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

class NonlinearFormIntegrator
{
protected:
	const QuadratureRule *qrule;

	////////////////////////////////////////////////////////////
public:
	NonlinearFormIntegrator(const QuadratureRule *qr = nullptr)
		: qrule(qr)
	{
		// nothing
	}

	virtual ~NonlinearFormIntegrator() {
		// nothing
	}

	////////////////////////////////////////////////////////////
public:

	void setQuadratureRule(const QuadratureRule *qr) { qrule = qr; }

	const QuadratureRule* getQuadratureRule() const { return qrule; }

	////////////////////////////////////////////////////////////
public:

	virtual void assembleElementVector(
		const FiniteElement &elem, 
		Transformation &trans,
		const VectorXd &v, 
		VectorXd &elvec) = 0;



};


////////////////////////////////////////////////////////////////////////////////

//
// Abstract class for bilinear form integrators.
//
class BilinearFormIntegrator : public NonlinearFormIntegrator
{
protected:

	typedef NonlinearFormIntegrator super_type;

	////////////////////////////////////////////////////////////
public:

	BilinearFormIntegrator(const QuadratureRule *qr = nullptr)
		: super_type(qr)
	{}

	virtual ~BilinearFormIntegrator()
	{}

	////////////////////////////////////////////////////////////
public:

	// Assemble element matrix
	// (not implemented)
	virtual void assembleElementMatrix(
		const FiniteElement &elem,
		Transformation &trans,
		MatrixXd &elmat);


	// Assemble element vector
	// (not implemented)
	// (extended from the super class)
	virtual void assembleElementVector(
		const FiniteElement &elem,
		Transformation &trans,
		const VectorXd &v,
		VectorXd &elvec);


	////////////////////////////////////////////////////////////
public:

	virtual const QuadratureRule* ensureQuadratureRule() const;

};



////////////////////////////////////////////////////////////////////////////////


//
// Integrator for Laplacian problem: 
// - Lap(u) = f
// The corresponding bilinear form is:
// a(u,v) = (grad u, grad v)
// 
class LaplacianIntegrator : public BilinearFormIntegrator
{
protected:

	////////////////////////////////////////////////////////////
public:

	LaplacianIntegrator() 
	{}

	////////////////////////////////////////////////////////////
public:

	// Assemble element matrix
	virtual void assembleElementMatrix(
		const FiniteElement &elem,
		Transformation &trans,
		MatrixXd &elmat) override;

	// Assemble element vector
	virtual void assembleElementVector(
		const FiniteElement &elem,
		Transformation &trans,
		const VectorXd &v,
		VectorXd &elvec) override;

};



////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

