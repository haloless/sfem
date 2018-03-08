#pragma once

#include "sfem.h"
#include "linalg.h"
#include "transform.h"
#include "quadrature.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


class ScalarCoefficient
{
protected:

	double time;

public:

	ScalarCoefficient() : time(0) {}

	virtual ~ScalarCoefficient() {}

public:

	void setTime(double t) { time = t; }
	double getTime() const { return time; }

public:

	// Evaluate
	virtual double eval(
		Transformation &trans,
		const QuadraturePoint &qp) = 0;

};

////////////////////////////////////////////////////////////////////////////////

class ConstantCoefficient : public ScalarCoefficient
{
protected:
	double value;

public:
	explicit ConstantCoefficient(double c = 0) : value(c)
	{
		// nothing
	}

	double getValue() const { return value; }

	void setValue(double c) { value = c; }


	// Evaluation returns constant value
	virtual double eval(Transformation &trans, const QuadraturePoint &qp) override {
		return value;
	}
};

////////////////////////////////////////////////////////////////////////////////

class PiecewiseConstantCoefficient : public ScalarCoefficient
{
protected:
	VectorXd values;

public:
	explicit PiecewiseConstantCoefficient(int num = 0) {
		if (num > 0) values.resize(num);
		values.setZero();
	}

	//
	VectorXd& getValues() { return values; }
	const VectorXd& getValues() const { return values; }

	// set all values to const
	void setValues(double c) { values.fill(c); }

	// Evaluation
	virtual double eval(Transformation &trans, const QuadraturePoint &qp) override {
		// TODO lookup values
		throw std::exception("not implemented: " __FUNCTION__);
		return 0;
	}
};


////////////////////////////////////////////////////////////////////////////////

//
// Templated coefficient that wraps function or functors
// 
template<typename T>
class TFunctionCoefficient : public ScalarCoefficient
{
protected:
	//typedef double (T)(const VectorXd &v);

	T fun;

public:

	TFunctionCoefficient(T t) : fun(t) {}

	// Directly call wrapped function
	double call(const VectorXd &x) {
		double val = fun(x);
		return val;
	}

	// Directly call wrapped function
	double operator()(const VectorXd &x) { return call(x); }

	// Evaluate function value
	virtual double eval(Transformation &trans, const QuadraturePoint &qp) override {

		VectorXd x;
		trans.transform(qp, x);

		double val = fun(x);

		return val;
	}
};

//
// Helper function 
// e.g., auto coef = MakeTFuncCoef(lambda)
//
template<typename T>
TFunctionCoefficient<T> MakeTFuncCoef(const T& fun)
{
	return TFunctionCoefficient<T>(fun);
}





////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////
