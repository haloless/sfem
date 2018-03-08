#pragma once

#include "sfem.h"
#include "linalg.h"

BEGIN_SFEM_NS;

struct QuadraturePoint
{
public:

	VectorXd p;
	double w;

	////////////////////////////////////////////////////////////

	QuadraturePoint() : p(sfem::SpaceDim()), w(0) 
	{}

	const VectorXd& position() const { return p; }
	double weight() const { return w; }
	int dim() const { return p.size(); }

	void reset() {
		p.setZero();
		w = 0;
	}

	void setDim(int ndim) {
		p.resize(ndim);
		reset();
	}

	void set(const int ndim, const double *pos, double wgt) {
		p.resize(ndim);
		for (int i = 0; i < ndim; i++) {
			p(i) = pos[i];
		}
		w = wgt;
	}
	
	void set1d(double x, double wgt) {
		const int nd = 1;
		std::array<double, nd> pos = { x };
		set(nd, pos.data(), wgt);
	}
	void set2d(double x, double y, double wgt) {
		const int nd = 2;
		std::array<double, nd> pos = { x, y };
		set(nd, pos.data(), wgt);
	}
	void set3d(double x, double y, double z, double wgt) {
		const int nd = 3;
		std::array<double, nd> pos = { x, y, z };
		set(nd, pos.data(), wgt);
	}

};

class QuadratureRule
{
protected:
	// quadrature points
	Array<QuadraturePoint> qpts;
	// order of quadrature
	int qorder;
	// dimension of quadrature
	int qdim;

public:
	QuadratureRule() : qpts(), qorder(0), qdim(0) {}

	QuadratureRule(int ngp, int order=-1)
		: qpts(ngp), qorder(order)
	{
		for (int i = 0; i < qpts.size(); i++) {
			qpts[i].reset();
		}
	}

	~QuadratureRule() {}

	////////////////////////////////////////////////////////////

	int getOrder() const { return qorder; }
	void setOrder(int order) { qorder = order; }

	int numPoints() const { return qpts.size(); }
	
	void setNumPoints(int ngp) { qpts.resize(ngp); }

	const Array<QuadraturePoint>& points() const { return qpts; }

	QuadraturePoint& getPoint(int i) { return qpts[i]; }
	const QuadraturePoint& getPoint(int i) const { return qpts[i]; }

	double calcWeightSum() const {
		double wsum = 0;
		for (QuadraturePoint const &qp : qpts) {
			wsum += qp.w;
		}
		return wsum;
	}

	// Get rule[i]
	QuadraturePoint& operator[](int i) { return qpts[i]; }
	const QuadraturePoint& operator[](int i) const { return qpts[i]; }

public:

	void printInfo(std::ostream &os) const;

public:

	// Global function to get proper Rule on Geometry with Order
	static const QuadratureRule* GetRule(int geom, int order);
	
	// 
	static const QuadratureRule* GetProperRule(int geom, int order);

};






END_SFEM_NS;



