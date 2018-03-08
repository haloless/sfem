#include "stdafx.h"

#include "quadrature.h"
#include "element.h"

#include <memory>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


void QuadratureRule::printInfo(std::ostream & os) const
{
	const int np = numPoints();
	const int order = getOrder();
	os << "Rule: np=" << np << "; order=" << order << std::endl;

	os << "sum(w)=" << calcWeightSum() << std::endl;

	for (int i = 0; i < np; i++) {
		const QuadraturePoint &pt = getPoint(i);
		os << "p=(" << pt.p.transpose() << "); w=" << pt.w << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////

// Important: 
// Remember this is defined on [0,1], not [-1,+1]!!!
// Also make sure the weights sum to unity!!!
class QuadratureFunction1D
{
public:
	virtual ~QuadratureFunction1D() {}

	// calculate rule for N points
	virtual void calcRule(int np, QuadratureRule &qr) = 0;
};

//
// Gauss-Legendre on 1D
class QuadratureFunctionGaussLegendre1D : public QuadratureFunction1D
{
	// QuadratureFunction1D ‚ð‰î‚µ‚ÄŒp³‚³‚ê‚Ü‚µ‚½
	virtual void calcRule(int np, QuadratureRule &qr) override {
		assert(np > 0);

		qr.setNumPoints(np);

		if (np == 1) {
			qr[0].set1d(0.5, 1.0);
		}
		else if (np == 2) {
			// +-sqrt(3)/3
			qr[0].set1d(0.21132486540518711775, 0.5);
			qr[1].set1d(0.78867513459481288225, 0.5);
		}
		else if (np == 3) {
			// +-sqrt(15)/5
			qr[0].set1d(0.11270166537925831148, 5.0 / 18.0);
			qr[1].set1d(0.5, 4.0 / 9.0);
			qr[2].set1d(0.88729833462074168852, 5.0 / 18.0);
		}
		else {
			std::cerr << __FUNCTION__ << ": unsupported np=" << np << std::endl;
			exit(1);
		}
	}
};




////////////////////////////////////////////////////////////////////////////////


struct RuleStore
{

	// pre-calculated quadrature rules
	Array<std::unique_ptr<QuadratureRule> > rules_triangle;
	Array<std::unique_ptr<QuadratureRule> > rules_quadrilateral;



	RuleStore() {
		initTriangleRules();
	}

	~RuleStore() {
		//std::cout << "delete store" << std::endl;
	}

	//
	// Guassian Quadrature use N points for (2*N-1) poly order
	// Therefore, (2*N) even order must be clipped to (2*N+1) odd order
	//
	int clipOrder(int order) const {
		if (order % 2 == 0) {
			return order + 1;
		}
		else {
			return order;
		}
	}

	void initTriangleRules() {
		const int max_order = 16;

		rules_triangle.resize(max_order);

		const int ndim = 2; // 2D triangle

		{ // 5th-order, 7-points
			const int order = 5;
			const int ngp = 7;

			constexpr double w1 = 0.225000000000000;
			constexpr double w2 = 0.132394152788506;
			constexpr double w3 = 0.125939180544827;
			constexpr double a1 = 0.333333333333333;
			constexpr double a2 = 0.059715871789770;
			constexpr double a3 = 0.797426985353087;
			constexpr double b2 = (1 - a2) / 2;
			constexpr double b3 = (1 - a3) / 2;


			// Gauss points
			constexpr double gp[ngp][3] = {
				a1,a1,1 - a1 - a1,
				a2,b2,1 - a2 - b2,
				b2,a2,1 - b2 - a2,
				b2,b2,1 - b2 - b2,
				a3,b3,1 - a3 - b3,
				b3,a3,1 - b3 - a3,
				b3,b3,1 - b3 - b3,
			};
			// Gauss weights, note they sum to 1, not 0.5
			constexpr double w[ngp] = { w1,w2,w2,w2,w3,w3,w3 };
			
			QuadratureRule *rule = new QuadratureRule(ngp, order);
			rules_triangle[order].reset(rule);

			for (int i = 0; i < ngp; i++) {
				rule->getPoint(i).set(ndim, gp[i], w[i] * 0.5);
			}

		}
	}

	void initQuadrilateralRules() {
		const int max_order = 16;
		rules_quadrilateral.resize(max_order);

		const int ndim = 2; // 2D square

		// TODO
	}



};

// singleton
static RuleStore all_rules;




const QuadratureRule* QuadratureRule::GetRule(int geom, int order)
{
	if (geom == BasicGeomType::Triangle) {
		const QuadratureRule *rule = all_rules.rules_triangle[order].get();
		if (!rule) {
			std::cerr << "Unsupported rule for triangle order=" << order << std::endl;
			return nullptr;
		}
		else {
			assert(rule->getOrder() == order);
			return rule;
		}
	}
	else {
		std::cerr << "Unsupported quadrature rule for geom=" << geom << std::endl;
		return nullptr;
	}
}



END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


