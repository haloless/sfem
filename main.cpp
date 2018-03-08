#include "stdafx.h"

#include "sfem.h"
#include "linalg.h"
#include "mesh.h"
#include "quadrature.h"
#include "fe.h"
#include "fespace.h"
#include "transform.h"
#include "linearform.h"
#include "bilinearform.h"

#include <cmath>
#include <iostream>
#include <functional>

#include <cstdarg>

static sfem::FiniteElement* selectFE(int elemType) {

	sfem::FiniteElement *fe = nullptr;

	if (elemType == sfem::ElementType::Tri) {
		fe = &sfem::TheTriangleFE;
		std::cout << "Linear triangle" << std::endl;
	}
	else if (elemType == sfem::ElementType::Tri2) {
		fe = &sfem::TheQuadraticTriangleFE;
		std::cout << "Quadratic triangle" << std::endl;
	}
	else {
		std::cerr << "Unsupported element type = " << elemType << std::endl;
		exit(1);
	}

	return fe;
}



static int main_test(int argc, char *argv[]) {
	if (0) {
		sfem::VectorXd a(2);
		sfem::VectorXd b(3);
		b.setRandom();

		a = b;
		std::cout << a << std::endl;
		std::cout << a.size() << std::endl;

		std::cout << std::endl;
	}

	if (0) {
		sfem::IsoparametricTransformation trans;

		// set to triangle
		trans.setSameAs(sfem::BasicGeomType::Triangle);

		//
		trans.getNodesMat() << 
			1, 2, 3,
			1, 0, 2;

		// for some quadrature rule
		const auto *rule = sfem::QuadratureRule::GetRule(sfem::BasicGeomType::Triangle, 5);
		assert(rule);

		for (const auto &qpoint : rule->points()) {
			std::cout << "Quadrature point = " << qpoint.p.transpose() << std::endl;
			trans.setPoint(qpoint);

			std::cout << "J = \n" << trans.Jacobian() << std::endl;
		}

		std::cout << trans.Jacobian().determinant() << std::endl;
		

		//

	}

	if (0) {
		int order = 1;
		int dim = 2;
		sfem::FEDescH1 desc(order, dim);

		std::cout << desc.getName() << std::endl;
	}

	if (1) {
		double a = 2.0;
		double b = 3.0;
		int cnt = 0;

		{
			auto coef = sfem::MakeTFuncCoef([&](const sfem::VectorXd &v) -> double {
				cnt += 1;
				double x = v(0);
				double y = v(1);
				double val = a*x + b*y;
				return val;
			});

			sfem::VectorXd x(2); x << 1, 2;
			std::cout << coef(x) << "<>" << (a*x(0) + b*x(1)) << std::endl;
			std::cout << cnt << std::endl;
		}
		{
			//std::function
			std::cout << std::hex << 15 << std::endl;
		}
	}

	return 0;
}



static int main_sfem(int argc, char *argv[])
{
	std::cout << "hoge" << std::endl;

	if (argc != 2) {
		std::cout << "Usage: sfem.exe <mesh>" << std::endl;
		return 1;
	}

	const char *meshfilename = argv[1];

	// 
	const int sdim = 2;
	sfem::SetSpaceDim(sdim);
	sfem::SetCoordSys(sfem::CoordType::Cartesian);

	//
	sfem::Mesh mesh;
	if (0) {
		mesh.readGmsh(meshfilename);
	} else {
		const int ncell[] = { 16, 16 };
		const double len[] = { 1.0, 1.0 };
		const double xlo[] = { 0.0, 0.0 };
		mesh.build2DSquare(ncell, len, xlo);
	}

	mesh.countVertices();
	mesh.printInfo(std::cout);

	// mesh.writeVtk("mesh.vtk");

	//sfem::LinearFEDescription fedesc;
	//sfem::FiniteElementDescription *fedesc = nullptr;
	//fedesc = new sfem::LinearFEDescription();

	//sfem::FiniteElementSpace *fespace = nullptr;
	//fespace = new sfem::FiniteElementSpace();
	//sfem::FiniteElementSpace fespace(&mesh, &fedesc, 1);


	//
	{
		const int ndof = mesh.numNodes();
		std::cout << "ndof=" << ndof << std::endl;

		sfem::VectorXd bvec(ndof);
		bvec.setZero();

		//sfem::MatrixXd amat(ndof, ndof);
		sfem::SparseMatrix amat(ndof, ndof);
		amat.setZero();

		auto fun_u = [](const sfem::VectorXd &p) -> double {
			double x = p(0);
			double y = p(1);
			double val = x*x + 2 * y*y + 1;
			//double val = 2 * x + y + 1;
			return val;
		};
		auto fun_f = [](const sfem::VectorXd &p) -> double {
			double val = -6.0;
			//double val = 0.0;
			return val;
		};

		double area = 0;
		double volume = 0;

		for (int i = 0; i < mesh.numElements(); i++) {
			const sfem::Element *elem = &mesh.getElement(i);
			//assert(elem->type() == sfem::ElementType::Tri);
			assert(elem->getBaseGeom() == sfem::BasicGeomType::Triangle);

			// select finite element
			const sfem::FiniteElement *fe = selectFE(elem->type());
			assert(fe);

			// collect node positions
			sfem::MatrixXd nodePos;
			mesh.getElementNodes(*elem, nodePos);

			// set isoparam transform
			sfem::IsoparametricTransformation trans;
			trans.setFE(fe);
			trans.getNodesMat() = nodePos;

			//
			sfem::VectorXd be(fe->numDof());
			be.setZero();
			sfem::MatrixXd ae(fe->numDof(), fe->numDof());
			ae.setZero();

			// 
			const int qorder = 5;
			const auto *qrule = sfem::QuadratureRule::GetRule(elem->getBaseGeom(), qorder);
			assert(qrule);

			for (const sfem::QuadraturePoint &qpt : qrule->points()) {

				trans.setPoint(qpt);

				sfem::VectorXd xpos;
				trans.transform(qpt, xpos);
				double r = xpos(0);

				sfem::MatrixXd Jmat = trans.Jacobian();
				double detJ = Jmat.determinant();

				area += detJ * qpt.w;
				volume += r * 2.0 * M_PI * detJ * qpt.w;

				sfem::VectorXd shape;
				fe->calcShape(qpt, shape);
				
				sfem::MatrixXd dshape;
				fe->calcDeriv(qpt, dshape);

				sfem::MatrixXd der = Jmat.inverse().transpose() * dshape.transpose();
				//std::cout << der << std::endl;
				ae += der.transpose() * der * detJ * qpt.w;

				double fval = fun_f(xpos);
				be += fval * detJ * qpt.w * shape;
			}

			for (int i = 0; i < elem->numNodes(); i++) {
				int ii = elem->getIndices()[i];
				bvec(ii) += be(i);
				
				for (int j = 0; j < elem->numNodes(); j++) {
					int jj = elem->getIndices()[j];
					//amat(ii, jj) += ae(i, j);
					amat.coeffRef(ii, jj) += ae(i, j);
				}
			}
		}

		std::cout << "area=" << area << std::endl;
		std::cout << "volume=" << volume << std::endl;

		std::cout << "BC" << std::endl;
		for (const sfem::Element* e : mesh.getBoundary()) {
			for (int i : e->getNodes()) {
				sfem::VectorXd xpos = mesh.getNodes()[i];
				double val = fun_u(xpos);
				
				// penalty method
				const double s = 1.0e9;
				bvec(i) += s * val;
				//amat(i, i) += s;
				amat.coeffRef(i, i) += s;

				//bvec(i) = val;
				//amat.row(i).setZero();
				//amat(i, i) = 1.0;
			}
		}

		std::cout << "solve" << std::endl;
		amat.makeCompressed();
		//sfem::VectorXd uvec = amat.householderQr().solve(bvec);
		Eigen::SparseLU<sfem::SparseMatrix> solver;
		solver.analyzePattern(amat);
		solver.factorize(amat);
		sfem::VectorXd uvec = solver.solve(bvec);

		// 
		sfem::VisDataSet vis("hoge", mesh);
		vis.addNodalScalarField("b", bvec);
		vis.addNodalScalarField("u", uvec);
		vis.writeVtk();
	}

	{
		const sfem::QuadratureRule *rule = nullptr;
		rule = sfem::QuadratureRule::GetRule(sfem::BasicGeomType::Triangle, 5);
		if (rule == nullptr) {
			std::cout << "no proper rule" << std::endl;
		}
		else {
			rule->printInfo(std::cout);
		}
	}


	return 0;
}


static int main_sfem2(int argc, char *argv[])
{

	// 
	const int sdim = 2;
	sfem::SetSpaceDim(sdim);
	sfem::SetCoordSys(sfem::CoordType::Cartesian);

	//
	sfem::Mesh mesh;
	if (argc == 2) {
		const char *meshfilename = argv[1];
		mesh.readGmsh(meshfilename);
	}
	else {
		const int ncell[] = { 16, 16 };
		const double len[] = { 1.0, 1.0 };
		const double xlo[] = { 0.0, 0.0 };
		mesh.build2DSquare(ncell, len, xlo);
	}

	mesh.countVertices();
	mesh.printInfo(std::cout);

	//
	{
		const int ndof = mesh.numNodes();
		std::cout << "ndof=" << ndof << std::endl;

		sfem::VectorXd bvec(ndof);
		bvec.setZero();

		//sfem::MatrixXd amat(ndof, ndof);
		sfem::SparseMatrix amat(ndof, ndof);
		amat.setZero();

		auto fun_u = sfem::MakeTFuncCoef([](const sfem::VectorXd &p) -> double {
			double x = p(0);
			double y = p(1);
			double val = x*x + 2 * y*y + 1;
			//double val = 2 * x + y + 1;
			return val;
		});

		auto fun_f = sfem::MakeTFuncCoef([](const sfem::VectorXd &p) -> double {
			double val = -6.0;
			//double val = 0.0;
			return val;
		});


		sfem::DomainLinearFormIntegrator binteg(fun_f, nullptr);
		sfem::LaplacianIntegrator ainteg;

		for (int i = 0; i < mesh.numElements(); i++) {
			const sfem::Element *elem = &mesh.getElement(i);
			assert(elem->getBaseGeom() == sfem::BasicGeomType::Triangle);

			// select finite element
			const sfem::FiniteElement *fe = selectFE(elem->type());
			assert(fe);

			// quadrature
			const int qorder = 5;
			const auto *qrule = sfem::QuadratureRule::GetRule(elem->getBaseGeom(), qorder);
			assert(qrule);

			// set isoparam transform
			sfem::IsoparametricTransformation trans;
			trans.setFE(fe);
			mesh.getElementNodes(*elem, trans.getNodesMat());

			// element vector
			sfem::VectorXd be(fe->numDof());
			be.setZero();
			// element matrix
			sfem::MatrixXd ae(fe->numDof(), fe->numDof());
			ae.setZero();

			//
			ainteg.setQuadratureRule(qrule);
			binteg.setQuadratureRule(qrule);

			//
			ainteg.assembleElementMatrix(*fe, trans, ae);

			//
			binteg.assembleElementVector(*fe, trans, be);

			// TODO global assemble
			for (int i = 0; i < elem->numNodes(); i++) {
				int ii = elem->getIndices()[i];
				bvec(ii) += be(i);

				for (int j = 0; j < elem->numNodes(); j++) {
					int jj = elem->getIndices()[j];
					amat.coeffRef(ii, jj) += ae(i, j);
				}
			}
		}

		std::cout << "BC" << std::endl;
		for (const sfem::Element* e : mesh.getBoundary()) {
			for (int i : e->getNodes()) {
				sfem::VectorXd xpos = mesh.getNodes()[i];
				double val = fun_u(xpos);

				// penalty method
				const double s = 1.0e9;
				bvec(i) += s * val;
				amat.coeffRef(i, i) += s;
			}
		}

		std::cout << "solve" << std::endl;
		amat.makeCompressed();
		//sfem::VectorXd uvec = amat.householderQr().solve(bvec);
		Eigen::SparseLU<sfem::SparseMatrix> solver;
		solver.analyzePattern(amat);
		solver.factorize(amat);
		sfem::VectorXd uvec = solver.solve(bvec);

		// 
		sfem::VisDataSet vis("hoge", mesh);
		vis.addNodalScalarField("b", bvec);
		vis.addNodalScalarField("u", uvec);
		vis.writeVtk();
	}

	{
		const sfem::QuadratureRule *rule = nullptr;
		rule = sfem::QuadratureRule::GetRule(sfem::BasicGeomType::Triangle, 5);
		if (rule == nullptr) {
			std::cout << "no proper rule" << std::endl;
		}
		else {
			rule->printInfo(std::cout);
		}
	}


	return 0;
}


static int myprintf(const char fmt[], ...) {
	va_list ap;
	va_start(ap, fmt);
	int ret = vfprintf(stdout, fmt, ap);
	va_end(ap);
	return ret;
}


int main(int argc, char *argv[])
{
	//main_sfem(argc, argv);
	//main_sfem2(argc, argv);

	//main_test(argc, argv);

	
	myprintf("argc = %d %d %d\n", argc, argc, argc);

	return 0;
}



