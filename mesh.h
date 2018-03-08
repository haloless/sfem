#pragma once

#include "sfem.h"
#include "linalg.h"
#include "element.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


class Mesh
{
protected:

	int nVert;
	//int nNode;
	//int nElem;

	// 
	Array<VectorXd> nodes;
	
	//
	Array<Element*> elements;
	
	//
	Array<Element*> boundary;


	////////////////////////////////////////////////////////////
public:
	Mesh();

	virtual ~Mesh();

	////////////////////////////////////////////////////////////
public:

	int numNodes() const { return nodes.size(); }
	int numElements() const { return elements.size(); }
	int numBoundaryElements() const { return boundary.size(); }

	Array<VectorXd>& getNodes() { return nodes; }
	Array<Element*>& getElements() { return elements; }
	Array<Element*>& getBoundary() { return boundary; }

	const Element& getElement(int i) const { return *elements[i]; }
	Element& getElement(int i) { return *elements[i]; }

	// CAUTION: vertices need count!!!
	int numVertices() const { return nVert; }
	int countVertices();

	// Get node positions (NDIM * NPOINT)
	void getElementNodes(int i, MatrixXd &pos) const;
	// Get node positions (NDIM * NPOINT)
	void getElementNodes(const Element &e, MatrixXd &pos) const;

	////////////////////////////////////////////////////////////
	// Mesh IO functions
public:
	 
	void printInfo(std::ostream &os);

	void build2DSquare(const int ncell[], const double len[], const double xlo[]);

	// Load GMSH-2 file
	int readGmsh(const char *filename);

	int writeVtk(const char *filename);


	////////////////////////////////////////////////////////////
private:

};

////////////////////////////////////////////////////////////////////////////////

class VisDataSet
{
protected:
	std::string title;

	Mesh *mesh;

	std::map<std::string, const VectorXd*> dataset;

public:
	VisDataSet(const std::string &t, Mesh &m)
		: title(t), mesh(&m)
	{}

	virtual ~VisDataSet() {}

public:

	VisDataSet& addNodalScalarField(const std::string& name, const VectorXd &data) {
		dataset[name] = &data;
		return *this;
	}

	void writeVtk();

};

////////////////////////////////////////////////////////////////////////////////



END_SFEM_NS;


