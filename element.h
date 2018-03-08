#pragma once

#include "sfem.h"
#include "geometry.h"
#include "fe.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;

enum ElementType
{
	Seg, Tri, Quad,

	Seg2, Tri2, Quad2,
};

////////////////////////////////////////////////////////////////////////////////

class Element
{
protected:
	// 
	int attribute;
	// BasicGeomType
	int baseGeom;

public:
	Element(int base=BasicGeomType::Undefined)
		: baseGeom(base), attribute(0)
	{}

	virtual ~Element() {}

public:

	// return type identifier
	virtual int type() const = 0;

	// nodes are the real points
	virtual int numNodes() const = 0;

	// vertices are only for corners
	virtual int numVertices() const = 0;

	// get internal buffer
	virtual int* getIndices() = 0;

	virtual const int* getIndices() const {
		return const_cast<Element*>(this)->getIndices();
	}

	//
	virtual void setIndices(const int *inds) {
		int n = numNodes();
		int *buf = getIndices();
		std::copy(inds, inds + n, buf);
	}


public:

	int getAttribute() const { return attribute; }
	void setAttribute(int a) { attribute = a; }

	int getBaseGeom() const { return baseGeom; }


	void getVertices(Array<int> &inds) const {
		int nv = numVertices();
		const int *buf = getIndices();
		inds.assign(buf, buf + nv);
	}

	void getNodes(Array<int> &inds) const {
		int nn = numNodes();
		const int *buf = getIndices();
		inds.assign(buf, buf + nn);
	}

	Array<int> getVertices() const {
		Array<int> inds;
		getVertices(inds);
		return inds;
	}

	Array<int> getNodes() const {
		Array<int> inds;
		getNodes(inds);
		return inds;
	}

	int operator[](int inode) const {
		return getIndices()[inode];
	}

};


////////////////////////////////////////////////////////////////////////////////

class SegmentElement : public Element
{
protected:
	int indices[2];
public:
	SegmentElement(int i0, int i1)
		: Element(BasicGeomType::Line)
	{
		indices[0] = i0;
		indices[1] = i1;
	}
	SegmentElement(const int *ind)
		: Element(BasicGeomType::Line)
	{
		indices[0] = ind[0];
		indices[1] = ind[1];
	}

	////////////////////////////////////////////////////////////
public:

	virtual int type() const override { return ElementType::Seg; }

	virtual int numNodes() const override { return 2; }

	virtual int numVertices() const override { return 2; }

	virtual int * getIndices() override { return indices; }

public:

	void setIndices(int i0, int i1) {
		indices[0] = i0;
		indices[1] = i1;
	}
};

////////////////////////////////////////////////////////////////////////////////

class TriangleElement : public Element
{
protected:
	int indices[3];

public:
	TriangleElement(int i0, int i1, int i2)
		: Element(BasicGeomType::Triangle)
	{
		indices[0] = i0;
		indices[1] = i1;
		indices[2] = i2;
	}
	TriangleElement(const int *ind)
		: Element(BasicGeomType::Triangle)
	{
		indices[0] = ind[0];
		indices[1] = ind[1];
		indices[2] = ind[2];
	}

	////////////////////////////////////////////////////////////
public:

	virtual int type() const override { return ElementType::Tri; }

	virtual int numNodes() const override { return 3; }

	virtual int numVertices() const override { return 3; }

	virtual int * getIndices() override { return indices; }

	////////////////////////////////////////////////////////////
public:
	void setIndices(int i0, int i1, int i2) {
		indices[0] = i0;
		indices[1] = i1;
		indices[2] = i2;
	}
};

////////////////////////////////////////////////////////////////////////////////

// Quadratic line
// Node config: i0 - i2 - i1
class QuadraticSegmentElement : public Element
{
protected:
	int indices[3];
public:
	QuadraticSegmentElement(int i0, int i1, int i2)
		: Element(BasicGeomType::Line)
	{
		indices[0] = i0;
		indices[1] = i1;
		indices[2] = i2;
	}
	QuadraticSegmentElement(const int *ind)
		: Element(BasicGeomType::Line)
	{
		indices[0] = ind[0];
		indices[1] = ind[1];
		indices[2] = ind[2];
	}

	////////////////////////////////////////////////////////////
public:

	virtual int type() const override { return ElementType::Seg2; }

	virtual int numNodes() const override { return 3; }

	virtual int numVertices() const override { return 2; }

	virtual int * getIndices() override { return indices; }

	////////////////////////////////////////////////////////////
public:
	void setIndices(int i0, int i1, int i2) {
		indices[0] = i0;
		indices[1] = i1;
		indices[2] = i2;
	}
};

////////////////////////////////////////////////////////////////////////////////

class QuadraticTriangleElement : public Element
{
protected:
	int indices[6];
public:
	QuadraticTriangleElement(int i0, int i1, int i2, int i3, int i4, int i5)
		: Element(BasicGeomType::Triangle)
	{
		indices[0] = i0;
		indices[1] = i1;
		indices[2] = i2;
		indices[3] = i3;
		indices[4] = i4;
		indices[5] = i5;
	}
	QuadraticTriangleElement(const int *ind)
		: Element(BasicGeomType::Triangle)
	{
		for (int i = 0; i < 6; i++) {
			indices[i] = ind[i];
		}
	}

	////////////////////////////////////////////////////////////
public:

	virtual int type() const override { return ElementType::Tri2; }

	virtual int numNodes() const override { return 6; }

	virtual int numVertices() const override { return 3; }

	virtual int * getIndices() override { return indices; }
};

////////////////////////////////////////////////////////////////////////////////
// 
// The corresponding finite elements
//

// linear & bilinear
extern Linear1DFiniteElement TheSegmentFE;
extern Linear2DFiniteElement TheTriangleFE;
extern Bilinear2DFiniteElement TheQuadrilateralFE;

// quadratic 
extern Quadratic2DFiniteElement TheQuadraticTriangleFE;






////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////




