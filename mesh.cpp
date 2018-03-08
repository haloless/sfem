#include "stdafx.h"

#include "mesh.h"

#include <set>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;

Mesh::Mesh()
{
	nVert = 0;

	// TODO
}

Mesh::~Mesh()
{

	// TODO
}

////////////////////////////////////////////////////////////////////////////////

int Mesh::countVertices()
{
	// use a set to count unique vertices
	std::set<int> vert_set;

	for (const Element *e : elements) {
		int n = e->numNodes();
		const int *inds = e->getIndices();
		vert_set.insert(inds, inds + n);
	}

	// cache vertex number
	nVert = vert_set.size();
	return nVert;
}

void Mesh::getElementNodes(int i, MatrixXd & pos) const
{
	getElementNodes(*elements[i], pos);
}

void Mesh::getElementNodes(const Element &e, MatrixXd & pos) const
{
	int nn = e.numNodes();
	const int *inds = e.getIndices();

	int nd = nodes[0].size();

	pos.resize(nd, nn);

	for (int k = 0; k < nn; k++) {
		pos.col(k) = nodes[inds[k]];
	}

}


END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


