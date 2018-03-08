#pragma once

#include "sfem.h"
#include "geometry.h"
#include "fe.h"
#include "fedesc.h"
#include "mesh.h"
#include "table.h"

#include <memory>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


class DofOrder
{
protected:
	int ndof;
	int ndim;

public:

	// empty
	DofOrder() : ndof(0), ndim(0) {}

	// given
	DofOrder(int ndof_, int ndim_)
		: ndof(ndof_), ndim(ndim_)
	{}

	////////////////////////////////////////////////////////////

	//
	int numDof() const { return ndof; }
	void setNumDof(int n) { ndof = n; }

	//
	int vecDim() const { return ndim; }
	void setVecDim(int v) { ndim = v; }

	// map 
	int map(int idof, int idim) const { return Map(ndof, ndim, idof, idim); }

	// map
	int operator()(int idof, int idim) const { return map(idof, idim); }

	//
	void mapDofs(const Array<int> &dofs, Array<int> &inds) const {
		MapDofs(ndof, ndim, dofs, inds);
	}


	////////////////////////////////////////////////////////////
public:
	// map (idof,idim) to index
	static int Map(int ndofs, int vdim, int idof, int idim);

	// 
	static void MapDofs(int ndofs, int vdim, const Array<int> &dofs, Array<int> &inds);

	//
	static void MapDofs(int ndofs, int vdim, const Array<int> &dofs, int idim, Array<int> &inds);
};


////////////////////////////////////////////////////////////////////////////////


// 
// Manage DOFs
// 
class FiniteElementSpace
{
protected:

	Mesh *mesh;

	const FiniteElementDescription *desc;

	// number of components
	const int ncomp;

	// number of DOFs
	// total number of unknowns is thus NDOFS*NCOMP
	int ndofs;

	//
	std::unique_ptr<Table> elem_dof;

	////////////////////////////////////////////////////////////
public:

	FiniteElementSpace(Mesh *m, const FiniteElementDescription *fedesc, int ncomp);

	virtual ~FiniteElementSpace();

	////////////////////////////////////////////////////////////
protected:

	void buildElementDofTable();


};








////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


