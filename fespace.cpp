#include "stdafx.h"

#include "fespace.h"

#include <iostream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


int DofOrder::Map(int ndofs, int vdim, int idof, int idim)
{
	// order as (dim,dof)
	int ind = idim + vdim * idof;
	return ind;
}

void DofOrder::MapDofs(int ndofs, int vdim, const Array<int>& dofs, Array<int>& inds)
{
	const int size = dofs.size();
	inds.resize(size*vdim);

	for (int idim = 0; idim < vdim; idim++) {
		for (int i = 0; i < size; i++) {
			inds[i + idim*size] = Map(ndofs, vdim, dofs[i], idim);
		}
	}
}

void DofOrder::MapDofs(int ndofs, int vdim, const Array<int>& dofs, int idim, Array<int>& inds)
{
	const int size = dofs.size();
	inds.resize(size);

	for (int i = 0; i < size; i++) {
		inds[i] = Map(ndofs, vdim, dofs[i], idim);
	}
}



////////////////////////////////////////////////////////////////////////////////


FiniteElementSpace::FiniteElementSpace(
	Mesh * m, const FiniteElementDescription * fedesc, int nc)
	: mesh(m), desc(fedesc), ncomp(nc)
{
	// TODO
}

FiniteElementSpace::~FiniteElementSpace()
{
	// TODO
}

////////////////////////////////////////////////////////////////////////////////

void FiniteElementSpace::buildElementDofTable()
{
	if (elem_dof) return;

	Table *tbl = new Table;

	const int nelem = mesh->numElements();

	for (int i = 0; i < nelem; i++) {

	}

	// TODO

	elem_dof.reset(tbl);
}


////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

