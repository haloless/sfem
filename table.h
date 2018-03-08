#pragma once

#include "sfem.h"
//#include "linalg.h"

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;

//
// A pair of integer indices
// Used to keep (i->j) connectivity
// 
struct IntPair
{
	int i, j;

	// default
	IntPair() : i(-1), j(-1) {}

	// 
	IntPair(int ii, int jj) : i(ii), j(jj) {}

	// equal test
	bool operator==(const IntPair &rhs) const {
		return (this->i == rhs.i && this->j == rhs.j);
	}

	// lesser test, allows sort
	bool operator<(const IntPair &rhs) const {
		if (this->i == rhs.i) {
			return this->j < rhs.j;
		}
		else {
			return this->i < rhs.i;
		}
	}

};


//
// CSR style connectivity table.
// 
//
class Table
{

protected:
	Array<int> iarr;
	Array<int> jarr;

	////////////////////////////////////////////////////////////
public:

	// empty table
	Table() {}

	// table with rows and allocate space for columns
	Table(int nrow, int ncol);

	// 
	virtual ~Table() {}

	////////////////////////////////////////////////////////////
public:

	void setSize(int nrow, int ncol);

	void setCapacity(int nrow, int nnz);

	// capacity, max number 
	int getCapacity() const;

	void makeI(int nrows);
	void makeJ();

	////////////////////////////////////////////////////////////
public:

	// number of rows
	int numRows() const { return iarr.size() - 1; }


	//
	bool hasRows() const { return numRows() > 0; }

	//
	int rowSize(int irow) const { return iarr[irow + 1] - iarr[irow]; }


	// get internal data buffer of a row
	const int* getRow(int i) const { return jarr.data() + iarr[i]; }
	int* getRow(int i) { return jarr.data() + iarr[i]; }

	// collect row data
	// the table must be packed before this
	void getRow(int i, Array<int> &row) const;

	// sort the indices in a row
	void sortRow(int i);
	// sort all rows
	void sortRows();

	// cound number of non-zero (NNZ)
	int countNonZeros() const;

	// get value at (irow,jcol)
	// If not found, return -1
	int operator()(int i, int j) const;

	////////////////////////////////////////////////////////////
public:

	// alloc array buffers only
	void alloc(int nrows, int nnz);

	// clear array buffers
	void clear();

	// pack up to mack compact in memory
	void pack();

	// Remember pairs must be sorted
	void buildFromPairs(int nrows, const Array<IntPair> &pairs);





};

class RowMajorTable
{
public:
	typedef Array<int> RowVec;

protected:

	Array<RowVec> tbl;

	////////////////////////////////////////////////////////////



};







END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////


