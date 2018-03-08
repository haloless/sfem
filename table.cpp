#include "stdafx.h"

#include "table.h"


////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

Table::Table(int nrow, int ncol)
{
	const int size = nrow;
	const int len = nrow * ncol;

	// NOTE the size+1
	iarr.resize(size + 1);
	//
	jarr.resize(len);

	// 
	for (int i = 0; i <= size; i++) {
		iarr[i] = ncol * i;
	}

	//
	std::fill(jarr.begin(), jarr.end(), -1);

}

////////////////////////////////////////////////////////////////////////////////

void Table::setSize(int nrow, int ncol)
{
	// 
	setCapacity(nrow, nrow*ncol);

	if (hasRows()) {
		// 
		for (int i = 0; i <= nrow; i++) {
			iarr[i] = ncol * i;
		}

		//
		std::fill(jarr.begin(), jarr.end(), -1);
	}
}

void Table::setCapacity(int nrow, int nnz)
{
	iarr.resize(nrow + 1);
	jarr.resize(nnz);

	// set I-array bounds
	if (nrow >= 0) {
		iarr[0] = 0;
		iarr[nrow] = nnz;
	}
}

int Table::getCapacity() const
{
	int size = numRows();
	return iarr[size];
}

void Table::makeI(int nrows)
{
	setCapacity(nrows, 0);

	for (int i = 0; i <= nrows; i++) {
		iarr[i] = 0;
	}
}

void Table::makeJ()
{
	int k = 0;

	const int size = numRows();
	
	for (int i = 0; i < size; i++) {
		int j = iarr[i];
		iarr[i] = k;
		k += j;
	}

	iarr[size] = k;

	jarr.resize(k);
}

////////////////////////////////////////////////////////////////////////////////



void Table::getRow(int i, Array<int>& row) const {
	int n = rowSize(i);
	const int *data = getRow(i);

	row.assign(data, data + n);
}

void Table::sortRow(int i)
{
	int *data = jarr.data();
	std::sort(data + iarr[i], data + iarr[i + 1]);
}

void Table::sortRows()
{
	const int nrow = numRows();
	for (int irow = 0; irow < nrow; irow++) {
		sortRow(irow);
	}
}

int Table::countNonZeros() const
{
	int nnz = 0;

	const int len = getCapacity();
	for (int k = 0; k < len; k++) {
		if (jarr[k] != -1) {
			nnz += 1;
		}
	}

	return nnz;
}

int Table::operator()(int i, int j) const
{
	if (i < 0 || i >= numRows()) return -1;

	// search in the row
	const int ibegin = iarr[i];
	const int iend = iarr[i + 1];
	for (int k = ibegin; k < iend; k++) {
		if (jarr[k] == j) {
			return k;
		}
		else if (jarr[k] == -1) {
			return - 1;
		}
	}

	// nothing found
	return -1;
}


////////////////////////////////////////////////////////////////////////////////

void Table::alloc(int nrows, int nnz)
{
	iarr.resize(nrows + 1);
	jarr.resize(nnz);
}

void Table::clear()
{
	iarr.clear();
	jarr.clear();
}

void Table::pack()
{
	const int size = numRows();
	const int nnz = countNonZeros();

	if (nnz != getCapacity()) {
		// non-zero elements do not match buffer capacity
		// need packing of data

		Array<int> jnew((size_t) nnz, -1);

		int icnt = 0;
		int idx = 0;
		for (int i = 0; i < size; i++) {
			for (int j = iarr[i]; j < iarr[i + 1]; j++) {
				if (jarr[j] == -1) break;
				// push to new J-array
				jnew[idx] = jarr[j];
				idx += 1;
			}

			iarr[i] = icnt;
			icnt = idx;
		}
		iarr[size] = nnz;

		jarr.swap(jnew);
	}
}

void Table::buildFromPairs(int nrows, const Array<IntPair>& pairs)
{
	clear();

	const int size = nrows;
	const int nnz = pairs.size();

	alloc(size, nnz);

	int cnt = 0;
	for (int i = 0; i <= size; i++) {
		iarr[i] = cnt;

		while (cnt < nnz && pairs[cnt].i == i) {
			jarr[cnt] = pairs[cnt].j;
			cnt += 1;
		}
	}

}


////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

