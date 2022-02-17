//----------------------------------------------------------------------------
//
//   Group       : Risk Methodology
//
//   Filename    : SparseMatrix.cpp
//
//   Description : Implementation for a sparse matrix object
//
//   Author      : Anatoly V Morosov
//
//   Date        : August 2005
//
//
//----------------------------------------------------------------------------

/* ===== SPARSE MATRIX ==== */
#include "edginc/config.hpp"
#include "edginc/SparseMatrix.hpp"
#include "edginc/Format.hpp"
#include <map>
#include <set>
/*========================== */
using namespace std;
DRLIB_BEGIN_NAMESPACE


void SparseDoubleMatrix::allocateMem()
{
	if (!nelem)
	{	
		elem = 0; 
		pos = 0; 
		return; 
	}

    elem = new double[nelem];
    pos = new long[nelem];

    if ((elem == NULL) || (pos == NULL))
        throw ModelException("SparseDoubleMatrix::allocateMem","out of memory");
	// Note: tmp is allocated when the leftMultiply is called for the first time
	// Corollary: if leftMultiply is never called - tmp is not allocated
}

void SparseDoubleMatrix::deallocateMem()
{
    delete [] elem; elem = 0;
    delete [] pos; pos = 0;
	delete [] tmp; tmp = 0;
}


/********************** Sparse Matrix Constructor ***************/
/* 
This function creates the matrix from 3 arrays: 
M[row_k,col_k] = d_k

The consistency requirement (not checked explicitly)is 
that the input arrays contain no duplicates, i.e. any 
{row, col} pair appears only once.

*/
SparseDoubleMatrix::SparseDoubleMatrix (
    int _ncols, 
    int _nrows, 
    const IntArray& col_in,
    const IntArray& row_in,
    const DoubleArray& d_in)
    :
nrows((short)_nrows), ncols((short)_ncols), nelem(0), 
pos(0), elem(0), sorted(false), tmp(0)
{
    size_t  i,j;

	const static string module = "SparseDoubleMatrix::SparseDoubleMatrix(3 arrays)";

    if (col_in.size() != row_in.size() || col_in.size() != d_in.size())
		throw ModelException(module, "3 arrays have different sizes");

    /* check matrix is not too big or too small */
    if ((nrows > SHRT_MAX) || (nrows <= 0))
		throw ModelException(module, "Sparse Matrix cannot handle matrices longer than 32767 rows or less than 1");

	if (row_in.empty())
		return; // empty, but that's ok


    for (i=0; i<size_t(d_in.size()); ++i)
    {
        if ((row_in[i] >= nrows || row_in[i] < 0) ||
            (col_in[i] >= ncols || col_in[i] < 0))
        { 
            nelem = 0;
			throw ModelException(module, "The (col, row) = ( "+Format::toString(col_in[i])
				+" , "+Format::toString(col_in[i])+" ) is outside of the declared ncols/nrows");
        }

        if (d_in[i] != 0)   // non-trivial element
            ++nelem;
    }

    allocateMem();

    for (i=0, j=0; i<size_t(d_in.size()); ++i)
    {
        if (d_in[i] != 0) // non-trivial element
        {
            elem[j] = d_in[i];
            pos[j]  = getPosFromIJ((short)row_in[i], (short)col_in[i]);
            j++;
        }
    }

    /* one should want to call sorting here or afterwards if Mij is ever needed 
    or if we need to assure the uniqueness of every {row, col} pair       */

    return /*SUCCESS*/; 
}

/********************** Sparse Matrix Constructor ***************/
SparseDoubleMatrix::SparseDoubleMatrix() 
:
nrows(0), ncols(0), nelem(0), 
pos(0), elem(0), sorted(true), tmp(0)
{}


/**
This function creates a sparse matrix object from DoubleMatrix 
*/

SparseDoubleMatrix::SparseDoubleMatrix(const DoubleMatrix& data) 
:
nrows((short)data.numRows()), ncols((short)data.numCols()), nelem(0), 
pos(0), elem(0), sorted(true), tmp(0)
{
    int i, j, k;
	const static string module = "SparseDoubleMatrix::SparseDoubleMatrix(const DoubleMatrix&)";

    /* check matrix is not too big */
    if ((nrows > (SHRT_MAX)) || (nrows <= 0))
		throw ModelException(module, "Sparse Matrix cannot handle matrices longer than 32767 rows or cols");


    for (i=0; i<nrows; i++)
        for (j=0; j<ncols; j++)
            if (data[j][i] != 0) /* if it is non-diagonal non-zero element */
                ++nelem;

    allocateMem();

    for (i=0, k=0; i<nrows; ++i)
        for (j=0; j<ncols; ++j)
            if (data[j][i] != 0)
            {
                elem[k] = data[j][i];
                pos[k]  = getPosFromIJ((short)i,(short)j);
                ++k;
            }

    return /*SUCCESS*/;
}


/********************** Sparse Matrix Constructor ***************/
/* 
This function creates a sparse matrix object from (double**) 
matrix 

*/


SparseDoubleMatrix::SparseDoubleMatrix (
                                            int _nrows,
                                            int _ncols,
                                            const double ** M)
                                            :
nrows((short)_nrows), ncols((short)_ncols), nelem(0), 
pos(0), elem(0), sorted(true), tmp(0)
{
    int i, j, k;
	const static string module = "SparseDoubleMatrix::SparseDoubleMatrix(int,int,double**)";


    /* check matrix is not too big */
    if ((nrows > (SHRT_MAX)) || (nrows <= 0))
		throw ModelException(module, "Sparse Matrix cannot handle matrices longer than 32767 rows or cols or less than 1");

    if (M == NULL)
		throw ModelException(module, "called with NULL data pointer");


    for (i=0; i<nrows; i++)
        for (j=0; j<ncols; j++)
            if (M[i][j] != 0) /* if it is non-diagonal non-zero element */
                ++nelem;

    allocateMem();

    for (i=0, k=0; i<nrows; i++)
        for (j=0; j<ncols; j++)
            if (M[i][j] != 0)
            {
                elem[k] = M[i][j];
                pos[k]  = getPosFromIJ((short)i,(short)j);
                k++;
            }

    return /*SUCCESS*/;
}


SparseDoubleMatrix& SparseDoubleMatrix::operator=(SparseCollection & sc)
{
    deallocateMem();
    nrows = sc.max_row+1;
    ncols = sc.max_col+1; 
    nelem = 0;
    pos = 0;
    elem = 0;

    sorted = true;
    size_t i,j;
    sc.condenseAsSum();
    for(i=0; i<sc.size(); ++i)
        if (sc.paired_data[i].second != 0.0)
            ++nelem;

    allocateMem();

    for(i=0, j=0; i<sc.size(); ++i)
        if (sc.paired_data[i].second != 0.0)
        {
            elem[j] = sc.paired_data[i].second;
            pos[j] = sc.paired_data[i].first;
            ++j;
        }
    return *this;
}


SparseDoubleMatrix::SparseDoubleMatrix(SparseCollection & sc):
    nrows(0), ncols(0), nelem(0), 
    pos(0), elem(0), sorted(true), tmp(0)
{
    *this = sc;
    return /*SUCCESS*/;
}



SparseDoubleMatrix::~SparseDoubleMatrix()
{
    deallocateMem();
}


SparseDoubleMatrix::SparseDoubleMatrix(const SparseDoubleMatrix& rhs)
: 
nrows(rhs.nrows), ncols(rhs.ncols), nelem(rhs.nelem), 
pos(0), elem(0), sorted(rhs.sorted), tmp(0)
{
    allocateMem();
    for(size_t i=0; i<nelem; ++i)
    {
        pos[i] = rhs.pos[i];
        elem[i] = rhs.elem[i];
    }
}

SparseDoubleMatrix& SparseDoubleMatrix::operator=(const SparseDoubleMatrix& rhs)
{
    if (this != &rhs)
    {
        nrows = rhs.nrows;
        ncols = rhs.ncols;
        if (nelem != rhs.nelem)
        {
            deallocateMem();
            nelem = rhs.nelem;
            allocateMem();
        }
        sorted = rhs.sorted;
        for(size_t i=0; i<nelem; ++i)
        {
            pos[i] = rhs.pos[i];
            elem[i] = rhs.elem[i];
        }
		tmp = 0;
    }
    return *this;
}


/***************** LeftMultiplyTo *******************************/ 
/* 
The most efficient operation on this matrix take M.V and put result into Vout. 

vector Vout has to be allocated before calling this routine and have same size as V

*/
void SparseDoubleMatrix::leftMultiplyTo(DoubleArray& V) const
{
    size_t i;

    if (ncols != V.size()) 
		throw ModelException("SparseDoubleMatrix::leftMultiplyTo(DoubleArray& V)", 
		"size of the vector is not the same as numCols in the Matrix");

	if (!tmp)
		tmp = new double[ncols];

    for(i=0; i<ncols; ++i)
	{
        tmp[i] = V[i];
		V[i] = 0;
	}

    for (i=0; i<nelem; ++i)
        V[getRowFromPos(pos[i])] += elem[i] * tmp[getColFromPos(pos[i])];

    return /*SUCCESS*/;

}

void SparseDoubleMatrix::leftMultiplyTo(double* V) const
{
    size_t i;

    if (V == NULL) 
		throw ModelException("SparseDoubleMatrix::leftMultiplyTo(double* V)", 
							 "called with a null pointer");

	if (!tmp)
		tmp = new double[ncols];

    for(i=0; i<ncols; ++i)
	{
		tmp[i] = V[i];
        V[i] = 0;
	}

    for (i=0; i<nelem; ++i)
        V[getRowFromPos(pos[i])] += elem[i] * tmp[getColFromPos(pos[i])];

    return /*SUCCESS*/;

}
/***************** LeftMultiplyTo (2D) *******************************/ 
/* 
Another very efficient operation on this matrix, to multiply M.Z and return back Zout 
where Z is defined as double** with size [M.ncols x N]

double **Zout must be allocated before calling and be the same size as Z  (i.e. [ncols x N])

for the SRM3 purpose, N is equal to the number of paths for regular dice and can be 
something else for advanced dice.
*/

void SparseDoubleMatrix::leftMultiplyTo(double **Z, size_t N) const
{
    size_t i, j;

    if ( Z==NULL ) 
		throw ModelException("SparseDoubleMatrix::leftMultiplyTo(double **Z, size_t N)", 
							 "called with a null pointer");

	if (!tmp)
		tmp = new double[ncols];

	for (j=0; j<N; ++j)
	{
		for(i=0; i<ncols; ++i)
		{
			tmp[i] = Z[i][j];
			Z[i][j] = 0;
		}
		for (i=0; i<nelem; ++i)
		{
			Z[getRowFromPos(pos[i])][j] += elem[i] * tmp[getColFromPos(pos[i])];
		}

	}

    return /*SUCCESS*/;
}

void SparseDoubleMatrix::leftMultiplyTo(DoubleMatrix &Z) const
{
	size_t i, j;

	if ( Z.numCols() != nrows)  // sic!
		throw ModelException("SparseDoubleMatrix::leftMultiplyTo(DoubleMatrix &Z)", 
							 "Z.numCols() != numRows in the Matrix");

	if (!tmp)
		tmp = new double[nrows];

	size_t nr = Z.numRows();
	for (j=0; j<nr; ++j)
	{
		for(i=0; i<ncols; ++i)
		{
			tmp[i] = Z[i][j];
			Z[i][j] = 0;
		}
		for (i=0; i<nelem; ++i)
		{
			Z[getColFromPos(pos[i])][j] += elem[i] * tmp[getRowFromPos(pos[i])];
		}

	}

	return /*SUCCESS*/;
}


// special extension (for GaussTerm:Dependence)
void SparseDoubleMatrix::leftMultiplyToRow(DoubleMatrix& Z, size_t row) const
{
	size_t i, j=row;

	if ( Z.numCols() != nrows)  // sic!
		throw ModelException("SparseDoubleMatrix::leftMultiplyTo(DoubleMatrix &Z, size_t row)", 
		"Z.numCols() != numRows in the Matrix");

	if (!tmp)
		tmp = new double[nrows];

	for(i=0; i<ncols; ++i)
	{
		tmp[i] = Z[i][j];
		Z[i][j] = 0;
	}
	for (i=0; i<nelem; ++i)
	{
		Z[getColFromPos(pos[i])][j] += elem[i] * tmp[getRowFromPos(pos[i])];
	}

	return /*SUCCESS*/;
}

// aux function
static int cmplong(const void *aa, const void *bb)
{
    return (*(const long*)aa - *(const long*)bb);
}

/**************************** getElem **********************************/ 
/* 
This function getElem(i,j) returns Mij element of the matrix.
It is not very fast if array is not sorted, and reasonably quick if it is.

If you expect to use it often, please consider to call M.Sort() beforehand.

Note, a sparse matrix created from (double**) matrix is created sorted.
A separate sorting is necessary only if a sparse matrix is created 
from 3 arrays (constructor #1).
*/
double SparseDoubleMatrix::getElem(int col, int row)
{
    size_t n = nelem;
    double result;

    if ((row > nrows) || (row < 0)) return 0;
    if ((col > ncols) || (col < 0)) return 0;

    long PosFromRowCol = getPosFromIJ((short)row, (short)col);

    if (!sorted)
            sort();

    //        p = (long*)bsearch(&PosFromRowCol, pos, n, sizeof(long), cmplong);
    pair<long *, long *> p = equal_range(pos, pos + n, PosFromRowCol);
    result = (p.second - p.first == 0) ? 0.0 : elem[p.first - pos]; 

    return result;

}/* getElem */


double SparseDoubleMatrix::GetMMTElem(size_t col, size_t row) const
{
	double val = 0;
	size_t i;
	if (!tmp)
		tmp = new double[nrows];
	fill(tmp, tmp+nrows, 0.0);

	for(i=0; i<nelem; ++i)
		if(getColFromPos(pos[i])==row)
			tmp[getRowFromPos(pos[i])] = elem[i];

	for(i=0; i<nelem; ++i)
		if(getColFromPos(pos[i])==col)
			val += tmp[getRowFromPos(pos[i])]*elem[i];

	return val;
}



/**************************** SortSparseMatrix *******************************/ 
/*
*/
void SparseDoubleMatrix::sort()
{
  const static string method = "SparseDoubleMatrix::sort";

    if (sorted) 
        return /*SUCCESS*/; // nothing to do here

    map<long, double> m;
    for(size_t i=0; i<nelem; ++i)
    {
        if (m.count(pos[i]))
        {
            char buff[256];
            sprintf(buff, "error: duplicates in the sparse matrix! ( row = %d, col = %d )", getRowFromPos(pos[i]), getColFromPos(pos[i]));
            throw ModelException(method, buff);
       }
        m[pos[i]] = elem[i];
    }

    size_t j=0;
    for(map<long, double>::const_iterator it=m.begin(); it!=m.end(); ++it, ++j)
    {
        pos[j] = it->first;
        elem[j] = it->second;
    }

    sorted = true;

    return /*SUCCESS*/;

} 

SparseDoubleMatrix SparseDoubleMatrix::symmToCorrel(double* squareErr, double eigenValueFloor) const 
{
    // for the time being, the approach here will be a "dense core approach"
    // which should work properly for any actual sparse correlation matrix
    // 

    // (1) get a mapping of all the rows containing anything except the diagonal 1


	if (nrows != ncols)
		throw ModelException("SparseDoubleMatrix::computeSquareRoot","The matrix is not square");

    set<unsigned short> dense_all;

    size_t i;
    for(i=0; i<nelem; ++i)
    {
        unsigned short row = getRowFromPos(pos[i]);
        unsigned short col = getColFromPos(pos[i]);
        double val = elem[i];

        if (val != delta(row, col))
        {
            dense_all.insert(row);
            dense_all.insert(col);
        }
    }


    if (dense_all.empty()) // if this was called with a diag-1 matrix
        return *this;

    vector<unsigned short>   dense2sparse; 
    dense2sparse.reserve(dense_all.size());

    map<unsigned short, size_t> sparse2dense;

    for (set<unsigned short>::iterator it = dense_all.begin(); it != dense_all.end(); ++it) 
    {
        sparse2dense[*it] = dense2sparse.size();
        dense2sparse.push_back(*it);
    }

    DoubleMatrix denseMtx(dense_all.size(), dense_all.size());
    SparseCollection sc;
    for(i=0; i<nelem; ++i)
    {
        unsigned short row = getRowFromPos(pos[i]);
        unsigned short col = getColFromPos(pos[i]);
        double val = elem[i];

        if (dense_all.count(row) && dense_all.count(col)) // it's in a dense part of the mtx, all right
        {
            denseMtx[sparse2dense[col]][sparse2dense[row]] = val;
        }
        else
        {
            sc.push_back(row, col, val);
        }
    }

    denseMtx = denseMtx.symmToCorrel(squareErr, eigenValueFloor);

    for(size_t dcol = 0; dcol < dense_all.size(); ++dcol)
    {
        unsigned short col = dense2sparse[dcol];
        for(size_t drow = 0; drow < dense_all.size(); ++drow)
            if (denseMtx[dcol][drow] != 0.0)
            {
                unsigned short row = dense2sparse[drow];
                sc.push_back(row, col, denseMtx[dcol][drow]);
            }
    }

    return SparseDoubleMatrix(sc);
}

/** computes the square root of a matrix aka Choleski / Cholesky etc. */
SparseDoubleMatrix SparseDoubleMatrix::computeSquareRoot() const
{
    // for the time being, the approach here will be a "dense core approach"
    // which should work properly for any actual sparse correlation matrix
    // 

    // (1) get a mapping of all the rows containing anything except the diagonal 1


	if (nrows != ncols)
		throw ModelException("SparseDoubleMatrix::computeSquareRoot","The matrix is not square");

    set<unsigned short> dense_all;

    size_t i;
    for(i=0; i<nelem; ++i)
    {
        unsigned short row = getRowFromPos(pos[i]);
        unsigned short col = getColFromPos(pos[i]);
        double val = elem[i];

        if (val != delta(row, col))
        {
            dense_all.insert(row);
            dense_all.insert(col);
        }
    }


    if (dense_all.empty()) // if this was called with a diag-1 matrix
        return *this;

    vector<unsigned short>   dense2sparse; 
    dense2sparse.reserve(dense_all.size());

    map<unsigned short, size_t> sparse2dense;

    for (set<unsigned short>::iterator it = dense_all.begin(); it != dense_all.end(); ++it) 
    {
        sparse2dense[*it] = dense2sparse.size();
        dense2sparse.push_back(*it);
    }

    DoubleMatrix denseMtx(dense_all.size(), dense_all.size());
    SparseCollection sc;
    for(i=0; i<nelem; ++i)
    {
        unsigned short row = getRowFromPos(pos[i]);
        unsigned short col = getColFromPos(pos[i]);
        double val = elem[i];

        if (dense_all.count(row) && dense_all.count(col)) // it's in a dense part of the mtx, all right
        {
            denseMtx[sparse2dense[col]][sparse2dense[row]] = val;
        }
        else
        {
            sc.push_back(row, col, val);
        }
    }

    denseMtx = denseMtx.computeSquareRoot();
    for(size_t dcol = 0; dcol < dense_all.size(); ++dcol)
    {
        unsigned short col = dense2sparse[dcol];
        for(size_t drow = 0; drow < dense_all.size(); ++drow)
            if (denseMtx[dcol][drow] != 0.0)
            {
                unsigned short row = dense2sparse[drow];
                sc.push_back(row, col, denseMtx[dcol][drow]);
            }
    }

    return SparseDoubleMatrix(sc);
}



void SparseDoubleMatrix::selfMultiplyToSmallmtx(const DoubleMatrix& M, unsigned short colOffset, unsigned short  rowOffset)
{
    SparseCollection sc;
    for(size_t n=0; n < nelem; ++n)
    {
       unsigned short j = getColFromPos(pos[n]);
       if (j<rowOffset || j >= rowOffset + M.numRows())
           continue;
       unsigned short i = getRowFromPos(pos[n]);

       for(unsigned short k = colOffset; k < M.numCols() + colOffset; ++k)
            if (M[k-colOffset][j-rowOffset] != delta(k, j))
                sc.push_back(i, k, elem[n]*(M[k-colOffset][j-rowOffset] - delta(k, j)));
    }
    *this += sc;
}

void SparseDoubleMatrix::smallmtxMultiplyToSelf(const DoubleMatrix& M, unsigned short  colOffset, unsigned short  rowOffset)
{
    SparseCollection sc;
    for(size_t n=0; n < nelem; ++n)
    {
        unsigned short j = getRowFromPos(pos[n]);
        if (j<colOffset || j >= colOffset + M.numCols())
            continue;
        unsigned short k = getColFromPos(pos[n]);

        for(unsigned short i = rowOffset; i < M.numRows() + rowOffset; ++i)
            if (M[j-colOffset][i-rowOffset] != delta(i, j))
                sc.push_back(i, k, elem[n]*(M[j-colOffset][i-rowOffset] - delta(i, j)));
    }
    *this += sc;
}

void SparseDoubleMatrix::print(ostream& out, int mode) const
{
    out<<endl;
    if (mode == 1) // print as full
    {
        DoubleMatrix x(numCols(), numRows());
        for(size_t k=0; k<nelem; ++k)
            x[getColFromPos(pos[k])][getRowFromPos(pos[k])] = elem[k];
        for(size_t i=0; i<numRows(); ++i)
        {
            for(size_t j=0; j<numCols(); ++j)
                out<<x[j][i]<<" \t";
            out<<endl;
        }
    }
    else if (mode == 0) // print as sparse
    {
        for(size_t k=0; k<nelem; ++k)
            out<<"[R: "<<getRowFromPos(pos[k])<<"  C: "<<getColFromPos(pos[k])<<" ]  :=  "<<elem[k]<<endl;
    }
}

CDoubleMatrixSP SparseDoubleMatrix::toDoubleMatrixSP() const
{
	CDoubleMatrixSP pmtx(new CDoubleMatrix( numCols(), numRows() ));
	CDoubleMatrix &mtx = *pmtx;
	for(size_t i=0; i<nelem; ++i)
	{
		mtx[getColFromPos(pos[i])][getRowFromPos(pos[i])] = elem[i];
	}
	return pmtx;
}



//private operation, needed for not-so-efficient xMultiplySelf implementation
// used word no-to- above to highlight that the matrix must be immutable for 
// all the efficient operations, and this operation is mutating the matrix.
SparseDoubleMatrix& SparseDoubleMatrix::operator+= (const SparseCollection& sc)
{
    SparseCollection scfull(sc);
    scfull.reserve(sc.size()+size());

    for(size_t k=0; k<nelem; ++k)
        scfull.push_back(getRowFromPos(pos[k]), getColFromPos(pos[k]), elem[k]);
    scfull.condenseAsSum(); // now condense sparse collection

    *this = scfull;

    return *this;

}


int  SparseDoubleMatrix::checkNoDiags() const 
{
	for(size_t i=0; i<nelem; ++i)
		if (getColFromPos(pos[i]) == getRowFromPos(pos[i]))
			return 1;
	return 0;
}





DRLIB_END_NAMESPACE

