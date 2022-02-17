/* ======================== SPARSE MATRIX ======================= */
/** Represents a Sparse m x n matrix where the elements are doubles. */


/******************** Implementing SparseMatrix class ************/
/* 
This implementation is likely to be the most memory and time 
efficient one for a immutable (once initialized and only 
read afterwards) SQUARE sparse matrix M which is likely to have 
non-zero numbers on a diagonal and is used to multiply by 
(double*)  vectors V, with UPDATING those as    

M.V -> V  [or as  V.M -> V ]

or by (double**) matrices U as

Mij Zjp -> Zip   [ or as  Zpi Mij -> Zpj ]

The matrix M is assumed to be not greater than 32768 x 32768

Number of non-zero elements is not limited.

Matrix is REQUIRED to be square.

Complimentary memeber method GetElem(i, j) returns a single 
element Mij, but it should be avoided if possible since 
it is not very efficient.
*/
    /************************************************************************
     *   Implementation Notes
     *   
     *   This implementation is good for a sparse matrix which 
     *   will be routinely multiplied on the left (  V' = S.V ) 
     *   or on the right ( V'' = V.S )
     *
     *   A sparse matrix is best used as an immutable object, or at the 
     *   very least must have an immutable structure. 
     *
     ************************************************************************/
#ifndef SPARSE_DOUBLE_MATRIX_HPP
#define SPARSE_DOUBLE_MATRIX_HPP

#include "edginc/config.hpp"
#include "edginc/DoubleMatrix.hpp"
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include "edginc/smartPtr.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE


class UTIL_DLL StaticPositionHelpers {
public:

/********* Static (local) helper functions ***********************/
/*
These three internal (auxiliary) functions are used to convert 
between (long) position and {(short) row, (short) col} pair

Two (short) coordinates in a sparse array are packed into 
a single (long) 'position'. In hexadecimal format it looks like:

Row:  0xRRRR
Col:  0xCCCC
Position:       0xRRRRCCCC 
*/

static inline long getPosFromIJ(unsigned short row, unsigned short col)
{
    return (((long)row)<<16) + (long)col;
}

static inline short getRowFromPos(long pos)
{
    return (short)(pos>>16);
}

static inline short getColFromPos(long pos)
{
    return (short)pos;
}

};

///  SparseCollection:
// what it is -- an unstructured bunch of elements, 
//  (col, row, value)
// with these properties:
//   (1) nobody cares if (col, row) are repeated
//   (2) new elements can be added at any time, but nothing can be deleted
//   (3) a sparse matrix can be quickly and conveniently built from sparse collection

class UTIL_DLL SparseCollection : private StaticPositionHelpers
{
public:

    SparseCollection() : min_col(0), max_col(0), min_row(0), max_row(0), is_condensed(true) {}
    SparseCollection(unsigned short nrows, unsigned short ncols) : 
    min_col(0), max_col(ncols-1), min_row(0), max_row(nrows-1), is_condensed(true) 
    {
        if (ncols==0 || nrows==0)
            throw ModelException("SparseCollection(nr, nc)","constructor with (nrows,ncols) cannot construct collection with 0 cols or 0 rows");
    }

    void reserve(size_t n)
    {
        paired_data.reserve(n);
    }

    size_t push_back(size_t _row, size_t _col, double elem)
    {
        unsigned short row = _row; 
        unsigned short col = _col;
        QLIB_VERIFY(size_t(row) == _row, "Input row is out of range");
        QLIB_VERIFY(size_t(col) == _col, "Input column is out of range");
        paired_data.push_back(std::make_pair(getPosFromIJ(row, col), elem));
        min_col = std::min(min_col, col);
        max_col = std::max(max_col, col);
        min_row = std::min(min_row, row);
        max_row = std::max(max_row, row);
        if(paired_data.size() > 1 && paired_data.rbegin()->first <= (++paired_data.rbegin())->first)
            is_condensed = false;

        return paired_data.size();
    }
    size_t condenseAsSum()
    {
        if (!is_condensed)
        {
            std::sort(paired_data.begin(), paired_data.end());
			size_t i,j;
            for(i=j=0; j<paired_data.size(); ++i)
            {
                paired_data[i] = paired_data[j];
                while ( ++j < paired_data.size() && paired_data[i].first == paired_data[j].first)
                    paired_data[i].second += paired_data[j].second;
            }
            paired_data.erase(paired_data.begin()+i, paired_data.end());
            is_condensed = true;
        }
        return paired_data.size();
    }
    size_t size() const { return paired_data.size();}
    void print(std::ostream& out)
    {
        for(size_t i=0; i<paired_data.size(); ++i)
            out<<getRowFromPos(paired_data[i].first)<<" "
                <<getColFromPos(paired_data[i].first)<<" "
                <<paired_data[i].second<<std::endl;
    }

	// calc X = a * I + b * X, whrere a and b are scalars, I is unit matrix and X is SparseCollection
	void linearTransform(double a, double b)
    {
    	if (max_row-min_row != max_col-min_col)
            throw ModelException("linearTransform(a,b)","The same number of rows and columns are expected.");

		for(size_t j=0; j<paired_data.size(); ++j){
			paired_data[j].second *= b;
        }
		if (!Maths::isZero(a)){
			for (unsigned short i=0; i < max_row-min_row+1; ++i){
				push_back(i + min_row, i + min_col, a);
			}
			is_condensed = false;
		}
		condenseAsSum();
    }

	void leftMultiplyTo(double* V) const
	{
		if (V == NULL) 
			throw ModelException("SparseDoubleMatrix::leftMultiplyTo(double* V)", 
								"called with a null pointer");

		int ncols = max_col - min_col + 1;
		double *tmp = new double[ncols];

		for(unsigned short i=0; i<ncols; ++i)
		{
			tmp[i] = V[i];
			V[i] = 0;
		}

		for (size_t i=0; i<size(); ++i)
			V[getRowFromPos(paired_data[i].first)] += paired_data[i].second * tmp[getColFromPos(paired_data[i].first)];

		return /*SUCCESS*/;
	}

	void convert( size_t i, Imsl_d_sparse_elem & sparse_elem ) const 
	{
		sparse_elem.row = getRowFromPos(paired_data[i].first);
		sparse_elem.col = getColFromPos(paired_data[i].first);
		sparse_elem.val = paired_data[i].second;
		return;
	}

	void convert( size_t i, int& row, int& col, double& val ) const 
	{
		row = getRowFromPos(paired_data[i].first);
		col = getColFromPos(paired_data[i].first);
		val = paired_data[i].second;
		return;
	}

	void empty()
	{
		paired_data.erase(paired_data.begin(), paired_data.end());
		is_condensed = true;
	}

	size_t nrows() const {
		return max_row - min_row + 1;
	}

	size_t ncols() const {
		return max_col - min_col + 1;
	}

private:
    friend class SparseDoubleMatrix;

    typedef std::vector<std::pair<long, double> > VectorPairedValues;
    VectorPairedValues paired_data;

    unsigned short min_col, max_col, min_row, max_row;
    bool is_condensed;

};

typedef refCountPtr<SparseCollection> SparseCollectionSP;

class UTIL_DLL SparseDoubleMatrix : private StaticPositionHelpers
{
public:
    /** Creates a sparse n x m matrix from 3 arrays (posX, posY, valXY) */
    SparseDoubleMatrix(int numCols, 
                        int numRows, 
                        const IntArray& colPos,
                        const IntArray& rowPos,
                        const DoubleArray& values);

    /** Creates a n x m matrix from a 2-dimensional data[col][row]      */
    SparseDoubleMatrix(int numCols, int numRows, const double** data);

    /** Creates a n x m matrix from a CDoubleMatrix     */
    SparseDoubleMatrix(const DoubleMatrix& data);

    /** Creates a n x m matrix from a sparse collection up to max element */
    SparseDoubleMatrix& operator=(SparseCollection & data);
    SparseDoubleMatrix(SparseCollection & data);

    SparseDoubleMatrix(const SparseDoubleMatrix& rhs);
    SparseDoubleMatrix& operator=(const SparseDoubleMatrix& rhs);

    ~SparseDoubleMatrix();

    /** Returns the number of columns in the matrix. Inline for performance */
    inline size_t numCols() const{
        return ncols;
    }

    /** Returns the number of rows in the matrix. Inline for performance */
    inline size_t numRows() const{
        return nrows;
    }

    inline size_t size() const{
        return nelem;
    }


	/** A transformation to CDoubleMatrix is possible while not recommended */
	CDoubleMatrixSP toDoubleMatrixSP() const;


    /** Self-modificators (useful when building the matrix )              */

    /** The functions below are used to L- and R- multiply a sparse matrix 
        to an almost identity matrix with SMALL (e.g. 3x3) adjustment 
        -- most typical in Get3A() correlation approach 

        The matrix to multiply self to will be properly sized and equal to
        Identity everywhere except the square superimposed by the provided 
        matrix M starting from colOffset and rowOffset. 
        
        For example: multiplying (3x3) matrix "*this" by using:

        this->selfMultiplyToSmallmtx(DoubleMatrix(1, 1, 3.5), 1, 2)  means that
                
                      [ 1   0   0 ]
        this = this * [ 0   1   0 ]
                      [ 0   3.5 1 ]
 
        this->smallmtxMultiplyToSelf(DoubleMatrix(2, 2, -2.1), 0, 0) means that

                [ -2.1  -2.1   0 ]
        this =  [ -2.1  -2.1   0 ] * this
                [  0     0     1 ]

        */
    void selfMultiplyToSmallmtx(const DoubleMatrix& M, unsigned short  colOffset, unsigned short  rowOffset);
    void smallmtxMultiplyToSelf(const DoubleMatrix& M, unsigned short  colOffset, unsigned short  rowOffset);



    /**  The most important use of the sparse matrix -- in-place multiplying a vector or a matrix
         or an array of doubles   
		 
		 Vi  := Mij Vj or
		 Zik := Mik Zjk
	
		Note that the temporary space necessary for this computation is allocated automatically once
	*/
    // with class types
	void leftMultiplyTo(DoubleArray& V) const;
	void leftMultiplyTo(DoubleMatrix& V) const;

	// with primitive types
    void leftMultiplyTo(double* V) const;
    void leftMultiplyTo(double **Z, size_t N) const;


    /** A special multiplication (by request) to a single row of the DoubleMatrix */
	void leftMultiplyToRow(DoubleMatrix& V, size_t k) const;

    /** computes the square root of a matrix aka Choleski / Cholesky etc. */
    SparseDoubleMatrix computeSquareRoot() const;

    /** input:  symmetric matrix 
        output: correlation matrix, average squared difference between symmetric and covariance matrices
        (condense the matrix and call the DoubleMatrix::symmToCorrel(...) ) */
    SparseDoubleMatrix symmToCorrel(double* squareErr, double eigenValueFloor) const;

	// get an element from  M.Mt, i.e. returns Eij = sum_k (Mik Mjk)
	// without computing the whole matrix. (restoring the correlations)
	double GetMMTElem(size_t col, size_t row) const;

    // For debugging purposes:
    //   mode can be :  0 for compressed format [pos_i pos_j val_ij]
    //                  1 for full matrix
    void print(std::ostream& out, int mode = 0) const;

	// For checking input beta matrix
	int  checkNoDiags() const;

	SparseDoubleMatrix();

    inline void erase()
    {
        deallocateMem();
        nrows = ncols = nelem = 0;
        sorted = true;
    }


private:

    unsigned short   nrows;  /* total nb of rows in full matrix */
    unsigned short   ncols;  /* total nb of cols in full matrix */ 
    size_t  nelem;  /* nb of non-trivial elements in sparse matrix list */
    long   *pos;    /* list of 'scrambled' positions of non-trivial elements */
    double *elem;   /* list of values of non-trivial elements */
    bool    sorted; /* whether the pos list is sorted or not */

	mutable double *tmp;    /* array to store values for left multiplication */ 

    void allocateMem();
    void deallocateMem();

    // should be avoided as it is very inefficient
    double getElem(int col, int row);

	// presorting makes GetElem work a tad faster
    void sort();

    SparseDoubleMatrix& operator+= (const SparseCollection& sc);

    inline double delta(size_t i, size_t j) const
    {
        return (i==j)?1.0:0.0;
    }
};

typedef refCountPtr<SparseDoubleMatrix> SparseDoubleMatrixSP;


DRLIB_END_NAMESPACE

#endif
