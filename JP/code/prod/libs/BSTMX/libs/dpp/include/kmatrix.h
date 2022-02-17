/***************************************************************
 * Module:	
 * Submodule:	
 * File:	kmatrix.h
 * Function:	Standard Definition and Include files.
 * Author:	A. Chou (modified C. Daher)
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kmatrix_H
#define	_kmatrix_H

#include "kstdinc.h"		// Exceptions, io, etc.

#if !defined (__SUNPRO_CC) || (__SUNPRO_CC != 0x0500) 

//--------------------------------------------------------------
/**
 *  Matrix class mimicking STL.
 */

template <class T> class KMatrix {
public:
	
	/** Default constructor. */
	KMatrix (int nRows=0, int nCols=0, const T**ptr=NULL)
	{
		memAlloc(nRows, nCols);
		if (ptr) {
			int	i, j;
			for (i=0; i<rows(); i++)
			for (j=0; j<cols(); j++) {
				(*this)[i][j] = ptr[i][j];
			}
		}
	}

	/** Copy constructor. */
	KMatrix (const KMatrix&);

	/** Destructor. */
	~KMatrix ()
		{ memFree(); }

	/** Copy operator. */
	KMatrix& operator=(const KMatrix&); 

	/** Set to constant operator. */
	KMatrix& operator=(const T&); 
	

	/** Resizes the matrix.  Values in the matrix are undefined. */
	void resize (int, int);

	/**
	 * Returns the size of the corresponding dimension.
	 * WARNING: C style indexing (i.e. 0=rows, 1=column).
	 */
	int size (int dim) const;


	/** Returns number of rows. */
	int rows() const { return mRows; }

	/** Returns number of rows. */
	int cols() const { return mCols; }


	/** Cast to T**. */
	operator T**();

	/** Cast to T**. */
	operator T**() const;

	/** Access element pointer style. */
	const T* operator[](int idxRow) const
		{ return mPtr[idxRow];}

	/** Access element pointer style. */
	T* operator[](int idxRow)
		{ return mPtr[idxRow];}



	/**
	 * Returns a 1D array
	 * when dim = 1 or 2, and loc tells you which row/col
	 * e.g. get_vector(0, 1) means set the dim=1 (row) to be 0
	 * In other words, get the first column.
	 */
	KVector(T) get_vector(int loc, int dim) const;


	/**
	 * Stores 1D vector.
	 * Same notation as get_vector.
	 */
	KMatrix& store_vector(const KVector(T)&, int loc, int dim);
	

	/** Output stream operator. */
	friend ostream &operator<<(ostream& s, const KMatrix& a) //
	{
		s << "KMatrix(" << a.mRows << "," << a.mCols << ")" << endl;
		
		for (int i = 0; i< a.mRows; i++) {
			for (int j = 0; j< a.mCols; j++) s << a.mPtr[i][j] << " ";
			s << endl;
		}
		s << endl;	
		return s;
	}
	
protected:
				/** Internally allocates memory. */
	void	memAlloc(int nRows, int nCols);
				/** Frees memory. */
	void	memFree();
				/** Number of rows. */
	int	mRows;
				/** Number of columns. */
	int	mCols;
				/** Pointer. */
	T**	mPtr;
};



//--------------------------------------------------------------
// Member function implementation
//



template <class T> inline 
KMatrix<T>::KMatrix(const KMatrix& a)
{
	memAlloc (a.mRows, a.mCols);
	
	for (int i = 0 ; i < mRows; i++)
		for (int j = 0; j < mCols; j++) mPtr[i][j] = a.mPtr[i][j];
}

template <class T> inline 
KMatrix<T>& KMatrix<T>::operator= (const KMatrix<T>& a)
{
	if (this == &a) return *this;
	
	resize(a.mRows, a.mCols);
	
	for (int i = 0 ; i < mRows; i++) {
		for (int j = 0; j < mCols; j++) 
			mPtr[i][j] = a.mPtr[i][j];
	}
	
	return *this;
}

template <class T> inline 
KMatrix<T>& KMatrix<T>::operator= (const T& a)
{
	for (int i = 0 ; i < mRows; i++) {
		for (int j = 0; j < mCols; j++) 
			mPtr[i][j] = a;
	}
	
	return *this;
}

template <class T> inline 
void KMatrix<T>::resize (int nRows, int nCols)
{
	if (nRows != mRows || nCols != mCols) {
		memFree();		
		memAlloc(nRows, nCols);
	}
}

template <class T> inline 
int KMatrix<T>::size(int dim) const
{
	if (dim == 0) return mRows;
	if (dim == 1) return mCols;
	else
		throw KFailure(
		"KMatrix<T>::size(int dim)"
		": illegal matrix dimension %d.\n",
			dim);
	return 0;
}

template <class T> inline 
KMatrix<T>::operator T**() const
{
	return mPtr;
}

template <class T> inline 
KMatrix<T>::operator T**() 
{
	return mPtr;
}

template <class T> inline 
KVector(T) KMatrix<T>::get_vector(int loc, int dim) const
{
	if (dim != 1 && dim != 2) 
	    throw KFailure("KMatrix<T>::get_vector:"
		" illegal dimension %d.\n", dim);

	int len;
	if (dim == 1) len = mCols;
	else len = mRows;

	KVector(T)	ans(len);
	for (int i=0; i< len; i++) {
		if (dim == 1) {
			ans[i] = mPtr[loc][i];
		}
		else {
			ans[i] = mPtr[i][loc];
		}
	}

	return ans;
}

template <class T> inline 
KMatrix<T>& KMatrix<T>::store_vector (const KVector(T)& a, int loc, int dim) 
{
	if (dim != 1 && dim != 2) 
	    throw KFailure("KMatrix<T>::store_vector:"
		" illegal dimension %d.\n", dim);

	int len;
	if (dim == 1) len = mCols;
	else len = mRows;

	if (len != a.size())
	    throw KFailure("KMatrix<T>::store_vector:"
		" wrong-size array (%d != %d).\n",
		len, a.size());

	for (int i = 0; i< len; i++) {
		if (dim == 1) {
			mPtr[loc][i] = a[i];
		}
		else {
			mPtr[i][loc] = a[i];
		}
	}

	return *this;
}


template <class T> inline 
void KMatrix<T>::memAlloc (int nRows, int nCols)
{
	mRows = nRows;
	mCols = nCols;
	
	if (nRows < 0 || nCols < 0)
	    throw KFailure("KMatrix<T>::store_vector:"
		"nRows=%d < 0 or nCols=%d < 0.\n", nRows, nCols);
	
	if (mRows == 0) {
	    mPtr = NULL;
	} else {
	    mPtr = new T* [mRows];
	    if (!mPtr) 
	   	throw KFailure("KMatrix<T>::store_vector: memory.\n");
		
	    for (int i=0; i< mRows; i++) {
		if (mCols == 0) {
		    mPtr[i] = NULL;
		} else {
		    mPtr[i]=new T [mCols]; 
		    if (!mPtr[i]) 
			throw KFailure("KMatrix<T>::store_vector: memory.\n");
	        }
	    }
	}
}

template <class T> inline 
void KMatrix<T>::memFree ()
{
	if (mPtr) 
	{
		for (int i=0; i<mRows; i++) {
			if (mPtr[i]) delete [] mPtr[i];
		}

		delete [] mPtr;
	}
}





//--------------------------------------------------------------


typedef	KVector(double)	KDVector;



//--------------------------------------------------------------
/**
 *   DOuble matrix class with matrix algebra operations.
 */

class KDMatrix : public KMatrix<double>  {
public:
	enum KDMatrixType
		{ IDENTITY };

	/**
	 * Default constructor.
	 */
	KDMatrix (int nRows=0, int nCols=0, const double**ptr=NULL)
		: KMatrix<double>(nRows, nCols, ptr) {};

	/**
	 * Diagonal constructor.
	 */
	KDMatrix(int nRows, KDMatrixType type);


	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KDMatrix& mat);


	/**
	 * Multiplies a matrix and a vector and returns the result vector.
	 */
friend	KDVector operator*(const KDMatrix &aMat, const KDVector &bVect);


	/**
	 * Multiplies two matrices and returns the result.
	 */
	KDMatrix operator*(const KDMatrix &bMat) const;


	/**
	 * Returns the direct sum of two matrices.
	 */
friend KDMatrix directsum(const KDMatrix& aMat, const KDMatrix &bMat);


	/**
	 * Computes and returns the inverse.
	 */
	KDMatrix inverse() const;

	/**
	 * Computes eigen values and vectors.
	 */
	void eigen(KDVector &eigVal, KDMatrix &eigVec) const;



	/**
	 * Computes a pseudo inverse of a matrix using SVD.
	 *
	 * Computes the pseudo inverse of a matrix
	 * obtained from an SVD. The arguments "regType", "alpha"
	 * and "param" determine the type and amount
	 * of regularization applied. Assume
	 * (w_0, ..., w_N) are the eigen values of the input matrix, then
	 * If "regType" is 0, no regulariztion is performed, i.e.
	 * lambda_i = 1/w_i if w_i < 0, 0 otherwise.
	 * If "regType" is 1, performs Tikhonov regulariztion, i.e.
	 * lambda_i = (w_i / alpha w^2 + w_i^2}$$
	 * where w=max|w_i|.
	 * If "regType" is 2, performs Landweber regulariztion, i.e.
	 * lambda_i = {1 - (1- (1 - param w_i^2)^(1/alpha) / w_i.
	 * The constant is dtermined by
	 * a = param / sup w_i^2.
	 * If "regType" is 3, performs cutoff regulariztion, i.e. 
	 * lambda_i = 1/w_i if w_i >= alpha w <> 0, otherwise 0.
	 */

	KDMatrix pseudoinverse(
		int regType = 0,	// (I) see below
		double alpha = 0e0,	// (I) see below
		double param = 0e0)	// (I) see below
 			const;


protected:
				/** Format for I/O. */
	int 	mFormat;
};




#endif












#endif

