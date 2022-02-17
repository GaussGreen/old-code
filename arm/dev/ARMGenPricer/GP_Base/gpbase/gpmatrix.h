/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gpmatrix.h
 *
 *  \brief a STL wrapped up matrix
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPINFRA_GPMATRIXTEMP_H
#define _INGPINFRA_GPMATRIXTEMP_H

#include "removeidentifiedwarning.h"
#include "port.h"
#include "env.h"
#include "ostringstream.h"
#include "gpvector.h"
#include "checkarg.h"
#include "gpvector.h"
#include "assignop.h"
#include "rootobject.h"

#include <vector>
CC_USING_NS(std,vector)
#include <string>
CC_USING_NS(std,string)
#include <iomanip>

/// kernel
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

#if defined(__ARM_MATRIX_NO_RANGE_CHECK)
	#define CHECKRANGE(i,j)  
#else
	#define CHECKRANGE(i,j) CheckRange(i,j);
#endif

template <typename T> class ARM_GP_T_Vector;
template <typename T> class ARM_GP_T_Tensor;

//////////////////////////////////////////////
/// \class ARM_GP_T_Matrix 
/// \brief template matrix class based on the STL vector ... 
///			- this wrapper around the STL vector allows us to
///				retrieve all the STL vector interface
///			- matrix in Column major (meaning that the position of row=i, col=j 
///				is i*itsColsNb+j
//////////////////////////////////////////////
template <typename T> class ARM_GP_T_Matrix : public ARM_RootObject
{
public:
	/// standard typedefs
	typedef T value_type;

#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef typename CC_NS(std,vector)<T>::iterator			iterator;
	typedef typename CC_NS(std,vector)<T>::const_iterator	const_iterator;   
#else
	typedef CC_NS(std,vector)<T>::iterator iterator;
	typedef CC_NS(std,vector)<T>::const_iterator const_iterator;
#endif 

    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef T const & const_reference;
    typedef size_t size_type;

	/// constructors
    inline explicit ARM_GP_T_Matrix(size_t row=0, size_t col=0, const T& initValue = T())
	:	ARM_RootObject(), itsRowsNb(row), itsColsNb(col), itsValues(row*col,initValue){ CC_ARM_SETNAME(ARM_GP_MATRIX);}
    inline explicit ARM_GP_T_Matrix(size_t row, size_t col, T* values)
	:	ARM_RootObject(), itsRowsNb(row), itsColsNb(col), itsValues(values,values+row*col){ CC_ARM_SETNAME(ARM_GP_MATRIX); }
    inline explicit ARM_GP_T_Matrix(size_t row, size_t col, const vector<T>& stl_vector )
	:	ARM_RootObject(), itsRowsNb(row), itsColsNb(col), itsValues(stl_vector )
	{ CC_NS(ARM_Check,CheckArgSize)( stl_vector, "stl_vector", row*col ); CC_ARM_SETNAME(ARM_GP_MATRIX); }
    inline explicit ARM_GP_T_Matrix(size_t row, size_t col, const ARM_GP_T_Vector<T>& gp_vector )
	:	ARM_RootObject(), itsRowsNb(row), itsColsNb(col), itsValues(gp_vector.GetValues() )
	{ CC_NS(ARM_Check,CheckArgSize)( gp_vector, "gp_vector", row*col ); CC_ARM_SETNAME(ARM_GP_MATRIX); }

	/// constructor to create a submatrix from an already constructred matrix
	/// Row = row-index of the start point
	/// Col = column-index of the start point
	/// N_rows = number of rows to add to (Row,Col)
	/// N_cols = number of cols to add to (Row,Col)
    inline explicit ARM_GP_T_Matrix( ARM_GP_T_Matrix<T> const& mat, size_t row, size_t col, size_t nRows, size_t nCols );

	/// conversion from matrix and tensor! (FIX FIX)
	inline explicit ARM_GP_T_Matrix(ARM_GP_T_Vector<T> const& matrix ) {CC_ARM_SETNAME(ARM_GP_MATRIX);};
	inline explicit ARM_GP_T_Matrix(ARM_GP_T_Tensor<T> const& tensor ) {CC_ARM_SETNAME(ARM_GP_MATRIX);};

	/// copy constructor
	inline ARM_GP_T_Matrix(ARM_GP_T_Matrix<T> const& rhs )
	:	ARM_RootObject(rhs),
		itsRowsNb(rhs.itsRowsNb),
		itsColsNb(rhs.itsColsNb), 
		itsValues( rhs.itsValues ) {};
	/// assignment operator
	ASSIGN_OPERATOR(ARM_GP_T_Matrix<T>)
	
	/// destructor
    virtual ~ARM_GP_T_Matrix(){};
	/// assignment to a constant of type T
	const ARM_GP_T_Matrix<T>& operator=( T const & a )
	{
		itsValues.assign(itsValues.size(),a);
		return *this;
	}

	/// traditional accessor!
	inline vector<T> GetValues() const { return itsValues;}
	inline void SetValues(const vector<T>& values)  { itsValues = values;}

    inline T operator()(size_t i,size_t j) const{ CHECKRANGE(i,j); return itsValues[i*itsColsNb+j];}
    inline T& operator()(size_t i,size_t j)		{ CHECKRANGE(i,j); return itsValues[i*itsColsNb+j];}
    inline T Elt(size_t i,size_t j) const		{ CHECKRANGE(i,j); return itsValues[i*itsColsNb+j];}
    inline T& Elt(size_t i,size_t j)			{ CHECKRANGE(i,j); return itsValues[i*itsColsNb+j];}

	/// matrix specific operation
    inline ARM_GP_T_Matrix<T>& transpose();
    inline bool IsSquared() const { return itsRowsNb == itsColsNb;}
	inline bool IsDiagonal() const;
    inline bool IsSymmetric() const;
    inline void CheckSquaredMatrix() const;
    inline void CheckSymmetricMatrix() const;
	inline void CheckCorrelMatrix() const;
	inline void Rot( double s, double tau, size_t i, size_t j, size_t k, size_t l );
	inline T trace() const;

    /// Be careful, this accessors return a copy  not pointor
    /// please don't forget to free memory
    inline ARM_GP_T_Vector<T>* GetColumn(size_t i) const;
    inline ARM_GP_T_Vector<T>* GetRow(size_t i) const;

	//for lazy people, it return a reference 
	inline ARM_GP_T_Vector<T> GetColumns(size_t i) const;
    inline ARM_GP_T_Vector<T> GetRows(size_t i) const;

	/// standard STD Vector function
	/// size and memory
	inline size_t rows() const	{ return itsRowsNb; }
	inline size_t cols() const	{ return itsColsNb; }
	inline size_t size() const	{ return itsColsNb * itsRowsNb; }
	inline bool empty() const { return itsValues.empty(); }
	inline void reserve(size_t rowsNb, size_t colsNb) { itsValues.reserve(rowsNb*colsNb); }
	inline void resize(size_t rowsNb,size_t colsNb){ itsValues.resize(rowsNb*colsNb); itsRowsNb=rowsNb; itsColsNb=colsNb;}
    inline size_t GetRowsNb() const {return itsRowsNb;}
    inline size_t GetColsNb() const {return itsColsNb;}
	inline void assign( iterator beg, iterator last) { itsValues.assign(beg,last); }

    ///added new line or row
    inline void push_backRow(const ARM_GP_T_Vector<double>& RowValues);
    inline void push_backColumn(const ARM_GP_T_Vector<double>& ColValues);

	/// Standard ARM Object support
	virtual ARM_Object* Clone() const { return new ARM_GP_T_Matrix<T>(*this);}
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// iterator support
	iterator begin() { return itsValues.begin(); }
	const_iterator begin() const { return itsValues.begin(); }
	iterator end() { return itsValues.end(); }
	const_iterator end() const { return itsValues.end(); }

	/// comparison operator
	inline bool operator==( ARM_GP_T_Matrix<T> const& rhs ) { return itsRowsNb == rhs.itsRowsNb && itsColsNb == rhs.itsColsNb && itsValues == rhs.itsValues; }
	inline bool operator!=( ARM_GP_T_Matrix<T> const& rhs ) { return itsRowsNb != rhs.itsRowsNb || itsColsNb != rhs.itsColsNb || itsValues != rhs.itsValues; }

	///  operation on a simple elem: addition, subtraction, multiplication, division
	ARM_GP_T_Matrix<T>& operator+=( const T&  val ){ return GeneralTOp(val, CC_NS(std,plus)<T>() ); }
	ARM_GP_T_Matrix<T>& operator-=( const T& val ) { return GeneralTOp(val, CC_NS(std,minus)<T>() ); }
	ARM_GP_T_Matrix<T> & operator*=( const T& val ){ return GeneralTOp(val, CC_NS(std,multiplies)<T>() ); }
	ARM_GP_T_Matrix<T> & operator/=( const T& val ){ return GeneralTOp(val, CC_NS(std,divides)<T>() ); }

	// Element-by-element: addition, subtraction, multiplication, division, 
	ARM_GP_T_Matrix<T>& operator+=( ARM_GP_T_Matrix<T> const& rhs ){ return GeneralMatrixOp(rhs, CC_NS(std,plus)<T>() ); }
	ARM_GP_T_Matrix<T>& operator-=( ARM_GP_T_Matrix<T> const& rhs ){ return GeneralMatrixOp(rhs, CC_NS(std,minus)<T>() ); }
	ARM_GP_T_Matrix<T>& operator*=( ARM_GP_T_Matrix<T> const& rhs );
	ARM_GP_T_Matrix<T>& operator/=( ARM_GP_T_Matrix<T> const& rhs ){ return GeneralMatrixOp(rhs, CC_NS(std,divides)<T>() ); }



private:
	size_t		itsRowsNb,
				itsColsNb;
    vector<T>	itsValues;

	void CheckRange(size_t i, size_t j) const
	{
#if !defined( __ARM_MATRIX_NO_RANGE_CHECK )
		if ( i< 0 || j<0 || i>itsRowsNb-1 || j>itsColsNb-1)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "range error in ARM_GP_T_Matrix");
#endif 
	}

	/// general purpose operator for T operator op
	template <typename TOp> ARM_GP_T_Matrix<T>&	GeneralTOp( const T& rhs, TOp& Op )
	{
		for( iterator iter=begin(), last=end(); iter!=last; ++iter )
			*iter = Op( *iter, rhs );
		return *this;
	}

	/// general purpose operator for T operator op on a matrix
	template <typename TOp> ARM_GP_T_Matrix<T>&	GeneralMatrixOp(ARM_GP_T_Matrix<T> const& rhs, TOp& Op )
	{
		if( rows()!=rhs.rows() || cols()!=rhs.cols()) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "range error in ARM_GP_T_Matrix::GeneralMatrixOp");

		for( iterator iter=begin(), iterator last=end(), const_iterator iterRhs=rhs.begin(); iter!=last; ++iter; ++iterRhs )
			*iter = Op( *iter, *iterRhs );
		return *this;
	}

};

/// standard comparison operator (symmetric version)
template <typename T> bool operator==( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs ){	return lhs == rhs; }
template <typename T> bool operator!=( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs ){	return lhs != rhs; }


/// standard multiplication *= operator
template <typename T>
ARM_GP_T_Matrix<T>& ARM_GP_T_Matrix<T>::operator*=( ARM_GP_T_Matrix<T> const& rhs )
{
	size_t rowsNb, colsNb, other;
	rowsNb = itsRowsNb;
	colsNb = rhs.cols();
	other = itsColsNb;

	if( itsColsNb!=rhs.rows() )  // Sanity Checks
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "size of matrix are incompatible");

	vector<T> results( rowsNb*colsNb,0 );

	for ( size_t i = 0 ; i < rowsNb ; i++ )
		for ( size_t j = 0 ; j < colsNb ; j++ )
			for ( size_t k = 0 ; k < other ; k++ )
				results[i*colsNb+j] += itsValues[i*itsColsNb+k]*rhs(k,j);

	itsRowsNb = rowsNb;
	itsColsNb = colsNb;
	itsValues.swap( results );
	return *this;
}

/// standard addition, subtraction, multiplication, division, 
template <typename T> ARM_GP_T_Matrix<T> operator+( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(lhs) += rhs; }
template <typename T> ARM_GP_T_Matrix<T> operator+( ARM_GP_T_Matrix<T> const& lhs, const T& val ) { return ARM_GP_T_Matrix<T>(lhs) += val; }
template <typename T> ARM_GP_T_Matrix<T> operator+( const T& val, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(rhs.size(),val) += rhs; }

template <typename T> ARM_GP_T_Matrix<T> operator-( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(lhs) -= rhs; }
template <typename T> ARM_GP_T_Matrix<T> operator-( ARM_GP_T_Matrix<T> const& lhs, const T& val ) { return ARM_GP_T_Matrix<T>(lhs) -= val; }
template <typename T> ARM_GP_T_Matrix<T> operator-( const T& val, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(rhs.size(),val) -= rhs; }

template <typename T> ARM_GP_T_Matrix<T> operator*( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs )  { return ARM_GP_T_Matrix<T>(lhs) *= rhs; }
template <typename T> ARM_GP_T_Matrix<T> operator*( ARM_GP_T_Matrix<T> const& lhs, const T& val ) { return ARM_GP_T_Matrix<T>(lhs) *= val; }
template <typename T> ARM_GP_T_Matrix<T> operator*( const T& val, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(rhs.size(),val) *= rhs; }

template <typename T> ARM_GP_T_Matrix<T> operator/( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(lhs) /= rhs; }
template <typename T> ARM_GP_T_Matrix<T> operator/( ARM_GP_T_Matrix<T> const& lhs, const T& val ) { return ARM_GP_T_Matrix<T>(lhs) /= val; }
template <typename T> ARM_GP_T_Matrix<T> operator/( const T& val, ARM_GP_T_Matrix<T> const& rhs ) { return ARM_GP_T_Matrix<T>(rhs.size(),val) /= rhs; }

/// standard multiplication by a Vector 
template <typename T> ARM_GP_T_Vector<T> * operator*( ARM_GP_T_Matrix<T> const& lhs, ARM_GP_T_Vector<T> const& rhs )
{
	size_t size_result = lhs.rows();
	size_t size = rhs.size();
	ARM_GP_T_Vector<T> * result = new ARM_GP_T_Vector<T>( size_result, 0.0 );

	if ( lhs.cols() != size )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "size of matrix and vector do not match");
	
	for ( size_t i = 0 ; i < size_result ; i++ )
		for ( size_t k = 0 ; k < size ; k++ )
			( *result )[i] += lhs(i,k)*rhs[k];

	return result;
}

/// left multiplication by a Vector 
template <typename T> ARM_GP_T_Vector<T> * operator*( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Matrix<T> const& rhs )
{
	size_t size_result = rhs.rows();
	size_t size = lhs.size();
	ARM_GP_T_Vector<T> * result = new ARM_GP_T_Vector<T>( size_result, 0.0 );

	if ( rhs.cols() != size )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "size of matrix and vector do not match");
	
	for ( size_t i = 0 ; i < size_result ; i++ )
		for ( size_t k = 0 ; k < size ; k++ )
			( *result )[i] += lhs[k]*rhs(k,i);

	return ( *result );
}

/// constructor from submatrix!
template <typename T> ARM_GP_T_Matrix<T>::ARM_GP_T_Matrix( ARM_GP_T_Matrix<T> const& mat, size_t row, size_t col, size_t nRows, size_t nCols )
:	itsRowsNb(nRows), itsColsNb(nCols), itsValues(nRows*nCols)
{
	if( mat.rows() < row+nRows)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "submatrix row too big");
	if( mat.cols() < col+nCols)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "submatrix col too big");
	for( size_t i=0; i<nRows; ++i )
		for( size_t j=0; j<nCols; ++j )
			itsValues[i*nCols+j] = mat(i+row,j+col);
	CC_ARM_SETNAME(ARM_GP_MATRIX);
}

///added new row
template <typename T> void ARM_GP_T_Matrix<T>::push_backRow(const ARM_GP_T_Vector<double>& RowValues)
{
	if( empty() )
	{
        itsRowsNb = 1;
        itsColsNb = RowValues.size();
        itsValues = RowValues.GetValues();
	}
	else
	{
        CC_NS(ARM_Check,CheckArgSize)( RowValues, "RowValues", itsColsNb);
		itsValues.reserve( size()+itsColsNb );
		for( size_t i=0; i<itsColsNb; ++i )
                itsValues.push_back(RowValues[i]);
		++itsRowsNb;
    }
};


///added new line
template <typename T> void ARM_GP_T_Matrix<T>::push_backColumn(const ARM_GP_T_Vector<double>& ColValues)
{
    CC_NS(ARM_Check,CheckArgSize)( ColValues, "ColValues", itsRowsNb);
	if( empty() )
	{
        itsColsNb = 1;
        itsRowsNb = ColValues.size();;
        itsValues = ColValues.GetValues();
	}
	else
	{
        CC_NS(ARM_Check,CheckArgSize)( ColValues, "ColValues", itsRowsNb);
		itsValues.resize( size()+itsRowsNb );
		for(int i=itsRowsNb-1; i>0; i--)
			for( int j=itsColsNb-1; j>=0; j-- )
				itsValues[i*(itsColsNb+1)+j]= itsValues[i*itsColsNb+j];
		for(i=itsRowsNb-1; i>=0; --i )
			itsValues[i*(itsColsNb+1)+itsColsNb]= ColValues[i];
		++itsColsNb;
	}
};


/// toString
template <typename T> string ARM_GP_T_Matrix<T>::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
	os << "GP_MATRIX[" << itsRowsNb << "," << itsColsNb << "]\n";

	for(size_t i=0; i<itsRowsNb; ++i)
	{
		for(size_t j=0; j<itsColsNb; ++j)
		{
			os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
				<< CC_NS(std,setw)(8) << Elt(i, j) << "\t";
		}
		os << CC_NS(std,endl);
	}
    return os.str();
}



/// transpose the matrix
template <typename T> ARM_GP_T_Matrix<T>& ARM_GP_T_Matrix<T>::transpose()
{
	vector<T> results( itsRowsNb*itsColsNb ); 
	for(size_t i=0;i<itsRowsNb;++i)
		for(size_t j=0;j<itsColsNb;++j)
			results[j*itsRowsNb+i] = itsValues[i*itsColsNb+j];
	CC_NS(std, swap)(itsRowsNb,itsColsNb);
	itsValues.swap( results );
	return *this;
}



/// test if a matrix is symmetric
template <typename T> bool ARM_GP_T_Matrix<T>::IsSymmetric() const
{ 
	if(!IsSquared()) 
		return false; 
	
	for(size_t i=0;i<itsRowsNb;++i)
		for(size_t j=0;j<i;++j)  
			if( fabs((*this)(i,j)-(*this)(j,i)) > K_NEW_DOUBLE_TOL ) 
				return false;
	return true; 
};



template <typename T> void ARM_GP_T_Matrix<T>::CheckSquaredMatrix() const
{
    if(!IsSquared())
    {
        char msg[255];
		sprintf( msg, "%s: itsRowsNb: %d while itsColsNb : %d", ARM_USERNAME.c_str(),itsRowsNb, itsColsNb);
        sprintf( msg, " matrix is not squared"); 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
    }
}


template <typename T> bool ARM_GP_T_Matrix<T>::IsDiagonal() const
{
	if(IsSquared())
    {
		for(size_t i=0;i<itsRowsNb;++i)
		{ 
			for(size_t j=0;j<i;++j)
			{
				if( fabs((*this)(i,j)) > K_NEW_DOUBLE_TOL ) 
					return false;
				if( fabs((*this)(j,i)) > K_NEW_DOUBLE_TOL ) 
					return false;
			}
		}
		return true;
	}
	else
		return false;
}


template <typename T> void ARM_GP_T_Matrix<T>::CheckCorrelMatrix() const
{
	CheckSquaredMatrix();
	CheckSymmetricMatrix();
	size_t i,j;

	/// check diagonal filled with one!
	for(i=0;i<itsRowsNb;++i)
	{ 
		if( fabs((*this)(i,i)-1.0) > K_NEW_DOUBLE_TOL ) 
		{
			char msg[255];
			sprintf( msg, "%s: Corraltion matrix: term (%d,%d) of diagonal = %d != 1.0", ARM_USERNAME.c_str(), i,i, (*this)(i,i) );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		}
	}

	/// check terms of the upper and lower triangular matrices between -1.0 and 1.0
	for(i=0;i<itsRowsNb;++i)
	{ 
		for(j=0;j<i;++j)  
		{
			if( (*this)(i,j) < -1.0-K_NEW_DOUBLE_TOL )
			{
				char msg[255]; sprintf( msg, "%s: elem(%d,%d) = %d < -1.0!", ARM_USERNAME.c_str(), i, j, (*this)(i,j) );
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
			}

			if( (*this)(i,j) > 1.0+K_NEW_DOUBLE_TOL )
			{
				char msg[255]; sprintf( msg, "%s: elem(%d,%d) = %d > 1.0!", ARM_USERNAME.c_str(), i, j, (*this)(i,j) );
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
			}
		}
	}
};


template <typename T> void ARM_GP_T_Matrix<T>::CheckSymmetricMatrix() const 
{
    if(!IsSymmetric())
    {
        char msg[255];
        sprintf( msg, " matrix is not symmetrical"); 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
    }
}


/// Rotate
template <typename T> void ARM_GP_T_Matrix<T>::Rot( double s, double tau,
	size_t i, size_t j, size_t k, size_t l )
{
    double g = Elt(i,j), h = Elt(k,l);
    Elt(i,j) = g - s*(h+g*tau); Elt(k,l) = h + s*(g-h*tau);
}


/// Warning no check to have symmetric matrix to avoid time penalty!
template <typename T> T ARM_GP_T_Matrix<T>::trace() const
{
	double tr=0.0;
	for(size_t i=0;i<itsRowsNb;++i)
		tr += Elt(i,i);
	return tr;
}


template <typename T> ARM_GP_T_Vector<T>* ARM_GP_T_Matrix<T>::GetColumn(size_t j) const
{ 
	ARM_GP_T_Vector<T>* result = new ARM_GP_T_Vector<T>(itsRowsNb);
	for( size_t i=0; i<itsRowsNb; ++i )
		result->Elt(i) = Elt(i,j);
	return result;
}

template <typename T> ARM_GP_T_Vector<T>* ARM_GP_T_Matrix<T>::GetRow(size_t i) const
{
    ARM_GP_T_Vector<T>* result = new ARM_GP_T_Vector<T>(itsColsNb);
	for( size_t j=0; j<itsColsNb; ++j )
		result->Elt(j) = Elt(i,j);
	return result;
}
template <typename T> ARM_GP_T_Vector<T> ARM_GP_T_Matrix<T>::GetColumns(size_t j) const
{ 
	ARM_GP_T_Vector<T> result(itsRowsNb);
	for( size_t i=0; i<itsRowsNb; ++i )
		result[i] = Elt(i,j);
	return ARM_GP_T_Vector<T>(result);
}

template <typename T> ARM_GP_T_Vector<T> ARM_GP_T_Matrix<T>::GetRows(size_t i) const
{
    ARM_GP_T_Vector<T> result(itsColsNb);
	for( size_t j=0; j<itsColsNb; ++j )
		result->Elt(j) = Elt(i,j);
	return ARM_GP_T_Vector<T>(result);
}


//FIXMEFRED : will do it later...

//template <> 
//bool ARM_GP_T_Matrix<bool>::operator()(size_t i,size_t j) const
//{ CHECKRANGE(i,j); bool item; item = itsValues[i*itsColsNb+j]; return item; }

//template <> 
//bool& ARM_GP_T_Matrix<bool>::operator()(size_t i,size_t j)		
//{ CHECKRANGE(i,j); bool item; item = itsValues[i*itsColsNb+j]; return item; }

//template <> 
//bool ARM_GP_T_Matrix<bool>::Elt(size_t i,size_t j) const		
//{ CHECKRANGE(i,j); bool item; item = itsValues[i*itsColsNb+j]; return item; }

//template <> 
//bool& ARM_GP_T_Matrix<bool>::Elt(size_t i,size_t j)			
//{ CHECKRANGE(i,j); bool item; item = itsValues[i*itsColsNb+j]; return item; }


/// undef MACRO
#undef CHECKRANGE

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

