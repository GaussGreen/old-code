/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gpmatrixtridiag.h
 *
 *  \brief a tridiagonal matrix implementation
 *	\author  A. Schauly
 *	\version 1.0
 *	\date September 2005
 */

#ifndef _INGPBASE_TRIDIAGONALMATRIX_H
#define _INGPBASE_TRIDIAGONALMATRIX_H

#include "env.h"
#include "port.h"
#include "gpmatrix.h"
#include "gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////
/// \class ARM_GP_T_TridiagMatrix
/// \brief template tridiagMatrix class
//////////////////////////////////////////////
template <typename T> class ARM_GP_T_TridiagMatrix : public ARM_RootObject
{
private:
	vector<T> itsDiagonal, itsUpperDiag, itsLowerDiag; /// The 3 diagonals of the matrix
	vector<T> itsNormalizationTerms, itsOtherCoeffs;  /// 2 vectors used for faster inversion for the thomas algortihm
	size_t itsSize;

public:
	/// Constructor/Destructors
	ARM_GP_T_TridiagMatrix<T>( const ARM_GP_T_TridiagMatrix<T>& rhs ) : ARM_RootObject(rhs), itsSize(rhs.itsSize), 
		itsDiagonal(rhs.itsDiagonal), itsUpperDiag(rhs.itsUpperDiag), itsLowerDiag(rhs.itsLowerDiag),itsNormalizationTerms(0), itsOtherCoeffs(0)  {}
	ARM_GP_T_TridiagMatrix<T>() : ARM_RootObject(), itsSize(0), itsDiagonal(0), itsUpperDiag(0), itsLowerDiag(0), itsNormalizationTerms(0), itsOtherCoeffs(0) {}
	ARM_GP_T_TridiagMatrix<T>(size_t size) : ARM_RootObject(), itsSize(size), itsDiagonal(size), itsUpperDiag(size-1), itsLowerDiag(size-1), itsNormalizationTerms(0), itsOtherCoeffs(0) {}
	virtual ~ARM_GP_T_TridiagMatrix<T>() {}

	ARM_GP_T_TridiagMatrix<T>& CopyNoInverse( const ARM_GP_T_TridiagMatrix<T>& rhs );

	/// Standard ARM Object support
	virtual ARM_Object* Clone() const { return new ARM_GP_T_TridiagMatrix<T>(*this);}
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// Accessors
	inline void setUpperDiag( const vector<T>& vec ) { itsUpperDiag = vec; }
	inline void setLowerDiag( const vector<T>& vec ) { itsLowerDiag = vec; }
	inline void setDiagonal( const vector<T>& vec ) { itsDiagonal = vec; }

	const ARM_GP_T_TridiagMatrix<T>& operator=( T const & a );

	/// For Inversion
	/// Precomputation of Inverse Matrix
	void PreComputeInverse();
	ARM_GP_T_Vector<T>& MultiplyByInverseInplace( ARM_GP_T_Vector<T>& rhs ) const;
	ARM_GP_T_Vector<T> MultiplyByInverse( ARM_GP_T_Vector<T>& rhs ) const;

	/// Mutliplication By a Vector
	ARM_GP_T_Vector<T> operator*( const ARM_GP_T_Vector<T>& rhs ) const;
	ARM_GP_T_Vector<T>& operator*=( ARM_GP_T_Vector<T>& rhs ) const;

	/// Operations with matrixes
	ARM_GP_T_TridiagMatrix<T> operator*( const ARM_GP_T_TridiagMatrix<T>& rhs );
	ARM_GP_T_TridiagMatrix<T>& operator*=( const ARM_GP_T_TridiagMatrix<T>& rhs );
	ARM_GP_T_TridiagMatrix<T> operator+( const ARM_GP_T_TridiagMatrix<T>& rhs );
	ARM_GP_T_TridiagMatrix<T>& operator+=( const ARM_GP_T_TridiagMatrix<T>& rhs );
	ARM_GP_T_TridiagMatrix<T> operator-( const ARM_GP_T_TridiagMatrix<T>& rhs );
	ARM_GP_T_TridiagMatrix<T>& operator-=( const ARM_GP_T_TridiagMatrix<T>& rhs );

};

template <typename T> 
ARM_GP_T_TridiagMatrix<T>& ARM_GP_T_TridiagMatrix<T>::CopyNoInverse( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	if( itsSize != rhs.itsSize )
	{
		itsSize = rhs.itsSize;
		itsDiagonal = rhs.itsDiagonal;
		itsUpperDiag = rhs.itsUpperDiag;
		itsLowerDiag = rhs.itsLowerDiag;
	}
	else
	{
		vector<T>::iterator iterThis, iterRhs, iterEnd;

		iterThis = itsDiagonal.begin();
		iterRhs = rhs.itsDiagonal.begin();
		iterRhsEnd = rhs.itsDiagonal.end();

		for( itsDiagonal.begin() ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
			(*iterThis) = (*iterRhs);

		iterThis = itsUpperDiag.begin();
		iterRhs = rhs.itsUpperDiag.begin();
		iterRhsEnd = rhs.itsUpperDiag.end();

		for( ; iterRhs !=!iterRhsEnd; ++iter, ++iterRhs )
			(*iterThis) = (*iterRhs);

		iterThis = itsLowerDiag.begin();
		iterRhs = rhs.itsLowerDiag.begin();
		iterRhsEnd = rhs.itsLowerDiag.end();

		for( ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
			(*iterThis) = (*iterRhs);
	}

	return (*this);
}

template <typename T> 
const ARM_GP_T_TridiagMatrix<T>& ARM_GP_T_TridiagMatrix<T>::operator=( T const & a )
{
	itsDiagonal.assign(itsDiagonal.size(),a);
	itsUpperDiag.assign(itsUpperDiag.size(),a);
	itsLowerDiag.assign(itsLowerDiag.size(),a);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////:
/// PreComputeInverse
//////////////////////////////////////////////////////////////////////////////////////////////////////:

template <typename T>
void ARM_GP_T_TridiagMatrix<T>::PreComputeInverse()
{
	if( itsNormalizationTerms.size() != itsSize )
		itsNormalizationTerms = vector<T>( itsSize );

	if( itsOtherCoeff.size() != itsSize-1 )
		itsOtherCoeff = vector<T>( itsSize-1 );

	size_t i;
	T Normalisation = itsDiagonal[0];
	T OtherCoeff = itsUpperDiag[0] / Normalisation;

	itsNormalizationTerms[0] = Normalisation;
	itsOtherCoeff[0] = OtherCoeff;

	for( i=1 ; i<itsSize-1 ; ++i )
	{
		Normalisation = itsDiagonal[i]-itsLowerDiag[i-1]*OtherCoeff;
		OtherCoeff = itsUpperDiag[i]/Normalisation;
		itsNormalizationTerms[i] = Normalisation;
		itsOtherCoeff[i] = OtherCoeff;
	}

	Normalisation = itsDiagonal[itsSize-1]-itsLowerDiag[itsSize-2]*OtherCoeff;
	itsNormalizationTerms[itsSize-1] = Normalisation;
}

template <typename T>
ARM_GP_T_Vector<T>& ARM_GP_T_TridiagMatrix<T>::MultiplyByInverseInplace( ARM_GP_T_Vector<T>& rhs ) const
{

#if defined(__GP_STRICT_VALIDATION)
	if( itsNormalizationTerms.size() != itsSize ||
		itsOtherCoeff.size() != itsSize-1 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "Inverse not precomputed!");
#endif

	size_t i;
	T CurrentState = vec[0]/itsNormalizationTerms[0];

	vec[0] = CurrentState;

	for( i=1 ; i< itsSize ; ++i)
	{
		CurrentState = (vec[i]-LowerTerm[i-1]*CurrentState)/itsNormalizationTerms[i];
		vec[i] = CurrentState;
	}

	for( i = itsSize-2 ; i>0 ; --i )
	{
		CurrentState = vec[i]-itsOtherCoeff[i]*CurrentState;
		vec[i] = CurrentState;
	}

	vec[0] = vec[0]-itsOtherCoeff[0]*CurrentState;
}

template <typename T>
ARM_GP_T_Vector<T> ARM_GP_T_TridiagMatrix<T>::MultiplyByInverse( ARM_GP_T_Vector<T>& rhs ) const
{
	ARM_GP_T_Vector<T> returnedVector( rhs );
	MultiplyByInverseInplace( returnedVector );
	return returnedVector;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////:
//// Basic Multiplication with vectors
//////////////////////////////////////////////////////////////////////////////////////////////////////:

template <typename T>
ARM_GP_T_Vector<T> ARM_GP_T_TridiagMatrix<T>::operator*( const ARM_GP_T_Vector<T>& rhs ) const
{
	ARM_GP_T_Vector<T> vec( rhs );
	return 	(*this)*=vec;
}

template <typename T>
ARM_GP_T_Vector<T>& ARM_GP_T_TridiagMatrix<T>::operator*=( ARM_GP_T_Vector<T>& rhs ) const
{
	size_t i = 0;
	T N,NM1,NP1;
	N=rhs[0];
	NP1=rhs[1];
	NM1=0;

	rhs[0] = itsDiagonal[0]*N + itsUpperDiag[0] * NP1;
	NM1=N;
	N=NP1;

	for( i=1; i<itsSize-1; ++i)
	{
		NP1 = rhs[i+1];
		rhs[i] = itsLowerDiag[i-1]*NM1 + itsDiagonal[i]*N + itsUpperDiag[i] * NP1;
		NM1 = N;
		N=NP1;
	}
	rhs[itsSize-1] = itsLowerDiag[itsSize-2]*NM1 + itsDiagonal[itsSize-1]*N;

	return rhs;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////:
//// Basic Operations with matrixes
//////////////////////////////////////////////////////////////////////////////////////////////////////:

template <typename T>
ARM_GP_T_TridiagMatrix<T> ARM_GP_T_TridiagMatrix<T>::operator*( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "unimplented *");
}


template <typename T>
ARM_GP_T_TridiagMatrix<T>& ARM_GP_T_TridiagMatrix<T>::operator*=( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "unimplented *=");
}

template <typename T>
ARM_GP_T_TridiagMatrix<T>& ARM_GP_T_TridiagMatrix<T>::operator+=( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	vector<T>::iterator iterThis, iterRhs, iterEnd;

	iterThis = itsDiagonal.begin();
	iterRhs = rhs.itsDiagonal.begin();
	iterRhsEnd = rhs.itsDiagonal.end();

	for( itsDiagonal.begin() ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
		(*iterThis) += (*iterRhs);

	iterThis = itsUpperDiag.begin();
	iterRhs = rhs.itsUpperDiag.begin();
	iterRhsEnd = rhs.itsUpperDiag.end();

	for( ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
		(*iterThis) += (*iterRhs);

	iterThis = itsLowerDiag.begin();
	iterRhs = rhs.itsLowerDiag.begin();
	iterRhsEnd = rhs.itsLowerDiag.end();

	for( ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
		(*iterThis) += (*iterRhs);

	return (*this);

}

template <typename T>
ARM_GP_T_TridiagMatrix<T>& ARM_GP_T_TridiagMatrix<T>::operator-=( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	vector<T>::iterator iterThis, iterRhs, iterEnd;

	iterThis = itsDiagonal.begin();
	iterRhs = rhs.itsDiagonal.begin();
	iterRhsEnd = rhs.itsDiagonal.end();

	for( itsDiagonal.begin() ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
		(*iterThis) -= (*iterRhs);

	iterThis = itsUpperDiag.begin();
	iterRhs = rhs.itsUpperDiag.begin();
	iterRhsEnd = rhs.itsUpperDiag.end();

	for( ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
		(*iterThis) -= (*iterRhs);

	iterThis = itsLowerDiag.begin();
	iterRhs = rhs.itsLowerDiag.begin();
	iterRhsEnd = rhs.itsLowerDiag.end();

	for( ; iterRhs != iterRhsEnd; ++iter, ++iterRhs )
		(*iterThis) -= (*iterRhs);

	return (*this);

}

template <typename T>
ARM_GP_T_TridiagMatrix<T> ARM_GP_T_TridiagMatrix<T>::operator+( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	ARM_GP_T_TridiagMatrix<T> returnedMatrix(rhs.itsSize);
	returnedMatrix.CopyNoInverse( rhs );
	return returnedMatrix+=rhs;

}

template <typename T>
ARM_GP_T_TridiagMatrix<T> ARM_GP_T_TridiagMatrix<T>::operator-( const ARM_GP_T_TridiagMatrix<T>& rhs )
{
	ARM_GP_T_TridiagMatrix<T> returnedMatrix(rhs.itsSize);
	returnedMatrix.CopyNoInverse( rhs );
	return returnedMatrix-=rhs;

}


/// toString
template <typename T> string ARM_GP_T_TridiagMatrix<T>::toString(const string& indent, const string& nextIndent) const
{
	size_t i;
    CC_Ostringstream os;
	os << "GP_TRIDIAG_MATRIX[" << itsSize << "]\n";

	os << "UpperTerms: \n";
	for(i=0; i<itsSize-1; ++i)
	{
			os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
				<< CC_NS(std,setw)(8) << itsUpperDiag[i] << "\t";
	}

	os << "\n Diagonal Terms : \n";

	for(i=0; i<itsSize; ++i)
	{
			os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
				<< CC_NS(std,setw)(8) << itsDiagonal[i] << "\t";
	}

	os << "\n LowerTerms : \n";

	for(i=0; i<itsSize-1; ++i)
	{
			os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
				<< CC_NS(std,setw)(8) << itsLowerDiag[i] << "\t";
	}
	os << "\n";

    return os.str();
}

CC_END_NAMESPACE()

#endif