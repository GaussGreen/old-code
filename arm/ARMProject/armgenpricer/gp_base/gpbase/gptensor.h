/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gptensor.h
 *
 *  \brief tensor class
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPINFRA_GPTENSOR_H
#define _INGPINFRA_GPTENSOR_H

#include "port.h"
#include "env.h"
#include "rootobject.h"

CC_BEGIN_NAMESPACE( ARM )

#if defined(__ARM_LINALG_NO_RANGE_CHECK)
	#define CHECKRANGE(i)  
#else
	#define CHECKRANGE(i) CheckRange(i);
#endif


/// forward declaration
template <typename T> class ARM_GP_T_Vector;
template <typename T> class ARM_GP_T_Tensor;

//////////////////////////////////////////////
/// \class ARM_GP_T_Tensor
/// \brief template tensor class based on the STL vector ... 
///			- this wrapper around the STL vector allows us to
///				retrieve all the STL vector interface
///			- tensor in Column major (meaning that the position of row=i, col=j 
///				is i*itsColsNb+j
//////////////////////////////////////////////

template <typename T> class ARM_GP_T_Tensor: public ARM_RootObject
{
public:
	/// standard typedefs
	typedef T value_type;
	typedef CC_NS(std,vector)<T>::iterator iterator;
	typedef CC_NS(std,vector)<T>::const_iterator const_iterator;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef T const & const_reference;
    typedef size_t size_type;
	typedef ARM_GP_T_Vector<size_t> index;

	/// constructors
    inline explicit ARM_GP_T_Tensor(const index& dimension, T initValue = T())
	:	ARM_RootObject(), itsDimension(dimension), itsValues(sizeFromIndex(),initValue){ SetName(ARM_GP_TENSOR);}
   	inline explicit ARM_GP_T_Tensor(const index& dimension, T* values)
	:	ARM_RootObject(), itsDimension(dimension), itsValues(values,values+sizeFromIndex()){ SetName(ARM_GP_TENSOR); }

    inline explicit ARM_GP_T_Tensor(const index& dimension, const vector<T>& vecValues )
	:	ARM_RootObject(), itsDimension(dimension), itsValues(vecValues )
	{ CC_NS(ARM_Check,CheckArgSize)( vecValues, "vecValues", sizeFromIndex() ); SetName(ARM_GP_TENSOR); }

	/// constructor to create a subtensor from an already constructred tensor
    inline explicit ARM_GP_T_Tensor( ARM_GP_T_Tensor<T> const& tensor, const index& lowerBoundIndex, const index& upperBoundIndex );

	/// conversion from vector and matrix!
	inline explicit ARM_GP_T_Tensor(ARM_GP_T_Vector<T> const& vector)
	: itsDimension(1,vector.size()), itsValues( vector.begin(), vector.end()){ SetName(ARM_GP_TENSOR); }
	
	inline explicit ARM_GP_T_Tensor(ARM_GP_T_Matrix<T> const& matrix )
	: itsDimension(2,), itsValues( matrix.begin(), matrix.end())
	{ 
		itsDimension[0] = matrix.rows();
		itsDimension[1] = matrix.columns();
		SetName(ARM_GP_TENSOR); 
	}

	/// copy constructor
	inline ARM_GP_T_Tensor(ARM_GP_T_Tensor<T> const& rhs )
	:	ARM_RootObject(rhs), itsDimension(rhs.itsDimension), itsValues( rhs.itsValues ) {};
	
	/// assignment operator
	inline ARM_GP_T_Tensor<T>& operator=( const ARM_GP_T_Tensor<T>& rhs )
	{	
		if( this != &rhs )
		{
			ARM_RootObject::operator =(rhs);
			itsDimension= rhs.itsDimension;
			itsValues	= rhs.itsValues;
		}
		return *this;
	}

	/// destructor
    virtual ~ARM_GP_T_Tensor(){};
	/// assignment to a constant of type T
	const ARM_GP_T_Tensor<T>& operator=( T const & a )
	{
		itsValues.assign(itsValues.size(),a);
		return *this;
	}

	/// traditional accessor!
    inline T operator()(const index& i) const	{ CHECKRANGE(i); return itsValues[elemOffset(i)];}
    inline T& operator()(const index& i)		{ CHECKRANGE(i); return itsValues[elemOffset(i)];}
    inline T Elt(const index& i) const			{ CHECKRANGE(i); return itsValues[elemOffset(i)];}
    inline T& Elt(const index& i)				{ CHECKRANGE(i); return itsValues[elemOffset(i)];}

    /// Be careful, this accessors return a copy  not pointor
    /// please don't forget to free memory
    inline ARM_GP_T_Vector<T>* GetVector(size_t dimensionNb) const;
    inline ARM_GP_T_Vector<T>* GetTransVector(size_t dimensionNb) const;

	/// standard STD Vector function
	/// size and memory
	inline size_t size() const{ return itsValues.size(); }
	inline size_t sizeFromIndex() const;

	inline bool empty() const { return itsValues.empty(); }
	inline void reserve(size_t size) { itsValues.reserve(size); }
    inline size_t GetDimNb(size_t i) const {return itsDimension[i];}

	/// Standard ARM Object support
	virtual ARM_Object* Clone() const { return new ARM_GP_T_Tensor<T>(*this);}
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// iterator support
	iterator begin() { return itsValues.begin(); }
	const_iterator begin() const { return itsValues.begin(); }
	iterator end() { return itsValues.end(); }
	const_iterator end() const { return itsValues.end(); }

	/// comparison operator
	inline bool operator==( ARM_GP_T_Tensor<T> const& rhs ) { return itsDimension== rhs.itsDimension && itsValues == rhs.itsValues; }
	inline bool operator!=( ARM_GP_T_Tensor<T> const& rhs ) { return itsDimension== rhs.itsDimension || itsValues != rhs.itsValues; }

	///  operation on a simple elem: addition, subtraction, multiplication, division
	ARM_GP_T_Tensor<T>& operator+=( const T&  val ){ return GeneralTOp(val, CC_NS(std,plus)<T>() ); }
	ARM_GP_T_Tensor<T>& operator-=( const T& val ) { return GeneralTOp(val, CC_NS(std,minus)<T>() ); }
	ARM_GP_T_Tensor<T> & operator*=( const T& val ){ return GeneralTOp(val, CC_NS(std,multiplies)<T>() ); }
	ARM_GP_T_Tensor<T> & operator/=( const T& val ){ return GeneralTOp(val, CC_NS(std,divides)<T>() ); }

	// Element-by-element: addition, subtraction, multiplication, division, 
	ARM_GP_T_Tensor<T>& operator+=( ARM_GP_T_Tensor<T> const& rhs ){ return GeneralTensorOp(rhs, CC_NS(std,plus)<T>() ); }
	ARM_GP_T_Tensor<T>& operator-=( ARM_GP_T_Tensor<T> const& rhs ){ return GeneralTensorOp(rhs, CC_NS(std,minus)<T>() ); }
	ARM_GP_T_Tensor<T>& operator*=( ARM_GP_T_Tensor<T> const& rhs ){ return GeneralTensorOp(rhs, CC_NS(std,multiplies)<T>() ); }
	ARM_GP_T_Tensor<T>& operator/=( ARM_GP_T_Tensor<T> const& rhs ){ return GeneralTensorOp(rhs, CC_NS(std,divides)<T>() ); }

private:
	index		itsDimension;
    vector<T>	itsValues;

	/// function to compute the position of an index
	/// goes backward to get the dimension!
	
	void CheckRange(const index& position) const
	{
#if !defined( __ARM_LINALG_NO_RANGE_CHECK )
		if( itsDimension.size() != position.size() )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "position and dimension is of different size!");

		for( size_t i=0; i<itsDimension.size(); ++i)
			if( position[i] >= itsDimension[i] )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "out of range!");
#endif 
	}

	/// general purpose operator for T operator op
	template <typename TOp> ARM_GP_T_Tensor<T>&	GeneralTOp( const T& rhs, TOp& Op )
	{
		for( iterator iter=begin(), end=end(); iter!=end; ++iter )
			*iter = Op( *iter, rhs );
		return *this;
	}

	/// general purpose operator for T operator op on a tensor
	template <typename TOp> ARM_GP_T_Tensor<T>&	GeneralTensorOp(ARM_GP_T_Tensor<T> const& rhs, TOp& Op )
	{
		if( itsDimension()!=rhs.itsDimension ) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "the two tensors have incompatible dimensions!");

		for( iterator iter=begin(), iterator end=end(), const_iterator iterRhs=rhs.begin(); iter!=end; ++iter; ++iterRhs )
			*iter = Op( *iter, *iterRhs );
		return *this;
	}

	/// function to return the size_t from an intvector
	inline size_t elemOffset(const index& position ) const
	{
		typedef index::const_iterator iciter;
		size_t multiplier = 1, result = 0;
	
		for( iciter posElem=position.end()-1, iciter dim = itsDimension.end()-1; posElem>=position.begin(); --posElem, --dim, multiplier *= *dim )
			result += *posElem * multiplier;
		return result;
	}
};


/// standard comparison operator (symmetric version)
template <typename T> bool operator==( ARM_GP_T_Tensor<T> const& lhs, ARM_GP_T_Tensor<T> const& rhs ){	return lhs == rhs; }
template <typename T> bool operator!=( ARM_GP_T_Tensor<T> const& lhs, ARM_GP_T_Tensor<T> const& rhs ){	return lhs != rhs; }

/// standard addition, subtraction, multiplication, division, 
template <typename T> ARM_GP_T_Tensor<T> operator+( ARM_GP_T_Tensor<T> const& lhs, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(lhs) += rhs; }
template <typename T> ARM_GP_T_Tensor<T> operator+( ARM_GP_T_Tensor<T> const& lhs, const T& val ) { return ARM_GP_T_Tensor<T>(lhs) += val; }
template <typename T> ARM_GP_T_Tensor<T> operator+( const T& val, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(rhs.size(),val) += rhs; }

template <typename T> ARM_GP_T_Tensor<T> operator-( ARM_GP_T_Tensor<T> const& lhs, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(lhs) -= rhs; }
template <typename T> ARM_GP_T_Tensor<T> operator-( ARM_GP_T_Tensor<T> const& lhs, const T& val ) { return ARM_GP_T_Tensor<T>(lhs) -= val; }
template <typename T> ARM_GP_T_Tensor<T> operator-( const T& val, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(rhs.size(),val) -= rhs; }

template <typename T> ARM_GP_T_Tensor<T> operator*( ARM_GP_T_Tensor<T> const& lhs, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(lhs) *= rhs; }
template <typename T> ARM_GP_T_Tensor<T> operator*( ARM_GP_T_Tensor<T> const& lhs, const T& val ) { return ARM_GP_T_Tensor<T>(lhs) *= val; }
template <typename T> ARM_GP_T_Tensor<T> operator*( const T& val, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(rhs.size(),val) *= rhs; }

template <typename T> ARM_GP_T_Tensor<T> operator/( ARM_GP_T_Tensor<T> const& lhs, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(lhs) /= rhs; }
template <typename T> ARM_GP_T_Tensor<T> operator/( ARM_GP_T_Tensor<T> const& lhs, const T& val ) { return ARM_GP_T_Tensor<T>(lhs) /= val; }
template <typename T> ARM_GP_T_Tensor<T> operator/( const T& val, ARM_GP_T_Tensor<T> const& rhs ) { return ARM_GP_T_Tensor<T>(rhs.size(),val) /= rhs; }


/// constructor from subtensor!
template <typename T> ARM_GP_T_Tensor<T>::ARM_GP_T_Tensor( ARM_GP_T_Tensor<T> const& tensor, const index& lowerBoundIndex, const index& upperBoundIndex )
:	ARM_RootObject(), itsDimension(), itsValues()
{
	if( lowerBoundIndex.size() <  upperBoundIndex.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "lower bound and upper bound index do  not have the same size!");

	for(size_t i=0, i<lowerBoundIndex.size(); ++ i)
		if( lowerBoundIndex[i] >= upperBoundIndex[i] )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "lower bound is above upper bound index!");

	itsDimension	= index(lowerBound.size())
	index ii		= lowerBound;

	for( size_t i=0, size totalIndex = 0; i<lowerBoundIndex.size(); ++i )
	{
		for( size_t j=0, size_t k=0; j<upperBoundIndex[i]-lowerBoundIndex[i]; ++j, ++totalIndex, ++ii[i] )
			itsValues[totalIndex] = tensor(ii)
		itsDimension[i]=j;
	}
	
	SetName(ARM_GP_TENSOR);
}



/// toString
template <typename T> string ARM_GP_T_Tensor<T>::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
	os << "GP_TENSOR\n";
	os << "dimension :" << itsDimension.toString() << "\n";
	os << "values    :\n";

	for(size_t i=0, totalIndex=0; i<itsDimension.size(); ++i)
	{
		for(size_t j=0; j<itsDimension[i]; ++j, ++totalIndex)
		{
			os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
				<< CC_NS(std,setw)(8) << itsValues[totalIndex] << "\t";
		}
		os << CC_NS(std,endl);
	}
    return os.str();
}


/// size computed from the index
template <typename T> size_t ARM_GP_T_Tensor<T>::sizeFromIndex() const
{
	index::const_iterator 
		dim    = itsDimension.begin(),
		dimEnd = itsDimension.end();
	
	size_t size = 1;
	for( ; dim!= dimEnd; ++dim )
		size *= *dim;
	return Size;
}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

