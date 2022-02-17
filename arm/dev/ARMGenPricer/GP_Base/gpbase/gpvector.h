/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gpvector.h
 *
 *  \brief a STL wrapped up vector
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPINFRA_GPVECTOR_H
#define _INGPINFRA_GPVECTOR_H

#include "removeidentifiedwarning.h"
#include "port.h"
#include "env.h"
#include "ostringstream.h"
#include "rootobject.h"

#include <vector>
CC_USING_NS(std,vector)
#include <string>
CC_USING_NS(std,string)
#include <algorithm>
#include <ostream>
using std::ostream;

#include <iomanip>

/// kernel
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

#if defined(__ARM_VECTOR_NO_RANGE_CHECK)
	#define CHECKRANGE(i)  
#else
	#define CHECKRANGE(i) CheckRange(i);
#endif

/// forward declaration
template <typename T> class ARM_GP_T_Matrix;
template <typename T> class ARM_GP_T_Tensor;


///	Super macro for opeartion
///		-UNARYOP

#define UNARYOP_ONEVAL(op)												\
	iterator iter;											\
	for( iter=itsValues.begin(); iter!=itsValues.end(); ++iter )		\
		*iter op val;													\
	return *this;

#define UNARYOP_VECTOR(op)																\
	CC_NS(ARM_Check,CheckSameArgSize)( *this, rhs, "lhs", "rhs" );						\
	iterator iter;																		\
	const_iterator rhsIter;																\
	for( iter=begin(), rhsIter=rhs.begin(); iter!=itsValues.end(); ++iter, ++rhsIter )	\
		*iter op *rhsIter;																\
	return *this;

//////////////////////////////////////////////
/// \class ARM_GP_T_Vector 
/// \brief template vector class based on the STL vector ... 
///			this wrapper around the STL vector allows us to
///			retrieve all the STL vector interface
//////////////////////////////////////////////
template <typename T> class ARM_GP_T_Vector : public ARM_RootObject
{
public:
	/// standard typedefs
	typedef T value_type;
    typedef ARM_GP_T_Vector<T>  GP_T_Vector;

// FIXMEFRED: mig.vc8 (21/05/2007 10:54:04): typedef typename
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef typename CC_NS(std,vector)<T>::iterator				iterator;
	typedef typename CC_NS(std,vector)<T>::const_iterator		const_iterator;   
#else					// Visual C++ 6
	typedef CC_NS(std,vector)<T>::iterator iterator;
	typedef CC_NS(std,vector)<T>::const_iterator const_iterator;
#endif 

    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef T const & const_reference;
    typedef size_t size_type;

	/// constructors
    inline explicit ARM_GP_T_Vector (size_t size=0, const T& initValue = T())
	:	ARM_RootObject(), itsValues(size,initValue){ SetName(ARM_GP_VECTOR);}
    inline explicit ARM_GP_T_Vector (size_t size, T* values)
	:	ARM_RootObject(), itsValues(values,values+size){ SetName(ARM_GP_VECTOR); }
    inline explicit ARM_GP_T_Vector(T* begin, T* end)
	:	ARM_RootObject(), itsValues(begin,end){ SetName(ARM_GP_VECTOR); }

// FIXMEFRED: mig.vc8 (28/05/2007 11:08:38): added to help migration
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	inline explicit ARM_GP_T_Vector (size_t size, const_iterator values)
	:	ARM_RootObject(), itsValues(&(*values),&(*values)+size){ SetName(ARM_GP_VECTOR); }
	inline explicit ARM_GP_T_Vector(const_iterator begin, const_iterator end)
	:	ARM_RootObject(), itsValues(&(*begin),&(*end)){ SetName(ARM_GP_VECTOR); }
    inline explicit ARM_GP_T_Vector (size_t size, iterator values)
	:	ARM_RootObject(), itsValues(&(*values),&(*values)+size){ SetName(ARM_GP_VECTOR); }
	inline explicit ARM_GP_T_Vector(iterator begin, iterator end)
	:	ARM_RootObject(), itsValues(&(*begin),&(*end)){ SetName(ARM_GP_VECTOR); }
#endif

	/// conversion from a std vectors
	inline explicit ARM_GP_T_Vector( const vector<T>& stdVec )
	:	ARM_RootObject(), itsValues(stdVec){ SetName(ARM_GP_VECTOR); }

	/// conversion from matrix and tensor! FIX FIX FIX (to be done)
	inline explicit ARM_GP_T_Vector(ARM_GP_T_Matrix<T> const& matrix ) {SetName(ARM_GP_VECTOR);};
	inline explicit ARM_GP_T_Vector(ARM_GP_T_Tensor<T> const& tensor ) {SetName(ARM_GP_VECTOR);};

	/// copy constructor
	inline ARM_GP_T_Vector(ARM_GP_T_Vector<T> const& rhs )
	:	ARM_RootObject(rhs), itsValues( rhs.itsValues ) {};
	/// assignment operator
	inline ARM_GP_T_Vector<T>& operator=( const ARM_GP_T_Vector<T>& rhs )
	{	
		if( this != &rhs )
		{
			ARM_RootObject::operator =(rhs);
			itsValues = rhs.itsValues;
		}
		return *this;
	}

	/// destructor
    virtual ~ARM_GP_T_Vector(){};
	/// assignment to a constant of type T
	const ARM_GP_T_Vector<T>& operator=( T const & a )
	{
		itsValues.assign(itsValues.size(),a);
		return *this;
	}

	/// traditional accessor!
    inline const T& operator()(size_t i) const	{ CHECKRANGE(i); return itsValues[i];}
    inline T& operator()(size_t i)				{ CHECKRANGE(i); return itsValues[i];}
    inline const T& operator[](size_t i) const	{ CHECKRANGE(i); return itsValues[i];}
    inline T& operator[](size_t i)				{ CHECKRANGE(i); return itsValues[i];}
    inline const T& Elt(size_t i) const			{ CHECKRANGE(i); return itsValues[i];}
    inline T& Elt(size_t i)						{ CHECKRANGE(i); return itsValues[i];}

	/// standard STD Vector function
	/// size and memory
	inline size_t size() const	{ return itsValues.size(); }
	inline bool empty() const { return itsValues.empty(); }
	inline void reserve(size_t size) { itsValues.reserve(size); }

	/// assignment
	inline void assign(size_t size, T value= T()) { itsValues.assign(size,value); }
	inline void assign(T* begin,T* end ) { itsValues.assign(begin,end); }
	inline void swap(ARM_GP_T_Vector<T>& rhs ) {itsValues.swap(rhs.itsValues);}
	inline void swap(const ARM_GP_T_Vector<T>& rhs ) {itsValues.swap(rhs.itsValues);}

	/// insertion deletion...
	inline iterator insert( iterator pos, const T& elem ) { return itsValues.insert(pos,elem);}
	inline void insert( iterator pos, size_t n, const T& elem) { itsValues.insert(pos,n,elem);}
	inline void insert( iterator pos, iterator begin, iterator end ) { itsValues.insert(pos,begin,end);}
	inline void push_back( const T& elem ) { itsValues.push_back(elem);}
	inline void pop_back() { itsValues.pop_back();}
    inline void fill(const GP_T_Vector& T_Vector) {for(int i = 0; i < T_Vector.size(); ++i) itsValues.push_back(T_Vector[i]);}
	inline iterator erase(iterator pos ) { return itsValues.erase(pos); }
	inline iterator erase(iterator begin, iterator end) { return itsValues.erase(begin,end);}
	inline void resize(size_t size, const T& value = T()) {	itsValues.resize(size,value);}
	inline void clear() { itsValues.clear(); }

	/// use of the algorithm
	inline ARM_GP_T_Vector<T>& sort() { CC_NS(std,sort)(begin(),end()); return *this; }
	inline iterator find( const T& val) { return CC_NS(std,find)( begin(),end(),val); }
	inline bool contain( const T& val) { return  (find(val) != end () ); }
	inline ARM_GP_T_Vector<T>& unique() { iterator pos = CC_NS(std,unique)(begin(),end()); itsValues.resize(pos-begin()); return *this;}

	/// Standard ARM Object support
	virtual ARM_Object* Clone() const { return new ARM_GP_T_Vector<T>(*this);}
    virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LGVEC";}


	/// iterator support
	iterator begin() { return itsValues.begin(); }
	const_iterator begin() const { return itsValues.begin(); }
	iterator end() { return itsValues.end(); }
	const_iterator end() const { return itsValues.end(); }

	/// comparison operator
	inline bool operator==( ARM_GP_T_Vector<T> const& rhs ) const { return itsValues == rhs.itsValues;}
	inline bool operator!=( ARM_GP_T_Vector<T> const& rhs ) const { return itsValues != rhs.itsValues;}
	inline bool operator< ( ARM_GP_T_Vector<T> const& rhs ) const { return itsValues < rhs.itsValues;	}
	inline bool operator> ( ARM_GP_T_Vector<T> const& rhs ) const { return itsValues > rhs.itsValues;	}
	inline bool operator<=( ARM_GP_T_Vector<T> const& rhs ) const { return itsValues <= rhs.itsValues;}
	inline bool operator>=( ARM_GP_T_Vector<T> const& rhs ) const { return itsValues >= rhs.itsValues;}

	///  operation on a simple elem: addition, subtraction, multiplication, division
	ARM_GP_T_Vector<T>& operator+=( const T& val ){ UNARYOP_ONEVAL(+=) }
	ARM_GP_T_Vector<T>& operator-=( const T& val ){ UNARYOP_ONEVAL(-=) }
	ARM_GP_T_Vector<T>& operator*=( const T& val ){ UNARYOP_ONEVAL(*=) }
	ARM_GP_T_Vector<T>& operator/=( const T& val ){ UNARYOP_ONEVAL(/=) }

	/// Element-by-element: addition, subtraction, multiplication, division, 
	ARM_GP_T_Vector<T>& operator+=( ARM_GP_T_Vector<T> const& rhs ){ UNARYOP_VECTOR(+=) }
	ARM_GP_T_Vector<T>& operator-=( ARM_GP_T_Vector<T> const& rhs ){ UNARYOP_VECTOR(-=) }
	ARM_GP_T_Vector<T>& operator*=( ARM_GP_T_Vector<T> const& rhs ){ UNARYOP_VECTOR(*=) }
	ARM_GP_T_Vector<T>& operator/=( ARM_GP_T_Vector<T> const& rhs ){ UNARYOP_VECTOR(/=) }

    inline double sum() const{ double summ=0.0;for( int i =0; i<itsValues.size(); ++i )	summ+=itsValues[i]; return summ;}
    ///Accessors
    inline const vector<T>& GetValues() const {return itsValues;}

private:
    vector<T> itsValues;
	void CheckRange(size_t i) const
	{
#if !defined( __ARM_VECTOR_NO_RANGE_CHECK )
		if(i<0 || i>itsValues.size()-1)
		{
			char msg[255];
			sprintf( msg, "Range error in ARM_GP_T_Vector: index %d  not in the range 0 %d ", i, itsValues.size()-1 );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, msg );
		}
#endif
	}
};


/// standard comparison operator (symmetric version)
template <typename T> bool operator==( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ){	return lhs.operator==(rhs); }
template <typename T> bool operator!=( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ){	return lhs.operator!=(rhs); }
template <typename T> bool operator<( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ){	return lhs.operator<(rhs); }
template <typename T> bool operator>( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ){	return lhs.operator>(rhs); }
template <typename T> bool operator<=( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ){	return lhs.operator<=(rhs); }
template <typename T> bool operator>=( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ){	return lhs.operator>=(rhs); }
template <typename T> void swap( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ) { lhs.swap(rhs);}

/// standard addition, subtraction, multiplication, division, 
template <typename T> ARM_GP_T_Vector<T> operator+( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(lhs) += rhs; }
template <typename T> ARM_GP_T_Vector<T> operator+( ARM_GP_T_Vector<T> const& lhs, const T& val ) { return ARM_GP_T_Vector<T>(lhs) += val; }
template <typename T> ARM_GP_T_Vector<T> operator+( const T& val, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(rhs.size(),val) += rhs; }

template <typename T> ARM_GP_T_Vector<T> operator-( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(lhs) -= rhs; }
template <typename T> ARM_GP_T_Vector<T> operator-( ARM_GP_T_Vector<T> const& lhs, const T& val ) { return ARM_GP_T_Vector<T>(lhs) -= val; }
template <typename T> ARM_GP_T_Vector<T> operator-( const T& val, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(rhs.size(),val) -= rhs; }

template <typename T> ARM_GP_T_Vector<T> operator*( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(lhs) *= rhs; }
template <typename T> ARM_GP_T_Vector<T> operator*( ARM_GP_T_Vector<T> const& lhs, const T& val ) { return ARM_GP_T_Vector<T>(lhs) *= val; }
template <typename T> ARM_GP_T_Vector<T> operator*( const T& val, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(rhs.size(),val) *= rhs; }

template <typename T> ARM_GP_T_Vector<T> operator/( ARM_GP_T_Vector<T> const& lhs, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(lhs) /= rhs; }
template <typename T> ARM_GP_T_Vector<T> operator/( ARM_GP_T_Vector<T> const& lhs, const T& val ) { return ARM_GP_T_Vector<T>(lhs) /= val; }
template <typename T> ARM_GP_T_Vector<T> operator/( const T& val, ARM_GP_T_Vector<T> const& rhs ) { return ARM_GP_T_Vector<T>(rhs.size(),val) /= rhs; }

/// toString
template <typename T> string ARM_GP_T_Vector<T>::toString(const string& indent, const string& nextIndent) const
{ 
    CC_Ostringstream os;
	os << "ARM_GP_T_VECTOR[" << itsValues.size() << "]={";

	if(itsValues.size())
	{
		int i=0;
		for(; i<itsValues.size()-1; ++i)
			os	<< CC_NS(std,setw)(20)<< CC_NS(std,scientific) << CC_NS(std,setprecision)(10) << itsValues[i] << ", ";
		os << CC_NS(std,setw)(20)<< CC_NS(std,scientific) << CC_NS(std,setprecision)(10) << itsValues[i];
	}
	os << "}\n";
    return os.str();
}


template <typename T> ostream& operator<< ( ostream& os, const ARM_GP_T_Vector<T>& vec )
{
	os << vec.toString();
	return os;
}

/// undef MACRO
#undef UNARYOP_ONEVAL
#undef UNARYOP_VECTOR
#undef CHECKRANGE


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
