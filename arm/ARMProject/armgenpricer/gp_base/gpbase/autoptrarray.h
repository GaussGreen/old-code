/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * %comment%
 *
 * $Log: autoptrarray.h,v $
 * Revision 1.1  2003/10/08 16:45:06  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file autoptrarray.h
 *
 *  \brief class to handle array in an auto_ptr way
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _INGPBASE_AUTOPTRARRAY_H
#define _INGPBASE_AUTOPTRARRAY_H

/// \class auto_ptr_array
/// this is the analogy of the auto_ptr but for the array case
/// \brief this has to be used only with arrays!
/// if this is not the case, use auto_ptr instead!!!


/// use our macro for namespace
#include "port.h"

CC_BEGIN_NAMESPACE( std )

/// auxiliary type to enable copies and assignments (now global)
/// to provide easy access
template<class Y> 
struct auto_ptr_array_ref 
{
	/// refer to something new Y[n] but there is no difference between this and new Y!
	/// the type is Y* indeed 
	Y* yp; 
	auto_ptr_array_ref (Y* rhs): yp(rhs) {}
};


/// \class auto_ptr_array
/// \brief the real class that implements ownership mechanism
/// on array of pointor
/// to be compliant with the ISO C++ standard
/// function do not throw hence the nothrow signature (throw() )
template<class T> 
class auto_ptr_array 
{
private:
	T* ap;    ///refers to the actual owned object (if any)
public:
	/// like the auto_ptr
	typedef T element_type;
	
	///constructor has to be explicit like in the auto_ptr case
	explicit auto_ptr_array (T* ptr = 0) throw()
		: ap(ptr) {}
	
	///copy constructors (with implicit conversion)
	///- note: nonconstant parameter
	auto_ptr_array (auto_ptr_array& rhs) throw()
		: ap(rhs.release()) {}

	
	////VC6 does not compile this
	/// hence commented out but should really be there
	/// the template parameter enables very large conversion
	////template<class Y>
	///	auto_ptr_array (auto_ptr_array<Y>& rhs) throw()
	///	: ap(rhs.release()) {}
	

	///assignments (with implicit conversion)
	///- note: nonconstant parameter
	auto_ptr_array& operator= (auto_ptr_array& rhs) throw() 
	{
		reset(rhs.release());
		return *this;
	}
	

	/// VC6 does not compile this
	/// same issue as above

	/// template<class Y>
	///	auto_ptr_array& operator= (auto_ptr_array<Y>& rhs) throw() 
	///{
	///	reset(rhs.release());
	///	return *this;
	///}
	
	/// destructor
	/// it automatically deletes the object
	/// the difference with the auto_ptr lies here
	/// auto_ptr should use delete
	/// while auto_ptr_array uses delete[]!
	~auto_ptr_array() throw() 
	{
		delete[] ap;
	}
	
	///////////////////////////////
	/// value access functions
	///////////////////////////////
	T* get() const throw() 
	{
		return ap;
	}
	
	////the const accessor safe version
	const T& operator[]( int i ) const throw() 
	{
		return ap[i];
	}

	/// this accessor is dangerous!!
	/// as it provides write mode access
	/// handle with care
	T& operator[]( int i ) throw() 
	{
		return ap[i];
	}
	
	/// release ownership
	/// makes our pointor to NULL
	/// and return the pointor
	T* release() throw() 
	{
		T* tmp = ap;
		ap = NULL;
		return tmp;
	}
	
	/// reset value
	/// reset deletes the pointor (with delete[] because of the
	/// array nature
	/// and reinitialises the pointor to the new value!
	void reset (T* ptr=0) throw() 
	{
		if (ap != ptr) 
		{
			delete[] ap;
			ap = ptr;
		}
	}
	
	/// special conversions with auxiliary type to enable copies and assignments
	auto_ptr_array(auto_ptr_array_ref<T> rhs) throw()
		: ap(rhs.yp) {}

	auto_ptr_array& operator= (auto_ptr_array_ref<T> rhs) throw() 
	{  ///new
		reset(rhs.yp);
		return *this;
	}
	
	
	
	/// Note that this should normally be part
	/// template<class Y>
	/// operator auto_ptr_array_ref<Y>() throw() 
	///{
	///	return auto_ptr_array_ref<Y>(release());
	///}
	
	///template<class Y>
	///	operator auto_ptr_array<Y>() throw() 
	///{
	///	return auto_ptr_array<Y>(release());
	///}
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
