/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: countedptr.h,v $
 * Revision 1.2  2003/13/08 16:47:44  ebenhamou
 * Reformating
 *
 * Revision 1.1  2003/10/08 16:47:44  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file countedptr.h
 *
 *  \brief class to provide reference counted pointor
 *	
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _INGPBASE_COUNTEDPTR_H
#define _INGPBASE_COUNTEDPTR_H


///////////////////////////////////////////////////////////
/// \class ARM_Counted_Ptr
/// this is a simple reference counted templated smart pointor
/// \brief Because Loki is not so much supported by VC6
/// I prefered this simple implementation!!!
/// we may want to move to Loki implementation when using
/// smarter compiler
///////////////////////////////////////////////////////////


/// use our macro for namespace
#include "port.h"
/// this header comes firts as it includes some preprocessor constants!
#include "removeidentifiedwarning.h"
#include "env.h"

#include <ostream>
using std::ostream;

CC_BEGIN_NAMESPACE( ARM )



///////////////////////////////////////////////////////////
/// \class for counted reference semantics
/// deletes the object to which it refers when the last ARM_CountedPtr
/// that refers to it is destroyed
///////////////////////////////////////////////////////////
 template <class T> class ARM_CountedPtr;
    
    class ARM_CountedPtrCopier {
      public:
        template <class T, class U>
        static void copy(const ARM_CountedPtr<T>& from, const ARM_CountedPtr<U>& to) {
            if (from.itsCount != to.itsCount) {
                // if to was the last reference to its object...
                //to.dispose() ;
				if (--(*to.itsCount) == 0) {
                    // ...delete the latter and the counter
                    if (to.itsPtr != 0 )
                        delete to.itsPtr;
                    delete to.itsCount;
                }
                // cast to the new type - the resulting pointer will
                // be null if the two types are not compatible
                to.itsPtr  = dynamic_cast<U*>(from.itsPtr);
                to.itsCount    = from.itsCount;
                (*to.itsCount)++;
            }
        }
    };

template <class T>
class ARM_CountedPtr 
{
	friend class ARM_CountedPtrCopier;
private:
	mutable T* itsPtr;        /// pointer to the value
	mutable long* itsCount;   /// shared number of owners
	
	/// dispose is THE function that is responsible
	/// for deleting if necessary
	inline void dispose() 
	{
#if defined(__SET_PTR_TO_NULLL)
		/// reseting pointor to null is better to get crashes...but it should not be necessary
		if (--*itsCount == 0)
		{
			delete itsCount;
			itsCount = 0;
			if (itsPtr != 0) {
				delete itsPtr;
				itsPtr = 0;
			}
		}
#else
		if (--*itsCount == 0)
		{
			delete itsCount;
			delete itsPtr;
		}
#endif
	}

public:
	/// initialize pointer with existing pointer
	/// - requires that the pointer p is a return value of new
	explicit ARM_CountedPtr (T* p=0)
		: itsPtr(p), itsCount(new long(1)) 
	{}
	
	/// copy pointer (one more owner)

	template <class U>
	inline ARM_CountedPtr (const ARM_CountedPtr<U>& p) throw()
		: itsPtr(NULL), itsCount(new long(1)) 
	{
		ARM_CountedPtrCopier::copy(p,*this);
		//++*itsCount;
	}


	inline ARM_CountedPtr (const ARM_CountedPtr<T>& p) throw()
		: itsPtr(p.itsPtr), itsCount(p.itsCount) 
	{
		//ARM_CountedPtrCopier::copy(p,*this);
		++*itsCount;
	}	

	/// destructor (delete value if this was the last owner)
	inline ~ARM_CountedPtr () throw() 
	{
		dispose();
	}
	
	/// assignment (unshare old and share new value)
	template <class U>
	inline ARM_CountedPtr<T>& operator= (const ARM_CountedPtr<U>& p) throw() 
	{

		ARM_CountedPtrCopier::copy(p,*this);
        return *this;
		/*		
		if (this != &p) 
		{
			dispose();
			itsPtr = p.itsPtr;
			itsCount = p.itsCount;
			++*itsCount;
		}
		return *this;
		*/
	}
	inline ARM_CountedPtr<T>& operator= (const ARM_CountedPtr<T>& p) throw() 
	{	
		if (this != &p) 
		{
			dispose();
			itsPtr = p.itsPtr;
			itsCount = p.itsCount;
			++*itsCount;
		}
		return *this;	
	}
	/// accessors
	/// everything is const to avoid overwritting
	/// access the value to which the pointer refers
	/// operator*
	inline T& operator*() const throw() 
    {
		return *itsPtr;
	}
	
	/// operator ->
	inline T* operator->() const throw() 
    {
		return itsPtr;
	}

    inline bool IsNull(void) const
    {
        return( itsPtr == NULL );
    }
};

/// some comparison functions
/// equal comparison
template <class T>
    inline 	bool operator==( const ARM_CountedPtr<T>& lhs, const ARM_CountedPtr<T>& rhs )
{
	/// use operator() to avoid deferencing
	return lhs.operator->() == rhs.operator->();
}

template <class T>
    inline 	bool operator==( const ARM_CountedPtr<T>& lhs, T* rhs )
{
	/// use operator() to avoid deferencing
	return lhs.operator->() == rhs;
}


/// not equal comparions function
template <class T>
	inline 	bool operator!=( const ARM_CountedPtr<T>& lhs, const ARM_CountedPtr<T>& rhs )
{	return !operator==( lhs, rhs ); }


template <class T>
	inline 	bool operator!=( const ARM_CountedPtr<T>& lhs, T* rhs )
{	return !operator==( lhs, rhs ); }


/// dump operator
template <class T> ostream& operator<< ( ostream& os, const ARM_CountedPtr<T>& rhs )
{
	os << (*rhs);
	return os;
}


CC_END_NAMESPACE()


#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
