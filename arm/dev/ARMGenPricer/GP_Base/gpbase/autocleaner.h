/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: autocleaner.h,v $
 * Revision 1.1  2003/10/24 16:45:06  ebenhamou
 * Initial revision
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file autocleaner.h
 *
 *  \brief In an ideal world with only smart pointor
 *		autocleaner should not exist! But for legacy reason,
 *		autocleaner allow user to hold pointor and make object
 *		exception safe... The standard autocleaner has a subset
 *		of functionalities of an autopointor.
 *		The main interest of autocleaner is therefore more in:
 *			- its vector version autocleanervector and 
 *			- its array version autocleanerarray
 * 		
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

/*----------------------------------------------------------------------------*/

#ifndef _INGPBASE_AUTOCLEANER_H
#define _INGPBASE_AUTOCLEANER_H


/// use our macro for namespace
#include "port.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )



/// \class ARM_AutoCleaner
/// \brief functionality is to
///		hold an object and delete
///		it when going out of the scope
///		unless being told to release it

template <typename T>
	class ARM_AutoCleaner
{
public:
	ARM_AutoCleaner( T* t):itsPtr(t), released(false) {}
	~ARM_AutoCleaner()
	{
		if(!released)
		{
			delete itsPtr;
			itsPtr = NULL;
		}
	}
	void Release(){ released=true;}

	void setPtr(T* ptr) {
		itsPtr = ptr;
		released = false;
	}
private:
	T* itsPtr;
	bool released;
};



/// \class ARM_AutoCleanerArray
/// \brief functionality is to
///		hold an object and delete
///		it when going out of the scope
///		unless being told to release it

template <typename T>
	class ARM_AutoCleanerArray
{
public:
	ARM_AutoCleanerArray( T* t ):itsPtr(t), released(false) {}
	~ARM_AutoCleanerArray() {
		if(!released)
		{
			delete[] itsPtr;
			itsPtr = NULL;
		}
	}
	void Release(){ released=true;}
private:
	T* itsPtr;
	bool released;
};


/// \class ARM_AutoCleanerVector
/// \brief note that we gives a pointor
///		to the vector to avoid expensive
///		copy!

template <typename T>
	class ARM_AutoCleanerVector
{
public:
	ARM_AutoCleanerVector( const vector<T*>& t) : itsVecPtr(t), released(false) {}
	~ARM_AutoCleanerVector()
	{
		if(!released)
		{
			for( int i=0; i<itsVecPtr.size(); ++i )
				delete itsVecPtr[i];
			itsVecPtr.resize(0);
		}
	}

	T*& operator[](size_t i) { return itsVecPtr[i]; }

	void Release(){ released=true;}
private:
	vector<T*> itsVecPtr;
	bool released;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
