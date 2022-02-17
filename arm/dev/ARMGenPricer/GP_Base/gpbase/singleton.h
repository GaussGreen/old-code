/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file singleton.h
 *
 *  \brief files to give design patterns template
 *	we remove the original design of singleton holder based on 
 *	the standard design pattern as it forced us to explicitly use
 *	a clean up method at end time... rather
 *	we use the power of extern to define a global scope variable
 *	whose destruction is automatically taken care by the compiler at
 *	the end of the program...
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date January 2004
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPBASE_SINGLETON_H
#define _INGPBASE_SINGLETON_H

#include "port.h"

CC_BEGIN_NAMESPACE( ARM )



template <typename T> class ARM_SingletonHolder
{
public:
	ARM_SingletonHolder()
	:	itspInstance(new T )
	{}
	
	T* Instance()
	{
		return itspInstance;
	}

	~ARM_SingletonHolder()
	{ 
		delete itspInstance; 
#if defined(__SET_PTR_TO_NULLL)
		itspInstance=NULL; 
#endif
	}

private:
	/// to avoid copy and assignment
	ARM_SingletonHolder(const ARM_SingletonHolder<T>& );
	ARM_SingletonHolder& operator=(const ARM_SingletonHolder<T>& );

	/// pointor to the implementation of the singleton
	T* itspInstance;
};



////////////////////////////////////////////////////////////////////////////////////
/// \class ARM_ClassicalSingletonHolder
/// \brief this is the classical singleton object design... compared to ARM_SingletonHolder
/// it forces to call CleanUp at the end!
////////////////////////////////////////////////////////////////////////////////////

template <typename T> class ARM_ClassicalSingletonHolder
{
public:
	static T* Instance(){return itspInstance;}
	static void CleanUp()
	{
		delete itspInstance;
#if defined(__SET_PTR_TO_NULLL)
		itspInstance=NULL; 
#endif
	}
private:
	/// to avoid creation and copy
	ARM_ClassicalSingletonHolder();
	ARM_ClassicalSingletonHolder(const ARM_ClassicalSingletonHolder<T>& );

	/// static pointor to make initialisation and deletion a solved problem
	static T* itspInstance;
};

template <typename T> T* ARM_ClassicalSingletonHolder<T>::itspInstance = new T;

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
