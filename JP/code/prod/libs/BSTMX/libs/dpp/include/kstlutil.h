/***************************************************************
 * Module:	PenGuin
 * File:	kstlutil.h
 * Function:	Standard Definition and Include files.
 * Author:	Christian Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kstlutil_H
#define	_kstlutil_H

#include "kstdinc.h"


/**
 * A class to cast vectors a C style arrays (pointers).
 */

template <class T> class KCArray {
public:
	/** Constructor. */
	KCArray(KVector(T)& v)
	{
		int	i, n = v.size();
		if (n) {
			mPtr = new T[n];
			for (i=0; i<n; i++)
				mPtr[i] = v[i];
		}
	}

	/** Destructor. */
	~KCArray()
	{
		delete [] mPtr;
	}

	/** Cast operator. */
	operator T*()
	{
		return mPtr;
	}



private:
	T*	mPtr;
};



//---------------------------------------------------------------
//
//

/**
 * Creates a C style array (using new) from an STL vector.
 */

template <class T> T*
DppNewFromVector(KVector(T) &v)
{
	T*	ptr = NULL;
	int	i, n = v.size();
	if (n) {
		ptr = new T[n];
		if (ptr == NULL) return(NULL);
		for (i=0; i<n; i++)
			ptr[i] = v[i];
	}
	return(ptr);
}

/**
 * Creates an STL vector form a C style array.
 */

template <class T> KVector(T)
DppVectorFromArray(int n, T* array)
{
	int	i;
	KVector(T)	v;
	for (i=0; i<n; i++)
		v.insert(v.end(), array[i]);
	return(v);
}


/**
 * Creates an STL vector from a Magnet Array.
 */

template <class T> KVector(T)
DppVectorFromCMLIBArray(const Array<T> &a)
{
	int     	i;
	KVector(T)	v;
	for (i=0; i<a.size(); i++)
		v.insert(v.end(), a[i]);
	return(v);
}


/**
 * Create TDate vector from Magnet Date array
 */
KVector(TDate)
DppDateArrayToTDateVector(const Array<Date> &a);



/**
 * Creates a string vector form a char-block style array.
 */

KVector(String) DppStrVectorFromCharBlock(int n, char* array);


/**
 * Adds an element to an STL vector.
 */

template <class T> KVector(T)
DppAddToVector(KVector(T) &v, const T& x, int removeDouble = FALSE)
{
	if (removeDouble) {
		int	i, n = v.size();
		if (x == v[i]) return(v);
	}
	v.insert(v.end(), x);

	return(v);
}





#endif


