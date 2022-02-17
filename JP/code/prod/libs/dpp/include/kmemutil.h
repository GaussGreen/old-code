/***************************************************************
 * Module:	PenGuin
 * File:	kmemutil.h
 * Function:	Standard Definition and Include files.
 * Author:	Christian Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kmemutil_H
#define	_kmemutil_H

#include "kstdinc.h"


//--------------------------------------------------------------

/**
 * A template class for garbage collection.
 * Used to allocate dynamically objects (using the member function
 * New) without having to keep track of the returned pointer
 * to delete.
 * The pointers can be freed either one by one by calling the Free(T*)
 * method (and ONLY this way, i.e. do not call directly delete) or
 * all at once by calling Free().
 * Caution: the destructor ~KGCollect frees all pointers allocated,
 * so the pointers allocated will not be available
 * outside the scope of the KGCollect object.
 */

template <class T> class KGCollect {
public:
	/**
	 * Destructor.
	 */
	~KGCollect()
	{
		Free();
	}

	/**
	 * Allocates and returns a new object (default constructor).
	 */
	T* New()
	{
		T *ptr = new T;
		if (ptr == NULL)
			throw KFailure("KGCollect::new: failed.\n");
		mPtr.insert(mPtr.end(), ptr);
		return(ptr);
	}

	/**
	 * Deletes a specified object pointer.
	 */
	void Free(T *ptr)
	{
		KGCTable::iterator p = mPtr.find(ptr);
		if (p == mPtr.end()) {
			throw KFailure("KGCollect::Free: pointer %p "
				"not allocated with GC.\n", ptr);
		}
		delete ptr;
		mPtr.erase(p);
	}

	/**
	 * Delete all memory.
	 */
	void Free()
	{
		for (KGCTable::iterator p = mPtr.begin();
		     p != mPtr.end(); ++p) 
		{
			delete *p;
		}
		mPtr.clear();
	}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KGCollect& object)
	{
		int	idx = 0;
		os << format("KGCollect:\n") << endl;

		for (KGCTable::iterator p = object.mPtr.begin();
		     p != object.mPtr.end();
		     ++p) {
			os << format(" %4d  %p", idx++, *p) << " `" <<
				*(*p) << "' " << endl;
		}
		return(os);
	}


private:
	/*
	 * Just a set of pointers.
	 */
	typedef	KVector(T*)	KGCTable;

	KGCTable	mPtr;
};



//--------------------------------------------------------------

/**
 * Template function to copy an array of length n of type T.
 * The new copied array must be freed using the delete operator.
 */

template <class T> T*
DppNewArrayCopy(int n, const T* array)
{
	if (n == 0) return(NULL);

	T	*v = new T[n];
	if (v == NULL)
		throw KFailure("new failed (size %d).\n", n);

        for (int i=0; i<n; i++) v[i] = array[i];

        return(v);
}
 




#endif


