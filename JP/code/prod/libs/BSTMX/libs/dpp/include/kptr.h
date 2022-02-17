/***************************************************************
 * Module:	DR Library C++ Utilities
 * Submodule:	
 * File:	kptr.h
 * Function:	Standard Definition and Include files.
 * Author:	Andrew Chou (revised c. Daher).
 * Revision:	$Header$
 ***************************************************************/
#ifndef _kptr_H
#define _kptr_H

#include "kstdinc.h"		// Exceptions, io, etc.



//--------------------------------------------------------------
/**
 * Reference counted pointer class.
 *
 * Performs reference counting on heap-allocated objects.
 * For each object created, we count the number of pointers to that object.
 * When this number goes to zero, the object is automatically deleted.
 * Think of it as simple garbage collecting.
 * This is NOT a true garbage collector.  There
 * are two flaws in this scheme that you should be aware of.
 * 
 * The general syntax looks like <BR>
 *	<TT> typedef KPtr<Object> ObjectPtr; </TT><BR>
 *	<TT> ObjectPtr x = new Object; </TT><BR>
 * 
 * Note: operator= only works between KPtrs and adjusts
 * the count(s) as necessary.

 * <ol>
 * <li> Circular references -- pointers that point to each other
 *      will not get destroyed.
 * <li> Inheritance -- pointers to classes cannot be up or down cast
 *      while preserving the reference counting features.
 *      We really need to count the total number of pointers
 *      pointing to an object (regardless from which point
 *      in the hierachy chain it comes from).
 * </ul>
 * 
 * Otherwise, you can use reference counted pointers,
 * EXACTLY like normal pointers, but you never
 * need to explicitly delete anything.<br>
 * 
 * Example:
 * 
 * <dir> <pre>
 *   @@typedef KPtr<Object> ObjectPtr; 
 *   @@ObjectPtr p = new Object(args ...);
 *   @@p->objectMethod ();
 *   @@(*p).objectVariable;
 * </pre></dir>
 *
 * Note that Sun-CC4.2 does not support templated operator->,
 * so you must use * to dereference.
 */

template <class T>
class KPtr {
public:
	/**
	 * Constructs from a pointer
	 */
	KPtr (T* ptr = NULL);

	/**
	 * Checks whether the pointer is NULL.
	 */
	bool null() const; 

	/**
	 * Destructor.
	 */
	~KPtr();

	/**
	 * Copy constructor.
	 */
	KPtr (const KPtr& ptr);

	/**
	 * Copy operator.
	 */
	KPtr& operator=(const KPtr& rhs);

	/**
	 * Emulates ordinary pointer operator->.
	 */
	T* operator->() const; 

	/**
	 * Emulates ordinary pointer operator*.
	 */
	T& operator*() const;  

	/**
	 * Print.
	 */
	friend ostream& operator<<(ostream& s, const KPtr& a)
	{s << a.m_ptrValue->m_ptr; return s;}

protected:
	struct KPtrValue {
		T*	m_ptr;
		int	m_refCount;

		KPtrValue(T* ptr) : m_refCount(1), m_ptr(ptr) {}

		~KPtrValue() {if (m_ptr != NULL) delete m_ptr;}
	};

	KPtrValue	*m_ptrValue;
};

template <class T> inline
KPtr<T>::KPtr (T* ptr) : m_ptrValue(new KPtrValue(ptr)) 
{
}

template <class T> inline
KPtr<T>::KPtr(const KPtr &ptr)
{
	m_ptrValue = ptr.m_ptrValue;
	++m_ptrValue->m_refCount;
}

template <class T> inline
KPtr<T>& KPtr<T>::operator=(const KPtr<T> &rhs)
{
	if (m_ptrValue == rhs.m_ptrValue) {
		return *this;
	}
	
	if (--m_ptrValue->m_refCount == 0) {
		delete m_ptrValue;
	}

	m_ptrValue = rhs.m_ptrValue;
	++m_ptrValue->m_refCount;
	return *this;
}

template <class T> inline
KPtr<T>::~KPtr ()
{
	if (--m_ptrValue->m_refCount == 0) {
		delete m_ptrValue;
	}
}

template <class T> inline
bool KPtr<T>::null() const 
{
	return (m_ptrValue->m_ptr == NULL);
}

template <class T> inline
T* KPtr<T>::operator->() const 
{
	return m_ptrValue->m_ptr;
}

template <class T> inline
T& KPtr<T>::operator*() const 
{
	return *(m_ptrValue->m_ptr);
}



#endif	//_kptr_H
