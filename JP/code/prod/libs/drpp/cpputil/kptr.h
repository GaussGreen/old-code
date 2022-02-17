#ifndef	__KPtr__h
/**@#-*/
#define __KPtr__h
/**@#+*/

#include "kplatdep.h"
#include IOSTREAM_H

//k	This is my reference counted pointer class.  Reference counting means
//k	that I count the number of pointers to an object and when this
//k	count goes to zero, I automatically delete this object.  It only makes
//k	sense for this object to be heap allocated.  The general syntax looks like
//k
//k	typedef KPtr<Object>	ObjectPtr;
//k	ObjectPtr x = new Object;
//k
//k	The standard pointer operations * and -> are supported.
//k
//k	Note: operator= only works between KPtrs and adjusts the count(s) as necessary.
//k
//k You can not assign ObjectPtr = Object*.  The reason is that ObjectPtr will
//k destroy the object when its count goes to zero.  Therefore, we must
//k make sure that only ObjectPtr's (no Object*'s) point to the Object.
//k
//k In general, you will mostly want to create pointers to the base and 
//k  		typedef KPtr<base> BasePtr;
//k			BasePtr x = new derived;
//k	will be all you need.  
//k	
//k	If you really need to get the C pointer, there is a way to get it 
//k (without changing the source, but using very ugly syntax) and 
//k then you can apply up and down casts.
//k
//k Limitations of reference counted pointers:
//k
//k		1)	Circular references -- pointers that point to each other will not
//k			get destroyed
//k		2)	Inheritance -- pointers to classes cannot be up or down cast while
//k			preserving the reference counting features.  We really need to count the
//k			total number of pointers pointing to an object (regardless from which point
//k			in the hierachy chain it comes from).


/**  Reference counted pointer class.

  Performs reference counting on heap-allocated objects.  For each object created, we count
  the number of pointers to that object.  When this number goes to zero, the object is automatically
  deleted.  Think of it as simple garbage collecting.  This is NOT a true garbage collector.  There
  are two flaws in this scheme that you should be aware of.

  <ol>
  <li> Circular references -- pointers that point to each other will not get destroyed.
  <li> Inheritance -- pointers to classes cannot be up or down cast while preserving 
   the reference counting features.  We really need to count the
   total number of pointers pointing to an object (regardless from which point
   in the hierachy chain it comes from).
  </ul>

  Otherwise, you can use reference counted pointers, EXACTLY like normal pointers, but you never
  need to explicitly delete anything.
  <br>
  Example:
  <dir> <pre>
  @@typedef KPtr<Object> ObjectPtr; 
  @@ObjectPtr p = new Object(args ...);
  @@p->objectMethod ();
  @@(*p).objectVariable;
  </pre></dir>

  Note that Sun-CC4.2 does not support templated operator->, so you must use * to dereference.
*/

template <class T>
class KPtr {
public:
	/** Constructs from a pointer */
	KPtr (T* ptr = NULL); //
	/** Checks whether the pointer is NULL */
	bool null() const; 

	/**@#-*/
	~KPtr(); //
	KPtr (const KPtr& ptr); //
	KPtr& operator=(const KPtr& rhs); //
	/**@#+*/

	/** Emulates ordinary pointer operator -> */
	T* operator->() const; 
	/** Emulates ordinary pointer operator * */
	T& operator*() const;  

	friend ostream& operator<<(ostream& s, const KPtr& a)
	{s << a.m_ptrValue->m_ptr; return s;}

protected:
	struct KPtrValue {
		T* m_ptr;
		int m_refCount;
		KPtrValue(T* ptr): m_refCount(1), m_ptr(ptr) {}
		~KPtrValue() {if (m_ptr != NULL) delete m_ptr;}
	};

	KPtrValue* m_ptrValue;
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

/**@#-*/

template <class T>
class DRPtr {
public:
	DRPtr (T* ptr = NULL); //
	~DRPtr(); //

	/** Checks whether the pointer is NULL */
	bool null() const; 

	DRPtr (const DRPtr& ptr); //
	DRPtr& operator=(const DRPtr& rhs); //

	/** Emulates ordinary pointer operator -> */
	T* operator->() const; 
	/** Emulates ordinary pointer operator * */
	T& operator*() const;  

	friend ostream& operator<<(ostream& s, const DRPtr& a)
	{s << a.m_ptrValue->m_ptr; return s;}

protected:
	struct DRPtrValue {
		T* m_ptr;
		int m_refCount;
		DRPtrValue(T* ptr): m_refCount(1), m_ptr(ptr) {}
		~DRPtrValue() {if (m_ptr != NULL) delete m_ptr;}
	};

	DRPtrValue* m_ptrValue;
};

template <class T> inline
DRPtr<T>::DRPtr (T* ptr) : m_ptrValue(new DRPtrValue(ptr)) 
{
}

template <class T> inline
DRPtr<T>::DRPtr(const DRPtr &ptr)
{
	m_ptrValue = ptr.m_ptrValue;
	++m_ptrValue->m_refCount;
}

template <class T> inline
DRPtr<T>& DRPtr<T>::operator=(const DRPtr<T> &rhs)
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
DRPtr<T>::~DRPtr ()
{
	if (--m_ptrValue->m_refCount == 0) {
		delete m_ptrValue;
	}
}

template <class T> inline
bool DRPtr<T>::null() const 
{
	return (m_ptrValue->m_ptr == NULL);
}

template <class T> inline
T* DRPtr<T>::operator->() const 
{
	return m_ptrValue->m_ptr;
}

template <class T> inline
T& DRPtr<T>::operator*() const 
{
	return *(m_ptrValue->m_ptr);
}

/**@#+*/

#endif

