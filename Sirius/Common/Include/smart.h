//	smart.h
//
//	This is an extension to the classical reference counting
//  object and pointer implementation. My intention here is
//  to support a relationship between a const and non-const
//  smart pointer to mirror the dumb pointer relationship.
//
//	It is possible to assign a non-const smart pointer to a
//	const smart pointer but not the other way round. The trick
//  here is to use an inheritance structure to implement this
//  behaviour and a union of a const and non-const dumb pointer
//  in the base class.
//
//	NEVER OVERRIDE THE ADDRESS OF OPERATOR OR THE STL CONTAINERS
//  WILL BREAK. USING SOME KIND OF ADAPTER CLASS WOULD OBFUSCATE
//  THE ANALYTICS.
//
//	Authors:							David Cuin
//										Richard See
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _SMART_H
#define _SMART_H

/////////////////////////////////////////////////////////////////////////////
//	RCObject
//
//	Reference counting object.
//
class RCObject {
public:
	void AddReference()
	{
		++m_nReferenceCount;
	}
	void RemoveReference()
	{
		if (--m_nReferenceCount == 0) delete this;
	}
	bool IsShared() const
	{
		return m_nReferenceCount > 1;
	}
protected:
	RCObject() : m_nReferenceCount(0) {}
	RCObject(const RCObject& t) : m_nReferenceCount(0) {}
	const RCObject& operator=(const RCObject&)
	{
		return *this;
	}
	virtual ~RCObject() {};

private:
	int									m_nReferenceCount;	
};


/////////////////////////////////////////////////////////////////////////////
//	RCConstPtr
//
//	Const-smart pointer class. This is the base for the non-const class
//
template<class T>
class RCConstPtr {
public:
	RCConstPtr(const T* realPtr = 0) : m_p_const(realPtr), m_dumb(false)
	{
		init();
	}	
	
	RCConstPtr(const RCConstPtr<T>& t) : m_p(t.m_p), m_dumb(t.m_dumb)
	{
		init();
	}
	
	~RCConstPtr()
	{
		if (m_p && !m_dumb) m_p->RemoveReference();
	}

	const RCConstPtr& operator=(const RCConstPtr& t)
	{
		assign(t);
		return *this;
	}
	
	const T* operator->() const
	{
		return m_p_const;
	}
	
	const T& operator*() const
	{
		return *m_p_const;
	}
	
	// Null tester - NEVER convert to void* or bool or you weaken == comparisons
	bool operator!(void) const
	{
		return m_p_const ? false : true;
	}
	
	// Equality operator
	bool operator==(const RCConstPtr& t) const
	{
		return m_p_const == t.m_p_const ? true : false;
	}

	// returns a dumb form of the pointer
	RCConstPtr dumb(void)
	{
		RCConstPtr h = *this;
		h._dumb();
		return h;
	}
	
protected:
	union
	{
		T*								m_p;
		const T*						m_p_const;
	};
	bool								m_dumb;
		
	void init()
	{
		if (!m_p_const) return;		
		if (!m_dumb) m_p->AddReference();
	}

	void assign(const RCConstPtr& t)
	{
		if (m_p != t.m_p){
			T *p_old = m_p;
			bool m_dumb_old = m_dumb;
			m_p = t.m_p;
			m_dumb = t.m_dumb;
			init();
			if (p_old && !m_dumb_old) p_old->RemoveReference();
		}
	}

	// declare as dumb pointer
	void _dumb(void)
	{
		if (!m_dumb){
			if (m_p){
				bool bShared = m_p->IsShared();
				m_p->RemoveReference();
				// If bShared = false then this means that RemoveReference has deleted the object
				if (!bShared) m_p = NULL;
			}
			m_dumb = true;
		}		
	}
};


/////////////////////////////////////////////////////////////////////////////
//	RCPtr
//
//	Smart pointer class. This is the child of the the const class
//
template<class T>
class RCPtr: public RCConstPtr<T>
{
public:	
	RCPtr(T* realPtr = 0) : RCConstPtr<T>(realPtr)
	{
	}
	
	RCPtr(const RCPtr<T>& t) : RCConstPtr<T>(t)
	{		
	}

	RCPtr& operator=(const RCPtr& t)
	{
		assign(t);		
		return *this;
	}
	
	T* operator->() const
	{
		return m_p;
	}

	T& operator*() const
	{
		return *m_p;
	}

	// returns a dumb form of the pointer
	RCPtr dumb(void)
	{
		RCPtr h = *this;
		h._dumb();
		return h;		
	}
};


#endif
