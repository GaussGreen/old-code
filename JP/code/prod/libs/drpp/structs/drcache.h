// DRCache.h: interface for the DRCache class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRCACHE_H__D4BA5D0F_3D2F_11D2_97D3_00C04FD8EB9A__INCLUDED_)
#define AFX_DRCACHE_H__D4BA5D0F_3D2F_11D2_97D3_00C04FD8EB9A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drstring.h"
#include <list>
#include <map>

//	Classes that support a DRCache.  
//
//	These classes form the base for the cache.  To use the cache, you need
//	to inherit from DRCacheInput and DRCacheItem.
//
//	In the following, I will use caching a tree as an example.  
//
//	The basic idea is that we have a set of inputs (eg critical rates,
//	critical dates, rate environment, etc) that generates an item
//	(eg complete timeline with all the nodes, edges, probabilities and
//	rates).
//
//	In our code, we are constantly using items.  To save time, we would
//	like to cache items, so that we don't have to constantly recreate
//	identical or near identical items.
//
//  We assume that that inputs use up very little memory, but the item is large.
//
//	Our cache does the following:
//
//		1)  For a given input, it checks whether any cached input
//			satisfies any existing input.  For example, an input
//			that describe a tree that is a subtree of a cached
//			tree can simply use the cached tree.
//
//			a)	If it finds a suitable cached item, it will
//				return this item.
//			b)	Otherwise, it will create a new item and put
//				it in cache along with the input description.
//		2)	Maintains a count of how much memory is the cache.  It
//			only keeps track of the item size, since the input is
//			assumed to be small.
//			If we exceed a designated capacity, we will clean the
//			cache of all UNREFERENCED items.  Any non-malicous
//			code will be safe from dangling references.
//
//	Our implementation uses four classes.
//	1)	DRCache		maintains the cache
//	2)	DRCacheInput	input base class
//	3)	DRCacheItem	item base class
//	4)	DRCacheItemPtr	item ptr class (performs reference counting)
//
//	We store the cache in a global variable called:		theCache
//
//	To use:
//		1)	Create an input and an item class that inherit from
//			the DRCacheInput DRCacheItem, respectively.
//		2)	Instantiate the pure virtual functions (DRCacheItem
//			has 1, DRCacheInput has 4)
//		3)	Make sure the constructor of the input class passes
//			the input type to the constructor of DRCacheInput.
//		4)	Whenever you need a cache item, write code that
//			looks like Input input;
//			DRCacheItemPtr a = theCache.get(input);
//			Item* b = (Item*) a.c_ptr();
//
//			Now, you can use b just like an Item*.  To properly
//			maintain the count, you must keep a and b in the
//			same scope.  Usually, a and b are both part of the
//			same class, so this is automatic.
//		
//			Note: this syntax and usage is regrettable, but I
//			can't template operator casts and doing everything
//			via RTTI didn't seem right.

//k
class DRCacheItem 
{
public:
	virtual ~DRCacheItem(){} //
	virtual int bytes() = 0;	// Tells how big the item is
};

//k	Performs reference counting of the cache item
class DRCacheItemPtr
{
public:
	DRCacheItemPtr () : m_ptr(NULL) {} //
	DRCacheItemPtr (DRCacheItem*); //
	DRCacheItemPtr (const DRCacheItemPtr&); //

	DRCacheItemPtr& operator=(const DRCacheItemPtr&); //
	DRCacheItem* c_ptr() {return m_ptr;} // Returns the c-pointer for upward casting

	~DRCacheItemPtr(); //
protected:
	DRCacheItem* m_ptr;
};

//k Base class for the input
class DRCacheInput 
{
public:
	DRCacheInput (DRString type) : m_type (type) {} // String description of type
	virtual ~DRCacheInput(){} //

	bool good (DRCacheInput& a); // Checks whether the input in a is suitable for
	// the input in *this

	DRString type() {return m_type;} // Returns the type

	friend ostream& operator<< (ostream& s, const DRCacheInput& a)
	{a.Print(s); return s;} //
protected:
	friend class DRCache;
	virtual DRCacheItem* create () = 0;
	virtual DRCacheInput* copy () =0;
	virtual void Print(ostream& s) const = 0;
	virtual bool ok_to_use (DRCacheInput&) = 0;
	DRString m_type;
};

//k
class DRCache  
{
public:
	DRCache(int capacity = 70000000); //
	~DRCache(); //
	void set_capacity(int capacity); //

	DRCacheItemPtr get (DRCacheInput&); // Gets a DRCacheItemPtr that satisfies the input
	void clean (); // Cleans any unreferenced memory

	int size(); // Returns current size
	int capacity(); // Returns current capacity

	friend ostream& operator<< (ostream&, const DRCache&); //
protected:
	friend class DRCacheItemPtr;

	struct CacheElem {
		CacheElem (DRCacheInput* input, DRCacheItem* item) 	 
			: m_item(item), m_input(input), m_count(0) {}
		DRCacheItem* m_item;
		DRCacheInput* m_input;
		int m_count;
		bool operator==(const CacheElem& a) const
		{return (m_item == a.m_item && m_input == a.m_input);}
	};

	typedef map <DRCacheItem*, CacheElem*, less<DRCacheItem*>, MYALLOC(CacheElem*) > ItemElemMap;
	typedef list <CacheElem, MYALLOC (CacheElem) > ElemList;
	typedef map <DRString, ElemList, less<DRString>, MYALLOC(ElemList)> TypeElemMap;

	TypeElemMap m_elemMap;
	ItemElemMap m_itemMap;
	int m_capacity;
	int m_size;

	DRCacheItemPtr Add (DRCacheInput&);

	DRCacheItemPtr Use (CacheElem&);

	void CheckCapacity();

	void IncreaseCount (DRCacheItem* item) 
	{
		(m_itemMap[item]->m_count)++;
	}

	void DecreaseCount (DRCacheItem* item) 
	{
		(m_itemMap[item]->m_count)--;
	}
};

extern DRCache theCache; //v Global cache

inline
DRCacheItemPtr::DRCacheItemPtr (DRCacheItem* p):m_ptr(p) 
{if (p) theCache.IncreaseCount(p);}

inline
DRCacheItemPtr::DRCacheItemPtr (const DRCacheItemPtr& a)
{
	m_ptr = a.m_ptr;
	if (m_ptr) theCache.IncreaseCount(m_ptr);
}

inline
DRCacheItemPtr& DRCacheItemPtr::operator=(const DRCacheItemPtr& a)
{
	m_ptr = a.m_ptr;
	if (m_ptr) theCache.IncreaseCount(m_ptr);
	return *this;
}

inline
DRCacheItemPtr::~DRCacheItemPtr()
{
	if (m_ptr)
		theCache.DecreaseCount(m_ptr);
}



#endif // !defined(AFX_DRCACHE_H__D4BA5D0F_3D2F_11D2_97D3_00C04FD8EB9A__INCLUDED_)
