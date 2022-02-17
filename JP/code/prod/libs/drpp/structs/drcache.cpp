// DRCache.cpp: implementation of the drcache class.
//
//////////////////////////////////////////////////////////////////////

#include "drcache.h"
#include "drexception.h"

DRCache theCache;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


bool DRCacheInput::good (DRCacheInput& a)
{
	if (a.m_type != m_type) return false;
	return ok_to_use (a);
}


DRCacheItemPtr DRCache::get (DRCacheInput& a)
{
	TypeElemMap::iterator iter = m_elemMap.find(a.type());
	if (iter != m_elemMap.end()) {
		ElemList& elist = (*iter).second;
		ElemList::iterator eIter = elist.begin();
		while (eIter != elist.end()) {
			if (a.good(*((*eIter).m_input)))
				return Use (*eIter);
			eIter++;
		}
	}
	return Add(a);
}

DRCacheItemPtr DRCache::Add (DRCacheInput& a) 
{
	CheckCapacity();
	DRCacheItem* item = a.create();
	DRCacheInput* input = a.copy();
	CacheElem newItem (input, item);
	
	ElemList& elist = m_elemMap[a.type()];
	elist.push_back(newItem);
	
	m_itemMap[item] = &elist.back();

	m_size += item->bytes();
	
//	cout << "Adding " << item << endl;
//	cout << input->m_type << " " << *input << endl;
	return DRCacheItemPtr(item);
}

DRCacheItemPtr DRCache::Use (CacheElem& a) 
{
//	cout << "Using " << a.m_item << endl;
	return DRCacheItemPtr(a.m_item);
}

void DRCache::CheckCapacity()
{
	if (m_size > m_capacity) {
		clean ();
		if (m_size > m_capacity)
			cout << "Warning Failed to Reduce Cache Below Capacity.  Try to eliminate references." << endl;
	}
}

DRCache::DRCache(int capacity) : m_capacity(capacity), m_size(0) 
{}

void DRCache::set_capacity(int capacity) 
{
	m_capacity = capacity; 
	CheckCapacity();
}

int DRCache::size()
{
	return m_size;
}

int DRCache::capacity()
{
	return m_capacity;
}

void DRCache::clean ()
{
	TypeElemMap::iterator mapIter = m_elemMap.begin();
	for (; mapIter != m_elemMap.end(); mapIter++) {
		ElemList& elist = (*mapIter).second;
		ElemList::iterator iter = elist.begin();
		while (iter != elist.end()) {
			if ((*iter).m_count == 0) {
//				cout << "Deleting " << (*iter).m_item << endl;
				m_size -= (*iter).m_item->bytes();
				m_itemMap.erase((*iter).m_item);
				delete (*iter).m_item;
				delete (*iter).m_input;
				iter = elist.erase(iter);
			}
			else {
				iter++;
			}
		}
	}
}

DRCache::~DRCache()
{
	TypeElemMap::iterator mapIter = m_elemMap.begin();
	for (; mapIter != m_elemMap.end(); mapIter++) {
		ElemList& elist = (*mapIter).second;
		ElemList::iterator iter = elist.begin();
		for (; iter != elist.end(); iter++) {
			delete (*iter).m_item;
			delete (*iter).m_input;
		}
	}
}

ostream& operator<< (ostream& s, const DRCache& a)
{
	s << "Size: " << a.m_size << "\t";
	s << "Capacity: " << a.m_capacity << endl;

	DRCache::TypeElemMap::const_iterator mapIter = a.m_elemMap.begin();
	for (; mapIter != a.m_elemMap.end(); mapIter++) {
		s << "Type: " << (*mapIter).first << endl;
		const DRCache::ElemList& elist = (*mapIter).second;
		DRCache::ElemList::const_iterator iter = elist.begin();
		for (; iter != elist.end(); iter++) {
			s << "\tInput: " << *((*iter).m_input) << endl;
			s << "\tItemPtr: " << (*iter).m_item << endl;
			s << "\tItemSize: " << (*iter).m_item->bytes() << endl;
			s << "\tCount: " << (*iter).m_count << endl;
		}
	}
	
	return s;
}
