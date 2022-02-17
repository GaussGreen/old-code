// KCache.cpp: implementation of the drcache class.
//
//////////////////////////////////////////////////////////////////////

#include "kcache.h"
#include "kexception.h"

KCache theCache;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


bool KCacheInput::good (KCacheInput& a)
{
	if (a.m_type != m_type) return false;
	return ok_to_use (a);
}


KCacheItemPtr KCache::get (KCacheInput& a)
{
	TypeElemMap::iterator iter = m_elemMap.find(a.type());
	if (iter != m_elemMap.end()) {
		ElemList& elist = (*iter).second;
		ElemList::iterator eIter = elist.begin();
		while (eIter != elist.end()) {
			if (a.good(*((*eIter).m_input)))
				return use (*eIter);
			eIter++;
		}
	}
	return add(a);
}

KCacheItemPtr KCache::add (KCacheInput& a) 
{
	checkCapacity();
	KCacheItem* item = a.create();
	KCacheInput* input = a.copy();
	CacheElem newItem (input, item);
	
	ElemList& elist = m_elemMap[a.type()];
	elist.push_back(newItem);
	
	m_itemMap[item] = &elist.back();

	m_size += item->bytes();
	
//	cout << "Adding " << item << endl;
//	cout << input->m_type << " " << *input << endl;
	return KCacheItemPtr(item);
}

KCacheItemPtr KCache::use (CacheElem& a) 
{
//	cout << "Using " << a.m_item << endl;
	return KCacheItemPtr(a.m_item);
}

void KCache::checkCapacity()
{
	if (m_size > m_capacity) {
		clean ();
		if (m_size > m_capacity)
			cout << "Warning Failed to Reduce Cache Below Capacity.  Try to eliminate references." << endl;
	}
}

KCache::KCache(int capacity) : m_capacity(capacity), m_size(0) 
{}

void KCache::setCapacity(int capacity) 
{
	m_capacity = capacity; 
	checkCapacity();
}

int KCache::size()
{
	return m_size;
}

int KCache::capacity()
{
	return m_capacity;
}

void KCache::clean ()
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

KCache::~KCache()
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

ostream& operator<< (ostream& s, const KCache& a)
{
	s << "Size: " << a.m_size << "\t";
	s << "Capacity: " << a.m_capacity << endl;

	KCache::TypeElemMap::const_iterator mapIter = a.m_elemMap.begin();
	for (; mapIter != a.m_elemMap.end(); mapIter++) {
		s << "Type: " << (*mapIter).first << endl;
		const KCache::ElemList& elist = (*mapIter).second;
		KCache::ElemList::const_iterator iter = elist.begin();
		for (; iter != elist.end(); iter++) {
			s << "\tInput: " << *((*iter).m_input) << endl;
			s << "\tItemPtr: " << (*iter).m_item << endl;
			s << "\tItemSize: " << (*iter).m_item->bytes() << endl;
			s << "\tCount: " << (*iter).m_count << endl;
		}
	}
	
	return s;
}
