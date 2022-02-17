//	objectmanager.cpp: implementation of the CObjectManager class.
//
//	Author :		   David Cuin
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "guidex.h"

GuidEx::GuidEx(void)
{
	m_h = m_l = 0Ui64;
}
GuidEx::GuidEx(const GUID& guid)
{
	void* p = (void*)&guid;
	m_h = *((unsigned __int64*)p + 0);
	m_l = *((unsigned __int64*)p + 1);
}
bool GuidEx::operator<(const GuidEx& ge) const
{
	if (m_l < ge.m_l) return true;
	if (m_l > ge.m_l) return false;
	if (m_h < ge.m_h) return true;
	return false;	
}
bool GuidEx::operator>(const GuidEx& ge) const
{
	if (m_l > ge.m_l) return true;
	if (m_l < ge.m_l) return false;
	if (m_h > ge.m_h) return true;
	return false;
}
GuidEx::operator CLSID(void) const
{
	ATLASSERT(sizeof(CLSID) == 2 * sizeof(unsigned __int64));	// Failure here is catastrophic - this class will have to be rewritten.
	CLSID clsid;
	void* p = (void*)&clsid;
	*((unsigned __int64*)p + 0) = m_h;
	*((unsigned __int64*)p + 1) = m_l;
	return clsid;
}
GuidEx::operator std::string(void) const
{
	OLECHAR		szGuid[40];	
	::StringFromGUID2(*this, szGuid, 40);
	return estring(szGuid);
}
const GuidEx& GuidEx::operator++(void)							// This is the prefix operator.
{			
	if (m_l == std::numeric_limits<unsigned __int64>::max()){
		if (m_h == std::numeric_limits<unsigned __int64>::max()){
			// We have used up all the available GUIDs. This will never happen.
			ATLASSERT(false);
			throw "The Universe has ended. Why is this program still running?";
		}
		m_l = 0Ui64;
		m_h++;
	} else {
		m_l++;
	}
	return *this;
}
