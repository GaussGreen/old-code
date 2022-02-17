//	guidex.h :	This wraps a guid to implement < and > operators so GUIDs can be
//				used in STL containers.
//
//	Author :	David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifndef _GUIDEX_H
#define _GUIDEX_H

class GuidEx
{
public:	
	GuidEx(void);
	GuidEx(const GUID& guid);
	bool operator<(const GuidEx& ge) const;
	bool operator>(const GuidEx& ge) const;
	operator CLSID(void) const;
	operator std::string(void) const;				// Don't ever define an estring cast since estring casting operators may cause clashes.
	const GuidEx& operator++(void);					// This is the prefix operator.
	
protected:
	unsigned __int64					m_l;
	unsigned __int64					m_h;
};

													
#endif
