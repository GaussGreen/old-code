//	CFastStream implementation
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "faststream.h"

/*static*/ const std::string			CFastStream::s_szEncoding = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";

CFastStream::CFastStream(void)
{	
	m_nWatermark = 0x1000000;
	m_sz = new char[m_nWatermark];
	m_nSize = 0;
}

CFastStream::CFastStream(const std::string& sz)
{
	m_nWatermark = 0x1000000;
	m_sz = new char[m_nWatermark];
	m_nSize = 0;
	this->operator<<(sz);
}

CFastStream::~CFastStream(void)
{
	delete m_sz;
}

CFastStream& CFastStream::operator<<(const std::string& s)
{	
	size_t								nSize = s.size();		// does not include any null terminator

	if (!nSize) return *this;
	if (nSize + m_nSize >= m_nWatermark){						// '>=' allows for sufficient memory for the null terminator
		// we need more memory
		m_nWatermark = (m_nWatermark + nSize) * 2;
		char* sz = new char[m_nWatermark];
		::memcpy(sz, m_sz, m_nSize);
		delete m_sz;
		m_sz = sz;
	}

	::memcpy(m_sz + m_nSize, s.data(), nSize);
	m_nSize += nSize;
	::memset(m_sz + m_nSize, NULL, 1);							// add the termination character
	return *this;
}

CFastStream& CFastStream::operator<<(long n)
{
	return this->operator<<(estring(n));	
}

CFastStream::operator char*(void) const
{
	return m_sz;
}

/*static*/ void CFastStream::FromXml(estring* psz, bool bLeaveAmpersand)
{
	psz->ReplaceStrInStr("&gt;", ">");
	psz->ReplaceStrInStr("&lt;", "<");	
	psz->ReplaceStrInStr("&quot;", "\"");	
	psz->ReplaceStrInStr("&apos;", "'");
	if (!bLeaveAmpersand) psz->ReplaceStrInStr("&amp;", "&");			// We must do this one last
}

/*static*/ HRESULT CFastStream::WriteScalar(const std::string& szPrefix, const CComVariant& vIn, xmlstreamer& ssOut)
{
	HRESULT								hr;
	CComVariant							vElement;

	if (vIn.vt == VT_DATE){
		// we always remap dates to prevent problems arising due to machines with different date conventions
		if (vElement.ChangeType(VT_R8, &vIn)) ATLASSERT(false);
		if (vElement.ChangeType(VT_BSTR)) ATLASSERT(false);
	} else if (hr = vElement.ChangeType(VT_BSTR, &vIn)){
		ssOut << szPrefix << "\tChange of variable failed\r\n";
		return hr;
	}			
	estring szElement(vElement.bstrVal);	
	if (!::strncmp(szElement.data(), s_szEncoding.data(), s_szEncoding.size())){
		// This is already Xml. Therefore, copy anything that follows the encoding string verbatim.
		ssOut << szPrefix << "\t" << szElement.mid(s_szEncoding.size()) << "\r\n";
	} else { 
		CFastStream::ToXml(&szElement, false);
		ssOut << szPrefix << "\t" << szElement << "\r\n";
	}
	return S_OK;
}

/*static*/ void CFastStream::ToXml(estring* psz, bool bLeaveAmpersand)
{
	if (!bLeaveAmpersand) psz->ReplaceStrInStr("&", "&amp;");	// We must do this one first!
	psz->ReplaceStrInStr(">", "&gt;");
	psz->ReplaceStrInStr("<", "&lt;");			
	psz->ReplaceStrInStr("\"", "&quot;");
	psz->ReplaceStrInStr("'", "&apos;");
}
