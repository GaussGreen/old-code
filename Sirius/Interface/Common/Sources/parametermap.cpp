//	parametermap.cpp:	implementation of CParameterMap, CStringResource and estring classes
//
//	author:				David Cuin
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "lmcons.h"
#include "parametermap.h"
#undef max
#undef min

#ifdef FIND_CELL_HELPER
#error macro 'FIND_CELL_HELPER' already defined
#endif


#define FIND_CELL_HELPER		GetValue(nRow, nCol, &sz);		\
								if (bIgnoreCaseAndSpaces){		\
									sz.StripWhiteSpace();		\
									sz.lc();					\
								}								\
								if (sz == szMatch){				\
									if (pnRow) *pnRow = nRow;	\
									if (pnCol) *pnCol = nCol;	\
									return S_OK;				\
								}

/*static*/ CEnumMap::enum_map			CEnumMap::s_map;
/*static*/ const VARTYPE				CParameterMap::vt = VT_USERDEFINED;
/*static*/ const IDispatch*				CParameterMap::pdispVal = NULL;

/*static*/ const bool					CParameterMap::s_avIsVariantTypeNumeric[13] = 
{
	false,	/*VT_EMPTY*/
	false,	/*VT_NULL*/
	true,	/*VT_I2*/
	true,	/*VT_I4*/
	true,	/*VT_R4*/
	true,	/*VT_R8*/
	true,	/*VT_CY*/
	true,	/*VT_DATE*/
	false,  /*VT_BSTR*/
	false,  /*VT_DISPATCH*/
	true,   /*VT_ERROR*/
	true,	/*VT_BOOL*/
	false  /*VT_VARIANT*/
};

/*static*/ const VARTYPE				CParameterMap::s_avVariantMap[13][13] = 
	/*These are the standard variable types (VT_UI1 = 17 is missing for brevity)
		VT_EMPTY	= 0,
		VT_NULL		= 1,
		VT_I2		= 2,
		VT_I4		= 3,
		VT_R4		= 4,
		VT_R8		= 5,
		VT_CY		= 6,
		VT_DATE		= 7,
		VT_BSTR		= 8,
		VT_DISPATCH	= 9,
		VT_ERROR	= 10,
		VT_BOOL		= 11,
		VT_VARIANT	= 12,*/
{
	{VT_EMPTY,		VT_NULL,		VT_I2,		VT_I4,		VT_R4,		VT_R8,		VT_CY,		VT_DATE,	VT_BSTR,	VT_DISPATCH,	VT_ERROR,	VT_BOOL,	VT_VARIANT},
	{VT_NULL,		VT_NULL,		VT_I2,		VT_I4,		VT_R4,		VT_R8,		VT_CY,		VT_DATE,	VT_BSTR,	VT_DISPATCH,	VT_ERROR,	VT_BOOL,	VT_VARIANT},
	{VT_I2,			VT_I2,			VT_I2,		VT_I4,		VT_R4,		VT_R8,		VT_CY,		VT_DATE,	VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_I2,		VT_VARIANT},
	{VT_I4,			VT_I4,			VT_I4,		VT_I4,		VT_R8,		VT_R8,		VT_CY,		VT_R8,		VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_I4,		VT_VARIANT},
	{VT_R4,			VT_R4,			VT_R4,		VT_R8,		VT_R4,		VT_R8,		VT_R8,		VT_R8,		VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_R4,		VT_VARIANT},
	{VT_R8,			VT_R8,			VT_R8,		VT_R8,		VT_R8,		VT_R8,		VT_R8,		VT_R8,		VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_R8,		VT_VARIANT},	
	{VT_CY,			VT_CY,			VT_CY,		VT_CY,		VT_R8,		VT_R8,		VT_CY,		VT_R8,		VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_CY,		VT_VARIANT},
	{VT_DATE,		VT_DATE,		VT_DATE,	VT_R8,		VT_R8,		VT_R8,		VT_R8,		VT_DATE,	VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_R8,		VT_VARIANT},
	{VT_BSTR,		VT_BSTR,		VT_BSTR,	VT_BSTR,	VT_BSTR,	VT_BSTR,	VT_BSTR,	VT_BSTR,	VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_BSTR,	VT_VARIANT},
	{VT_DISPATCH,	VT_DISPATCH,	VT_VARIANT,	VT_VARIANT, VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_DISPATCH,	VT_VARIANT, VT_VARIANT,	VT_VARIANT},
	{VT_ERROR,		VT_ERROR,		VT_VARIANT,	VT_VARIANT, VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_DISPATCH,	VT_ERROR,	VT_VARIANT,	VT_VARIANT},
	{VT_BOOL,		VT_BOOL,		VT_I2,		VT_I4,		VT_R4,		VT_R8,		VT_CY,		VT_R8,		VT_BSTR,	VT_VARIANT,		VT_VARIANT,	VT_BOOL,	VT_VARIANT},
	{VT_VARIANT,	VT_VARIANT,		VT_VARIANT,	VT_VARIANT, VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_VARIANT,	VT_DISPATCH,	VT_VARIANT, VT_VARIANT,	VT_VARIANT},
};

/*static*/ const long					CParameterMap::s_avVariantSizes[8] = 
	/*These are the sizes of the variant data element NOT the variant itself*/
{
	0L,	/*VT_EMPTY*/
	0L,	/*VT_NULL*/
	2L,	/*VT_I2*/
	4L,	/*VT_I4*/
	4L,	/*VT_R4*/
	8L,	/*VT_R8*/
	8L,	/*VT_CY*/
	8L	/*VT_DATE*/
};


//////////////////////////////////////////////////////////////////////
// CStringResource implementation
//////////////////////////////////////////////////////////////////////
CStringResource::CStringResource(UINT nResourceID)
{
	m_nResourceID = nResourceID;
}
CStringResource::operator CComVariant(void)
{
	CComBSTR s;
	s.LoadString(m_nResourceID);
	return CComVariant(s);	
}
CComBSTR CStringResource::bstr(void) const
{
	CComBSTR s;
	s.LoadString(m_nResourceID);
	return s;
}
std::string CStringResource::str(void) const
{
	estring sz;
	sz.LoadString(m_nResourceID);
	return sz;
}
CStringResource::operator std::string(void)
{
	return str();
}
int	CStringResource::CompareNoCase(const std::string& sz) const
{
	return estring::CompareNoCase(str(), sz);
}
int CStringResource::CompareNoCase(const CComBSTR& s) const
{
	return estring::CompareNoCase(str(), estring(s));
}
int CStringResource::CompareNoCaseAndSpace(const std::string& sz) const
{
	return estring::CompareNoCaseAndSpace(str(), sz);
}


//////////////////////////////////////////////////////////////////////
// estring implementation
//////////////////////////////////////////////////////////////////////
estring::estring()
{
}
estring::estring(double f)
{
	std::stringstream ss;
	ss << f;
	assign(ss.str());
}
estring::estring(long n)
{
	std::stringstream ss;
	ss << n;
	assign(ss.str());
}
estring::estring(int n)
{
	std::stringstream ss;
	ss << n;
	assign(ss.str());
}
estring::estring(unsigned int n)
{
	std::stringstream ss;
	ss << n;
	assign(ss.str());
}
estring::estring(char* sz)
{
	assign(sz);
}
estring::estring(const std::stringstream& ss)
{
	assign(ss.str());
}
estring::estring(const BSTR& b)
{	
	if (b) assign((LPCTSTR)_bstr_t(b));
}
estring::estring(const _bstr_t& b)
{
	assign((LPCTSTR)b);
}
estring::estring(const CLSID& clsid)
{
	LPOLESTR							psz_unicode = NULL;	
	LPSTR								psz_ansi = NULL;
		
	::StringFromCLSID(clsid, &psz_unicode);
	psz_ansi = UnicodeToAnsi(psz_unicode);	
	assign(psz_ansi);
	delete psz_ansi;
	::CoTaskMemFree(psz_unicode);
}
estring::estring(std::string::size_type n, char c)
{
	resize(n, c);
}
estring::estring(const CComVariant& vIn)
{	
	CComVariant v;
	if (v.ChangeType(VT_BSTR, &vIn)) return;
	assign((char*)_bstr_t(v.bstrVal));	
}

estring::estring(const _variant_t& vIn)
{
	CComVariant v;
	if (v.ChangeType(VT_BSTR, &vIn)) return;
	assign((char*)_bstr_t(v.bstrVal));
}

estring::estring(const VARIANT& vIn)
{
	CComVariant v;
	if (v.ChangeType(VT_BSTR, &vIn)) return;
	assign((char*)_bstr_t(v.bstrVal));
}

estring::estring(const std::string& sz)
{
	assign(sz);
}
estring::estring(const _com_error& e)
{
	try {
		assign(e.Description());
	} catch (...){
		assign(e.ErrorMessage());
	}
}
estring::~estring()
{
}
estring::operator CComBSTR(void) const
{
	return GetBSTR();
}
estring::operator CComVariant(void) const
{
	return CComVariant(GetBSTR());
}
bool estring::operator!=(CComVariant v) const
{
	if (v.ChangeType(VT_BSTR)) return false;
	estring sz(v.bstrVal);
	return strcmp(c_str(), sz.c_str()) ? true : false;
}
bool estring::operator==(UINT nResourceID) const
{
	estring								sz;
	sz.LoadString(nResourceID);
	return CompareNoCaseAndSpace(sz) ? false : true;
}
estring& estring::operator+=(const std::string& sz)
{
	return (estring&)((std::string*)this)->operator+=(sz);
}
estring& estring::operator+=(double f)
{
	return operator+=(estring(f));
}
estring& estring::operator+=(const BSTR& s)
{
	return operator+=(estring(s));	
}
//
//	convers an ANSI string to Unicode
//	
/*static*/ OLECHAR* estring::AnsiToUnicode(const std::string& szAnsi)
{
	return AnsiToUnicode(szAnsi.c_str());
}
/*static*/ OLECHAR* estring::AnsiToUnicode(LPCSTR lpAnsi)
{
   ULONG      ulCnt;
   OLECHAR   *lpUnicode;

   // see how big the unicode string will be
   ulCnt = MultiByteToWideChar(CP_ACP, MB_PRECOMPOSED, lpAnsi, -1, NULL, 0) + 1;

   // allocate a buffer for the unicode string
   lpUnicode = new OLECHAR [ulCnt];

   MultiByteToWideChar(CP_ACP, 0, lpAnsi, -1, lpUnicode, ulCnt);

   return (lpUnicode);
}
//
//	standard string comparison on a case-insensitive basis
//
int estring::CompareNoCase(estring sz) const
{
	estring								szThis(*this);
	szThis.lc();
	sz.lc();
	return strcmp(szThis.c_str(), sz.c_str());
}
/*static*/ int estring::CompareNoCase(std::string sz1, std::string sz2)
{
	estring::lc(&sz1);
	estring::lc(&sz2);
	return strcmp(sz1.c_str(), sz2.c_str());	
}
int estring::CompareNoCase(CComBSTR& s) const
{
	return CompareNoCase(estring(s));
}
/*static*/ int estring::CompareNoCase(const std::string& sz, CComBSTR& s)
{
	return CompareNoCase(sz, estring(s));
}
/*static*/ int estring::CompareNoCaseAndSpace(std::string sz1, std::string sz2)
{
	StripWhiteSpace(&sz1);
	StripWhiteSpace(&sz2);	
	lc(&sz1);
	lc(&sz2);
	return strcmp(sz1.c_str(), sz2.c_str());
}
int estring::CompareNoCaseAndSpace(const CComVariant& v, size_t nCharactersToCompare) const
{
	estring								sz(v);
	estring								szThis(*this);
	StripWhiteSpace(&sz);
	StripWhiteSpace(&szThis);
	szThis.lc();
	sz.lc();
	return strncmp(szThis.c_str(), sz.c_str(), nCharactersToCompare);
}
//
//	standard string comparison on a case and whitespace insensitive basis
//
int estring::CompareNoCaseAndSpace(const std::string& sz) const
{
	return CompareNoCaseAndSpace(*this, sz);
}
//
//	* and ? wildcard comparison
//
/*static*/ bool estring::CompareWildcard(const std::string& sz, const std::string& szWildcard)
{	
	return WildcardMatch(szWildcard.data(), sz.data());
}
bool estring::CompareWildcard(const std::string& szWildcard) const
{	
	return WildcardMatch(szWildcard.data(), data());
}
//
//	Enhanced version of substr. Returns the substring of szIn between the first
//	instance of szFrom and the first instance of szTo after the first instance
//	of szFrom.
//
/*static*/ std::string estring::esubstr(const std::string& szIn, std::string szFrom, std::string szTo)
{
	std::string::size_type nFrom;
	std::string::size_type nTo;

	if ((nFrom = szIn.find(szFrom, 0)) == std::string::npos) return "";
	nFrom += szFrom.length();
	if ((nTo = szIn.find(szTo, nFrom)) == std::string::npos) return "";
	return szIn.substr(nFrom, nTo - nFrom);
}
std::string	estring::esubstr(std::string szFrom, std::string szTo) const
{
	return esubstr(*this, szFrom, szTo);	
}
//
//	performs a case-insensitive find on the object string
//
/*static*/ std::string::size_type estring::findnocase(const std::string* psz, estring sz, std::string::size_type pos/* = 0*/)
{
	estring szThis(*psz);
	szThis.lc();
	sz.lc();
	return szThis.find(sz, pos);
}
std::string::size_type estring::findnocase(estring sz, std::string::size_type pos/* = 0*/) const
{
	estring szThis(*this);
	szThis.lc();
	sz.lc();
	return szThis.find(sz, pos);
}
//
//	Returns a BSTR representation of the string.
//
/*static*/ CComBSTR estring::GetBSTR(const std::string& sz)
// Memory leak fixed by Matt Harris
{
	CComBSTR b(sz.c_str());
	return b;
}
CComBSTR estring::GetBSTR(void) const
{
	return GetBSTR(*this);
}
HRESULT estring::GetBSTR(BSTR* ps) const
{
	return GetBSTR(*this, ps);
}
/*static*/ HRESULT estring::GetBSTR(const std::string& sz, BSTR* ps)
{
	return CComBSTR(sz.c_str()).CopyTo(ps);	
}
/*static*/ estring estring::GetNewGuid(void)
{
	GUID guid;
	OLECHAR		szGuid[40];
	if (FAILED(::CoCreateGuid(&guid))) throw "Unhandled excpetion in estring::GetNewGuid";
	::StringFromGUID2(guid, szGuid, 40);
	return estring(szGuid);
}
//
//	returns a variant form of the object
//
/*static*/ HRESULT estring::GetValue(const std::string& sz, VARIANT* pResult)
{
	return CComVariant(sz.c_str()).Detach(pResult);
}
HRESULT estring::GetValue(VARIANT* pResult) const
{
	return CComVariant(c_str()).Detach(pResult);
}
CComVariant estring::GetValue(void) const
{
	return CComVariant(c_str());
}
/*static*/ CComVariant estring::GetValue(const std::string& sz)
{
	return CComVariant(sz.c_str());
}
//
//	Return true if the object represents a Boolean. Some overloads return the
//	actual Boolean value.
//
/*static*/ bool estring::isbool(const std::string& sz, bool* pb)
{
	CComVariant v(estring::GetValue(sz));
	if (v.ChangeType(VT_BOOL)) return false;
	if (pb) *pb = v.boolVal ? true : false;
	return true;
}
bool estring::isbool(bool* pb) const
{
	return isbool(*this, pb);
}
//
//	Return true if the object represents a date. Some overloads return the
//	actual date value.
//
bool estring::isdate(long* pn) const
{
	return estring::isdate(*this, pn);
}
/*static*/ bool estring::isdate(const std::string& sz, long* pn)
{
	CComVariant v(estring::GetValue(sz));
	if (v.ChangeType(VT_DATE)) return false;
	if (v.ChangeType(VT_I4)) ATLASSERT(false);
	if (pn) *pn = v.lVal;
	return true;
}
//
//	Return true if the object represents a double. Some overloads return the
//	actual double value.
//
/*static*/ bool estring::isdouble(const std::string& sz, double* pf)
{
	char* endp = NULL;
	*pf = strtod(sz.c_str(), &endp);
	return *endp == NULL;
}
/*static*/ bool estring::isdouble(const std::string& sz)
{
	double f;
	return isdouble(sz, &f);
}
bool estring::isdouble(double* pf) const
{
	return isdouble(*this, pf);
}
bool estring::isdouble(void) const
{
	return isdouble(*this);
}
//
//	Returns true if the object represents a long integer. Some overloads
//	return the actual long value.
//
/*static*/ bool estring::islong(const std::string& sz, long* pn)
{
	char* endp = NULL;
	*pn = ::strtol(sz.c_str(), &endp, 10);
	return *endp == NULL;
}
/*static*/ bool estring::islong(const std::string& sz)
{
	long n;
	return islong(sz, &n);
}
bool estring::islong(long* pn) const
{
	return islong(*this, pn);
}
bool estring::islong(void) const
{
	return islong(*this);
}
//
//	Returns true if the object represents an unsigned short integer.
//
/*static*/ bool estring::isunsignedshort(const std::string& sz, unsigned short* pn)
{
	char* endp = NULL;
	long n = ::strtol(sz.c_str(), &endp, 10);
	if (*endp) return false;
	if (n < 0 || n > std::numeric_limits<unsigned short>::max()) return false;
	*pn = (unsigned short)n;
	return true;
}
//
//	Convert as string to lower case
//
/*static*/ void estring::lc(std::string* psz)
{
	std::transform(psz->begin(), psz->end(), psz->begin(), ::tolower);
}
void estring::lc(void)
{    
	lc(this);	
}
//
//	returns the first t characters of the string or the first characters
//  up to (and not including) a given string
//
/*static*/ estring estring::left(const std::string* pszIn, size_type t)
{
	return pszIn->substr(0, t);
}
estring estring::left(size_type t) const
{
	return this->substr(0, t);
}
/*static*/ estring estring::left(const std::string* pszIn, const std::string& sz)
{
	size_type t = pszIn->find(sz);
	return left(pszIn, t == npos ? pszIn->size() : t);
}
estring estring::left(const std::string& sz) const
{
	size_type t = this->find(sz);	
	return left(t == npos ? size() : t);
}
//
//	loads a string from the string resource files
//
HRESULT estring::LoadString(UINT nResourceID)
{
	int									nBufferUsed;
	int									nBufferLength = 1;
	
	do {
		nBufferLength *= 2;
		LPTSTR psz = new TCHAR[nBufferLength];
		nBufferUsed = ::LoadString(_pModule->m_hInstResource, nResourceID, psz, nBufferLength);
		if (nBufferUsed != nBufferLength - 1) assign(psz);
		delete psz;
	} while (nBufferUsed == nBufferLength - 1);
	return S_OK;
}
//
//	Removes the leading white space from a string.
//
/*static*/ void estring::ltrim(std::string* psz)
{
	size_t nFirst = psz->find_first_not_of(" \n\r\t\f\v");
	if (nFirst == psz->npos){
		psz->erase();
	} else {
		psz->assign(psz->substr(nFirst, psz->npos - nFirst));
	}
}
void estring::ltrim(void)
{
	ltrim(this);
}
//	returns the characters between (and not including) the first occurrences of two input strings
/*static*/ estring estring::mid(const std::string* psz, const std::string& szL, const std::string& szR, bool bIgnoreCase /*= false*/)
{
	std::string::size_type				posL;
	std::string::size_type				posR;

	if (bIgnoreCase){
		if ((posL = estring::findnocase(psz, szL)) == std::string::npos) return estring();
	} else {
		if ((posL = psz->find(szL)) == std::string::npos) return estring();
	}
	posL += szL.size();
	if (bIgnoreCase){	
		if ((posR = estring::findnocase(psz, szR, posL)) == std::string::npos) return estring();		
	} else {
		if ((posR = psz->find(szR, posL)) == std::string::npos) return estring();
	}
	return psz->substr(posL, posR - posL);
}
estring estring::mid(const std::string& szL, const std::string& szR, bool bIgnoreCase /*=false*/) const
{
	return mid(this, szL, szR, bIgnoreCase);
}
/*static*/ estring estring::mid(const std::string* psz, const std::string& szL)
{
	std::string::size_type				posL;
	if ((posL = psz->find(szL)) == std::string::npos) return estring();
	posL += szL.size();
	return psz->substr(posL);
}
estring estring::mid(const std::string& szL) const
{
	return mid(this, szL);	
}
estring estring::mid(size_type t) const
{
	return this->substr(t);
}

/*static*/ void estring::ReplaceWhiteSpaceWithSpace(std::string* psz)
{
	std::string::size_type st;
	while ((st = psz->find_first_of("\n\r\t\f\v")) != npos){
		(*psz)[st] = ' ';
	}
}
void estring::ReplaceWhiteSpaceWithSpace(void)
{
	ReplaceWhiteSpaceWithSpace(this);
}

//
//	Replace all occurences of one string with another. Returns the number of substitutions made.
//
/*static*/ long estring::ReplaceStrInStr(std::string* psz, std::string szFrom, std::string szTo, std::string::size_type nPos/* = 0*/)
//	psz - string to adjust
//	szFrom - The substring which should be replaced.
//	szTo - The replacement substring.
//	nPos - The starting position to make the replacement from.
{
	long nSubstitutions = 0L;
	size_t b = nPos;
    for (;;){
		b = psz->find(szFrom, b);
		if (b == npos) break;
		psz->replace(b, szFrom.size(), szTo);
		nSubstitutions++;
		b += szTo.size();
    }  
	return nSubstitutions;
}
long estring::ReplaceStrInStr(std::string szFrom, std::string szTo, std::string::size_type nPos/* =0*/)
{    
	return ReplaceStrInStr(this, szFrom, szTo, nPos);    
}
/*static*/ estring estring::right(const std::string* pszIn, const std::string& sz)
{
	size_type							t;
	
	if ((t = pszIn->rfind(sz)) == std::string::npos) return *pszIn;
	return pszIn->substr(t + sz.size());
}
estring estring::right(const std::string& sz) const
{
	return right(this, sz);
}
/*static*/ estring estring::right(const std::string* pszIn, size_type t)
{
	if (t >= pszIn->size()){
		return *pszIn;
	} else {
		return pszIn->substr(pszIn->size() - t);
	}
}
estring estring::right(size_type t) const
{
	return right(this, t);	
}

//
//	Removes the trailing white space from a string.
//
/*static*/ void estring::rtrim(std::string* psz)
{
	size_t nLast = psz->find_last_not_of(" \n\r\t\f\v");
	
	if (nLast == psz->npos){
		psz->erase();
	} else {		
		psz->assign(psz->substr(0, nLast + 1));
	}
}
void estring::rtrim(void)
{
	rtrim(this);
}
//
//	sets the object from various inputs
//
/*static*/ HRESULT estring::Set(const BSTR& b, std::string* psz)
{
	psz->assign((LPCTSTR)_bstr_t(b));
	return S_OK;
}
HRESULT estring::Set(const BSTR& b)
{	
	assign((LPCTSTR)_bstr_t(b));
	return S_OK;
}
HRESULT estring::Set(LPCOLESTR lpsz)
{
	assign((LPCTSTR)_bstr_t(lpsz));
	return S_OK;
}
HRESULT estring::Set(const std::stringstream& ss)
{
	assign(ss.str());
	return S_OK;
}
void estring::SetError(HRESULT hr)
{	             
	TCHAR* pTemp = NULL;
    int nLen = ::FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_IGNORE_INSERTS | FORMAT_MESSAGE_FROM_SYSTEM,
							   NULL,
							   hr,
							   MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
							   (LPTSTR)&pTemp,
							   1,
							   NULL);
	assign(pTemp);
	trim();
	if (right(1) == ".") assign(left(size() - 1));
	::LocalFree(pTemp);
}
//
//	Sets a vector of strings from this object that are the substrings 
//	separated by a given delimit string.
//
//	We return the size of the vector.
//
/*static*/ long estring::Split(std::string szIn, std::string szDelimit, std::vector<estring>* pVector)
{
	return estring::Split(szIn, szDelimit, (std::vector<std::string>*)pVector);
}
/*static*/ long estring::Split(std::string szIn, std::string szDelimit, std::vector<std::string>* pVector)
{	   
	pVector->clear();

    std::string::size_type offset = 0;
    while ( ( offset = szIn.find(szDelimit) ) != std::string::npos ){
		std::string temp = szIn.substr( 0, offset  );
		if ( temp.size() > 0 ){
			pVector->push_back( temp );
		}
		szIn = szIn.substr( offset+szDelimit.length() );
    }

    // Make sure we get the last segment
    if ( szIn.size() > 0 ){
		pVector->push_back(szIn);
    }
	
	return pVector->size();
}
long estring::Split(std::string szDelimit, std::vector<std::string>* pVector)
{
	return Split(*this, szDelimit, pVector);
}
//
//  Sets a parameter map matrix of strings from this object based on say tab
//  and carriage return string delimiters.
//
void estring::Split(const std::string& szColumnDelimit, const std::string& szRowDelimit, CParameterMap* ppm)
{
	std::vector<std::string>			aszLines;
				
	Split(szRowDelimit, &aszLines);
	ppm->SetRows(aszLines.size());
	for (long nRow = 0; nRow < aszLines.size(); nRow++){
		std::vector<std::string>		aszLine;
		estring::Split(aszLines[nRow], szColumnDelimit, &aszLine);
		if (ppm->GetCols() < aszLine.size()) ppm->SetColumns(aszLine.size());
		for (long nCol = 0; nCol < aszLine.size(); nCol++){
			ppm->SetValue(nRow, nCol, aszLine[nCol]);
		}
	}
}
//
//	Returns A and B from an input string of the form A (B)
//	B will be blank if not of this form.
//	The function returns the correctly spelt form of the input.
//
/*static*/ estring estring::SplitBracket(const std::string& szIn, std::string* pszA, std::string* pszB)
//	pszA - returned, nullable
//	pszB - returned, nullable
{
	std::string szA, szB;
	
	szA = left(&szIn, "(");
	trim(&szA);
	if (!szA.size()){
		szA.assign(szIn);
		trim(&szA);
		szB.clear();
	} else {
		szB = estring::mid(&szIn, "(", ")");
		trim(&szB);
	}
	
	if (pszA) pszA->assign(szA);
	if (pszB) pszB->assign(szB);
	
	if (szB.size()){		
		return szA + " (" + szB + ")";
	} else {
		return szA;
	}
}
//
//	Remove all whitespace characters from a string
//
/*static*/ void estring::StripWhiteSpace(std::string* psz)
{
	ReplaceStrInStr(psz, " ", "", 0);
	ReplaceStrInStr(psz, "\n", "", 0);
	ReplaceStrInStr(psz, "\r", "", 0);
	ReplaceStrInStr(psz, "\t", "", 0);
	ReplaceStrInStr(psz, "\f", "", 0);
	ReplaceStrInStr(psz, "\v", "", 0);
}
void estring::StripWhiteSpace(void)
{        
	StripWhiteSpace(this);
}
//
//	Convert a string to upper case
//
/*static*/ void estring::uc(std::string* psz)
{
	std::transform(psz->begin(), psz->end(), psz->begin(), ::toupper);
}
void estring::uc(void)
{        
	uc(this);
}
//
//	Removes the leading and trailing white space from a string.
//

/*static*/ CComBSTR estring::trim(const CComBSTR& s)
{
	estring sz(s);
	sz.trim();
	return sz;
}
/*static*/ void estring::trim(std::string* psz)
{
	size_t nFirst = psz->find_first_not_of(" \n\r\t\f\v");
	size_t nLast = psz->find_last_not_of(" \n\r\t\f\v");	
	if (nFirst == psz->npos){
		psz->erase();
	} else {
		psz->assign(psz->substr(nFirst, nLast + 1 - nFirst));
	}	
}
void estring::trim(void)
{
	trim(this);
}
//
//	convers a Unicode string to ANSI
//
/*static*/ LPSTR estring::UnicodeToAnsi(LPOLESTR wStr)
{
	ULONG len = WStrlen(wStr);
	LPSTR str = new char[len + 1];
	if (str) sprintf(str,"%ws", wStr);	
	return str;
}

//	Performs the reverse operation to SplitBracket
/*static*/ estring estring::UnSplitBracket(const std::string& szA, const std::string& szB)
{
	if (szA.size() && szB.size()){
		return szA + " (" + szB + ")";
	} else if (szA.size()){
		return szA;
	} else if (szB.size()){
		return szB;
	} else {
		return estring();
	}
}

/*static*/ bool estring::WildcardScan(const char*& szWildcard, const char*& sz)
{
	// First remove the '?' and '*'
	for (szWildcard++; *sz != '\0' && (*szWildcard == '?' || *szWildcard == '*'); szWildcard++){
		if (*szWildcard == '?') sz++;
	}
	while ( *szWildcard == '*'){
		szWildcard++;
	}
	
	// If sz is empty and szWildcard has more characters or szWildcard is empty then return.
	if (*sz == '\0' && *szWildcard != '\0') return false;
	if (*sz == '\0' && *szWildcard == '\0') return true;
	
			
	// Search the substring
	else
	{
		const char* wdsCopy = szWildcard;
		const char* szCopy = sz;
		bool  bYes     = true;
		do  {
			if (!WildcardMatch(szWildcard, sz)){
				szCopy++;
			}
			szWildcard = wdsCopy;
			sz		  = szCopy;
			while ((*szWildcard != *sz) && (*sz != '\0')) sz++;
			wdsCopy = szWildcard;
			szCopy = sz;
		} while ((*sz != '\0') ? !WildcardMatch(szWildcard, sz) : (bYes = false) != false);

		if (*sz == '\0' && *szWildcard == '\0')	return true;

		return bYes;
	}
}
/*static*/ bool estring::WildcardMatch(const char* szWildcard, const char* sz)
{
	bool bYes = true;

	// Iterate and delete '?' and '*' one by one
	while (*szWildcard != '\0' && bYes && *sz != '\0'){
		if (*szWildcard == '?'){
			sz++;
		} else if (*szWildcard == '*'){		
			bYes = WildcardScan(szWildcard, sz);
			szWildcard--;
		} else {
			bYes = (*szWildcard == *sz);
			sz++;
		}
		szWildcard++;
	}
	while (*szWildcard == '*' && bYes) szWildcard++;
	return bYes && *sz == '\0' && *szWildcard == '\0';
}
//
//	returns the length in characters of a Unicode string
//
/*static*/ ULONG estring::WStrlen(PWCHAR str)
{
	ULONG result = 0;
	while (*str++ != UNICODE_NULL) {
		result++;
	}
	return result;
}


//////////////////////////////////////////////////////////////////////
// CParameterMap implementation
//
// Construction/Destruction/Operators
//////////////////////////////////////////////////////////////////////
CParameterMap::CParameterMap() : m_bAutoGrow(false)
{
	m_nRows = m_nCols = 0;
	m_av = NULL;
}
CParameterMap::CParameterMap(const CParameterMap& pm) : m_bAutoGrow(false)
{
	m_nRows = m_nCols  = 0;
	m_av = NULL;
	if (SetValue(pm)) ATLASSERT(false);
}
CParameterMap::CParameterMap(const std::string& s) : m_bAutoGrow(false)
{
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetValue(s);
}
CParameterMap::CParameterMap(const BSTR& b) : m_bAutoGrow(false)
{
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetValue(b);
}
CParameterMap::CParameterMap(double f) : m_bAutoGrow(false)
{
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetValue(f);
}
CParameterMap::CParameterMap(const std::stringstream& ss) : m_bAutoGrow(false)
{
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetValue(ss.str());
}
CParameterMap::CParameterMap(int nRows, int nCols) : m_bAutoGrow(false)
{
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetSize(nRows, nCols);
}
CParameterMap::CParameterMap(const CMatrix& m) : m_bAutoGrow(false)
{	
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetValue(m);
}
CParameterMap::CParameterMap(const CVector& m) : m_bAutoGrow(false)
{	
	m_nRows = m_nCols = 0;
	m_av = NULL;	
	SetValue(m);
}
CParameterMap::~CParameterMap()
{
	Clear();
}
const CParameterMap& CParameterMap::operator=(const CParameterMap& pm)
{
	if (SetValue(pm)) ATLASSERT(false);
	return *this;	
}


//////////////////////////////////////////////////////////////////////////////	
//	AddToEnd
//
//	Adds the input parameter map to the bottom of the current object
//
const CParameterMap& CParameterMap::AddToEnd(const CParameterMap& pm)
{
	long								nRowsOld = GetRows();
	
	if (pm.IsBlank(-1, -1)) return *this;	
	// allocate the extra space in this object
	SetRows(m_nRows + pm.m_nRows);
	SetColumns(__max(m_nCols, pm.m_nCols));
	// copy pm to it
	for (long nRow = nRowsOld; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < pm.m_nCols; nCol++){
			::VariantCopy(m_av + nRow * m_nCols + nCol, pm.m_av + (nRow - nRowsOld) * pm.m_nCols + nCol);
		}
	}
	return *this;
}
const CParameterMap& CParameterMap::AddToEnd(const std::string& sz)
{
	CParameterMap						pm;
	pm.SetValue(sz);
	return AddToEnd(pm);
}
const CParameterMap& CParameterMap::AddToEnd(const std::stringstream& ss)
{
	CParameterMap						pm;
	pm.SetValue(ss);
	return AddToEnd(pm);
}
//
//	Adds two input parameter maps; a parameter map is returned with the first
//	parameter map above the second.
//
/*static*/ const CParameterMap& CParameterMap::AddToEnd(const CParameterMap& pmTop, const CParameterMap& pmBottom)
{	
	return CParameterMap(pmTop).AddToEnd(pmBottom);
}


//////////////////////////////////////////////////////////////////////////////	
//	AddToRHS
//
//	adds the input parameter map to the right of the current object
//
const CParameterMap& CParameterMap::AddToRHS(const CParameterMap& pm)
{
	long								nColsOld = GetCols();
	
	if (pm.IsBlank(-1, -1)) return *this;
	// allocate the extra space in this object	
	SetRows(__max(m_nRows, pm.m_nRows));
	SetColumns(m_nCols + pm.m_nCols);	
	// copy pm to it
	for (long nRow = 0; nRow < pm.m_nRows; nRow++){
		for (long nCol = nColsOld; nCol < m_nCols; nCol++){
			::VariantCopy(m_av + nRow * m_nCols + nCol, pm.m_av + nRow * pm.m_nCols + nCol - nColsOld);
		}
	}
	return *this;
}


//////////////////////////////////////////////////////////////////////////////
//	ArrayToVector
//
//	Transforms an input variant to a std::vector of parameter map objects and/or
//	a std::vector of variants.
//	If vIn is not a 1D safe array then return a vector with just one element in it.
//
/*static*/ HRESULT CParameterMap::ArrayToVector(const VARIANT& vIn, std::vector<CComVariant>* pvvOut, std::vector<CParameterMap>* pvpmOut, bool bDecompose/* = true*/)
//	vIn - input variant
//	pvvOut - (returned, nullable) array of variants corresponding to the decomposed vIn
//	pvpmOut - (returned, nullable) array of parameter maps corresponding to the decomposed vIn
{	
	HRESULT								hr;
	long								nDimensions;
	SAFEARRAY*							psa;							// safe array extracted from the variant		
	CComVariant							vElement;						// safe array element
	CParameterMap						pm;								// element to insert into vpmOut
			
	if (!pvpmOut && !pvvOut) return S_OK;	// no point in doing anything
	nDimensions = GetArrayDimensions(vIn);
	if (pvvOut) pvvOut->clear();
	if (pvpmOut) pvpmOut->clear();
	
	if (nDimensions == 1){				
		// vIn is a 1D safearray
		long							nl, nh;							// boundaries of vector (obvious notation)
		long							nIndex;							// safe array index
		if (vIn.vt & VT_BYREF){
			psa = *(vIn.pparray);
		} else {
			psa = vIn.parray;
		}					
		if (hr = ::SafeArrayGetLBound(psa, 1, &nl)) return hr;
		if (hr = ::SafeArrayGetUBound(psa, 1, &nh)) return hr;
		for (nIndex = nl; nIndex <= nh; nIndex++){				
			if (hr = SafeArrayElementToVariant(psa, &nIndex, &vElement)) return hr;
			if (pvvOut) pvvOut->push_back(vElement);
			if (pvpmOut){		
				if (hr = pm.SetValue(vElement)) return hr;
				pvpmOut->push_back(pm);
			}
		}		
	} else if (nDimensions == 2 && bDecompose){
		// vIn is a 2D safearray. However it could have only 1 row or 1 column; in which case it's a vector.
		long							ncl, nch, nrl, nrh;				// boundaries of matrix (column low, column high, row low and row high)
		long							anIndex[2];
		if (vIn.vt & VT_BYREF){
			psa = *(vIn.pparray);
		} else {
			psa = vIn.parray;
		}		
		if (hr = ::SafeArrayGetLBound(psa, 1, &nrl)) return hr;
		if (hr = ::SafeArrayGetUBound(psa, 1, &nrh)) return hr;
		if (hr = ::SafeArrayGetLBound(psa, 2, &ncl)) return hr;
		if (hr = ::SafeArrayGetUBound(psa, 2, &nch)) return hr;
		if (ncl == nch){
			// column vector
			nDimensions = 1;			
			anIndex[1] = ncl;
			for (anIndex[0] = nrl; anIndex[0] <= nrh; anIndex[0]++){
				if (hr = SafeArrayElementToVariant(psa, anIndex, &vElement)) return hr;
				if (pvvOut) pvvOut->push_back(vElement);
				if (pvpmOut){
					if (hr = pm.SetValue(vElement)) return hr;
					pvpmOut->push_back(pm);
				}
			}
		} else if (nrl == nrh){
			// row vector
			nDimensions = 1;			
			anIndex[0] = nrl;
			for (anIndex[1] = ncl; anIndex[1] <= nch; anIndex[1]++){
				if (hr = SafeArrayElementToVariant(psa, anIndex, &vElement)) return hr;
				if (pvvOut) pvvOut->push_back(vElement);
				if (pvpmOut){
					if (hr = pm.SetValue(vElement)) return hr;
					pvpmOut->push_back(pm);
				}
			}
		}
	}
	if (nDimensions != 1){	
		// pass through
		if (pvvOut) pvvOut->push_back(vIn);
		if (pvpmOut){				
			if (hr = pm.SetValue(vIn)) return hr;
			pvpmOut->push_back(pm);
		}
	}
	
	if (bDecompose){
		// decompose stuff
		if (hr = ArrayToVectorDecompose(pvpmOut)) return hr;
		if (hr = ArrayToVectorDecompose(pvvOut)) return hr;
	}
	return S_OK;
}
//	If pvpm contains only one element, this function checks to see if it is a vector and if so,
//  reconstitutes pvpm as a vector of these elements.
/*static*/ HRESULT CParameterMap::ArrayToVectorDecompose(std::vector<CComVariant>* pvv)
{
	if (!pvv) return S_OK;
	if (pvv->size() == 1){
		CParameterMap pm;
		pm.SetValue((*pvv)[0]);
		if (pm.IsVector()) return pm.GetValue(pvv);			
	}
	return S_OK;	
}
/*static*/ HRESULT CParameterMap::ArrayToVectorDecompose(std::vector<CParameterMap>* pvpm)
{
	if (!pvpm) return S_OK;
	if (pvpm->size() == 1 && (*pvpm)[0].IsVector()) return (*pvpm)[0].GetValue(pvpm);
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	Attach
//
//	Attaches the specified source object to this object and clears the source
//
HRESULT CParameterMap::Attach(CParameterMap* pSrc)
{
	Clear();
	m_nRows = pSrc->m_nRows;
	m_nCols = pSrc->m_nCols;
	m_av = pSrc->m_av;

	// clear the source object without freeing memory
	pSrc->m_nRows = 0;
	pSrc->m_nCols = 0;
	pSrc->m_av = NULL;
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	Clear
//
//	clears member variables
//
void CParameterMap::Clear()
{	
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){			
			::VariantClear(m_av + nRow * m_nCols + nCol);
		}
	}
	delete m_av;
	m_av = NULL;
	m_nRows = m_nCols = 0;
}


///////////////////////////////////////////////////////////////////////////////
//	CLSIDFromProgID
//
//	Wraps ::CLSIDFromProgID
//
/*static*/ HRESULT CParameterMap::CLSIDFromProgID(const CComBSTR& s, CLSID& clsid)
{
	return ::CLSIDFromProgID(s, &clsid);
}


///////////////////////////////////////////////////////////////////////////////
//	CLSIDToObjectName
//
//	Gets the name of the object from the CLSID value. To do this, we use
//	ProgIDFromCLSID which typically returns a string in the form a.b.c, then
//	we extract b.
//
/*static*/ HRESULT CParameterMap::CLSIDToObjectName(const CLSID& CLSID, std::string* psz)
{
	LPOLESTR							s = NULL;
	LPSTR								sz;
	HRESULT								hr;
	
	if (hr = ::ProgIDFromCLSID(CLSID, &s) != S_OK){
		// most likely due to an unregistered DLL		
		return hr;
	}
	sz = estring::UnicodeToAnsi(s);	
	::CoTaskMemFree(s);
	psz->assign(sz);
	delete sz;
	// this is of the form a.b.c, we only need b
	std::string::size_type nPos;
	nPos = psz->find('.');
	if (nPos != std::string::npos){
		psz->assign(psz->substr(nPos + 1, psz->length() - nPos));
		nPos = psz->find('.');
		if (nPos != std::string::npos){
			psz->assign(psz->substr(0, nPos));
		}
	}
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	CopyToClipboard
//
//	Prints the object to the clipboard. This is useful for debugging.
//
HRESULT CParameterMap::CopyToClipboard(void) const
{	
	long								nRow;							// row under consideration
	long								nCol;							// column under consideration		
	HANDLE								h;								// global handle to pass to the clipboard API
	int									nLength(0);						// length of the string
	int									nPos(0);						// position in the global memory string
	std::string							sz;								// string extracted from the matrix
	HRESULT								hr;
	
	// get nLength	
	for (nRow = 0; nRow < m_nRows; nRow++){
		if (nRow) nLength += 2;											// add space for carriage return and line feed
		for (nCol = 0; nCol < m_nCols; nCol++){
			if (hr = GetValue(nRow, nCol, &sz)) return hr;									
			nLength += (sz.size() + (nCol ? 1 : 0));				// add an extra space for a TAB if necessary
		}		
	}
	
	// allocate memory for h and copy the string data into it
	h = ::GlobalAlloc(GMEM_DDESHARE, nLength + 1);				// add space for null termination character
	::memset(h, NULL, nLength + 1);
	for (nRow = 0; nRow < m_nRows; nRow++){
		if (nRow){			
			strcpy((char*)h + nPos, "\x0d\x0a");
			nPos += 2;
		}
		for (nCol = 0; nCol < m_nCols; nCol++){
			if (nCol){
				strcpy((char*)h + nPos, "\t");
				nPos += 1;
			}
			GetValue(nRow, nCol, &sz);			
			strcpy((char*)h + nPos, sz.c_str());
			nPos += sz.size();			
		}		
	}	
	
	// attach h to the clipboard (conceptually, this transfers ownership of h to the clipboard)	
	hr = S_OK;
	if (!::OpenClipboard(NULL)) return E_FAIL;
	if (!::EmptyClipboard() && !hr) hr = E_FAIL;
	if (!::SetClipboardData(CF_TEXT, h) && !hr) hr = E_FAIL;	
	if (!::CloseClipboard() && !hr) hr = E_FAIL;
	return hr;
}

//////////////////////////////////////////////////////////////////////////////
//	CreateBlankVariant
//
//	Creates and returns a blank 2D variant safe-array with nRows and nCols.
//
/*static*/ void CParameterMap::CreateBlankVariant(long nRows, long nCols, VARIANT* pv)
{
	SAFEARRAY*							psa;
	SAFEARRAYBOUND						sab[2];							// safe array size specification vector
	CComVariant							v;
	
	// create the safe array
	sab[0].lLbound = 1;
	sab[0].cElements = nRows;
	sab[1].lLbound = 1;
	sab[1].cElements = nCols;	
	if (!(psa = ::SafeArrayCreate(VT_VARIANT, 2, sab))) throw "Failure to create safearray"; // ToDo - better error - use GetLastError?
					
	// ToDo - newer versions of ATL have a CComVariant constructor / assignment operator that takes a SAFEARRAY pointer
	v.vt = VT_ARRAY | VT_VARIANT;
	v.parray = psa;
	v.Detach(pv);
}


///////////////////////////////////////////////////////////////////////////////
//	DateToString
//
//	maps a date value to a string
//
/*static*/ HRESULT CParameterMap::DateToString(DATE date, std::string* psz)
{
	CComVariant							v;
	HRESULT								hr;
	
	v.vt = VT_DATE;
	v.date = date;
	if (hr = v.ChangeType(VT_BSTR)) return hr;
	psz->assign((char*)_bstr_t(v.bstrVal));
	return S_OK;
}


/////////////////////////////////////////////////////////////////////////////
//	DisplayError
//
//	display any stored error object in an appropriate message box
//
/*static*/ void CParameterMap::DisplayError(const std::string& szMessage, UINT nType)
{	
	int									nLength = GetWindowTextLength(::GetActiveWindow());
	LPTSTR								lpString = new TCHAR[nLength + 1];
	std::string							szTitle;
	estring								sz(szMessage);
	
	if (::GetWindowText(::GetActiveWindow(), lpString, nLength + 1)){
		szTitle.assign(lpString);
	} else {
		szTitle = "Sirius";
	}	
	delete lpString;		
	sz.ReplaceStrInStr("\x95", "");		// remove bullet point
	if (sz.size() && sz[sz.size() - 1] == '\n') sz.resize(sz.size() - 1);
	if (sz.size() && sz[sz.size() - 1] != '.') sz.assign(sz + ".");
	::MessageBox(::GetActiveWindow(), sz.c_str(), szTitle.c_str(), MB_OK | nType);
}
/*static*/ void CParameterMap::DisplayError(const _bstr_t& s, UINT nType)
{
	DisplayError(estring(s), nType);
}
/*static*/ void CParameterMap::DisplayError(UINT nType)
{
	CComVariant							v;	
	ErrorHandler(CComPtr<IDispatch>(), S_OK, &v);
	DisplayError(estring(v), nType);
}
/*static*/ void CParameterMap::DisplayError(UINT nResourceID, UINT nType)
{
	DisplayError(CStringResource(nResourceID), nType);
}


///////////////////////////////////////////////////////////////////////////////
//	ExcelDateToJulianDate
//
//	Maps an Excel-style date to a Julian date. 1904 date system is NOT supported.
//	Note we don't want parameter map to depend on excelinterface.h which is why
//  this function belongs here.
//
/*static*/ unsigned long CParameterMap::ExcelDateToJulianDate(long nExcelDate)
{
	if (nExcelDate <= 60){
		// Microsoft thinks that 1900 was a leap year! Day 60 is 29-Feb-1900.
		return nExcelDate + 2415019 + 1;		
	} else {	
		return nExcelDate + 2415019;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	ExcelGetFormula
//
//	Returns the formula text in the (nRow, nCol) value of an input range.
//	Note we don't want parameter map to depend on excelinterface.h which is why
//  this function belongs here.
//
/*static*/ HRESULT CParameterMap::ExcelGetFormula(CComDispatchDriverEx& ddRange, long nRow, long nCol, std::string* pszOut)
{	
	CComVariant							vNewRange;						// value of "Range.Cells"
	CComVariant							vFormula;						// value of "Range.Formula"
	HRESULT								hr;
	CComVariant							params[2];

	params[0].vt = VT_I4;
	params[0].lVal = nCol + 1;		// convert to 1-base and remember that parameters are passed in reverse order
	params[1].vt = VT_I4;
	params[1].lVal = nRow + 1;
	if (hr = ddRange.GetPropertyByName(L"Cells", params, 2, &vNewRange)) return hr;
	if (hr = CComDispatchDriverEx(vNewRange.pdispVal).GetPropertyByName(L"Formula", &vFormula)) return hr;
	pszOut->assign(estring(vFormula.bstrVal));
	return S_OK;	
}


///////////////////////////////////////////////////////////////////////////////
//	ExcelRangeToRowCol
//
//	Returns the top, left, bottom and right absolute values of an input
//	Excel range.
//	Note we don't want parameter map to depend on excelinterface.h which is why
//  this function belongs here.
//
/*static*/ HRESULT CParameterMap::ExcelRangeToRowCol(CComDispatchDriverEx& ddRange, long* pnTopRow, long* pnLeftCol, long* pnBottomRow, long* pnRightCol)
//	ddRange - wraps an input Excel range object
//	pnTopRow, pnLeftCol, pnBottomRow, pnRightCol  - (returned, nullable)..	
{
	HRESULT								hr;
	CComVariant							v;
	long								nTopRow;
	long								nLeftCol;

	if (hr = ddRange.GetPropertyByName(L"Row", &v)) return hr;
	nTopRow = v.lVal;
	if (hr = ddRange.GetPropertyByName(L"Column", &v)) return hr;
	nLeftCol = v.lVal;
		
	if (pnTopRow) *pnTopRow = nTopRow;
	if (pnLeftCol) *pnLeftCol = nLeftCol;
	if (pnBottomRow){
		if (hr = ddRange.GetPropertyByName(L"Rows", &v)) return hr;
		if (hr = CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Count", &v)) return hr;
		*pnBottomRow = nTopRow + v.lVal - 1;
	}
	if (pnRightCol){
		if (hr = ddRange.GetPropertyByName(L"Columns", &v)) return hr;
		if (hr = CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Count", &v)) return hr;
		*pnRightCol = nLeftCol + v.lVal - 1;
	}
	return S_OK;		
}


///////////////////////////////////////////////////////////////////////////////
//	Extract
//
//	Returns a paramter map which represents this object but containing only
//	those rows headed by specified string resource.
//
//	The unextracted rows are (optionally) returned into another input
//	parameter map
//
HRESULT CParameterMap::Extract(UINT nResourceIDHeader, bool bRemoveHeaderColumn, CParameterMap* ppmExtracted, CParameterMap* ppmRemaining/* = NULL*/) const
//	nResourceIDHeader - string resource to identify which rows we keep
//	bRemoveHeaderColumn - true if we remove the zeroth column at the end of the operation
//	ppmExtracted - parameter map returned containing the extracted rows
//	ppmRemaining - unextracted rows (bRemoveHeaderColumn has no effect on this parameter map)
{
	std::string							szHeader = CStringResource(nResourceIDHeader);
	estring								sz;

	if (ppmExtracted->SetValue(*this) || (ppmRemaining && ppmRemaining->SetValue(*this))){
		ATLASSERT(false);
		return E_FAIL;
	}	
	for (long nRow = m_nRows - 1; nRow >= 0; nRow--){		
		if (GetValue(nRow, 0, &sz)) continue;
		if (sz.CompareNoCaseAndSpace(szHeader)){
			// remove this row
			ppmExtracted->RemoveRow(nRow);
		} else if (ppmRemaining){
			// remove this row from the remaining map
			ppmRemaining->RemoveRow(nRow);
		}
	}
	if (bRemoveHeaderColumn){
		return ppmExtracted->RemoveColumn(0);
	} else {
		return S_OK;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	FindCell
//
//	Searches for a cell in the parameter object with a given value.
//
HRESULT CParameterMap::FindCell(std::string szMatch, bool bIgnoreCaseAndSpaces, find_cell_search fcs, long* pnRow, long* pnCol) const
//	szMatch - String to search for.
//	bIgnoreCaseAndSpaces - True if we perform the search on a case and space insensitive bases.
//	fcs - This is used to speed up the search.
//  pnRow, pnCol - cell coordinates (returned).
{
	long								nRow;
	long								nCol;
	estring								sz;

	if (bIgnoreCaseAndSpaces){
		estring::StripWhiteSpace(&szMatch);
		estring::lc(&szMatch);		
	}
	switch (fcs){
	case LeftColumnOnly:
		nCol = 0;
		for (nRow = 0; nRow < m_nRows; nRow++){
			FIND_CELL_HELPER
		}
		return E_FAIL;
	case TopRowOnly:
		nRow = 0;
		for (nCol = 0; nCol < m_nCols; nCol++){
			FIND_CELL_HELPER
		}
		return E_FAIL;
	case TopLeftByColumns:
		for (nRow = 0; nRow < m_nRows; nRow++){
			for (nCol = 0; nCol < m_nCols; nCol++){
				FIND_CELL_HELPER								
			}
		}
		return E_FAIL;
	case TopLeftByRows:
		for (nCol = 0; nCol < m_nCols; nCol++){
			for (nRow = 0; nRow < m_nRows; nRow++){
				FIND_CELL_HELPER
			}
		}
		return E_FAIL;
	default:
		ATLASSERT(false);
		return E_FAIL;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	GetArrayDimensions
//
//	Returns the number of dimensions in an input variant if it's an array, or
//	zero if it is a scalar.
//
//	We also, optionally, return a vector of the sizes of each dimension. (It's
//	the caller's responsibility to free this memory.)
//
/*static*/ short int CParameterMap::GetArrayDimensions(const CComVariant& v, long** ppsize /*= NULL*/)
//	v - input variant
//	ppsize - if not NULL, then this is set to the sizes of each dimension
{			
	SAFEARRAY*							pArray;							// array part of v
	HRESULT								hRes;							// result of SafeArrayGetDim
		
	if (!(v.vt & VT_ARRAY)) return 0;
	if (v.vt & VT_BYREF){
		pArray = *(v.pparray);
	} else {
		pArray = v.parray;
	}
	hRes = ::SafeArrayGetDim(pArray);
	if (hRes <= 0) return 0;
	if (ppsize){
		ULONG	nDim;		// safe array dimension under consideration
		LONG	nLower;		// upper bound of dimension nDim
		LONG	nUpper;		// lower bound of dimension nDim
		*ppsize = new long[hRes];
		for (nDim = 1; nDim <= hRes; nDim++){		
			::SafeArrayGetLBound(pArray, nDim, &nLower);
			::SafeArrayGetUBound(pArray, nDim, &nUpper);
			(*ppsize)[nDim - 1] = nUpper - nLower + 1;
		}
	}
	return (short int)hRes;
}


///////////////////////////////////////////////////////////////////////////////
//	GetAutoGrow
//
//	return the autogrow value
//
bool CParameterMap::GetAutoGrow(void) const
{
	return m_bAutoGrow;
}


///////////////////////////////////////////////////////////////////////////////
//	GetCols
//
//	returns the number of columns in the parameter map
//
long CParameterMap::GetCols(void) const
{
	return m_nCols;
}


///////////////////////////////////////////////////////////////////////////////
//	GetColumn
//
//	Sets another parameter map as a column of this map.
//
HRESULT CParameterMap::GetColumn(long nCol, CParameterMap* ppm, bool bRemoveBlanks/* = true*/) const
//	nCol - column to obtain
//	ppm - object set
{
	HRESULT								hr;
	if (nCol < 0 || nCol >= m_nCols) return E_FAIL;

	ppm->Clear();
	ppm->m_nRows = m_nRows;
	ppm->m_nCols = 1;
	ppm->m_av = new VARIANT[ppm->m_nRows * ppm->m_nCols];								// ToDo - memory exception
	
	// populate the map
	for (long nRow = 0; nRow < m_nRows; nRow++){				
		::VariantInit(ppm->m_av + nRow * ppm->m_nCols + 0);				// ToDo - out of bounds exception
		if (hr = ::VariantCopy(ppm->m_av + nRow * ppm->m_nCols + 0, m_av + nRow * m_nCols + nCol)) return hr;
	}	
	if (bRemoveBlanks){
		ppm->RemoveBlanksAtEnd();
		ppm->RemoveBlanksAtStart();
	}
	return S_OK;
}
//
//  Returns a 1D variant safearray with elements of a given input type
//  based from a column of the parameter map.
//
void CParameterMap::GetColumn(long nCol, VARTYPE vt, VARIANT* pVal) const
{
	SAFEARRAY*							psa;
	long								nIndex;							// element in the variant array
	SAFEARRAYBOUND						sab;							// safe array size specification vector
	
	// create the safe array
	sab.lLbound = 1;
	sab.cElements = m_nRows;
	if (!(psa = ::SafeArrayCreate(VT_VARIANT, 1, &sab))){
		throw "Could not create safe array";
	}
				
	// populate the safe array
	for (nIndex = 1; nIndex <= m_nRows; nIndex++){
		CComVariant vElement;
		if (vElement.ChangeType(vt, GetConstElementPtr(nIndex - 1, nCol))){
			::SafeArrayDestroy(psa);
			throw "Invalid element at row " + estring(nIndex) + " column " + estring(nCol + 1);
		}
		if (::SafeArrayPutElement(psa, &nIndex, &vElement)){
			::SafeArrayDestroy(psa);
			throw "Could not write element at row " + estring(nIndex);
		}
	}			
	// ToDo - newer versions of ATL have a CComVariant constructor / assignment operator that takes a SAFEARRAY pointer
	CComVariant v;	
	v.vt = VT_ARRAY | VT_VARIANT;
	v.parray = psa;
	if (v.Detach(pVal)){
		throw "Error encountered when detaching the variant";
	}
}


///////////////////////////////////////////////////////////////////////////////
//	GetConstElementPtr
//
//	Returns a constant pointer to the memory location of an element of the
//  parameter map.
//
const CComVariant* CParameterMap::GetConstElementPtr(long nRow, long nCol) const
{
	ATLASSERT(nRow >= 0 && nRow < m_nRows && nCol >= 0 && nCol < m_nCols);
	return (CComVariant*)(m_av + nRow * m_nCols + nCol);
}


///////////////////////////////////////////////////////////////////////////////
//	GetDistinct
//
//	Returns a 1 column parameter map that represents the distinct elements
//  in an input map.
//
void CParameterMap::GetDistinct(CParameterMap* ppm) const
{
	std::map<std::string, bool>			map;
	std::string							sz;
	
	ppm->Clear();
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (GetValue(nRow, nCol, &sz)) throw "Invalid value at R" + estring(nRow + 1) + "C" + estring(nCol + 1);
			map[sz] = true;
		}
	}	
	ppm->SetValue(map);
	ppm->SetColumns(1);
}


///////////////////////////////////////////////////////////////////////////////
//	GetElementPtr
//
//	returns a pointer to the memory location of an element of the parameter map
//
CComVariant* CParameterMap::GetElementPtr(long nRow, long nCol)
{
	ATLASSERT(nRow >= 0 && nRow < m_nRows && nCol >= 0 && nCol < m_nCols);
	return (CComVariant*)(m_av + nRow * m_nCols + nCol);
}


///////////////////////////////////////////////////////////////////////////////
//	GetElementRowCol
//
//	Performs the reverse operation of GetElementPtr
//
HRESULT CParameterMap::GetElementRowCol(const CComVariant* pvElement, long* pnRow, long* pnCol) const
//	pnRow, pnCol - returned (nullable)
{
	unsigned long						nArrayAddress = reinterpret_cast<unsigned long>(m_av);
	unsigned long						nElementAddress = reinterpret_cast<unsigned long>((VARIANT*)pvElement);

	if (nElementAddress < nArrayAddress || nElementAddress > nArrayAddress + m_nRows * m_nCols * sizeof(VARIANT)) return E_FAIL;	
	unsigned long						nDistance = nElementAddress - nArrayAddress;
	ATLASSERT(!(nDistance % sizeof(VARIANT)));	
	long								nElement = nDistance / sizeof(VARIANT);

	if (pnRow){
		*pnRow = nElement / m_nCols;
		ATLASSERT(*pnRow >= 0 && *pnRow < m_nRows);
	}
	if (pnCol){
		*pnCol = nElement % m_nCols;
		ATLASSERT(*pnCol >= 0 && *pnCol < m_nCols);
	}

#	ifdef _DEBUG
	if (pnRow && pnCol){
		CParameterMap*	ppm = const_cast<CParameterMap*>(this);
		CComVariant*	pv = ppm->GetElementPtr(*pnRow, *pnCol);
		ATLASSERT(pv == pvElement);
	}
#	endif
	
	return S_OK;
}

///////////////////////////////////////////////////////////////////////////////
//	GetNonBlankRows
//
//	returns the number of rows in the parameter map that are not blank
//
long CParameterMap::GetNonBlankRows(void) const
{
	long								nRet(0);
	
	for (long nRow = 0; nRow < m_nRows; nRow++){
		if (!IsBlank(nRow)) nRet++;
	}		
	return nRet;
}


////////////////////////////////////////////////////////////////////////////
//	GetObject
//
//	If the parameter map holds a single variant pointing to an object then
//	we return this.
//
HRESULT CParameterMap::GetObject(CComVariant* pvObject) const
{
	if (!m_av) return E_FAIL;
	if (m_nRows != 1 || m_nCols != 1 || m_av->vt != VT_DISPATCH) return E_FAIL;
	return pvObject->Copy(m_av);
}


////////////////////////////////////////////////////////////////////////////
//	GetObjectIID
//
//	returns the interface ID of an object associated with a given IDispatch
//  pointer
//
/*static*/ HRESULT CParameterMap::GetObjectIID(CComDispatchDriverEx& dd, IID* piidOut)
{
	return GetObjectIID(dd.p, piidOut);
}
/*static*/ HRESULT CParameterMap::GetObjectIID(CComPtr<IDispatch> sp, IID* piidOut)
//	sp - input object
//	piidOut - interface ID of pdisp
{
	CComPtr<ITypeInfo>					pti;	
	HRESULT								hr;		

	if (!sp) return E_FAIL;
	if (hr = sp->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;		
	TYPEATTR*	pTypeAttr;
	if (!(hr = pti->GetTypeAttr(&pTypeAttr))){
		*piidOut = pTypeAttr->guid;
		pti->ReleaseTypeAttr(pTypeAttr);
	}
	return hr;
}


////////////////////////////////////////////////////////////////////////////
//	GetObjectName
//
//	returns the name of the object associated with a given IDispatch pointer
//	
/*static*/ HRESULT CParameterMap::GetObjectName(CComDispatchDriverEx& dd, std::string* pszOut)
{
	return GetObjectName(CComPtr<IDispatch>(dd.p), pszOut);
}


///////////////////////////////////////////////////////////////////////////////
//	GetRow
//
//	sets another parameter map as a row of this map
//
HRESULT CParameterMap::GetRow(long nRow, CParameterMap* ppm) const
//	nRow - column to obtain
//	ppm - object set
{
	HRESULT								hr;
	if (nRow < 0 || nRow >= m_nRows) return E_FAIL;

	ppm->Clear();
	ppm->m_nRows = 1;
	ppm->m_nCols = m_nCols;
	ppm->m_av = new VARIANT[ppm->m_nRows * ppm->m_nCols];								// ToDo - memory exception
	
	// populate the map
	for (long nCol = 0; nCol < m_nCols; nCol++){				
		::VariantInit(ppm->m_av + 0 * ppm->m_nCols + nCol);				// ToDo - out of bounds exception
		if (hr = ::VariantCopy(ppm->m_av + 0 * ppm->m_nCols + nCol, m_av + nRow * m_nCols + nCol)) return hr;
	}	
	ppm->RemoveBlanksAtRHS();
	ppm->RemoveBlanksAtLHS();
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	GetRows
//
//	returns the number of rows in the parameter map
//
long CParameterMap::GetRows(void) const
{
	return m_nRows;
}


///////////////////////////////////////////////////////////////////////////////
//	GetString
//
//	Attempts to set a string type variable from the object. This only works for
//	a scalar value.
//
HRESULT CParameterMap::GetString(std::string* psz) const
{	
	if (!m_av) return E_FAIL;
	
	HRESULT			hr;	
	CComVariant		v(m_av[0]);
	
	if (!IsScalar()) return E_FAIL;		
	if (hr = v.ChangeType(VT_BSTR)) return hr;
	psz->assign((char*)_bstr_t(v.bstrVal));
	return S_OK;
}
HRESULT CParameterMap::GetString(BSTR* ps) const
{
	std::string							sz;
	HRESULT								hr;

	if (hr = GetString(&sz)) return hr;	
	return CComBSTR(sz.c_str()).CopyTo(ps);	
}
CComBSTR CParameterMap::GetString(void) const
{
	CComBSTR							s;
	GetString(&s);
	return s;
}


///////////////////////////////////////////////////////////////////////////////
//	GetUser
//
//	Returns a std::string denoting the user name of the current thread.
//	This is the name of the user currently logged onto the system.
//
/*static*/std::string CParameterMap::GetUser(void)
{
	TCHAR								lpsz[UNLEN + 1];				// user name
	DWORD								nLen;							// length of the name buffer
	
	nLen = UNLEN + 1;
	::GetUserName(lpsz, &nLen);
	estring sz(lpsz);
	sz.trim();	
	return sz;
}


///////////////////////////////////////////////////////////////////////////////
//	GetUnionVariantType
//
//	Returns the 'smallest' type of variant that could be used to store the
//  data in two variants of given types.
//
/*static*/ VARTYPE CParameterMap::GetUnionVariantType(VARTYPE vt1, VARTYPE vt2)
{
	if (vt1 < 0 || vt1 > 12 || vt2 < 0 || vt2 > 12) throw "Unknown or unsupported variant type in CResults::GetUnionVariantType";
	return s_avVariantMap[vt1][vt2];
}


///////////////////////////////////////////////////////////////////////////////
//	GetValue
//
//	maps the object to a 2D VB-style variant
//
HRESULT CParameterMap::_GetValue(VARIANT* pv) const
{
	// ToDo - we shouldn't need this function. We need it at present because
	// there is a contradiction in the way that IDispatch pointers are handled
	// in put values.
	if (m_nRows == 1 && m_nCols == 1){
		return ::VariantCopy(pv, m_av);
	} else {
		return GetValue(pv);
	}
}
HRESULT CParameterMap::GetValue(VARIANT* pv) const
{
	::VariantInit(pv);	
	if (IsBlank()){		
		// null
		return CComVariant().Detach(pv);
	} else {
		SAFEARRAY*							psa;
		long								nIndex[2];						// element in the variant array
		SAFEARRAYBOUND						sab[2];							// safe array size specification vector
		CComVariant							v;
		HRESULT								hr;
		
		// create the safe array
		sab[0].lLbound = 1;
		sab[0].cElements = m_nRows;
		sab[1].lLbound = 1;
		sab[1].cElements = m_nCols;	
		if (!(psa = ::SafeArrayCreate(VT_VARIANT, 2, sab))) return E_FAIL;	// ToDo - better error
					
		// populate the safe array
		for (nIndex[0] = 1; nIndex[0] <= m_nRows; nIndex[0]++){
			for (nIndex[1] = 1; nIndex[1] <= m_nCols; nIndex[1]++){						
				CComVariant vElement;						
				if (hr = vElement.Copy(m_av + (nIndex[0] - 1) * m_nCols + (nIndex[1] - 1))){
					::SafeArrayDestroy(psa);
					return hr;
				}
				if (hr = ::SafeArrayPutElement(psa, nIndex, &vElement)){				
					::SafeArrayDestroy(psa);
					return hr;
				}
			}
		}			
		// ToDo - newer versions of ATL have a CComVariant constructor / assignment operator that takes a SAFEARRAY pointer
		v.vt = VT_ARRAY | VT_VARIANT;
		v.parray = psa;
		return v.Detach(pv);
	}
}
CComVariant CParameterMap::GetValue(void) const
{	
	CComVariant	vRet;		
	if (GetValue(&vRet)) ATLASSERT(false);
	return vRet;	
}
//
//	returns a pointer to the element at (nRow, nCol) or NULL if invalid
//
VARIANT* CParameterMap::GetValue(long nRow, long nCol) const
{
	if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols) return NULL;
	return (m_av + nRow * m_nCols + nCol);
}
//
//	if nRow <= -1 and nCol <= -1 we return the whole matrix
//	if nRow <= -1 and nCol >= 0 we return the (nCol)th column as a vector
//	if nRow >= 0  and nCol <= -1 we return the (nRow)th row as a vector
//	if nRow >= 0 and nCol >= 0 we return a particular element
//
HRESULT CParameterMap::GetValue(long nRow, long nCol, VARIANT* pv) const
{
	if (nRow < 0 && nCol < 0){
		return GetValue(pv);
	} else if (nRow < 0){
		// return the (nCol)th column as a vector
		if (nCol >= m_nCols) return E_FAIL;
		// ToDo - finish
		ATLASSERT(false);
		return E_FAIL;		
	} else if (nCol < 0){
		// return the (nRow)th row as a vector
		if (nRow >= m_nRows) return E_FAIL;
		// ToDo - finish
		ATLASSERT(false);
		return E_FAIL;
	}
	if (nRow >= m_nRows || nCol >= m_nCols) return E_FAIL;
	HRESULT			hr;
	CComVariant		v;
	if (hr = v.Copy(m_av + nRow * m_nCols + nCol)) return hr;
	return v.Detach(pv);
}
//
//	return an element of the matrix
//
HRESULT CParameterMap::GetValue(long nRow, long nCol, bool* pb) const
{
	HRESULT								hr;
	CComVariant							v;
	if (nRow < 0 || nCol < 0) return E_FAIL;
	if (hr = GetValue(nRow, nCol, &v)) return hr;
	if (v.ChangeType(VT_BOOL)) return E_FAIL;		
	*pb = (v.boolVal == VARIANT_TRUE ? true : false);
	return S_OK;
}
HRESULT CParameterMap::GetValue(long nRow, long nCol, unsigned short* pn) const
{
	HRESULT								hr;
	long								n;
	if (hr = GetValue(nRow, nCol, &n)) return hr;
	if (n < 0 || n > std::numeric_limits<unsigned short>::max()) return E_FAIL;		// Valid because longs can be bigger than unsigned shorts.
	*pn = (unsigned short)n;
	return S_OK;
}
HRESULT	CParameterMap::GetValue(long nRow, long nCol, long* pn) const
{
	HRESULT								hr;
	CComVariant							v;
	if (nRow < 0 || nCol < 0) return E_FAIL;
	if (hr = GetValue(nRow, nCol, &v)) return hr;
	if (v.ChangeType(VT_I4)){
		// try changing to a date
		if (v.ChangeType(VT_DATE)) return E_FAIL;
		*pn = (long)v.date;
		return S_OK;
	}	
	*pn = v.lVal;
	return S_OK;
}
HRESULT CParameterMap::GetValue(long nRow, long nCol, double* pf) const
{
	HRESULT								hr;
	CComVariant							v;
	if (nRow < 0 || nCol < 0) return E_FAIL;
	if (hr = GetValue(nRow, nCol, &v)) return hr;
	if (v.ChangeType(VT_R8)){
		// try changing to a date
		if (v.ChangeType(VT_DATE)) return E_FAIL;
		*pf = v.date;
		return S_OK;
	}
	*pf = v.dblVal;
	return S_OK;
}
HRESULT	CParameterMap::GetValue(long nRow, long nCol, std::string* psz) const
{
	HRESULT								hr;
	CComVariant							v;
	if (nRow < 0 || nCol < 0) return E_FAIL;
	if (hr = GetValue(nRow, nCol, &v)) return hr;	
	if (v.ChangeType(VT_BSTR)) return E_FAIL;	
	psz->assign((char*)_bstr_t(v.bstrVal));
	return S_OK;
}
HRESULT CParameterMap::GetValue(long nRow, long nCol, BSTR* ps) const
{		
	std::string							sz;
	HRESULT								hr;

	if (hr = GetValue(nRow, nCol, &sz)) return hr;
	return CComBSTR(sz.c_str()).CopyTo(ps);	
}
HRESULT CParameterMap::GetValue(long nRow, long nCol, CParameterMap* ppm) const
{
	HRESULT								hr;
	CComVariant							v;

	if (hr = GetValue(nRow, nCol, &v)) return hr;
	return ppm->SetValue(v);
}
HRESULT CParameterMap::GetValue(long nRow, long nCol, CVector* pvector) const
{
	HRESULT								hr;
	CParameterMap						pm;

	if (hr = GetValue(nRow, nCol, &pm)) return hr;
	return pm.GetValue(pvector);
}
HRESULT CParameterMap::GetValue(double* pf) const
{
	if (!IsScalar()) return E_FAIL;
	return GetValue(0, 0, pf);
}
HRESULT CParameterMap::GetValue(long* pn) const
{
	if (!IsScalar()) return E_FAIL;
	return GetValue(0, 0, pn);
}
HRESULT CParameterMap::GetValue(bool* pb) const
{	
	if (!IsScalar()) return E_FAIL;
	return GetValue(0, 0, pb);
}
HRESULT CParameterMap::GetValue(std::string* psz) const
{	
	return GetString(psz);	
}
HRESULT CParameterMap::GetValue(BSTR* ps) const
{
	HRESULT	hr;
	
	ATLASSERT(false);	// ToDo - test
	if (!m_av || !IsScalar()) return E_FAIL;
	CComVariant v(m_av[0]);
	if (hr = v.ChangeType(VT_BSTR)) return hr;	
	CComBSTR s;
	s.Attach(v.bstrVal);	// ToDo - test that this does not take a deep copy - I expect v.bstrVal to be set to null after this.
	return s.CopyTo(ps);
	return S_OK;
}
HRESULT CParameterMap::GetValue(CParameterMap* ppm) const
{
	return ppm->SetValue(*this);
}
HRESULT CParameterMap::GetValue(CVector* pvector) const
{
	if (m_nCols == 1){
		pvector->resize(m_nRows);
		for (long nRow = 0; nRow < m_nRows; nRow++){		
			 GetValue(nRow, 0, &(*pvector)[(int)nRow]);
		}
		return S_OK;
	} else if (m_nRows == 1){			
		pvector->resize(m_nCols);
		for (long nCol = 0; nCol < m_nCols; nCol++){
			 GetValue(0, nCol, &(*pvector)[(int)nCol]);
		}
		return S_OK;
	}
	return E_FAIL;	
}
HRESULT CParameterMap::GetValue(CMatrix* pmatrix) const
{
	pmatrix->resize(m_nRows, m_nCols);		
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){			
			GetValue(nRow, nCol ,&(*pmatrix)[(int)nRow][(int)nCol]);
		}
	}
	return S_OK;
}
HRESULT CParameterMap::GetValue(std::stringstream& ssOut) const
{
	HRESULT								hr;		
	
	ssOut.str("");		
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			std::string sz;
			if (hr = GetValue(nRow, nCol, &sz)) return hr;
			ssOut << sz;
		}
	}
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	GetValueFromColumnHeader
//
//	This function finds the cell containing szHeader and returns the value of
//	the cell immediately below that row
//
HRESULT CParameterMap::GetValueFromColumnHeader(const std::string& szHeader, double* pf) const
{
	HRESULT								hr;
	long								nRow, nCol;
	
	if (hr = FindCell(szHeader, true, TopRowOnly, &nRow, &nCol)) return hr;
	return GetValue(nRow + 1, nCol, pf);
}


///////////////////////////////////////////////////////////////////////////////
//	GetValueFromRowHeader
//
//	This function finds the row starting with szHeader and returns the value of
//	the cell immediately to the RHS of that row
//
HRESULT CParameterMap::GetValueFromRowHeader(UINT nResourceIDHeader, bool bAllowBlank, std::string* pszData, long* pnRowUsed/* = NULL */) const
{
	return GetValueFromRowHeader(CStringResource(nResourceIDHeader), bAllowBlank, pszData, pnRowUsed);
}
HRESULT CParameterMap::GetValueFromRowHeader(const std::string& szHeader, bool bAllowBlank, std::string* pszData, long* pnRowUsed/* = NULL */) const
{
	HRESULT								hr;
	long								nRow, nCol;
		
	if (hr = FindCell(szHeader, true, LeftColumnOnly, &nRow, &nCol)) return hr;
	if (pnRowUsed) *pnRowUsed = nRow;
	if (hr = GetValue(nRow, nCol + 1, pszData)) return hr;
	estring::trim(pszData);
	if (!bAllowBlank && !pszData->size()) return E_FAIL;
	return S_OK;
}
HRESULT	CParameterMap::GetValueFromRowHeader(UINT nResourceIDHeader, double* pf, long* pnRowUsed/* = NULL*/) const
{
	return GetValueFromRowHeader(CStringResource(nResourceIDHeader), pf, pnRowUsed);
}
HRESULT CParameterMap::GetValueFromRowHeader(const std::string& szHeader, double* pf, long* pnRowUsed/* = NULL*/) const
{
	HRESULT								hr;
	long								nRow, nCol;

	if (hr = FindCell(szHeader, true, LeftColumnOnly, &nRow, &nCol)) return hr;	
	if (pnRowUsed) *pnRowUsed = nRow;
	return GetValue(nRow, nCol + 1, pf);
}
HRESULT	CParameterMap::GetValueFromRowHeader(UINT nResourceIDHeader, CComVariant* pv, long* pnRowUsed/* = NULL*/) const
{
	return GetValueFromRowHeader(CStringResource(nResourceIDHeader), pv, pnRowUsed);
}
HRESULT CParameterMap::GetValueFromRowHeader(const std::string& szHeader, CComVariant* pv, long* pnRowUsed/* = NULL */) const
{
	HRESULT								hr;
	long								nRow, nCol;
	
	if (hr = FindCell(szHeader, true, LeftColumnOnly, &nRow, &nCol)) return hr;
	if (pnRowUsed) *pnRowUsed = nRow;
	return GetValue(nRow, nCol + 1, pv);
}


///////////////////////////////////////////////////////////////////////////////
//	GetVariantSize
//
//	Returns the number of bytes in the data field of the variant (Not the size
//  of the variant itself).
//
/*static*/ long	CParameterMap::GetVariantSize(const CComVariant& v)
{
	if (v.vt == VT_BSTR){
		// special case - return the length of the string
		long n = ::SysStringLen(v.bstrVal);
		return n;
	} else if (v.vt < 0 || v.vt > 7){
		throw "Unhandled variant type " + estring(v.vt) + " in CParameterMap::GetVariantSize";
	} else {
		return s_avVariantSizes[v.vt];
	}
	return 0;
}


/////////////////////////////////////////////////////////////////////////////
//	IsExcelFunctionWizardVisible
//
//	returns TRUE if the excel function wizard is visible
//
/*static*/ bool CParameterMap::IsExcelFunctionWizardVisible(void)
{
	HWND								hWnd;							// window handle	
	std::string							szClassName;
	
	// The strategy here is to find the first blank caption window of class "bosa_sdm_XL8".
	// If the function wizard is showing for the active process then this will be the
	// window found. Hence, we check the class of the parent window and check its activation.
	
	if (!(hWnd = ::FindWindow(_T("bosa_sdm_XL8"), _T("")))){
		// the function wizard associated with Excel version 8 is not showing - try version 9
		if (!(hWnd = ::FindWindow(_T("bosa_sdm_XL9"), _T("Function Arguments")))){
			return false;
		}
	}
	
	// now examine the parent window
	if (!(hWnd =::GetWindow(hWnd, GW_OWNER))) return false;			// window is not owned by another window

	// check the class name - should be the main excel frame class
	if (CParameterMap::WindowToClassName(hWnd, &szClassName)) return false;
	if (szClassName != "XLMAIN") return false;
	
	// compare with the active window (the active window can be either the wizard window itself or one
	// of the ref-edit controls)
	if (hWnd != ::GetActiveWindow() && hWnd != ::GetParent(::GetActiveWindow())) return false;

	// wizard is showing in the current process if this point is reached
	return true;	
}


///////////////////////////////////////////////////////////////////////////////
//	IsBlank
//
//	if (nRow < 0 and nCol < 0) return true iff the whole matrix is blank
//	if (nRow < 0 and nCol >=0) return true iff column nCol is blank
//	if (nRow >=0 and nCol < 0) return true iff row nRow is blank
//	if (nRow >=0 and nCol >=0) return true iff the element (nRow, nCol) is blank
//	
bool CParameterMap::IsBlank(long nRow/* = -1*/, long nCol/* = -1*/) const
{
	if (!m_nRows || !m_nCols) return true;
	if (nRow < 0){
		// return true iff the entire matrix is blank (for nCol < 0 case) or return true iff column nCol is blank
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (!IsBlank(nRow, nCol)) return false;
		}					
		return true;
	} else if (nCol < 0){
		// return true iff row nRow is blank
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (!IsBlank(nRow, nCol)) return false;
		}
		return true;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		// out of bounds returns true - we are blank at this cell
		return true;
	}

	switch (m_av[nRow * m_nCols + nCol].vt){
	case VT_EMPTY: case VT_ERROR:
		return true;
	case VT_BSTR:
		return !_bstr_t(m_av[nRow * m_nCols + nCol].bstrVal).length();
	default:
		return false;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	IsBlankOrDouble
//
//	Returns the value for CParameterMap::IsBlank || CParameterMap::IsDouble.
//
//	For non-zero nResourceID, an error object is set if an actual element tests false.
//
bool CParameterMap::IsBlankOrDouble(long nRow/* = -1*/, long nCol/* = -1*/, UINT nResourceID/* = 0*/) const
{
	if (!m_nRows || !m_nCols) return true;
	if (nRow < 0){
		// return true iff the entire matrix is blank or double (for nCol < 0 case) or return true iff column nCol is blank or double
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (!IsBlankOrDouble(nRow, nCol, nResourceID)) return false;
		}					
		return true;
	} else if (nCol < 0){
		// return true iff row nRow is blank
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (!IsBlankOrDouble(nRow, nCol, nResourceID)) return false;
		}
		return true;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		// out of bounds returns true - we are blank at this cell
		return true;
	}

	return (IsBlank(nRow, nCol) || IsDouble(nRow, nCol, nResourceID) /*IsBlank should be called before IsDouble or nResourceID might be missed*/);
}


///////////////////////////////////////////////////////////////////////////////
//	IsColumnVector
//
//	returns true if this object has one row and many columns
//
bool CParameterMap::IsColumnVector(void) const
{
	return (m_nCols == 1);
}


///////////////////////////////////////////////////////////////////////////////
//	IsDate
//
//	if (nRow < 0 and nCol < 0) return true iff the all the cells in the matrix are set to dates
//	if (nRow < 0 and nCol >=0) return true iff all the cells in column nCol are set to dates
//	if (nRow >=0 and nCol < 0) return true iff all the cells in row nRow are set to dates
//	if (nRow >=0 and nCol >=0) return true iff the element (nRow, nCol) is a date
//
//	For non-zero nResourceID, an error object is set if an actual element tests false.
//
bool CParameterMap::IsDate(long nRow/* = -1*/, long nCol/* = -1*/, UINT nResourceID/* = 0*/) const
{
	if (!m_nRows || !m_nCols) return false;
	if (nRow < 0){
		// return true iff all the cells in the matrix matrix are set to dates (for nCol < 0 case)
		// or return true iff all the cells in column nCol are set to dates
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (!IsDate(nRow, nCol, nResourceID)) return false;
		}					
		return true;
	} else if (nCol < 0){
		// return true iff all the cells in row nRow are set to dates
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (!IsDate(nRow, nCol, nResourceID)) return false;
		}
		return true;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		// out of bounds returns false - we do not have a date at this cell
		return false;
	}	
	return !CComVariant().ChangeType(VT_DATE, m_av + nRow * m_nCols + nCol);
}


///////////////////////////////////////////////////////////////////////////////
//	IsDouble
//
//	if (nRow < 0 and nCol < 0) return true iff the whole matrix is double
//	if (nRow < 0 and nCol >=0) return true iff column nCol is double
//	if (nRow >=0 and nCol < 0) return true iff row nRow is double
//	if (nRow >=0 and nCol >=0) return true iff the element (nRow, nCol) is double
//
//	For non-zero nResourceID, an error object is set if an actual element tests false.
//
bool CParameterMap::IsDouble(long nRow/* = -1*/, long nCol/* = -1*/, UINT nResourceID/* = 0*/) const
{
	if (!m_nRows || !m_nCols) return false;
	if (nRow < 0){
		// return true iff the entire matrix is double (for nCol < 0 case) or return true iff column nCol is double
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (!IsDouble(nRow, nCol, nResourceID)) return false;
		}					
		return true;
	} else if (nCol < 0){
		// return true iff row nRow is double
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (!IsDouble(nRow, nCol, nResourceID)) return false;
		}
		return true;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		// out of bounds returns false - we do not have a double at this cell
		return false;
	}

	switch (m_av[nRow * m_nCols + nCol].vt){
	case VT_I2: case VT_I4: case VT_R4: case VT_R8: case VT_CY: case VT_DATE: case VT_BOOL: case VT_DECIMAL: case VT_I1: case VT_UI1: case VT_UI2: case VT_UI4: case VT_INT: case VT_UINT:
		return true;
	case VT_BSTR:
		if (estring(m_av[nRow * m_nCols + nCol]).isdouble()) return true;
		// follow-through to default is intentional
	default:
		;// cases not considered VT_EMPTY VT_NULL VT_DISPATCH VT_ERROR VT_VARIANT VT_UNKNOWN VT_RECORD VT_ARRAY		
	}
	if (nResourceID) CParameterMap::ReturnErrorRV(nResourceID, m_av[nRow * m_nCols + nCol]);
	return false;
}


///////////////////////////////////////////////////////////////////////////////
//	IsLong
//
//	if (nRow < 0 and nCol < 0) return true iff the whole matrix is long
//	if (nRow < 0 and nCol >=0) return true iff column nCol is long
//	if (nRow >=0 and nCol < 0) return true iff row nRow is long
//	if (nRow >=0 and nCol >=0) return true iff the element (nRow, nCol) is long
//	
//	For non-zero nResourceID, an error object is set if an actual element tests false.
//
bool CParameterMap::IsLong(long nRow/* = -1*/, long nCol/* = -1*/, UINT nResourceID/* = 0*/) const
{		
	if (!m_nRows || !m_nCols) return false;
	if (nRow < 0){
		// return true iff the entire matrix is long (for nCol < 0 case) or return true iff column nCol is long
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (!IsLong(nRow, nCol, nResourceID)) return false;
		}					
		return true;
	} else if (nCol < 0){
		// return true iff row nRow is long
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (!IsLong(nRow, nCol, nResourceID)) return false;
		}
		return true;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		// out of bounds returns false - we do not have a long at this cell
		return false;
	}

	switch (m_av[nRow * m_nCols + nCol].vt){
	case VT_I2: case VT_I4: case VT_BOOL: case VT_I1: case VT_UI1: case VT_UI2: case VT_UI4: case VT_INT: case VT_UINT:
		return true;
	case VT_R4: case VT_R8: case VT_CY: case VT_DATE: case VT_DECIMAL:
		// it might be possible that the number has no decimal even if its data type allows it
		{			
			CComVariant v;
			if (!v.ChangeType(VT_BSTR, m_av + nRow * m_nCols + nCol) && (estring(v).islong())) return true;			
		}	
	case VT_BSTR:
		if (estring(m_av[nRow * m_nCols + nCol]).islong()) return true;
	default:
		; // cases not considered VT_EMPTY VT_NULL VT_DISPATCH VT_ERROR VT_VARIANT VT_UNKNOWN VT_RECORD VT_ARRAY		
	}
	if (nResourceID) CParameterMap::ReturnErrorRV(nResourceID, m_av[nRow * m_nCols + nCol]);
	return false;
}


///////////////////////////////////////////////////////////////////////////////
//	IsObject
//
//	Returns true if the parameter map contains exactly one element of
//	type VT_DISPATCH that contains a non-null pointer.
//
bool CParameterMap::IsObject(void) const
{
	if (m_nRows != 1 || m_nCols != 1) return false;
	return GetConstElementPtr(0, 0)->vt == VT_DISPATCH && GetConstElementPtr(0, 0)->pdispVal ? true : false;
}


///////////////////////////////////////////////////////////////////////////////
//	IsRowVector
//
//	returns true if this object has one column and many rows
//
bool CParameterMap::IsRowVector(void) const
{
	return (m_nRows == 1);
}


///////////////////////////////////////////////////////////////////////////////
//	IsScalar
//
//	Returns true if this object represents a scalar.
//
bool CParameterMap::IsScalar(void) const
{
	return (m_nRows == 1 && m_nCols == 1 && (GetArrayDimensions(m_av[0]) == 0));
}


///////////////////////////////////////////////////////////////////////////////
//	IsVariantTypeNumeric
//
//	Returns true if an input variant type is a number.
//
/*static*/ bool CParameterMap::IsVariantTypeNumeric(VARTYPE vt)
{
	if (vt < 0 || vt > 12) throw "Unhandled variant type " + estring(vt) + " in CParameterMap::IsVariantTypeNumeric";
	return s_avIsVariantTypeNumeric[vt];
}


///////////////////////////////////////////////////////////////////////////////
//	IsVector
//
//	return true if this object represents a row or column vector
//
bool CParameterMap::IsVector(void) const
{
	return (m_nRows == 1 || m_nCols == 1);
}


///////////////////////////////////////////////////////////////////////////////
//	InsertAt
//
//	inserts a parameter map matrix into this object starting at the specified
//	row and column
//
HRESULT CParameterMap::InsertAt(long nStartRow, long nStartCol, const CParameterMap& pm)
{
	HRESULT								hr;
	
	// ensure that this object is large enough to accommodate pm
	if (hr = SetRows(__max(m_nRows, nStartRow + pm.m_nRows))) return hr;
	if (hr = SetColumns(__max(m_nCols, nStartCol + pm.m_nCols))) return hr;

	// insert pm
	for (long nRow = 0; nRow < pm.m_nRows; nRow++){
		for (long nCol = 0; nCol < pm.m_nCols; nCol++){
			if (hr = ::VariantCopy(m_av + (nStartRow + nRow) * m_nCols + (nStartCol + nCol), pm.m_av + nRow * pm.m_nCols + nCol)) return hr;
		}
	}
	return S_OK;	
}


///////////////////////////////////////////////////////////////////////////////
//	InsertColumn
//
//	inserts a column at the specified position in the parameter map
//
HRESULT CParameterMap::InsertColumn(long nAt)
{
	if (nAt < 0){		
		return E_FAIL;
	} else if (nAt <= m_nCols){
		// insert a new column at position nAt
		VARIANT* av = new VARIANT[m_nRows * (m_nCols + 1)];
		for (long nRow = 0; nRow < m_nRows; nRow++){
			for (long nCol = 0; nCol < (m_nCols + 1); nCol++){
				if (nCol < nAt){
					// preserve
					memcpy(av + nRow * (m_nCols + 1) + nCol, m_av + nRow * m_nCols + nCol, sizeof(VARIANT));
					::memset(m_av + nRow * m_nCols + nCol, NULL, sizeof(VARIANT));
				} else if (nCol == nAt){
					// set the new variant to blank
					::VariantInit(av + nRow * (m_nCols + 1) + nCol);					
				} else /*nCol > nAt*/{
					// displace the element in this object by one place to the right
					memcpy(av + nRow * (m_nCols + 1) + nCol, m_av + nRow * m_nCols + nCol - 1, sizeof(VARIANT));
					::memset(m_av + nRow * m_nCols + nCol - 1, NULL, sizeof(VARIANT));
				}
			}
		}
		delete m_av;		
		m_av = av;
		m_nCols++;
		return S_OK;
	} else if (nAt == m_nCols + 1){
		// insert a new column at the end
		return SetColumns(m_nCols + 1);
	} else {		
		return E_FAIL;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	InsertRow
//
//	inserts a row at the specified position in the parameter map
//
HRESULT CParameterMap::InsertRow(long nAt)
{
	if (nAt < 0){
		return E_FAIL;
	} else if (nAt <= m_nRows){
		// insert a new row at position nAt (remember that nAt is zero-based)
		VARIANT* av = new VARIANT[(m_nRows + 1) * m_nCols];		
		// move the first nAt rows from m_av to av (note that the memset call passes ownership without messing up reference counting etc.)	
		::memcpy(av, m_av, m_nCols * nAt * sizeof(VARIANT));						
		::memset(m_av, NULL, m_nCols * nAt * sizeof(VARIANT));		
		// blank the next m_nCols cells in av
		::memset(av + m_nCols * nAt, NULL, m_nCols * sizeof(VARIANT));				
		 // copy the remaining rows (m_nRows - nAt)
		::memcpy(av + m_nCols * (nAt + 1), m_av + m_nCols * nAt, m_nCols * (m_nRows - nAt) * sizeof(VARIANT));
		::memset(m_av + m_nCols * nAt, NULL, m_nCols * (m_nRows - nAt) * sizeof(VARIANT));
		// reassign members of this object		
		delete m_av;		
		m_av = av; 
		m_nRows++;
		return S_OK;
	} else if (nAt == m_nRows + 1){
		// insert a new row at the end
		return SetRows(m_nRows + 1);
	} else {
		return E_FAIL;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	JulianDateToExcelDate
//
//	maps a Julian date to an Excel-style date. 1904 date system is NOT supported.
//
/*static*/ long CParameterMap::JulianDateToExcelDate(unsigned long nJulianDate)
{
	if (nJulianDate <= 2415019 + 60){
		// Microsoft thinks that 1900 was a leap year! Day 60 is 29-Feb-1900.
		return nJulianDate - 2415019 - 1;
	} else {		
		return nJulianDate - 2415019;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	Max
//
//	Returns the maximum number in a parameter map. Fails if no number is found.
//
HRESULT CParameterMap::Max(double *pVal) const
{
	HRESULT								hr = E_FAIL;
	
	*pVal = std::numeric_limits<long>::min();	
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			double f(0.0);
			if (!GetValue(nRow, nCol, &f)){
				hr = S_OK;
				*pVal = __max(*pVal, f);
			}
		}
	}
	return hr;
}


///////////////////////////////////////////////////////////////////////////////
//	Min
//
//	Returns the minimum number in a parameter map. Fails if no number is found.
//
HRESULT CParameterMap::Min(double *pVal) const
{
	HRESULT								hr = E_FAIL;
	
	*pVal = std::numeric_limits<long>::max();	
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			double f(0.0);
			if (!GetValue(nRow, nCol, &f)){
				hr = S_OK;
				*pVal = __min(*pVal, f);
			}
		}
	}
	return hr;
}


//////////////////////////////////////////////////////////////////////////////	
//	ProgIDFromCLSID
//
//	Wraps ::ProdIDFromCLSID
//
/*static*/ HRESULT CParameterMap::ProgIDFromCLSID(const CLSID& clsid, CComBSTR& s)
{
	HRESULT								hr;
	LPOLESTR							pb;
	
	if (hr = ::ProgIDFromCLSID(clsid, &pb)) return hr;
	s = pb;
	CoTaskMemFree(pb);
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveAt
//
//	Removes nElements from the array m_av starting at nStart
//
HRESULT CParameterMap::RemoveAt(long nStart, long nElements)
{
	long								nElement;	
	for (nElement = nStart; nElement < m_nRows * m_nCols - nElements; nElement++){
		// ToDo - can we do this without actually recopying the variant?
		::VariantClear(m_av + nElement);
		::VariantCopy(m_av + nElement, m_av + nElement + nElements);
	}	
	for (nElement = m_nRows * m_nCols - nElements; nElement < m_nRows * m_nCols; nElement++){
		::VariantClear(m_av + nElement);				
	}
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveBlankMargins
//
//	removes any blank edge around the parameter map
//
void CParameterMap::RemoveBlankMargins(bool bLeaveStart)
{
	RemoveBlanksAtEnd();
	RemoveBlanksAtLHS();
	RemoveBlanksAtRHS();
	if (!bLeaveStart) RemoveBlanksAtStart();
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveBlankRows
//
//	removes any blank rows from the parameter map
//
HRESULT CParameterMap::RemoveBlankRows(void)
{
	HRESULT								hr;	
	for (long nRow = m_nRows - 1; nRow >= 0; nRow--){
		if (IsBlank(nRow)){
			if (hr = RemoveRow(nRow)) return hr;
		}
	}
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveBlanks[...]
//
//	remove any blank rows / columns in various parts of the matrix
//
void CParameterMap::RemoveBlanksAtEnd(void)
{
	while (m_nRows && IsBlank(m_nRows - 1, -1)){
		RemoveRow(m_nRows - 1);
	}
}
void CParameterMap::RemoveBlanksAtLHS(void)
{
	while (m_nCols && IsBlank(-1, 0)){
		RemoveColumn(0);
	}
}
void CParameterMap::RemoveBlanksAtRHS(void)
{
	while (m_nCols && IsBlank(-1, m_nCols - 1)){
		RemoveColumn(m_nCols - 1);
	}
}
void CParameterMap::RemoveBlanksAtStart(void)
{
	while (m_nRows && IsBlank(0, -1)){
		RemoveRow(0);
	}
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveCalculationNumberFromHandle
//
//	Transforms a handle of the form a::b::c to a::b
//
/*static*/ HRESULT CParameterMap::RemoveCalculationNumberFromHandle(const CLSID& clsid, std::string* psz, bool* pbIsValidSpreadsheetHandle/* = NULL*/)
//	clsid - object associated with the handle
//  pbIsValidSpreadsheetHandle - if given we perform further checks to decide whether or not the handle looks like a valid spreadsheet handle.
{		
	std::vector<std::string>			vector;							// map of sz
	std::string							szObject;						// object name associated with clsid
	
	if (pbIsValidSpreadsheetHandle) *pbIsValidSpreadsheetHandle = false;

	// first check to see if the handle is of the form a::b::c
	estring::Split(*psz, "::", &vector);
	if (vector.size() != 2 && vector.size() != 3) return S_FALSE;
	
	// the first element should be the object name associated with clsid
	if (clsid != CLSID_NULL){
		if (CLSIDToObjectName(clsid, &szObject)) return S_FALSE;
		if (szObject != vector[0]) return S_FALSE;
	}
	// the third element (if present) should be numeric
	if (vector.size() == 3 && !estring::islong(vector[2])) return S_FALSE;
	
	// We set up the return value here just in case the caller doesn't
	// care about the validitiy of the second element in vector.
	psz->assign(vector[0] + "::" + vector[1]);

	// The second element should look like a cell reference.
	// There are various ways (of varying rigour) to detect this property.
	// The one I have gone for is that the element must contain []!$ in that
	// order. This should trap all but the most perverse user handles.	
	if (pbIsValidSpreadsheetHandle){
		std::string::size_type st;
		if ((st = vector[1].find('[')) != std::string::npos){
			if ((st = vector[1].find(']', st)) != std::string::npos){
				if ((st = vector[1].find('!', st)) != std::string::npos){
					if ((st = vector[1].find('$', st)) != std::string::npos){
						*pbIsValidSpreadsheetHandle = true;
					}
				}
			}
		}
	}
	return vector.size() != 3 ? S_FALSE : S_OK;
}
/*static*/ HRESULT CParameterMap::RemoveCalculationNumberFromHandle(const CLSID& clsid, CComBSTR& ps, bool* pbIsValidSpreadsheetHandle/* = NULL*/)
//  pbIsValidSpreadsheetHandle - if given we perform further checks to decide whether or not the handle looks like a valid spreadsheet handle.
{
	HRESULT								hr;
	estring								sz(ps);

	if (hr = RemoveCalculationNumberFromHandle(clsid, &sz, pbIsValidSpreadsheetHandle)) return hr;	
	return sz.GetBSTR(&ps);
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveColumn
//
//	remove a column from the matrix
//
HRESULT CParameterMap::RemoveColumn(long nCol)
//	nCol - column to remove
{
	long								nRow;							// row under consideration	
	if (nCol < 0 || nCol >= m_nCols) return E_FAIL;
	for (nRow = m_nRows - 1; nRow >= 0; nRow--){
		RemoveAt(nRow * m_nCols + nCol, 1);
	}	
	m_nCols--;
	return S_OK;	
}


///////////////////////////////////////////////////////////////////////////////
//	RemoveRow
//
//	remove a row from the matrix
//
HRESULT CParameterMap::RemoveRow(long nRow)
{
	if (nRow < 0 || nRow >= m_nRows) return E_FAIL;
	RemoveAt(nRow * m_nCols, m_nCols);
	m_nRows--;
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	ReplaceBlanksWithSpaces
//
//	Replaces all empty variants in the map with blank spaces.
//
void CParameterMap::ReplaceBlanksWithSpaces(void)
{
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (GetElementPtr(nRow, nCol)->vt == VT_EMPTY){
				SetValue(nRow, nCol, std::string());
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
//	Resize
//
//	Resizes the parameter map.
//
HRESULT CParameterMap::Resize(long nRows, long nColumns, CComVariant vFill/* = CComVariant()*/)
//	vFill - if the matrix is expanded, then vFill specifies the value to 
//          write to the new area.
{
	HRESULT								hr;
	long								nRowsOld = GetRows();
	long								nColumnsOld = GetCols();

	if (hr = SetRows(nRows)) return hr;
	if (hr = SetColumns(nColumns)) return hr;

	// fill in any new values with vFill if non-blank
	if (vFill.vt == VT_EMPTY) return S_OK;
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = nRow < nRowsOld ? nColumnsOld : 0; nCol < m_nCols; nCol++){
			SetValue(nRow, nCol, vFill);
		}
	}
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	ReturnError
//
//	creates and returns an error object that can be trapped by VBA
//
/*static*/ HRESULT CParameterMap::ReturnErrorS(const std::string& sz, const IID& iid/* =_IID_ISiriusApplication*/)
{
	CComPtr<ICreateErrorInfo>	spcei;
	CComPtr<IErrorInfo>			spei;
	estring						szNew(sz);

	if (szNew.size() && szNew[szNew.size() - 1] != '.' && szNew[szNew.size() - 1] != '!' && szNew[szNew.size() - 1] != '?') szNew += ".";	
	if (!CreateErrorInfo(&spcei)){
		if (spcei->SetDescription((CComBSTR)szNew)) ATLASSERT(false);
		if (spcei->SetGUID(iid)) ATLASSERT(false);
		if (!spcei->QueryInterface(IID_IErrorInfo, (void**)&spei)){							
			if (::SetErrorInfo(0, spei)) ATLASSERT(false);
		}		
#		ifdef using_sirius_application
			_Module.AddError(szNew);
#		endif
	} else {
		ATLASSERT(false);
	}
	return E_FAIL;				
}
/*static*/ HRESULT CParameterMap::ReturnErrorS(const std::stringstream& ss, const IID& iid/* =_IID_ISiriusApplication*/)
{
	return ReturnErrorS(ss.str(), iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorS(LPCTSTR lpsz, const IID& iid/* =_IID_ISiriusApplication*/)
{
	return ReturnErrorS(std::string(lpsz), iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorS(const CComBSTR& s, const IID& iid/* =_IID_ISiriusApplication*/)
{
	return ReturnErrorS(estring(s), iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorSS(std::string sz1, const std::string& sz2, const IID& iid/* =_IID_ISiriusApplication*/)
{	
	if (sz2.size()){
		sz1 += " '";
		sz1 += sz2;
		sz1 += "'";
	}
	return ReturnErrorS(sz1, iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorR(UINT nResourceID, const IID& iid/* =_IID_ISiriusApplication*/)
{
	CComBSTR					s;
	if (!s.LoadString(nResourceID)) ATLASSERT(false);
	return ReturnErrorS(s, iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorRN(UINT nResourceID, long n, const IID& iid/* =_IID_ISiriusApplication*/)
{
	return ReturnErrorRS(nResourceID, estring(n), iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorRP(UINT nResourceID, const CParameterMap& pmIn, const IID& iid/* =_IID_ISiriusApplication*/)
{
	std::string					sz;	
	pmIn.GetString(&sz);
	return ReturnErrorRS(nResourceID, sz, iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorRR(UINT nResourceID1, UINT nResourceID2, const IID& iid/* = _IID_ISiriusApplication*/)
{
	return ReturnErrorRS(nResourceID1, CStringResource(nResourceID2), iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorRS(UINT nResourceID, const std::string& sz, const IID& iid/* =_IID_ISiriusApplication*/)
{
	CComBSTR					s;
	if (!s.LoadString(nResourceID)) ATLASSERT(false);			
	if (sz.size()){		
		if (estring::isdouble(sz)){
			s += L" ";
			s += CComBSTR(sz.c_str());
		} else {				
			s += L" '";
			s += CComBSTR(sz.c_str());
			s += L"'";
		}
	}
	return ReturnErrorS(s, iid);	
}
/*static*/ HRESULT CParameterMap::ReturnErrorRS(UINT nResourceID, const CComBSTR& sIn, const IID& iid/* = _IID_ISiriusApplication*/)
{
	CComBSTR					s;
	if (!s.LoadString(nResourceID)) ATLASSERT(false);
	if (sIn.Length()){
		s += L" '";
		s += sIn;
		s += L"'";
	}
	return ReturnErrorS(s, iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorRSR(UINT nResourceID1, const std::string& sz, UINT nResourceID2, const IID& iid/* =_IID_ISiriusApplication*/)
{
	CComBSTR					s1, s2;
	if (!s1.LoadString(nResourceID1)) ATLASSERT(false);
	if (!s2.LoadString(nResourceID2)) ATLASSERT(false);
	if (sz.size()){
		s1 += L" '";
		s1 += CComBSTR(sz.c_str());
		s1 += L"'.";
	} else {
		s1 += L".";
	}
	s1 += L" ";
	s1 += s2;	
	return ReturnErrorS(s1, iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorSRS(const std::string& sz1, UINT nResourceID, const std::string& sz2, const IID& iid /*=_IID_ISiriusApplication*/)
//	nResourceID must denote a string of the form ...%s...%s...
{
	estring								szFormat;
	char*								szBuffer;
	HRESULT								hr;
	
	if (szFormat.LoadString(nResourceID)) ATLASSERT(false);
	szBuffer = new char[szFormat.size() + sz1.size() + sz2.size() + 1];		// this will be big enough
	sprintf(szBuffer, szFormat.c_str(), sz1.c_str(), sz2.c_str());	
	hr = CParameterMap::ReturnErrorS(szBuffer, iid);
	delete szBuffer;
	return hr;
}
/*static*/ HRESULT CParameterMap::ReturnErrorRV(UINT nResourceID, const CComVariant& vIn, const IID& iid/* =_IID_ISiriusApplication*/)
{
	std::string							sz;
	CParameterMap						pm;
	pm.SetValue(vIn);
	pm.GetString(&sz);
	return ReturnErrorRS(nResourceID, sz, iid);
}
/*static*/ HRESULT CParameterMap::ReturnErrorSV(const std::string& sz, const CComVariant& vIn, const IID& iid/* =_IID_ISiriusApplication*/)
{
	std::string							sz2;
	CParameterMap						pm;
	pm.SetValue(vIn);
	pm.GetString(&sz2);
	return ReturnErrorSS(sz, sz2, iid);
}


///////////////////////////////////////////////////////////////////////////////
//	Round
//
//	1) rounds a double to the nearest whole and returns it
//
//	2) rounds a double to a certain number of decimal places and returns it
//
/*static*/ int CParameterMap::Round(double f)
//	f - number to round
{
	int									n;								// rounded form of f

	n = (int)f;

	if ((abs)((double)n - f) >= 0.5){
		n += sgn(f);
	}
	
	return n;
}
/*static*/ double CParameterMap::Round(double f, unsigned short int nPrecision)
{
	double								fOut;							// return value
	double								fTemp = pow(10.0, (double)nPrecision);
	
	if (!f || !fTemp) return 0.0;
	
	fOut = f * fTemp;
	if (fOut >= 0.0){
		fOut = floor(fOut + 0.5);
	} else {
		fOut = ceil(fOut - 0.5);
	}
	fOut /= fTemp;

	return fOut;
}


///////////////////////////////////////////////////////////////////////////////
//  SafeArrayElementToVariant
//
//	maps an element of a safe array to a variant
//
/*static*/ HRESULT CParameterMap::SafeArrayElementToVariant(SAFEARRAY* psa, long* nIndex, VARIANT* pvRet)
//	psa - safe array pointer
//	nIndex - index vector of the element to obtain
//	nType - data type of the safe array element
//	pvRet - variant to return
{		
	VARTYPE								vt;
	HRESULT								hr;

	if (hr = ::SafeArrayGetVartype(psa, &vt)) return hr;
	switch (vt){
	case VT_UI1:
		{
			BYTE b;
			if (::SafeArrayGetElement(psa, nIndex, &b)) return E_FAIL;
			return CComVariant(b).Detach(pvRet);
		}
		break;
	case VT_I2:
		{
			short ns;
			if (::SafeArrayGetElement(psa, nIndex, &ns)) return E_FAIL;
			return CComVariant(ns).Detach(pvRet);
		}
		break;
	case VT_I4:
		{
			long n;
			if (::SafeArrayGetElement(psa, nIndex, &n)) return E_FAIL;
			return CComVariant(n, VT_I4).Detach(pvRet);
		}		
		break;
	case VT_R8:
		{
			double f;
			if (::SafeArrayGetElement(psa, nIndex, &f)) return E_FAIL;
			return CComVariant(f).Detach(pvRet);		
		}
		break;
	case VT_BOOL:
		{
			bool b;
			if (::SafeArrayGetElement(psa, nIndex, &b)) return E_FAIL;
			return CComVariant(b).Detach(pvRet);
		}
		break;
	case VT_ERROR:
		{
			SCODE sc;
			if (::SafeArrayGetElement(psa, nIndex, &sc)) return E_FAIL;
			return CComVariant(sc).Detach(pvRet);
		}
		break;
	case VT_CY:
		{
			CY cy;
			if (::SafeArrayGetElement(psa, nIndex, &cy)) return E_FAIL;
			return CComVariant(cy).Detach(pvRet);
		}
		break;		
	case VT_DATE:
		{
			DATE date;
			if (::SafeArrayGetElement(psa, nIndex, &date)) return E_FAIL;
			return CComVariant(date).Detach(pvRet);
		}
		break;
	case VT_BSTR:
		{
			CComBSTR s;
			if (::SafeArrayGetElement(psa, nIndex, &s)) return E_FAIL;
			return CComVariant(s).Detach(pvRet);
		}
		break;
	default:
		{
			// Memory leak fixed by Gbenga Olatoye - April 4, 2005.
			CComVariant vRet;
			if (::SafeArrayGetElement(psa, nIndex, &vRet)) return E_FAIL;
			return vRet.Detach(pvRet);
		}
	}	
}


///////////////////////////////////////////////////////////////////////////////
//	SetAutoGrow
//
//	Enable / disable autogrow feature.
//
void CParameterMap::SetAutoGrow(bool b)
{
	m_bAutoGrow = b;
}


///////////////////////////////////////////////////////////////////////////////
//	SetBlank
//
//	if (nRow < 0 and nCol < 0) we blank all the elements in the matrix
//	if (nRow < 0 and nCol >=0) we blank all the elements in column nCol
//	if (nRow >=0 and nCol < 0) we blank all the elements in row nRow
//	if (nRow >=0 and nCol >=0) we blank the element at nRow, nCol
//
HRESULT CParameterMap::SetBlank(long nRow, long nCol)
{
	HRESULT								hr;
	if (nRow < 0){
		ATLASSERT(false);		// ToDo - test
		// set all the rows blank for nCol
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (hr = SetBlank(nRow, nCol)) return false;
		}					
		return S_OK;
	} else if (nCol < 0){
		ATLASSERT(false);		// ToDo - test
		// set all the columns blank for nRow
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (hr = SetBlank(nRow, nCol)) return false;
		}
		return S_OK;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		ATLASSERT(false);		// ToDo - test
		// out of bounds
		return E_FAIL;		
	}

	return ::VariantClear(m_av + nRow * m_nCols + nCol);
}


///////////////////////////////////////////////////////////////////////////////
//	SetCellBlankIfValue
//
//	If the value at an input location equals the value of an input string
//	then we blank it.
//
HRESULT CParameterMap::SetCellBlankIfValue(long nRow, long nCol, std::string szMatch, bool bIgnoreCaseAndSpaces)
{
	HRESULT								hr;
	estring								szValue;

	if (bIgnoreCaseAndSpaces){
		estring::StripWhiteSpace(&szMatch);
		estring::lc(&szMatch);
	}
	if (hr = GetValue(nRow, nCol, &szValue)) return hr;
	if (bIgnoreCaseAndSpaces){
		szValue.StripWhiteSpace();
		szValue.lc();
	}
	if (szValue == szMatch){
		// blank this cell
		if (hr = SetBlank(nRow, nCol)) return hr;
		// if the cell is on the edge of the parameter map then check appropriate remove row / column conditions
		if (nRow == 0) RemoveBlanksAtStart();
		if (nRow == m_nRows - 1) RemoveBlanksAtEnd();
		if (nCol == 0) RemoveBlanksAtLHS();
		if (nCol == m_nCols - 1) RemoveBlanksAtRHS();
	}
	return S_OK;	
}


///////////////////////////////////////////////////////////////////////////////
//	SetColumns
//
//	force the matrix to have a certain number of columns
//
HRESULT CParameterMap::SetColumns(long nCols)
//	nCols - new number of columns in the matrix
{
	if (nCols < 0) return E_FAIL;
		
	if (m_nCols == nCols){
		// already set at this size	
		return S_OK;
	} else if (!m_nRows){
		// blank matrix - all we need to do is set the m_nCols parameter
		m_nCols = nCols;
		return S_OK;
	} else if (nCols < m_nCols){
		// delete the last (m_nCols - nCols) columns on the matrix		
		for (long nRow = m_nRows - 1; nRow >= 0; nRow--){
			RemoveAt(nRow * m_nCols + nCols, m_nCols - nCols);
		}				
		m_nCols = nCols;
		return S_OK;		
	} else {
		// insert new columns on the right hand side		
		VARIANT* av = new VARIANT[m_nRows * nCols];
		for (long nRow = 0; nRow < m_nRows; nRow++){
			for (long nCol = 0; nCol < nCols; nCol++){				
				if (nCol < m_nCols){
					memcpy(av + nRow * nCols + nCol, m_av + nRow * m_nCols + nCol, sizeof(VARIANT));
					::memset(m_av + nRow * m_nCols + nCol, NULL, sizeof(VARIANT));
				} else {
					::VariantInit(av + nRow * nCols + nCol);
				}
			}
		}
		delete m_av;		
		m_av = av;
		m_nCols = nCols;
		return S_OK;		
	}
}


///////////////////////////////////////////////////////////////////////////////
//	SetRows
//
//	force the matrix to have a certain number of rows
//
HRESULT CParameterMap::SetRows(long nRows)
//	nRows - new number of rows in the matrix
{
	if (nRows < 0) return E_FAIL;
		
	if (m_nRows == nRows){
		// already set at this size
		return S_OK;
	} else if (!m_nCols){
		// blank matrix - all we need to do is set the m_nRows parameter
		m_nRows = nRows;
		return S_OK;
	} else if (nRows < m_nRows){
		// delete the last (m_nRows - nRows) rows on the matrix		
		RemoveAt(nRows * m_nCols, (m_nRows - nRows) * m_nCols);		
		m_nRows = nRows;
		return S_OK;
	} else {
		// insert new rows at the end		
		VARIANT* av = new VARIANT[nRows * m_nCols];		
		for (long nRow = 0; nRow < nRows; nRow++){
			for (long nCol = 0; nCol < m_nCols; nCol++){				
				if (nRow < m_nRows){
					memcpy(av + nRow * m_nCols + nCol, m_av + nRow * m_nCols + nCol, sizeof(VARIANT));
					::memset(m_av + nRow * m_nCols + nCol, NULL, sizeof(VARIANT));
				} else {
					::VariantInit(av + nRow * m_nCols + nCol);
				}
			}
		}
		delete m_av;		
		m_av = av;
		m_nRows = nRows;
		return S_OK;
	}
}


///////////////////////////////////////////////////////////////////////////////
//	SetSize
//
//	Set from an input number of rows and columns
//
HRESULT CParameterMap::SetSize(long nRows, long nCols)
{
	Clear();
	m_nRows = nRows;
	m_nCols = nCols;
	m_av = new VARIANT[m_nRows * m_nCols];								// ToDo - memory exception
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			::VariantInit(m_av + nRow * m_nCols + nCol);				// ToDo - out of bounds exception			
			::VariantClear(m_av + nRow * m_nCols + nCol);
		}
	}
	return S_OK;
}


///////////////////////////////////////////////////////////////////////////////
//	SetValue
//
HRESULT CParameterMap::SetValue(const CParameterMap& pm)
{	
	HRESULT								hr;
	Clear();
	
	m_nRows = pm.m_nRows;
	m_nCols = pm.m_nCols;	
	m_av = new VARIANT[m_nRows * m_nCols];							// ToDo - memory exception
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			::VariantInit(m_av + nRow * m_nCols + nCol);
			if (hr = ::VariantCopy(m_av + nRow * m_nCols + nCol, pm.m_av + nRow * m_nCols + nCol)) return hr;
		}
	}
	return S_OK;
}
HRESULT CParameterMap::SetValue(double f)
{
	m_nRows = 1;
	m_nCols = 1;
	m_av = new VARIANT[1];
	::VariantInit(m_av);
	return CComVariant(f).Detach(m_av);
}
HRESULT CParameterMap::SetValue(long nRow, long nCol, const std::stringstream& ss)
{
	return SetValue(nRow, nCol, ss.str());
}
HRESULT CParameterMap::SetValue(const std::string& sz)
{	
	CComBSTR s(sz.c_str());
	return SetValue(s);
}
HRESULT CParameterMap::SetValue(const std::stringstream& ss)
{
	return SetValue(ss.str());
}
HRESULT CParameterMap::SetValue(const BSTR& s)	// ToDo this should be const CComBSTR
{
	CComVariant							v(s);
	Clear();
	m_nRows = 1;
	m_nCols = 1;
	m_av = new VARIANT[1];							// ToDo - memory exception
	::VariantInit(m_av);
	return v.Detach(m_av);
}
//
//	Set the matrix member variables from an input variant. This includes a
//	VB-style variant or an Excel range object.
//
HRESULT CParameterMap::SetValue(const CComVariant& v, bool bRemoveBlankMargins/* = true*/)
{		
	int									nDimensions;
	SAFEARRAY*							psa;						// safe array extracted from the variant
	long								nx1, nx2, ny1, ny2;			// boundaries of matrix (obvious notation)
	long								nIndex[2];					// element in the variant array
	HRESULT								hr;
		
	if (v.vt == (VT_VARIANT | VT_BYREF)){
		return SetValue(*v.pvarVal);
	} else if (v.vt == VT_DISPATCH){
		IID						iid;
		CComDispatchDriverEx	dd(v.pdispVal);
		CComVariant				vValue;								// range value
		long					nTop, nBottom, nLeft, nRight;		// coordinates of dd
		if (!v.pdispVal){
			Clear();
			return S_OK;
		}
		if (hr = GetObjectIID(dd, &iid)) return hr;		
		// We decompose an excel range into its variant value.
		if (iid == IID_ExcelRange){
			if (hr = dd.GetPropertyByName(L"Value", &vValue)) return hr;
			if (hr = SetValue(vValue, false/*don't remove the blank margins yet or we break the incomplete calculation detection code*/)) return hr;
			// Now check that every blank value in the parameter map also has a null formula.
			// If this is not the case then we have incomplete calculation.
			if (m_nRows == 1 && m_nCols == 1 && vValue.vt != VT_EMPTY){
				// no need for incomplete test calculation so do nothing
			} else {
				if (ExcelRangeToRowCol(dd, &nTop, &nLeft, &nBottom, &nRight)){
					ATLTRACE("Cannot test incomplete calculation - CParameterMap::ExcelRangeToRowCol failed\n");
				} else {
					for (long nRow = nTop; nRow <= nBottom; nRow++){
						for (long nCol = nLeft; nCol <= nRight; nCol++){
							std::string szFormula;
							std::string szValue;
							long nElement = (nRow - nTop) * m_nCols + nCol - nLeft;
							ATLASSERT(nElement >= 0 && nElement < m_nRows * m_nCols);
							VARTYPE vt = m_av[nElement].vt;
							if (vt == VT_ERROR) throw "Excel range contains at least one error value";
							if (hr = ExcelGetFormula(dd, nRow - nTop, nCol - nLeft, &szFormula)) return hr;
							if (szFormula.size() && vt == VT_EMPTY){
								// this implies incomplete calculation							
								return E_FAIL;
							}
						}
					}
				}
			}
			if (bRemoveBlankMargins) RemoveBlankMargins(true);
			return S_OK;
		}
	}
	
	Clear();
	if (!(nDimensions = GetArrayDimensions(v))){
		// scalar case
		m_nRows = 1;
		m_nCols = 1;
		m_av = new VARIANT[1];						// ToDo - memory exception
		::VariantInit(m_av);
		if (hr = CComVariant(v).Detach(m_av)) return hr;
	} else {		
		if (v.vt & VT_BYREF){
			psa = *(v.pparray);
		} else {
			psa = v.parray;
		}
		switch (nDimensions){
		case 1:
			// vector		
			SafeArrayGetLBound(psa, 1, &ny1);
			SafeArrayGetUBound(psa, 1, &ny2);		
			m_nCols = 1;
			m_nRows = ny2 - ny1 + 1;
			m_av = new VARIANT[m_nRows];			// ToDo - memory exception
			for (nIndex[0] = ny1; nIndex[0] <= ny2; nIndex[0]++){
				::VariantInit(m_av + (nIndex[0] - ny1));
				if (hr = SafeArrayElementToVariant(psa, nIndex, m_av + (nIndex[0] - ny1))) return hr;
			}
			break;
		case 2:
			// matrix	
			SafeArrayGetLBound(psa, 1, &ny1);
			SafeArrayGetUBound(psa, 1, &ny2);
			SafeArrayGetLBound(psa, 2, &nx1);
			SafeArrayGetUBound(psa, 2, &nx2);
			m_nRows = ny2 - ny1 + 1;
			m_nCols = nx2 - nx1 + 1;
			m_av = new VARIANT[m_nRows * m_nCols];		// ToDo - memory exception
			for (nIndex[0] = ny1; nIndex[0] <= ny2; nIndex[0]++){
				for (nIndex[1] = nx1; nIndex[1] <= nx2; nIndex[1]++){
					CComVariant vElement;
					if (hr = SafeArrayElementToVariant(psa, nIndex, &vElement)) return hr;
					::VariantInit(m_av + (nIndex[0] - ny1) * m_nCols + (nIndex[1] - nx1));
					if (hr = vElement.Detach(m_av + (nIndex[0] - ny1) * m_nCols + (nIndex[1] - nx1))) return hr;
				}		
			}
			break;
		default:
			return E_FAIL;
		}
	}
	if (bRemoveBlankMargins) RemoveBlankMargins(true);
	return S_OK;
}

HRESULT CParameterMap::SetValue(const CMatrix& m)
//	ppf - array to copy. If this is NULL then we allocate the memory in the parameter map and clear the variants
{	
	Clear();
	m_nRows = m.rows();
	m_nCols = m.cols();
	m_av = new VARIANT[m_nRows * m_nCols];								// ToDo - memory exception
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			::VariantInit(m_av + nRow * m_nCols + nCol);				// ToDo - out of bounds exception			
			if (m.rows()){
				(m_av + nRow * m_nCols + nCol)->vt = VT_R8;
				(m_av + nRow * m_nCols + nCol)->dblVal = m[nRow][nCol];			
			} else {
				::VariantClear(m_av + nRow * m_nCols + nCol);
			}
		}
	}
	return S_OK;
}

HRESULT CParameterMap::SetValue(const CVector& m)
//	m - array to copy. If this is NULL then we allocate the memory in the parameter map and clear the variants
{	
	Clear();
	m_nRows = m.getsize();
	m_nCols = 1;
	m_av = new VARIANT[m_nRows * m_nCols];								// ToDo - memory exception
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			::VariantInit(m_av + nRow * m_nCols + nCol);				// ToDo - out of bounds exception			
			if (m.getsize()){
				(m_av + nRow * m_nCols + nCol)->vt = VT_R8;
				(m_av + nRow * m_nCols + nCol)->dblVal = m[nRow];			
			} else {
				::VariantClear(m_av + nRow * m_nCols + nCol);
			}
		}
	}
	return S_OK;
}
HRESULT CParameterMap::SetValue(const iVector& iv)
{
	Clear();
	m_nRows = iv.getsize();
	m_nCols = 1;	
	m_av = new VARIANT[m_nRows * m_nCols];								// ToDo - memory exception
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			::VariantInit(m_av + nRow * m_nCols + nCol);				// ToDo - out of bounds exception
			if (iv.getsize()){
				(m_av + nRow * m_nCols + nCol)->vt = VT_I4;
				(m_av + nRow * m_nCols + nCol)->dblVal = (long)iv[nRow];
			} else {
				::VariantClear(m_av + nRow * m_nCols + nCol);
			}
		}
	}
	return S_OK;
}
//
//	Inserts a value into the matrix at a given row and column subject to:
//	if (nRow < 0 and nCol < 0) then every cell in the matrix is set to the value.
//	if (nRow < 0 and nCol >=0) then every cell in column nCol is set to the value.
//	if (nRow >=0 and nCol < 0) then every cell in row nRow is set to the value.
//
HRESULT	CParameterMap::SetValue(long nRow, long nCol, const CComVariant& v)
{
	HRESULT								hr;
	if (nRow < 0){
		// set all the rows in column nCol with value v
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (hr = SetValue(nRow, nCol, v)) return hr;
		}					
		return S_OK;
	} else if (nCol < 0){
		// set all the columns in row nRow with value v
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (hr = SetValue(nRow, nCol, v)) return hr;
		}
		return S_OK;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){
		// out of bounds
		if (!m_bAutoGrow || nRow < 0 || nCol < 0) return ReturnErrorS("Array index out of bounds");
		if (nRow >= m_nRows && (hr = SetRows(nRow + 1))) return hr;
		if (nCol >= m_nCols && (hr = SetColumns(nCol + 1))) return hr;
		// no return here is deliberate - we want to follow on to the next statement
	}		
	return CComVariant(v).Detach(m_av + nRow * m_nCols + nCol);	
}
//
//	Inserts a string value into the matrix
//
HRESULT	CParameterMap::SetValue(long nRow, long nCol, const std::string& sz)
{
	return SetValue(nRow, nCol, CComVariant(sz.c_str()));
}
//
//	Inserts a BSTR into the matrix
//
HRESULT CParameterMap::SetValue(long nRow, long nCol, const CComBSTR& s)
{
	return SetValue(nRow, nCol, CComVariant(s));
}
//
//	Inserts a double into the matrix
//
HRESULT CParameterMap::SetValue(long nRow, long nCol, double f)
{
	return SetValue(nRow, nCol, CComVariant(f));
}
//
//.Inserts a string denoted by an input resoure ID into the matrix
//
HRESULT CParameterMap::SetValue(long nRow, long nCol, const CStringResource& sr)
{
	return SetValue(nRow, nCol, sr.str());
}
//
//	This function is odd. The element at nRow, nCol is set to the
//	first element of vector. Subsequent elements are set in
//	subsequent columns starting at n.
//
HRESULT CParameterMap::SetValue(long nRow, long nCol, const CVector& v)
{
	HRESULT								hr;
	ATLASSERT(false);
	// ToDo - test
	SetColumns(__max(m_nCols, nCol + v.getsize() - 1));
	for (long n = 0; n < v.getsize(); n++){
		if (hr = SetValue(nRow, nCol + n, v[n])) return hr;
	}
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	SetValues
//
//	This function allows the setting of many cells in the parameter map from
//	a string input. A delimiter is also supplied to enable separation of the
//	input string. The variant type parameter vt specifies the conversion that
//	should be made to each extracted field before inserting it into the
//	paramter map object
//
HRESULT CParameterMap::SetValues(VARTYPE vt, const std::string& szDelimit, const std::string& szValues)
{
	std::vector<std::string>			vector;							// map of szValues
	HRESULT								hr;
	int									nElement = 0;					// elemen of the parameter map
	
	estring::Split(szValues, szDelimit, &vector);
	ATLASSERT(vector.size() <= m_nRows * m_nCols);
	for (std::vector<std::string>::iterator it = vector.begin(); it < vector.end(); it++){
		CComVariant v(it->c_str());
		if (hr = v.ChangeType(vt)) return hr;		
		if (hr = v.Detach(m_av + nElement++)) return hr;
	}
	return S_OK;	
}
//
//	This function allows the setup of a parameter map from a dimension input
//  and an argument list
//
HRESULT CParameterMap::SetValues(long nRows, long nCols, /*CComVariant*/...)
{
	HRESULT								hr;
	va_list								vl;
		
	if (hr = SetSize(nRows, nCols)) return hr;
	va_start(vl, nCols);
	for (long nRow = 0; nRow < nRows; nRow++){
		for (long nCol = 0; nCol < nCols; nCol++){			
			SetValue(nRow, nCol, va_arg(vl, CComVariant));
		}
	}	
	va_end(vl);
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	Split
//
//	Splits this parameter map into n distinct parameter map objects. We do
//	this by splitting each (string) element of this parameter map according to
//	an input delimiter. Then each element is inserted into a new parameter
//	map object.
//
//	We error if any element this parameter map has a different number of
//	delimited strings or is not convertible to a string
//
HRESULT CParameterMap::Split(const std::string& szDelimit, VARTYPE vt, int nParameterMaps, /*CParameterMap**/...)
//	szDelimit - delimiter used to break up each cell in this parameter map
//	vt - variant type to which we set each element in each new parameter map
//	nParameterMaps - number of parameter maps passed in
{
	HRESULT								hr = S_OK;		
	std::vector<std::string>			vector;
	CComVariant							vElement;									// element of this parameter map
	int									nElementsPerCell;							// number of elements per parameter map cell
	CParameterMap						ppm[16];									// output parameter maps
	int									nParameterMap;								// element of ppm
	va_list								vl;
	
	ATLASSERT(nParameterMaps <= 16);
	for (long nRow = 0; nRow < m_nRows && !hr; nRow++){
		for (long nCol = 0; nCol < m_nCols && !hr; nCol++){
			if (hr = vElement.ChangeType(VT_BSTR, m_av + nRow * m_nCols + nCol)) continue;
			estring sz(vElement.bstrVal);
			sz.Split(szDelimit, &vector);
			if (!nRow && !nCol){
				nElementsPerCell = vector.size();
				for (nParameterMap = 0; nParameterMap < __min(nParameterMaps, nElementsPerCell); nParameterMap++){
					ppm[nParameterMap].SetSize(m_nRows, m_nCols);
				}
			} else {
				if (vector.size() != nElementsPerCell){
					hr = E_FAIL;
					continue;
				}
			}
			for (nParameterMap = 0; nParameterMap < __min(nParameterMaps, nElementsPerCell) && !hr; nParameterMap++){
				if (hr = vElement.ChangeType(vt, &CComVariant(vector[nParameterMap].c_str()))) continue;
				if (hr = ppm[nParameterMap].SetValue(nRow, nCol, vElement)) continue;
			}
		}
	}

	if (!hr){
		// set the output parameter maps		
		va_start(vl, nParameterMaps);
		for (nParameterMap = 0; nParameterMap < __min(nParameterMaps, nElementsPerCell); nParameterMap++){
			va_arg(vl, CParameterMap*)->Attach(ppm + nParameterMap);
		}
		va_end(vl);
	}	
	return hr;
}


/////////////////////////////////////////////////////////////////////////////
//	StringToDate
//
//	performs the reverse operation of DateToString
//
/*static*/ HRESULT CParameterMap::StringToDate(const std::string& sz, DATE* pdate)
{
	HRESULT								hr;
	CComVariant							v = estring::GetValue(sz);

	if (hr = v.ChangeType(VT_DATE)){		
		// Before throwing in the towel, try converting to VT_R8 then to VT_DATE
		if (hr = v.ChangeType(VT_R8)){
			return hr;
		}
		if (hr = v.ChangeType(VT_DATE)){					
			return hr;
		}
	}
	*pdate = v.date;
	return S_OK;
}
/*static*/ HRESULT CParameterMap::StringToDate(const std::string& sz, long* pnDate)
{
	HRESULT								hr;
	DATE								date;

	if (hr = StringToDate(sz, &date)) return hr;
	*pnDate = (long)date;
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	Substitute
//
//	Replaces one string in a parameter map with another one.
//
HRESULT CParameterMap::Substitute(const std::string& szOld, const std::string& szNew, bool bPartString, bool bIgnoreCaseAndSpaces, long nRow, long nCol)
{
	HRESULT								hr;	

	if (!m_nRows || !m_nCols) return S_OK;
	if (nRow < 0){		
		for (long nRow = 0; nRow < m_nRows; nRow++){
			if (hr = Substitute(szOld, szNew, bPartString, bIgnoreCaseAndSpaces, nRow, nCol)) return hr;
		}					
		return S_OK;
	} else if (nCol < 0){		
		for (long nCol = 0; nCol < m_nCols; nCol++){
			if (hr = Substitute(szOld, szNew, bPartString, bIgnoreCaseAndSpaces, nRow, nCol)) return hr;
		}
		return S_OK;
	} else if (nRow < 0 || nRow >= m_nRows || nCol < 0 || nCol >= m_nCols){		
		return E_FAIL;
	}	
	if (GetElementPtr(nRow, nCol)->vt != VT_BSTR) return S_OK;
	
	estring sz;
	if (GetValue(nRow, nCol, &sz)) return S_OK;		// nothing to substitute
	if (bPartString){
		ATLASSERT(!bIgnoreCaseAndSpaces);
		if (sz.ReplaceStrInStr(szOld, szNew)){
			// At least one substitution was made.
			return SetValue(nRow, nCol, sz);
		}
	} else {		
		if (bIgnoreCaseAndSpaces && !sz.CompareNoCaseAndSpace(szOld) || !bIgnoreCaseAndSpaces && sz == szOld){
			return SetValue(nRow, nCol, szNew);
		}
	}
	return S_OK;
}


//////////////////////////////////////////////////////////////////////////////
//	ThrowComError[...]
//
//	throws a _com_error object
//
/*static*/ void CParameterMap::ThrowComErrorR(UINT nResourceID, const IID& iid/* =_IID_ISiriusApplication*/)
{
	CComBSTR					s;
	if (!s.LoadString(nResourceID)) ATLASSERT(false);
	ThrowComErrorS(s, iid);
}
/*static*/ void CParameterMap::ThrowComErrorRR(UINT nResourceID1, UINT nResourceID2, const IID& iid/* = _IID_ISiriusApplication*/)
{
	ThrowComErrorRS(nResourceID1, CStringResource(nResourceID2), iid);
}
/*static*/ void CParameterMap::ThrowComErrorRS(UINT nResourceID, const std::string& sz, const IID& iid/* =_IID_ISiriusApplication*/)
{						
	CComBSTR					s;
	if (!s.LoadString(nResourceID)) ATLASSERT(false);
	if (sz.size()){
		s += L" '";
		s += CComBSTR(sz.c_str());
		s += L"'.";
	}
	ThrowComErrorS(s, iid);
}
/*static*/ void CParameterMap::ThrowComErrorRV(UINT nResourceID, const CComVariant& vIn, const IID& iid/* =_IID_ISiriusApplication*/)
{
	ThrowComErrorRS(nResourceID, estring(vIn), iid);	
}
/*static*/ void CParameterMap::ThrowComErrorS(const CComBSTR& s, const IID& iid/* =_IID_ISiriusApplication*/)
{
	CComPtr<ICreateErrorInfo>	spcei;
	CComPtr<IErrorInfo>			spei;
							
	if (!CreateErrorInfo(&spcei)){
		if (spcei->SetDescription(s)) ATLASSERT(false);
		if (spcei->SetGUID(iid)) ATLASSERT(false);
		if (!spcei->QueryInterface(IID_IErrorInfo, (void**)&spei)){
			if (::SetErrorInfo(0, spei)) ATLASSERT(false);
		}
#		ifdef using_sirius_application
			_Module.AddError(estring(s));
#		endif
	} else {
		ATLASSERT(false);
	}	
	throw _com_error(E_FAIL, spei, true);
}


/////////////////////////////////////////////////////////////////////////////
//	ThrowStringError
//
//	Throws a string error from a populated error object.
//
/*static*/ void CParameterMap::ThrowStringError(CComPtr<IDispatch> spObject)
{
	CComPtr<IErrorInfo>			spErrorInfo;
	CComPtr<ISupportErrorInfo>	spSupportErrorInfo;	
	CComBSTR					sError;
			
	if (!::GetErrorInfo(0, &spErrorInfo)){
		if (!spObject){
			// set it to the application object
			spObject.CoCreateInstance(_CLSID_SiriusApplication);
		}		
		if (spObject && !spObject->QueryInterface(IID_ISupportErrorInfo, (void**)&spSupportErrorInfo)){
			if (!spErrorInfo->GetDescription(&sError)) throw estring(sError);
		}
	}
	
	throw "Unhandled exception - please report this to the development team";
}


//////////////////////////////////////////////////////////////////////////////
//	Transpose
//
//	sets the parameter map to its transpose
//
HRESULT CParameterMap::Transpose(void)
{	
	VARIANT*							av;

	av = new VARIANT[m_nCols * m_nRows];								// ToDo - memory exception		
	for (long nRow = 0; nRow < m_nRows; nRow++){
		for (long nCol = 0; nCol < m_nCols; nCol++){
			::memcpy(av + nCol * m_nRows + nRow, m_av + nRow * m_nCols + nCol, sizeof(VARIANT));
			::memset(m_av + nRow * m_nCols + nCol, NULL, sizeof(VARIANT));
		}
	}
	delete m_av;
	m_av = av;			
	std::swap(m_nRows, m_nCols);
	return S_OK;	
}


//////////////////////////////////////////////////////////////////////////////
//	VariableArgumentListToArray
//
//	Converts a variable argument list of variants into a 1D safearray. If there
//	is only 1 variable argument then we assign the output to this directly.
//
/*static*/ HRESULT CParameterMap::VariableArgumentListToArray(CComVariant* pvOut, long nArgs, va_list vl)
{
	HRESULT								hr;	
	CComVariant							vElement;						// one of the variable argument list values
	SAFEARRAY*							psa;							// this is used to form the variant array
	long								nIndex;							// element in the variant array
	SAFEARRAYBOUND						sab;							// safe array size specification vector
	long								nSize;							// current size of the array
	long								nLastNonBlank;					// last non-blank element in the variant array
	
	// create a suitable safe array
	sab.lLbound = 1;
	sab.cElements = 1;
	if (!(psa = ::SafeArrayCreate(VT_VARIANT, 1, &sab))) return E_FAIL;
	hr = S_OK;
	nIndex = 0;
	nLastNonBlank = 0;

	// read the variable argument values into the safe array	
	while (!hr && nArgs--){
		vElement.Copy(&va_arg(vl, VARIANT));
		nIndex++;
		if (vElement.vt == VT_ERROR){
			// empty the variant
			vElement.Clear();
		} else {
			// set the watermark
			nLastNonBlank = nIndex;
		}
		if (hr = ::SafeArrayGetUBound(psa, 1, &nSize)) break;
		if (nIndex > nSize){
			// we need more space in the safe array - double it
			sab.cElements = nSize * 2;
			if (hr = ::SafeArrayRedim(psa, &sab)) break;
		}
		if (hr = ::SafeArrayPutElement(psa, &nIndex, &vElement)) break;		
	}	
	if (hr){
		if (::SafeArrayDestroy(psa)) ATLASSERT(FALSE);
		return hr;
	}

	if (!nLastNonBlank){
		// no data were passed in
		if (::SafeArrayDestroy(psa)) ATLASSERT(FALSE);
		pvOut->Clear();
		return S_OK;
	} else if (nLastNonBlank == 1){
		// only one element - detach this element from the safe array		
		hr = ::SafeArrayGetElement(psa, &nLastNonBlank, &vElement);
		if (::SafeArrayDestroy(psa)) ATLASSERT(FALSE);
		return pvOut->Attach(&vElement);
	} else {
		sab.cElements = nLastNonBlank;
		hr = ::SafeArrayRedim(psa, &sab);
		// ToDo - newer versions of ATL have a CComVariant constructor / assignment operator that takes a SAFEARRAY pointer
		pvOut->Clear();
		pvOut->vt = VT_ARRAY | VT_VARIANT;
		pvOut->parray = psa;
		return hr;
	}
}
/*static*/ HRESULT CParameterMap::VariableArgumentListToArray(CComVariant* pvOut, long nArgs, /*VARIANT*/...)
//	pvOut - output safearray
{
	HRESULT								hr;
	va_list								vl;	
	va_start(vl, nArgs);
	hr = VariableArgumentListToArray(pvOut, nArgs, vl);
	va_end(vl);	
	return hr;	
}
/*static*/ HRESULT CParameterMap::VariableArgumentListToArray(VARIANT* pvOut, long nArgs, /*CParameterMap*/...)
//	pvOut - output safearray
//	nArgs - number of arguments
{
	HRESULT								hr;
	va_list								vl;
	SAFEARRAY*							psa;							// this is used to form the variant array
	long								nIndex;							// element in the variant array
	SAFEARRAYBOUND						sab;							// safe array size specification vector

	// create a suitable safe array
	sab.lLbound = 1;
	sab.cElements = nArgs;
	if (!(psa = ::SafeArrayCreate(VT_VARIANT, 1, &sab))) return E_FAIL;
	hr = S_OK;
	nIndex = 0;

	// read the variable argument values into the safe array
	va_start(vl, nArgs);
	while (!hr && nArgs--){
		CComVariant vElement;
		if (!(hr = va_arg(vl, CParameterMap).GetValue(&vElement))){
			nIndex++;
			hr = ::SafeArrayPutElement(psa, &nIndex, &vElement);
		}
	}
	va_end(vl);
	if (!nIndex){
		// no data were passed in
		if (::SafeArrayDestroy(psa)) ATLASSERT(FALSE);
		::VariantClear(pvOut);
		return S_OK;
	} else if (nIndex == 1){
		// only one element - detach this element from the safe array
		CComVariant vElement;
		hr = ::SafeArrayGetElement(psa, &nIndex, &vElement);
		if (::SafeArrayDestroy(psa)) ATLASSERT(FALSE);
		return vElement.Detach(pvOut);
	} else {
		::VariantClear(pvOut);
		pvOut->vt = VT_ARRAY | VT_VARIANT;
		pvOut->parray = psa;
		return S_OK;
	}
}


//////////////////////////////////////////////////////////////////////////////
//	VectorToArray
//
//	Maps a vector of objects to a 1D safe array where each element is the
//	variant associated with the corresponding parameter map value.
//
/*static*/ HRESULT CParameterMap::VectorToArray(std::vector<CComVariant>& vv, VARIANT* pvOut)
{
	if (!vv.size()){
		// null input		
		::VariantClear(pvOut);
		return S_OK;
	} else if (vv.size() == 1){
		// set pvOut directly to the input
		__asm int 3;
		return CComVariant(vv[0]).Detach(pvOut);
	} else {
		SAFEARRAY*		psa;
		SAFEARRAYBOUND	sab;
		long			nIndex = 0;		
		HRESULT			hr = S_OK;
				
		sab.lLbound = 1;
		sab.cElements = vv.size();
		if (!(psa = ::SafeArrayCreate(VT_VARIANT, 1, &sab))) return E_FAIL;
		for (std::vector<CComVariant>::iterator it = vv.begin(); it < vv.end() && !hr; it++){
			nIndex++;			
			hr = ::SafeArrayPutElement(psa, &nIndex, it);
		}
		if (hr){
			if (::SafeArrayDestroy(psa)) ATLASSERT(false);
			return hr;
		} else {
			::VariantClear(pvOut);
			pvOut->vt = VT_ARRAY | VT_VARIANT;
			pvOut->parray = psa;
			return S_OK;
		}
	}
}
/*static*/ HRESULT CParameterMap::VectorToArray(std::vector<CParameterMap*>& vppm, CComVariant* pvOut)
{
	if (!vppm.size()){
		// null input
		ATLASSERT(false);
		pvOut->Clear();
		return S_OK;
	} else if (vppm.size() == 1){
		// set pvOut directly to the only parameter map input
		return vppm[0]->GetValue(pvOut);
	} else {
		SAFEARRAY*		psa;
		SAFEARRAYBOUND	sab;
		long			nIndex = 0;
		CComVariant		vElement;
		HRESULT			hr = S_OK;
				
		sab.lLbound = 1;
		sab.cElements = vppm.size();
		if (!(psa = ::SafeArrayCreate(VT_VARIANT, 1, &sab))) return E_FAIL;
		for (std::vector<CParameterMap*>::const_iterator it = vppm.begin(); it < vppm.end() && !hr; it++){			
			if (!(hr = (*it)->GetValue(&vElement))){
				nIndex++;
				hr = ::SafeArrayPutElement(psa, &nIndex, &vElement);
			}
		}
		if (hr){
			if (::SafeArrayDestroy(psa)) ATLASSERT(false);
			return hr;
		} else {
			pvOut->Clear();
			pvOut->vt = VT_ARRAY | VT_VARIANT;
			pvOut->parray = psa;
			return S_OK;
		}
	}
}
/*static*/ HRESULT CParameterMap::VectorToArray(std::vector<CParameterMap>& vpm, VARIANT* pvOut)
{	
	if (!vpm.size()){
		// null input		
		::VariantClear(pvOut);
		return S_OK;
	} else if (vpm.size() == 1){
		// set pvOut directly to the only parameter map input
		return vpm[0].GetValue(pvOut);
	} else {
		SAFEARRAY*		psa;
		SAFEARRAYBOUND	sab;
		long			nIndex = 0;
		CComVariant		vElement;
		HRESULT			hr = S_OK;
				
		sab.lLbound = 1;
		sab.cElements = vpm.size();
		if (!(psa = ::SafeArrayCreate(VT_VARIANT, 1, &sab))) return E_FAIL;
		for (std::vector<CParameterMap>::const_iterator it = vpm.begin(); it < vpm.end() && !hr; it++){			
			if (!(hr = it->GetValue(&vElement))){
				nIndex++;
				hr = ::SafeArrayPutElement(psa, &nIndex, &vElement);
			}
		}
		if (hr){
			if (::SafeArrayDestroy(psa)) ATLASSERT(false);
			return hr;
		} else {
			::VariantClear(pvOut);	
			pvOut->vt = VT_ARRAY | VT_VARIANT;
			pvOut->parray = psa;
			return S_OK;
		}
	}	
}


/////////////////////////////////////////////////////////////////////////////
//	WindowToClassName
//
//	obtains the window class name associated with an input HWND
//
/*static*/ HRESULT CParameterMap::WindowToClassName(HWND hWnd, std::string* pszOut)
//	hWnd - window to search
{	
	int									nBufferUsed;
	int									nBufferLength = 512;
	
	do {
		nBufferLength *= 2;
		LPTSTR psz = new TCHAR[nBufferLength];
		nBufferUsed = ::GetClassName(hWnd, psz, nBufferLength);
		if (nBufferUsed != nBufferLength - 1) pszOut->assign(psz);
		delete psz;
	} while (nBufferUsed == nBufferLength - 1);	
	return S_OK;
}

#undef FIND_CELL_HELPER