//	parametermap.h: Interface for auto_ptr_ex
//							      CStringResource
//								  estring
//								  CParameterMap
//
//					Also defines several variable mapping functions
//					that exploit CParameterMap.
//
//	Author:		    David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifndef _PARAMETERMAP_H
#define _PARAMETERMAP_H

#pragma once
#include <cmatrix.h>					// for CVector and CMatrix classes
#include <set>
class CParameterMap;


//////////////////////////////////////////////////////////////////////
//	Private Macros
//
//	Macros prefixed with _ should not be used explicitly as their
//  parameters are likely to change.
//
#define _map_parameter_err(InputName, MappedName)												\
	{																							\
		CParameterMap			_pmValue;														\
		std::string				szValue;														\
		std::stringstream		ss;																\
		_pmValue.SetValue(InputName);															\
		_pmValue.GetString(&szValue);															\
		estring szName(#InputName);																\
		if (szName.find('[') != std::string::npos){												\
			szName.assign(estring(#MappedName));												\
		}																						\
		if (szValue.size()){																	\
			ss << "Parameter '" << szName << "' could not be set to value '" << szValue << "'";	\
		} else {																				\
			ss << "Invalid value for parameter '" << szName << "'";								\
		}																						\
		throw ss.str();																			\
	}


//////////////////////////////////////////////////////////////////////
//	[begin/end]_function 
//
//  Prefix and suffix all COM functions with these macros to implement
//	excel-wizard checks and exception handling.
//
#define begin_function														\
	{																		\
		try {

#define end_function														\
		} catch (LPCTSTR lpsz){												\
			return CParameterMap::ReturnErrorS(lpsz);						\
		} catch (const std::string& sz) {									\
			return CParameterMap::ReturnErrorS(sz);							\
		} catch (const _com_error &e){										\
			return CParameterMap::ReturnErrorS(estring(e));					\
		} catch (const std::bad_alloc& rcException) {						\
			return CParameterMap::ReturnErrorS(rcException.what());			\
		} catch (const std::exception& rcException) {						\
			return CParameterMap::ReturnErrorS(rcException.what());			\
		} catch (const CStringResource& sr){								\
			return CParameterMap::ReturnErrorS(sr.str());					\
		} catch (const CComBSTR& s){										\
			return CParameterMap::ReturnErrorS(estring(s));					\
		} catch (...){														\
			TCHAR szError[512];												\
			sprintf(szError, _T("Unhandled exception in %s (%d)! Please report to the development team."), _T(__FILE__), __LINE__); \
			return CParameterMap::ReturnErrorS(CComBSTR(szError));			\
		}																	\
	}																		\
	return S_OK;


//////////////////////////////////////////////////////////////////////
//	map_
//
//	Use these to map input variants to useful data types.
//	Data validation is performed and exceptions (catchable by a
//  begin_function / end_function block) can be thrown.
//
#define map_parameter(InputName, Type, MappedName)					\
	Type			MappedName;										\
	{																\
		CParameterMap		_pm;									\
		if (_pm.SetValue(InputName) || _pm.GetValue(&MappedName)){	\
			_map_parameter_err(InputName, MappedName);				\
		}															\
	}

#define map_parameter_no_err(InputName, Type, MappedName)			\
	Type			MappedName;										\
	{																\
		CParameterMap		_pm;									\
		_pm.SetValue(InputName);									\
		_pm.GetValue(&MappedName);									\
	}

#define map_optional_parameter(InputName, Type, MappedName, DefaultValue)				\
	Type		MappedName = DefaultValue;												\
	{																					\
		CParameterMap	_pm;															\
		if (_pm.SetValue(InputName) || (!_pm.IsBlank() && _pm.GetValue(&MappedName))){	\
			_map_parameter_err(InputName, MappedName);									\
		}																				\
	}

#define map_optional_parameter_no_err(InputName, Type, MappedName, DefaultValue)		\
	Type		MappedName = DefaultValue;												\
	{																					\
		CParameterMap	_pm;															\
		_pm.SetValue(InputName);														\
		_pm.GetValue(&MappedName);														\
	}

#define map_enum_parameter(InputName, Type, MappedName)							\
	Type			MappedName;													\
	{																			\
		if (CEnumMap::GetEnum(#Type, LIBID_Sirius, InputName, &MappedName)){	\
			_map_parameter_err(InputName, MappedName);							\
		}																		\
	}

#define map_optional_enum_parameter(InputName, Type, MappedName, DefaultValue)					\
	Type			MappedName = DefaultValue;													\
	{																							\
		if (CEnumMap::GetEnum(#Type, LIBID_Sirius, InputName, DefaultValue, &MappedName)){		\
			_map_parameter_err(InputName, MappedName);											\
		}																						\
	}

#define map_enum_vector_parameter(InputName, Type, MappedName)				\
	std::vector<Type>	MappedName;											\
	{																		\
		CParameterMap	_pm;												\
		if (_pm.SetValue(InputName) || !_pm.IsVector()){					\
			_map_parameter_err(InputName, MappedName);						\
		}																	\
		if (_pm.GetCols() != 1){											\
			_pm.Transpose();												\
		}																	\
		for (long nRow = 0; nRow < _pm.GetRows(); nRow++){					\
			std::string sz;													\
			Type		type;												\
			if (_pm.GetValue(nRow, 0, &sz)){								\
				_map_parameter_err(InputName, MappedName);					\
			}																\
			if (CEnumMap::GetEnum(#Type, LIBID_Sirius, sz, &type)){			\
				estring szName(#InputName);									\
				if (szName.find('[') != std::string::npos){					\
					szName.assign(estring(#MappedName));					\
				}															\
				throw "'" + sz + "' is an invalid value for an element of parameter '" + szName + "'";		\
			}																\
			MappedName.push_back(type);										\
		}																	\
	}

#define map_optional_enum_vector_parameter(InputName, Type, MappedName)			\
	std::vector<Type>	MappedName;												\
	{																			\
		CParameterMap	_pm;													\
		if (_pm.SetValue(InputName)){											\
			_map_parameter_err(InputName, MappedName);							\
		}																		\
		if (!_pm.IsBlank()){													\
			if (!_pm.IsVector()){												\
				_map_parameter_err(InputName, MappedName);						\
			}																	\
			if (_pm.GetCols() != 1){											\
				_pm.Transpose();												\
			}																	\
			for (long nRow = 0; nRow < _pm.GetRows(); nRow++){					\
				std::string sz;													\
				Type		type;												\
				if (_pm.GetValue(nRow, 0, &sz)){								\
					_map_parameter_err(InputName, MappedName);					\
				}																\
				if (CEnumMap::GetEnum(#Type, LIBID_Sirius, sz, &type)){			\
					estring szName(#InputName);									\
					if (szName.find('[') != std::string::npos){					\
						szName.assign(estring(#MappedName));					\
					}															\
					throw "'" + sz + "' is an invalid value for an element of parameter '" + szName + "'";		\
				}																\
				MappedName.push_back(type);										\
			}																	\
		}																		\
	}


//////////////////////////////////////////////////////////////////////
//	unmap_parameter
//
//  Performs the reverse operation to the map_parameter style
//  macros. The output variable is a parameter map object.
//
#define unmap_parameter(MappedName, OutputName)					\
	CParameterMap				OutputName;						\
	{															\
		HRESULT hr;												\
		if (hr = OutputName.SetValue(MappedName)) return hr;	\
	}

																				
//////////////////////////////////////////////////////////////////////
//	auto_ptr_ex
//
//	This is an implementation of std::auto_ptr that can be used
//  easily by CAdapt.
//
template<class _Ty>
class auto_ptr_ex
{
public:
	typedef _Ty element_type;
	explicit auto_ptr_ex(_Ty *_P = 0) _THROW0()
		: _Owns(_P != 0), _Ptr(_P) {}
	auto_ptr_ex(const auto_ptr_ex<_Ty>& _Y) _THROW0()
		: _Owns(_Y._Owns), _Ptr(_Y.release()) {}
	auto_ptr_ex<_Ty>& operator=(const auto_ptr_ex<_Ty>& _Y) _THROW0()
		{if (this != &_Y)
			{if (_Ptr != _Y.get())
				{if (_Owns)
					delete _Ptr;
				_Owns = _Y._Owns; }
			else if (_Y._Owns)
				_Owns = true;
			_Ptr = _Y.release(); }
		return (*this); }
	~auto_ptr_ex()
		{if (_Owns)
			delete _Ptr; }
	_Ty& operator*() const _THROW0()
		{return (*get()); }
	_Ty *operator->() const _THROW0()
		{return (get()); }
	_Ty *get() const _THROW0()
		{return (_Ptr); }
	_Ty *release() const _THROW0()
		{((auto_ptr_ex<_Ty> *)this)->_Owns = false;
		return (_Ptr); }
private:
	bool _Owns;
	_Ty *_Ptr;
};


//////////////////////////////////////////////////////////////////////
//	CStringResource
//
class CStringResource
{
public:
	CStringResource(UINT nResourceID);	
	operator							CComVariant(void);
	operator							std::string(void);
	std::string							str(void) const;
	CComBSTR							bstr(void) const;
	int									CompareNoCase(const std::string& sz) const;
	int									CompareNoCase(const CComBSTR& s) const;
	int									CompareNoCaseAndSpace(const std::string& sz) const;

protected:
	UINT								m_nResourceID;	
};


//////////////////////////////////////////////////////////////////////
//	estring
//
class estring : public std::string
{
public:
	estring();
	estring(long n);
	estring(int n);
	estring(unsigned int n);
	estring(double f);
	estring(std::string::size_type n, char c);
	estring(char* sz);
	estring(const _bstr_t& b);
	estring(const BSTR& b);	
	estring(const CLSID& clsid);
	estring(const VARIANT& vIn);
	estring(const CComVariant& vIn);
	estring(const _variant_t& vIn);
	estring(const std::string& sz);
	estring(const std::stringstream& ss);
	explicit estring(const _com_error& e);
	virtual								~estring();
	
	operator CComVariant(void) const;
	operator CComBSTR(void) const;
	bool operator!=(CComVariant v) const;
	bool operator==(UINT nResourceID) const;
	estring& operator+=(double f);
	estring& operator+=(const BSTR& s);
	estring& operator+=(const std::string& sz);
	
	static OLECHAR*						AnsiToUnicode(const std::string& szAnsi);
	static OLECHAR*						AnsiToUnicode(LPCSTR lpAnsi);
	int									CompareNoCase(estring sz) const;
	static int							CompareNoCase(std::string sz1, std::string sz2);
	static int							CompareNoCase(const std::string& sz, CComBSTR& s);
	static int							CompareNoCaseAndSpace(std::string sz1, std::string sz2);
	int									CompareNoCaseAndSpace(const std::string& sz) const;
	int									CompareNoCaseAndSpace(const CComVariant& v, size_t nCharactersToCompare) const;
	int									CompareNoCase(CComBSTR& s) const;
	static bool							CompareWildcard(const std::string& sz, const std::string& szWildcard);
	bool								CompareWildcard(const std::string& szWildcard) const;
	static std::string					esubstr(const std::string& szIn, std::string szFrom, std::string szTo);
	std::string							esubstr(std::string szFrom, std::string szTo) const;
	static std::string::size_type		findnocase(const std::string* psz, estring sz, std::string::size_type pos = 0);
	std::string::size_type				findnocase(estring sz, std::string::size_type pos = 0) const;
	static CComBSTR						GetBSTR(const std::string& sz);
	CComBSTR							GetBSTR(void) const;	
	static HRESULT						GetBSTR(const std::string& sz, BSTR* ps);
	HRESULT								GetBSTR(BSTR* ps) const;
	static estring						GetNewGuid(void);
	static HRESULT						GetValue(const std::string& sz, VARIANT* pResult);
	HRESULT								GetValue(VARIANT* pResult) const;
	CComVariant							GetValue(void) const;
	static CComVariant					GetValue(const std::string& sz);	
	static bool							isbool(const std::string& sz, bool* pb);
	bool								isbool(bool* pb) const;
	static bool							isdate(const std::string& sz, long* pn);
	bool								isdate(long* pn) const;
	static bool							isdouble(const std::string& sz, double* pf);
	static bool							isdouble(const std::string& sz);
	bool								isdouble(double* pf) const;
	bool								isdouble(void) const;		
	static bool							islong(const std::string& sz, long* pn);
	static bool							islong(const std::string& sz);
	bool								islong(long* pn) const;
	bool								islong(void) const;
	static bool							isunsignedshort(const std::string& sz, unsigned short* pn);
	static void							lc(std::string* psz);
	void								lc(void);	
	static estring						left(const std::string* psz, size_type t);
	estring								left(size_type t) const;
	static estring						left(const std::string* psz, const std::string& sz);
	estring								left(const std::string& sz) const;
	static estring						mid(const std::string* psz, const std::string& szL, const std::string& szR, bool bIgnoreCase = false);
	static estring						mid(const std::string* psz, const std::string& szL);		
	estring								mid(const std::string& szL, const std::string& szR, bool bIgnoreCase = false) const;
	estring								mid(const std::string& szL) const;
	estring								mid(size_type t) const;
	HRESULT								LoadString(UINT nResourceID);		
	static void							ltrim(std::string* psz);
	void								ltrim(void);	
	static void							rtrim(std::string* psz);
	void								rtrim(void);	
	void								ReplaceWhiteSpaceWithSpace(void);
	static void							ReplaceWhiteSpaceWithSpace(std::string* psz);
	static long							ReplaceStrInStr(std::string* psz, std::string szFrom, std::string szTo, std::string::size_type nPos = 0);	
	long								ReplaceStrInStr(std::string szFrom, std::string szTo, std::string::size_type nPos = 0);	
	static estring						right(const std::string* pszIn, const std::string& sz);
	estring								right(const std::string& sz) const;
	static estring						right(const std::string* pszIn, size_type t);
	estring								right(size_type t) const;	
	static HRESULT						Set(const BSTR& b, std::string* psz);
	HRESULT								Set(const BSTR& b);
	HRESULT								Set(LPCOLESTR lpsz);
	HRESULT								Set(const std::stringstream& ss);
	void								SetError(HRESULT hr);
	static long							Split(std::string szIn, std::string szDelimit, std::vector<estring>* pVector);
	static long							Split(std::string szIn, std::string szDelimit, std::vector<std::string>* pVector);
	long								Split(std::string szDelimit, std::vector<std::string>* pVector);	
	void								Split(const std::string& szColumnDelimit, const std::string& szRowDelimit, CParameterMap* ppm);
	static estring						SplitBracket(const std::string& szIn, std::string* pszA, std::string* pszB);
	static void							StripWhiteSpace(std::string* psz);
	void								StripWhiteSpace(void);
	static CComBSTR						trim(const CComBSTR& s);
	static void							trim(std::string* psz);
	void								trim(void);	
	static void							uc(std::string* psz);
	void								uc(void);				
	static LPSTR						UnicodeToAnsi(LPOLESTR wStr);
	static estring						UnSplitBracket(const std::string& szA, const std::string& szB);

protected:
	static bool							WildcardMatch(const char* szWildcards, const char* sz);
	static bool							WildcardScan(const char*& szWildcards, const char*& sz);
	static ULONG						WStrlen(PWCHAR str);
};


//////////////////////////////////////////////////////////////////////
//	CParameterMap
//
class CParameterMap  
{
public:
	// The following members allow the use of map_[...] macros which expect a variant
	// and access the vt and pdispVal members.
	// Cynics would argue that this is a kludge! (At least I've made them const statics and 
	// set vt to VT_USERDEFINED and pdispVal to NULL.)
	static const VARTYPE				vt;
	static const IDispatch*				pdispVal;
public:	
	// Construction (we never have a constructor from a variant since this can fail due to incomplete UDF calculation)
	CParameterMap();
	explicit CParameterMap(const std::string& s);
	explicit CParameterMap(const BSTR& b);
	CParameterMap(int nRows, int nCols);
	explicit CParameterMap(const CMatrix& m);
	explicit CParameterMap(const CVector& v);
	CParameterMap(const CParameterMap& pm);
	explicit CParameterMap(double f);
	explicit CParameterMap(const std::stringstream& ss);

	// Operators
	virtual ~CParameterMap(void);
	const CParameterMap& operator=(const CParameterMap& pm);
	
public:
	// GetValue		
	CComVariant							GetValue(void) const;
	HRESULT								GetValue(double* pf) const;	
	HRESULT								GetValue(long* pn) const;
	HRESULT								GetValue(bool* pb) const;
	HRESULT								_GetValue(VARIANT* pv) const;
	HRESULT								GetValue(VARIANT* pv) const;
	HRESULT								GetValue(std::string* psz) const;
	HRESULT								GetValue(BSTR* ps) const;
	HRESULT								GetValue(CParameterMap* ppm) const;		
	HRESULT								GetValue(CVector* pvector) const;
	HRESULT								GetValue(CMatrix* pmatrix) const;
	HRESULT								GetValue(std::stringstream& ssOut) const;
	HRESULT								GetValue(std::stringstream* pssOut) const {return GetValue(*pssOut);}	
	template<class I> HRESULT			GetValue(long nRow, long nCol, CComPtr<I>& sp)
	{
		HRESULT hr;
		CComVariant v;
		if (hr = GetValue(nRow, nCol, &v)) return hr;
		if (v.vt != VT_DISPATCH) return E_FAIL;
		if (!v.pdispVal){
			// Null values are OK, but failed dynamic casts (below) are not.
			sp = NULL;
			return S_OK;
		}		
		sp = dynamic_cast<I*>(v.pdispVal);		
		return sp ? S_OK : E_FAIL;
	}
	// This specialisation solves a potential no RTTI data problem.
	template<>	HRESULT					GetValue(long nRow, long nCol, CComPtr<IDispatch>& sp)
	{
		HRESULT hr;
		CComVariant v;
		if (hr = GetValue(nRow, nCol, &v)) return hr;
		if (v.vt != VT_DISPATCH) return E_FAIL;
		sp = v.pdispVal;
		return S_OK;		
	}
	template<class T> HRESULT			GetValue(GVector<T>* pv) const
	{
		HRESULT							hr;
		std::vector<T>					vector;
		if (hr = GetValue(&vector)) return hr;
		pv->resize(vector.size());
		for (int n = 0; n < vector.size(); n++){
			(*pv)[n] = vector[n];
		}
		return S_OK;
	}
	template<class T> HRESULT			GetValue(std::vector<T>* pvector) const
	{		
		HRESULT							hr;
		if (!m_nCols || !m_nRows){
			// blank matrix
			pvector->clear();
			return S_OK;
		} else if (m_nCols == 1){
			pvector->resize(m_nRows);
			for (long nRow = 0; nRow < m_nRows; nRow++){		
				if (hr = GetValue(nRow, 0, &(*pvector)[nRow])) return hr;
			}
			return S_OK;
		} else if (m_nRows == 1){
			// This is the potential transpose case. If we cannot set pvector in this form then we
			// effectively transpose the input and proceed as in the m_nCols == 1 case.			
			pvector->resize(1);
			if (!GetValue(&(*pvector)[0])) return S_OK;
			// Return the transpose if this point is reached.			
			pvector->resize(m_nCols);
			for (long nCol = 0; nCol < m_nCols; nCol++){
				if (hr = GetValue(0, nCol, &(*pvector)[nCol])) return hr;
			}
			return S_OK;
		} else {
			// assume each row is an element of pvector of type T			
			pvector->resize(m_nRows);
			for (long nRow = 0; nRow < m_nRows; nRow++){
				CParameterMap pm;
				if (hr = GetRow(nRow, &pm)) return hr;				
				if (hr = pm.GetValue(&(*pvector)[nRow])) return hr;			
			}
			return S_OK;
		}
		return E_FAIL;		
	}
	template<class T> HRESULT			GetValue(long nCol, std::vector<T>& vector) const
	{
		HRESULT hr;		
		vector.clear();
		for (long nRow = 0; nRow < m_nRows; nRow++){
			T value;
			if (hr = GetValue(nRow, nCol, &value)) return hr;
			vector.push_back(value);
		}
		return S_OK;
	}
	// if nested templates were supported then we could write template<template<class Key, class Value> class Map<Key, Value> > HRESULT GetValue(long nColKey, long nColValue, Map<Key, Value>& map) const
	template<class Key, class Value> HRESULT GetValue(long nColKey, long nColValue, std::map<Key, Value>& map) const
	// nColValue - set to -1 if we are not interested in obtaining the value from the parameter map (useful for the 1 column parameter map case)
	{
		HRESULT hr;
		map.clear();

		for (long nRow = 0; nRow < m_nRows; nRow++){
			Key key;
			Value value;
			if (hr = GetValue(nRow, nColKey, &key)) return hr;
			if (nColValue >= 0 && (hr = GetValue(nRow, nColValue, &value))) return hr;
			map[key] = value;
		}
		return S_OK;
	}
	template<class Key, class Value> HRESULT GetValue(long nColKey, long nColValue, std::hash_map<Key, Value>& map) const
	// nColValue - set to -1 if we are not interested in obtaining the value from the parameter map (useful for the 1 column parameter map case)
	{
		HRESULT hr;
		map.clear();

		for (long nRow = 0; nRow < m_nRows; nRow++){
			Key key;
			Value value;
			if (hr = GetValue(nRow, nColKey, &key)) return hr;
			if (nColValue >= 0 && (hr = GetValue(nRow, nColValue, &value))) return hr;
			map[key] = value;
		}
		return S_OK;
	}
	// Again, if we had nested templates then the following two functions could be templatised into one.
	template<class Key, class Value1, class Value2> HRESULT GetValue(long nColKey, long nColValue1, long nColValue2, std::map<Key, Value1>& map1, std::map<Key, Value2>& map2) const
	{
		HRESULT hr;
		map1.clear();
		map2.clear();
		for (long nRow = 0; nRow < m_nRows; nRow++){
			Key key;
			if (hr = GetValue(nRow, nColKey, &key)) return hr;
			if (!IsBlank(nRow, nColValue1)){
				Value1 value;
				if (hr = GetValue(nRow, nColValue1, &value)) return hr;
				map1[key] = value;
			}
			if (!IsBlank(nRow, nColValue2)){
				Value2 value;
				if (hr = GetValue(nRow, nColValue2, &value)) return hr;
				map2[key] = value;
			}
		}
		return S_OK;
	}	
	template<class Key, class Value1, class Value2> HRESULT GetValue(long nColKey, long nColValue1, long nColValue2, std::hash_map<Key, Value1>& map1, std::hash_map<Key, Value2>& map2) const
	{
		HRESULT hr;
		map1.clear();
		map2.clear();
		for (long nRow = 0; nRow < m_nRows; nRow++){
			Key key;
			if (hr = GetValue(nRow, nColKey, &key)) return hr;
			if (!IsBlank(nRow, nColValue1)){
				Value1 value;
				if (hr = GetValue(nRow, nColValue1, &value)) return hr;
				map1[key] = value;
			}
			if (!IsBlank(nRow, nColValue2)){
				Value2 value;
				if (hr = GetValue(nRow, nColValue2, &value)) return hr;
				map2[key] = value;
			}
		}
		return S_OK;
	}		
	template<class Key, class Value> HRESULT GetValue(long nColKey, long nColValue, std::vector<Key>& vectorKey, std::vector<Value>& vectorValue) const
	{
		HRESULT hr;
		vectorKey.clear();
		vectorValue.clear();
		for (long nRow = 0; nRow < m_nRows; nRow++){
			Key key;
			Value value;
			if (hr = GetValue(nRow, nColKey, &key)) return hr;
			if (hr = GetValue(nRow, nColValue, &value)) return hr;
			vectorKey.push_back(key);
			vectorValue.push_back(value);
		}
		return S_OK;
	}

	VARIANT*							GetValue(long nRow, long nCol) const;	
	HRESULT								GetValue(long nRow, long nCol, bool* pb) const;
	HRESULT								GetValue(long nRow, long nCol, unsigned short* pn) const;
	HRESULT								GetValue(long nRow, long nCol, VARIANT* pv) const;
	HRESULT								GetValue(long nRow, long nCol, long* pn) const;	
	HRESULT								GetValue(long nRow, long nCol, double* pf) const;	
	HRESULT								GetValue(long nRow, long nCol, std::string* psz) const;
	HRESULT								GetValue(long nRow, long nCol, BSTR* ps) const;	
	HRESULT								GetValue(long nRow, long nCol, CParameterMap* ppm) const;
	HRESULT								GetValue(long nRow, long nCol, CVector* pvector) const;

public:
	// SetValue	
	HRESULT								SetValue(const CParameterMap& pm);
	HRESULT								SetValue(double f);
	HRESULT								SetValue(const BSTR& s);
	HRESULT								SetValue(const std::string& sz);
	HRESULT								SetValue(const std::stringstream& ss);
	HRESULT								SetValue(const CComVariant& v, bool bRemoveBlankMargins = true);
	HRESULT								SetValue(const CMatrix& m);
	HRESULT								SetValue(const CVector& v);
	HRESULT								SetSize(long nRows, long nCols);
	HRESULT								SetValue(const iVector& iv);
	template<class T> HRESULT			SetValue(const std::set<T>& set)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetSize(set.size(), 1)) return hr;
		for (std::set<T>::const_iterator it = set.begin(); it != set.end(); it++){
			if (hr = SetValue(nRow++, 0, *it)) return hr;
		}
		return S_OK;
	}
	
	template<class T> HRESULT			SetValue(const std::vector<T>& v)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetSize(v.size(), 1)) return hr;
		for (std::vector<T>::const_iterator it = v.begin(); it < v.end(); it++){
			if (hr = SetValue(nRow++, 0, *it)) return hr;
		}
		return S_OK;
	}

	// ToDo - if nested templates were supported then we could combine the following two functions
	template<class Key, class Value> HRESULT SetValue(const std::map<Key, Value>& map)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetSize(map.size(), 2)) return hr;
		for (std::map<Key, Value>::const_iterator it = map.begin(); it != map.end(); it++){
			if (hr = SetValue(nRow, 0, it->first)) return hr;
			if (hr = SetValue(nRow, 1, it->second)) return hr;
			nRow++;
		}
		return S_OK;
	}
	template<class Key, class Value> HRESULT SetValue(const std::multimap<Key, Value>& map)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetSize(map.size(), 2)) return hr;
		for (std::map<Key, Value>::const_iterator it = map.begin(); it != map.end(); it++){
			if (hr = SetValue(nRow, 0, it->first)) return hr;
			if (hr = SetValue(nRow, 1, it->second)) return hr;
			nRow++;
		}
		return S_OK;
	}
	// Again, if we had nested templates then the following two functions could be templatised into one.
	template<class Key, class Value1, class Value2> HRESULT SetValue(const std::map<Key, Value1>& map1, const std::map<Key, Value2>& map2)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetSize(map1.size() + map2.size(), map2.size() ? 3 : 2)) return hr;			// maximum possible size
		std::map<Key, Value1>::const_iterator it1 = map1.begin();
		std::map<Key, Value2>::const_iterator it2 = map2.begin();
		while (it1 != map1.end() || it2 != map2.end()){
			if (it2 == map2.end() || (it1 != map1.end() && it1->first < it2->first)){
				if (hr = SetValue(nRow, 0, it1->first)) return hr;
				if (hr = SetValue(nRow, 1, it1->second)) return hr;
				it1++;
				nRow++;
			} else if (it1 == map1.end() || (it2 != map2.end() && it2->first < it1->first)){
				if (hr = SetValue(nRow, 0, it2->first)) return hr;
				if (hr = SetValue(nRow, 2, it2->second)) return hr;
				it2++;
				nRow++;
			} else if (it1->first == it2->first){
				if (hr = SetValue(nRow, 0, it1->first)) return hr;
				if (hr = SetValue(nRow, 1, it1->second)) return hr;
				if (hr = SetValue(nRow, 2, it2->second)) return hr;
				it1++;
				it2++;
				nRow++;
			} else {
				ATLASSERT(false);
			}
		}		
		return SetRows(nRow);
	}
	template<class Key, class Value1, class Value2> HRESULT SetValue(const std::hash_map<Key, Value1>& map1, const std::hash_map<Key, Value2>& map2)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetValue(map1.size() + map2.size(), map2.size() ? 3 : 2)) return hr;			// maximum possible size
		std::map<Key, Value1>::const_iterator it1 = map1.begin();
		std::map<Key, Value2>::const_iterator it2 = map2.begin();
		while (it1 != map1.end() || it2 != map2.end()){
			if (it2 == map2.end() || (it1 != map1.end() && it1->first < it2->first)){
				if (hr = SetValue(nRow, 0, it1->first)) return hr;
				if (hr = SetValue(nRow, 1, it1->second)) return hr;
				it1++;
				nRow++;
			} else if (it1 == map1.end() || (it2 != map2.end() && it2->first < it1->first)){
				if (hr = SetValue(nRow, 0, it2->first)) return hr;
				if (hr = SetValue(nRow, 2, it2->second)) return hr;
				it2++;
				nRow++;
			} else if (it1->first == it2->first){
				if (hr = SetValue(nRow, 0, it1->first)) return hr;
				if (hr = SetValue(nRow, 1, it1->second)) return hr;
				if (hr = SetValue(nRow, 2, it2->second)) return hr;
				it1++;
				it2++;
				nRow++;
			} else {
				ATLASSERT(false);
			}
		}		
		return SetRows(nRow);
	}
	HRESULT								SetValue(long nRow, long nCol, const CComVariant& v);
	HRESULT								SetValue(long nRow, long nCol, const std::stringstream& ss);
	HRESULT								SetValue(long nRow, long nCol, const std::string& sz);
	HRESULT								SetValue(long nRow, long nCol, const CComBSTR& s);
	HRESULT								SetValue(long nRow, long nCol, double f);
	HRESULT								SetValue(long nRow, long nCol, const CStringResource& sr);	
	template <class I> HRESULT			SetValue(long nRow, long nCol, const CAdapt<CComQIPtr<I> >& sp)
	{
		CComVariant						v(sp.m_T);
		return SetValue(nRow, nCol, v);
	}
	template <class I> HRESULT			SetValue(long nRow, long nCol, CComPtr<I> sp)
	{		
		CComVariant						v(sp);
		return SetValue(nRow, nCol, v);
	}
	HRESULT								SetValue(long nRow, long nCol, const CVector& v);	
	HRESULT								SetValues(long nRows, long nCols, /*CComVariant*/...);

public:		
	enum find_cell_search{
		TopLeftByColumns = 1,
		TopLeftByRows = 2,
		LeftColumnOnly = 3,
		TopRowOnly = 4
	};
	
	template<class I, class MlEq> static long AnalyticToCom(const RCPtr<MlEq> h, CComPtr<I> sp)
	{			
		void*							p = &*h;
		long							n = reinterpret_cast<long>(p);		
		return sp->put_AnalyticObject(n) ? 0 : n;		// i.e. return zero if failed
	}

	template<class I, class MlEq> static long ComToAnalytic(CComPtr<I> sp, RCPtr<MlEq>& h)
	{
		long							n = 0;
		if (!sp){
			h = NULL;
			return 0;
		}
		if (!sp->get_AnalyticObject(&n)) h = reinterpret_cast<MlEq*>(n);
		return n;										// i.e. returns zero if failed
	}
	template<class I, class MlEq> static long ComToAnalytic(CComQIPtr<I> sp, RCPtr<MlEq>& h)
	{
		long							n = 0;
		if (!sp){
			h = NULL;
			return 0;
		}
		if (!sp->get_AnalyticObject(&n)) h = reinterpret_cast<MlEq*>(n);
		return n;										// i.e. returns zero if failed
	}
	//
	//	provides the same function as an error handler in VBA
	//
	template <class ISingular> static HRESULT ErrorHandler(CComPtr<ISingular> spObject, HRESULT hrError, std::string* psz)
	//  hrError - returns the text associated with this error if no richer information is available
	{
		CComVariant							v;
		ErrorHandler(spObject, hrError, &v);
		psz->assign(estring(v));
		return S_OK;
	}
	template <class ISingular> static HRESULT ErrorHandler(CComPtr<ISingular> spObject, HRESULT hrError, VARIANT* pResult)
	//  hrError - returns the text associated with this error if no richer information is available
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
				if (!spErrorInfo->GetDescription(&sError)){
					CComVariant(sError).Detach(pResult);
					return S_OK;
				}
			}
		}

		if (hrError != S_OK){		
			estring szError;		
			szError.SetError(hrError);
			if (szError.size()) return szError.GetValue(pResult);
		}
		
		return CComVariant("An error occurred. Sirius cannot return any more details.").Detach(pResult);
	}
	
	// ToDo - Join these GetObjectName functions once nested templates are supported - template <class I, template <class> class T>
	template <class I>
	static HRESULT GetObjectName(CComQIPtr<I> spObject, std::string* pszOut)
	{
		CComPtr<ITypeInfo>					pti;
		HRESULT								hr;
		CComBSTR							sName;
		
		if (!spObject) return E_FAIL;
		if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;		
		if (hr = pti->GetDocumentation(-1, &sName, NULL, NULL, NULL)) return hr;
			
		estring								sz(sName);
		
		if (!sz.size()) return E_FAIL;
		if (sz.c_str()[0] == 'I'){		
			*pszOut = sz.substr(1);		// remove the prefixed 'I'
		} else {
			pszOut->assign(sz);
		}
		return S_OK;	
	}
	
	template <class I>
	static HRESULT GetObjectName(CComPtr<I> spObject, std::string* pszOut)
	{
		CComPtr<ITypeInfo>					pti;
		HRESULT								hr;
		CComBSTR							sName;
		
		if (!spObject) return E_FAIL;
		if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;		
		if (hr = pti->GetDocumentation(-1, &sName, NULL, NULL, NULL)) return hr;
			
		estring								sz(sName);

		if (!sz.size()) return E_FAIL;
		if (sz.c_str()[0] == 'I'){		
			*pszOut = sz.substr(1);		// remove the prefixed 'I'
		} else {
			pszOut->assign(sz);
		}
		return S_OK;	
	}
	
	template <class T> static GVector<T> StdVectorToGVector(const std::vector<T>& vIn)
	{
		GVector<T>	vOut(vIn.size());
		for (int nElement = 0; nElement < vIn.size(); nElement++){
			vOut[nElement] = vIn[nElement];
		}
		return vOut;
	}

	template <class Key, class Value>	static void MergeMap(std::map<Key, Value>* pmap1, const std::map<Key, Value>& map2)
	{
		for (std::map<Key, Value>::const_iterator it2 = map2.begin(); it2 != map2.end(); it2++){
			std::map<Key, Value>::const_iterator it1 = pmap1->find(it2->first);
			if (it1 == pmap1->end()){
				(*pmap1)[it2->first] = it2->second;
			}
		}
	}
	
	static const CParameterMap&			AddToEnd(const CParameterMap& pmTop, const CParameterMap& pmBottom);		
	static HRESULT						ArrayToVector(const VARIANT& vIn, std::vector<CComVariant>* pvvOut, std::vector<CParameterMap>* pvpmOut, bool bDecompose = true);
	static HRESULT						CLSIDFromProgID(const CComBSTR& s, CLSID& clsid);
	static HRESULT						CLSIDToObjectName(const CLSID& CLSID, std::string* psz);
	static void							CreateBlankVariant(long nRows, long nCols, VARIANT* pv);
	static HRESULT						DateToString(DATE date, std::string* psz);
	static void							DisplayError(UINT nType);
	static void							DisplayError(const _bstr_t& s, UINT nType);
	static void							DisplayError(const std::string& szMessage, UINT nType);
	static void							DisplayError(UINT nResourceID, UINT nType);	
	static unsigned long				ExcelDateToJulianDate(long nExcelDate);
	static HRESULT						ExcelGetFormula(CComDispatchDriverEx& ddRange, long nRow, long nCol, std::string* pszOut);
	static HRESULT						ExcelRangeToRowCol(CComDispatchDriverEx& ddRange, long* pnTopRow, long* pnLeftCol, long* pnBottomRow, long* pnRightCol);
	static short int					GetArrayDimensions(const CComVariant& v, long** ppsize = NULL);
	static HRESULT						GetObjectIID(CComDispatchDriverEx& dd, IID* piidOut);
	static HRESULT						GetObjectIID(CComPtr<IDispatch> sp, IID* piidOut);
	static HRESULT						GetObjectName(CComDispatchDriverEx& dd, std::string* pszOut);	
	static std::string					GetUser(void);		
	static VARTYPE						GetUnionVariantType(VARTYPE vt1, VARTYPE vt2);
	static long							GetVariantSize(const CComVariant& v);
	static bool							IsExcelFunctionWizardVisible(void);
	static bool							IsVariantTypeNumeric(VARTYPE vt);
	static long							JulianDateToExcelDate(unsigned long nJulianDate);
	static HRESULT						ProgIDFromCLSID(const CLSID& clsid, CComBSTR& s);	
	static HRESULT						RemoveCalculationNumberFromHandle(const CLSID& clsid, std::string* psz, bool* pbIsValidSpreadsheetHandle = NULL);	
	static HRESULT						RemoveCalculationNumberFromHandle(const CLSID& clsid, CComBSTR& ps, bool* pbIsValidSpreadsheetHandle = NULL);
	static HRESULT						ReturnErrorS(const std::string& sz, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorS(const std::stringstream& ss, const IID& iid =_IID_ISiriusApplication);
	static HRESULT						ReturnErrorS(LPCTSTR lpsz, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorS(const CComBSTR& s, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorSRS(const std::string& sz1, UINT nResourceID, const std::string& sz2, const IID& iid =_IID_ISiriusApplication);
	static HRESULT						ReturnErrorSS(std::string sz1, const std::string& sz2, const IID& iid =_IID_ISiriusApplication);
	static HRESULT						ReturnErrorR(UINT nResourceID, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorRN(UINT nResourceID, long n, const IID& iid =_IID_ISiriusApplication);
	static HRESULT						ReturnErrorRP(UINT nResourceID, const CParameterMap& pmIn, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorRR(UINT nResourceID1, UINT nResourceID2, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorRS(UINT nResourceID, const std::string& sz, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorRS(UINT nResourceID, const CComBSTR& sIn, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorRSR(UINT nResourceID1, const std::string& sz, UINT nResourceID2, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorRV(UINT nResourceID, const CComVariant& vIn, const IID& iid = _IID_ISiriusApplication);
	static HRESULT						ReturnErrorSV(const std::string& sz, const CComVariant& vIn, const IID& iid =_IID_ISiriusApplication);
	static int							Round(double f);
	static double						Round(double f, unsigned short int nPrecision);
	static HRESULT						SafeArrayElementToVariant(SAFEARRAY* psa, long* nIndex, VARIANT* pvRet);		
	static HRESULT						StringToDate(const std::string& sz, DATE* pdate);
	static HRESULT						StringToDate(const std::string& sz, long* pnDate);
	static void							ThrowComErrorR(UINT nResourceID, const IID& iid = _IID_ISiriusApplication);
	static void							ThrowComErrorRR(UINT nResourceID1, UINT nResourceID2, const IID& iid = _IID_ISiriusApplication);
	static void							ThrowComErrorRS(UINT nResourceID, const std::string& sz, const IID& iid = _IID_ISiriusApplication);
	static void							ThrowComErrorRV(UINT nResourceID, const CComVariant& vIn, const IID& iid = _IID_ISiriusApplication);
	static void							ThrowComErrorS(const CComBSTR& s, const IID& iid =_IID_ISiriusApplication);
	static void							ThrowStringError(CComPtr<IDispatch> spObject);
	static HRESULT						VariableArgumentListToArray(CComVariant* pvOut, long nArgs, va_list vl);
	static HRESULT						VariableArgumentListToArray(CComVariant* pvOut, long nArgs, /*VARIANT*/...);
	static HRESULT						VariableArgumentListToArray(VARIANT* pvOut, long nArgs, /*CParameterMap*/...);
	static HRESULT						VectorToArray(std::vector<CComVariant>& vv, VARIANT* pvOut);
	static HRESULT						VectorToArray(std::vector<CParameterMap*>& vppm, CComVariant* pvOut);
	static HRESULT						VectorToArray(std::vector<CParameterMap>& vpm, VARIANT* pvOut);
	static HRESULT						WindowToClassName(HWND hWnd, std::string* pszOut);
							
	const CParameterMap&				AddToEnd(const CParameterMap& pm);
	const CParameterMap&				AddToEnd(const std::string& sz);
	const CParameterMap&				AddToEnd(const std::stringstream& ss);
	const CParameterMap&				AddToRHS(const CParameterMap& pm);	
	HRESULT								Attach(CParameterMap* pSrc);
	void								Clear(void);
	HRESULT								CopyToClipboard(void) const;
	HRESULT								Extract(UINT nResourceIDHeader, bool bRemoveHeaderColumn, CParameterMap* ppmExtracted, CParameterMap* ppmRemaining = NULL) const;
	HRESULT								FindCell(std::string szMatch, bool bIgnoreCaseAndSpaces, find_cell_search fcs, long* pnRow, long* pnCol) const;	
	bool								GetAutoGrow(void) const;
	long								GetCols(void) const;
	HRESULT								GetColumn(long nCol, CParameterMap* ppm, bool bRemoveBlanks = true) const;
	void								GetColumn(long nCol, VARTYPE vt, VARIANT* pVal) const;
	const CComVariant*					GetConstElementPtr(long nRow, long nCol) const;
	void								GetDistinct(CParameterMap* ppm) const;
	CComVariant*						GetElementPtr(long nRow, long nCol);
	HRESULT								GetElementRowCol(const CComVariant* pvElement, long* pnRow, long* pnCol) const;
	long								GetNonBlankRows(void) const;
	HRESULT								GetRow(long nRow, CParameterMap* ppm) const;
	long								GetRows(void) const;	
	HRESULT								GetString(std::string* psz) const;
	HRESULT								GetString(BSTR* ps) const;
	CComBSTR							GetString(void) const;
	HRESULT								GetObject(CComVariant* pvObject) const;	
	HRESULT								GetValueFromColumnHeader(const std::string& szHeader, double* pf) const;	
	HRESULT								GetValueFromRowHeader(UINT nResourceIDHeader, bool bAllowBlank, std::string* pszData, long* pnRowUsed = NULL) const;
	HRESULT								GetValueFromRowHeader(const std::string& szHeader, bool bAllowBlank, std::string* pszData, long* pnRowUsed = NULL) const;
	HRESULT								GetValueFromRowHeader(UINT nResourceIDHeader, double* pf, long* pnRowUsed = NULL) const;
	HRESULT								GetValueFromRowHeader(const std::string& szHeader, double* pf, long* pnRowUsed = NULL) const;
	HRESULT								GetValueFromRowHeader(UINT nResourceIDHeader, CComVariant* pv, long* pnRowUsed = NULL) const;
	HRESULT								GetValueFromRowHeader(const std::string& szHeader, CComVariant* pv, long* pnRowUsed = NULL) const;	
	HRESULT								InsertAt(long nStartRow, long nStartCol, const CParameterMap& pm);
	HRESULT								InsertColumn(long nAt);
	HRESULT								InsertRow(long nAt);
	bool								IsBlank(long nRow = -1, long nCol = -1) const;
	bool								IsBlankOrDouble(long nRow/* = -1*/, long nCol/* = -1*/, UINT nResourceID/* = 0*/) const;
	bool								IsColumnVector(void) const;
	bool								IsDate(long nRow = -1, long nCol = -1, UINT nResourceID = 0) const;
	bool								IsDouble(long nRow = -1, long nCol = -1, UINT nResourceID = 0) const;
	bool								IsLong(long nRow = -1, long nCol = -1, UINT nResourceID = 0) const;
	bool								IsObject(void) const;
	bool								IsRowVector(void) const;
	bool								IsScalar(void) const;
	bool								IsVector(void) const;	
	HRESULT								Max(double *pVal) const;
	HRESULT								Min(double *pVal) const;	
	HRESULT								RemoveBlankRows(void);
	void								RemoveBlanksAtEnd(void);
	HRESULT								RemoveColumn(long nCol);
	HRESULT								RemoveRow(long nRow);
	void								ReplaceBlanksWithSpaces(void);
	HRESULT								Resize(long nRows, long nColumns, CComVariant vFill = CComVariant());
	void								SetAutoGrow(bool b);
	HRESULT								SetBlank(long nRow, long nCol);
	HRESULT								SetCellBlankIfValue(long nRow, long nCol, std::string szMatch, bool bIgnoreCaseAndSpaces);
	HRESULT								SetColumns(long nCols);	
	HRESULT								SetRows(long nRows);
	template<class I> HRESULT			SetToObjectPointer(CComPtr<I> spObject)
	{
		HRESULT							hr;
		if (hr = SetSize(1, 1)) return hr;
		return CComVariant(spObject.p).Detach(m_av);
	}		
	HRESULT								SetValues(VARTYPE vt, const std::string& szDelimit, const std::string& szValues);
	HRESULT								Split(const std::string& szDelimit, VARTYPE vt, int nParameterMaps, /*CParameterMap**/...);
	HRESULT								Substitute(const std::string& szOld, const std::string& szNew, bool bPartString, bool bIgnoreCaseAndSpaces, long nRow, long nCol);
	HRESULT								Transpose(void);

protected:	
	template<class T> HRESULT			SetValueScalar(const std::vector<T>& vector)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetValue(vector.size(), 1)) return hr;
		for (std::vector<T>::const_iterator it = vector.begin(); it < vector.end(); it++){
			if (hr = SetValue(nRow++, 0, *it)) return hr;
		}
		return S_OK;
	}
	template<class T> HRESULT			SetValueVector(const std::vector<T>& vector)
	{
		HRESULT			hr;
		long			nRow = 0;
		if (hr = SetValue(vector.size(), 1)) return hr;
		for (std::vector<T>::const_iterator it = vector.begin(); it < vector.end(); it++){
			if (hr = SetValue(nRow++, 0, *it)) return hr;
		}
		return S_OK;
	}
	
	
	static HRESULT						ArrayToVectorDecompose(std::vector<CParameterMap>* pvpm);
	static HRESULT						ArrayToVectorDecompose(std::vector<CComVariant>* pvv);
					
	HRESULT								RemoveAt(long nStart, long nElements);
	void								RemoveBlankMargins(bool bLeaveStart);
	void								RemoveBlanksAtLHS(void);
	void								RemoveBlanksAtRHS(void);
	void								RemoveBlanksAtStart(void);	
	template <class T> static int		sgn(T t) {return t < (T)0 ? -1 : +1;}

	long								m_nRows;						// number of rows in the matrix
	long								m_nCols;						// number of columns in the matrix		
	VARIANT*							m_av;							// representation of the matrix elements
	bool								m_bAutoGrow;					// true if we grow the parameter map to an appropriate size if we attempt to insert an element out of bounds
	static const VARTYPE				s_avVariantMap[13][13];			// used by GetUnionVariantType()
	static const long					s_avVariantSizes[8];			// used by GetVariantSize()
	static const bool					s_avIsVariantTypeNumeric[13];	// used by IsVariantTypeNumeric()
};


//////////////////////////////////////////////////////////////////////
//	CEnumMap
//
class CEnumMap
{
typedef std::map<std::string, CAdapt<auto_ptr_ex<CEnumMap> > >	enum_map;

public:				
	// returns an enumerator value form of a string
	template<class T, class V>
	static HRESULT GetEnum(const std::string& szEnumeratorName, GUID libid, const V& v, T* enumerator)
	{
		HRESULT							hr;
		enum_map::const_iterator		it;
		CParameterMap					pm;

		if (hr = Find(szEnumeratorName, libid, it)) return hr;
		if (hr = pm.SetValue(v)) return hr;
		return it->second.m_T->ParameterMapToEnum(pm, enumerator);
	}
	template<class T, class V>
	static HRESULT GetEnum(CComPtr<ITypeInfo> spti, GUID libid, const V& v, T* enumerator)
	{
		HRESULT							hr;
		CComBSTR						sEnumeratorName;

		if (hr = spti->GetDocumentation(-1, &sEnumeratorName, 0,0,0)) return hr;
		return GetEnum(estring(sEnumeratorName), libid, v, enumerator);
	}
	// this overload supports a default value
	template<class T, class V>
	static HRESULT GetEnum(const std::string& szEnumeratorName, GUID libid, const V& v, T defaultvalue, T* enumerator)
	{
		HRESULT							hr;
		enum_map::const_iterator		it;
		CParameterMap					pm;

		if (hr = Find(szEnumeratorName, libid, it)) return hr;
		if (hr = pm.SetValue(v)) return hr;
		if (pm.IsBlank()){
			*enumerator = defaultvalue;
			return S_OK;
		}
		return it->second.m_T->ParameterMapToEnum(pm, enumerator);
	}
	// this overload supports a default value
	template<class T, class V>
	static T GetEnum(const std::string& szEnumeratorName, GUID libid, const V& v, const T& defaultvalue)
	{
		T								t = defaultvalue;
		GetEnum(szEnumeratorName, libid, v, &t);
		return t;
	}

	// returns a comma-separated list of enumerator values
	static HRESULT GetEnumList(const std::string& szEnumeratorName, GUID libid, std::string* psz)
	{
		HRESULT							hr;
		enum_map::const_iterator		it;
		if (hr = Find(szEnumeratorName, libid, it)) return hr;
		psz->assign(it->second.m_T->m_ss.str());
		return S_OK;
	}
	static HRESULT GetEnumList(CComPtr<ITypeInfo> spti, GUID libid, std::string* psz)
	{
		HRESULT							hr;
		CComBSTR						sEnumeratorName;		

		if (hr = spti->GetDocumentation(-1, &sEnumeratorName, 0,0,0)) return hr;
		return GetEnumList(estring(sEnumeratorName.m_str), libid, psz);
	}

	// returns a string vector of enumerator strings
	static HRESULT GetEnumList(const std::string& szEnumeratorName, GUID libid, std::vector<std::string>* pv)
	{
		HRESULT							hr;
		enum_map::const_iterator		itEnumMap;
		if (hr = Find(szEnumeratorName, libid, itEnumMap)) return hr;
		
		pv->clear();
		for (std::map<std::string, long>::const_iterator it = itEnumMap->second.m_T->m_mapStringEnum.begin(); it != itEnumMap->second.m_T->m_mapStringEnum.end(); it++){
			pv->push_back(it->first);
		}
		return S_OK;
	}
				
	// returns the string form of an enumerator value
	template<class T>
	static HRESULT GetString(const std::string& szEnumeratorName, GUID libid, T enumerator, std::string* psz)
	{
		HRESULT							hr;
		enum_map::const_iterator		it;

		if (hr = Find(szEnumeratorName, libid, it)) return hr;
		std::map<long, std::string>::const_iterator itEnumString = it->second.m_T->m_mapEnumString.find(enumerator);
		if (itEnumString == it->second.m_T->m_mapEnumString.end()) return E_FAIL;
		psz->assign(itEnumString->second);
		return S_OK;		
	}
	template<class T>
	static HRESULT GetString(CComPtr<ITypeInfo> spti, GUID libid, T enumerator, std::string* psz)
	{							
		HRESULT							hr;
		CComBSTR						sEnumeratorName;		

		if (hr = spti->GetDocumentation(-1, &sEnumeratorName, 0,0,0)) return hr;
		return GetString(estring(sEnumeratorName), libid, enumerator, psz);
	}
	template<class T>
	static estring GetString(const std::string& szEnumeratorName, GUID libid, T enumerator)
	{
		std::string						sz;
		GetString(szEnumeratorName, libid, enumerator, &sz);
		return sz;
	}		
	
	virtual								~CEnumMap(void)
	{
		ATLTRACE("Calling CEnumMap destructor\n");
	}
	
protected:
	CEnumMap(void)
	{
		ATLTRACE("Calling CEnumMap constructor\n");
	}

	static HRESULT Find(const std::string& szEnumeratorName, GUID libid, enum_map::const_iterator& it)
	{
		HRESULT							hr;		
		if ((it = s_map.find(szEnumeratorName)) == s_map.end()){
			// set the contribution to s_map				
			ATLTRACE("Setting up the %s enumerator\n", szEnumeratorName.c_str());
			CComPtr<ITypeLib>			sptl;
			CComPtr<ITypeInfo>			spti;
			OLECHAR*					sName;
			unsigned short				nFind = 1;
			MEMBERID					memid = -1;	

			if (hr = ::LoadRegTypeLib(libid, 1.0, 0, LOCALE_USER_DEFAULT, &sptl)) return hr;
			sName = estring::AnsiToUnicode(szEnumeratorName);
			HRESULT hrHash = LHashValOfNameSys(SYS_WIN32, LOCALE_USER_DEFAULT, sName);
			sptl->FindName(sName, hrHash, &spti, &memid, &nFind);
			delete sName;
			if (!spti) return CParameterMap::ReturnErrorS("Cannot find enumerator '" + szEnumeratorName + "'");
			auto_ptr_ex<CEnumMap> ap(new CEnumMap);									
			hr = ap->Set(spti);			
			s_map[szEnumeratorName] = ap;			// remember that this transfers ownership to the auto_ptr_ex
			it = s_map.find(szEnumeratorName);
			if (it == s_map.end()){
				ATLASSERT(false);
				return E_FAIL;
			}					
		}
		return S_OK;
	}

	template<class T>
	HRESULT ParameterMapToEnum(const CParameterMap& pm, T* penumerator) const
	{
		HRESULT											hr;
		std::map<std::string, long>::const_iterator		it;		
		estring											sz;
		long											n;

		if (pm.IsLong()){
			if (hr = pm.GetValue(&n)) return hr;
		} else {					
			if (hr = pm.GetString(&sz)) return hr;				
			sz.lc();
			sz.StripWhiteSpace();
			// Specific conversions - ToDo - make a little more elegant!
			if (sz == "30/360"){
				sz = "thirty360";
			} else if (sz == "mid"){
				sz = "middle";
			}
			sz.ReplaceStrInStr("/", "");
			if ((it = m_mapLCaseStringEnum.find(sz)) == m_mapLCaseStringEnum.end()) return E_FAIL;
			n = it->second;
		}

		if (m_mapEnumString.find(n) == m_mapEnumString.end()) return E_FAIL;		// the long cannot be converted to an enumerator type
		*penumerator = (T)n;
		return S_OK;		
	}

	HRESULT Set(CComPtr<ITypeInfo> spTypeInfo)
	{					
		TYPEATTR*						ta = NULL;
		HRESULT							hr;

		// clear all return structures
		m_ss.str("");
		m_mapStringEnum.clear();
		m_mapLCaseStringEnum.clear();
		m_mapEnumString.clear();
		
		if (hr = spTypeInfo->GetTypeAttr(&ta)) return hr;
		if (ta->typekind == TKIND_ENUM){
			for (UINT n = ta->cVars ; n > 0; n--){
				// We count backwards so that we support many enumerators with the
				// same value. The line m_mapEnumString[v.lVal] = sz will either
				// insert or replace. The default enumerator name for each similar
				// enumerator value is therefore the first one reached in enums.h.
				VARDESC*	vardesc;
				if (!spTypeInfo->GetVarDesc(n - 1, &vardesc)){
					CComBSTR		s;
					CComVariant		v;																					
					spTypeInfo->GetDocumentation(vardesc->memid, &s, 0, 0, 0);
					if (v.ChangeType(VT_I4, vardesc->lpvarValue)) ATLASSERT(false);
					estring			sz(s.m_str);				
					m_mapEnumString[v.lVal] = sz;
					m_mapStringEnum[sz] = v.lVal;
					if (n != ta->cVars) m_ss << ',';
					m_ss << sz;
					sz.lc();
					sz.StripWhiteSpace();
					m_mapLCaseStringEnum[sz] = v.lVal;					
					spTypeInfo->ReleaseVarDesc(vardesc);
				}
			}					
		} else {
			ATLASSERT(false);
			hr = E_NOTIMPL;
		}
		spTypeInfo->ReleaseTypeAttr(ta);
		return hr;
	}
		
	static enum_map						s_map;							// holds all the discovered enumerator maps
	std::map<long, std::string>			m_mapEnumString;				// map of enumerator values to strings
	std::map<std::string, long>			m_mapLCaseStringEnum;			// map of (case and space insensitive) string values to enumerator values
	std::map<std::string, long>			m_mapStringEnum;				// map of (case-sensitive) string values to enumerator values
	std::stringstream					m_ss;							// list of comma-separated values
};


#endif