//	siriusinterfacebase.h:	Include this in all Sirius dll's.
//
//	Author:					David Cuin
//
////////////////////////////////////////////////////////////////////////

#ifndef _SIRIUSINTERFACEBASE_H
#define _SIRIUSINTERFACEBASE_H

#define using_sirius_application											// Must be above #include <siriusbase.h>
#include <siriusbase.h>
#define g_pApplication					_Module._GetSiriusApplication()		// For just-in-time singleton object creation implementation.


/////////////////////////////////////////////////////////////////////////////
//	declare_member_...
//
//	These macros associate a COM level property with Get and Put functions
//	of an underlying analytic object.
//
//	Use this to associate an analytic-level member variable with the
//  COM interface. Note that the specialised macros should be used
//  for strings, Booleans and objects.
//
#define declare_member_variable(type, property)		\
	STDMETHOD(get_##property)(type* pVal)			\
	{												\
		*pVal = m_h->Get##property();				\
		return S_OK;								\
	}												\
	STDMETHOD(put_##property)(type newVal)			\
	{												\
		m_h->Put##property(newVal);					\
		return S_OK;								\
	}

//	Special case of declare_member_variable for string properties.
#define declare_member_string(property)						\
	STDMETHOD(get_##property)(BSTR* pVal)					\
	{														\
		return estring::GetBSTR(m_h->Get##property(), pVal);\
	}														\
	STDMETHOD(put_##property)(BSTR newVal)					\
	{														\
		estring sz(newVal);									\
		sz.trim(&sz);										\
		m_h->Put##property(sz);								\
		return S_OK;										\
	}

//	Special case of declare_member_variable for Boolean properties.
#define declare_member_bool(property)						\
	STDMETHOD(get_##property)(YesNoEnum* pVal)				\
	{														\
		*pVal = m_h->Get##property() ? Yes : No;			\
		return S_OK;										\
	}														\
	STDMETHOD(put_##property)(YesNoEnum newVal)				\
	{														\
		m_h->Put##property(newVal == Yes ? true : false);	\
		return S_OK;										\
	}

//	Special case of declare_member_variable for Object properties.
//  Note that the Type is IType on the interface side and MlEqType
//  on the analytics side. ToDo - remove the COM defined variable
//  and use the map_com_to_analytic and map_analytic_to_com
//  macros instead.
#define declare_member_object(type, property)		\
	STDMETHOD(get_##property)(I##type** pVal)		\
	{												\
		return __sp##property.CopyTo(pVal);			\
	}												\
	STDMETHOD(put_##property)(I##type* newVal)		\
	{												\
		__sp##property = newVal;					\
		MlEq##type##Handle h;						\
		if (newVal){								\
			long n;									\
			__sp##property->get_AnalyticObject(&n);	\
			h = reinterpret_cast<MlEq##type##*>(n);	\
		}											\
		m_h->Put##property(h);						\
		return S_OK;								\
	}												\
	CComPtr<I##type>				__sp##property;

//	Allows the association of a sirius collection (e.g. IAssets) with an
//	analytic-level std::vector of singular smart pointer handles.
#define declare_member_collection(Collection, Singular, Member, Property)			\
public:																				\
	STDMETHOD(get_##Property)(I##Collection** pVal)									\
	{																				\
		return _sp##Member.CopyTo(pVal);											\
	}																				\
	STDMETHOD(put_##Property)(I##Collection* newVal)								\
	{																				\
		_sp##Member = newVal;														\
		if (newVal){																\
			long			nCount = 0;												\
			HRESULT			hr;														\
			if (hr = _sp##Member->get_Count(&nCount)) return hr;					\
			for (long n = 1; n <= nCount; n++){										\
				CComPtr<I##Singular> sp;											\
				if (hr = _sp##Member->get_Item(CComVariant(n), &sp)) return hr;		\
				map_com_to_analytic(sp, Singular, h);								\
				Member.push_back(h);												\
			}																		\
		} else {																	\
			Member.clear();															\
		}																			\
		return S_OK;																\
	}																				\
protected:																			\
	CComPtr<I##Collection>					_sp##Member;							\
	std::vector<MlEq##Singular##Handle>		Member;


/////////////////////////////////////////////////////////////////////////////
//	begin_calculate / end calculate
//
//  These are used in product type Evaluate() functions to detect which values
//  should be calculated.
//
#define begin_calculate(ResultPtr, Parameter)						\
{																	\
	HRESULT _hr;													\
	VARIANT_BOOL _b;												\
	if (_hr = (ResultPtr)->Calculate(Parameter, &_b)) return _hr;	\
	if (_b){														\

#define end_calculate												\
		}															\
	}


/////////////////////////////////////////////////////////////////////////////
//	map_..._com_to_analytic / map_analytic_to_com / associate_analytic_object
//
//  These macros are used to get the analytic smart pointer associated with
//  a com smart pointer. They only work if the com class, say, IType has
//  a member variable m_h of type MlEqTypeHandle. Furthermore, you need to 
//  include associate_analytic_object in the interface class.
//
#define map_com_to_analytic(Com, Type, Analytic)					\
	MlEq##Type##Handle	Analytic;									\
	CParameterMap::ComToAnalytic(Com, Analytic);

#define map_bare_com_to_analytic(Com, Type, Analytic)				\
	CComPtr<I##Type> _sp_##Com = Com;								\
	map_com_to_analytic(_sp_##Com, Type, Analytic);

#define map_analytic_to_com(Analytic, Type, Com)					\
	CComPtr<I##Type> Com;											\
	Com.CoCreateInstance(CLSID_##Type);								\
	CParameterMap::AnalyticToCom(Analytic, Com);

#define associate_analytic_object(CObject, MlEqObject)							\
	MlEqObject##Handle				m_h;										\
	STDMETHOD(get_AnalyticObject)(/*[out, retval]*/ long* pVal)					\
	{																			\
		CObject* pCObject = dynamic_cast<CObject*>(this);						\
		if (!pCObject) return E_POINTER;										\
		void* pAnalyticObject = &*pCObject->m_h;								\
		if (!pAnalyticObject) return E_POINTER;									\
		*pVal = reinterpret_cast<long>(pAnalyticObject);						\
		return S_OK;															\
	}																			\
	STDMETHOD(put_AnalyticObject)(/*[in]*/ long newVal)							\
	{																			\
		MlEqObject* pAnalyticObject = reinterpret_cast<MlEqObject*>(newVal);	\
		m_h = pAnalyticObject;													\
		return S_OK;															\
	}


/////////////////////////////////////////////////////////////////////////////
//	declare_serialisable / implement_serialisable
//
//	Use these to implement the Xml read property (which is useful for object
//  serialisation).
//
#define declare_serialisable										\
	STDMETHOD(get_Xml)(/*[out, retval]*/ BSTR* pVal);

#define implement_serialisable(Singular)							\
	STDMETHODIMP C##Singular::get_Xml(BSTR* pVal)					\
	{																\
		HRESULT hr;													\
		xmlstreamer			ssXML;									\
		I##Singular*		p = dynamic_cast<I##Singular*>(this);	\
		hr = CXmlStreamer::GetXML(CComVariant(p), ssXML);			\
		return estring(ssXML).GetBSTR(pVal);						\
	}


/////////////////////////////////////////////////////////////////////////////
//	map_[optional_]object_parameter
//
//  These macros are used to map an input object handle (or key) to an
//  analytic smart pointer. Usage restrictions and interface requirements
//  are as for map_com_to_analytic and map_analytic_to_com.
//
//	The _optional_ macro returns a null pointer in MappedName if the input
//  name is null. The non-optional version throws a string exception in this
//  case.
//
#define map_optional_object_parameter(InputName, Type, MappedName)												\
	CComQIPtr<I##Type>					_s##MappedName;															\
	MlEq##Type##Handle					MappedName = NULL;														\
	{																											\
		if ((InputName).vt == VT_DISPATCH){																		\
			if ((InputName).pdispVal){																			\
				try {																							\
					_s##MappedName = dynamic_cast<I##Type*>(const_cast<IDispatch*>((InputName).pdispVal));		\
				} catch (...){																					\
				}																								\
			}																									\
		}																										\
		if (!_s##MappedName){																					\
			CParameterMap					pmInput;															\
			pmInput.SetValue(InputName);																		\
			if (pmInput.IsObject()){																			\
				_s##MappedName = dynamic_cast<I##Type*>(pmInput.GetElementPtr(0, 0)->pdispVal);					\
			}																									\
		}																										\
		if (!_s##MappedName){																					\
			CComPtr<IDispatch>				spDispatch;															\
			estring							sz;																	\
			try {																								\
				map_optional_parameter((InputName), estring, MappedName, "");									\
				sz.assign(MappedName);																			\
			} catch (const std::string& szError){																\
				/*the input name could be the data for an object*/												\
				map_optional_parameter((InputName), CComVariant, MappedName, CComVariant());					\
				if (_Module.GetSiriusApplication()->GetObject(MappedName, L#Type, (IDispatch**)&spDispatch)) throw szError;	\
				_s##MappedName = spDispatch;																	\
				if (!_s##MappedName) throw szError;																\
			}																									\
			if (!_s##MappedName){																				\
				sz.trim();																						\
				if (sz.size()){																					\
					if (_Module.GetSiriusApplication()->GetObject(sz, L#Type, (IDispatch**)&spDispatch)) CParameterMap::ThrowStringError(_Module.GetSiriusApplication().p);					\
					_s##MappedName = spDispatch;																																			\
					if (!_s##MappedName) throw "No object of the correct type is associated with handle '" #MappedName "'";																	\
				}																								\
			}																									\
		}																										\
		if (_s##MappedName){																					\
			long		n;																						\
			_s##MappedName->get_AnalyticObject(&n);																\
			MappedName = reinterpret_cast<MlEq##Type##*>(n);													\
		}																										\
	}

#define map_object_parameter(InputName, Type, MappedName)														\
	CComQIPtr<I##Type>					_s##MappedName;															\
	MlEq##Type##Handle					MappedName;																\
	{																											\
		if ((InputName).vt == VT_DISPATCH){																		\
			try {																								\
				_s##MappedName = dynamic_cast<I##Type*>(const_cast<IDispatch*>((InputName).pdispVal));			\
			} catch (...){																						\
			}																									\
		}																										\
		if (!_s##MappedName){																					\
			CComPtr<IDispatch>				spDispatch;															\
			map_parameter(InputName, CComVariant, v);															\
			if (_Module.GetSiriusApplication()->GetObject(v, L#Type, (IDispatch**)&spDispatch)) CParameterMap::ThrowStringError(_Module.GetSiriusApplication().p);	\
			_s##MappedName = spDispatch;																															\
			if (!_s##MappedName) throw "No object of the correct type is associated with handle '" + estring(v) + "'";												\
		}																										\
		long							n;																		\
		_s##MappedName->get_AnalyticObject(&n);																	\
		MappedName = reinterpret_cast<MlEq##Type##*>(n);														\
	}


/////////////////////////////////////////////////////////////////////////////
//	map_[optional_]com_object_parameter
//
//  These macros are used to map an input object handle (or key) to an
//  COM smart pointer. Usage restrictions and interface requirements
//  are as for map_com_to_analytic and map_analytic_to_com.
//
//	The _optional_ macro returns a null pointer in MappedName if the input
//  name is null. The non-optional version throws a string exception in this
//  case.
//
#define map_optional_com_object_parameter(InputName, Type, MappedName)		\
	CComPtr<I##Type>					MappedName;							\
	{																		\
		CComPtr<IDispatch>				spDispatch;							\
		estring							sz;									\
		{																	\
			map_optional_parameter(InputName, estring, MappedName, "");		\
			sz.assign(MappedName);											\
		}																	\
		sz.trim();															\
		if (sz.size()){														\
			if (_Module.GetSiriusApplication()->GetObject(sz, L#Type, (IDispatch**)&spDispatch)) CParameterMap::ThrowStringError(_Module.GetSiriusApplication().p);		\
			MappedName = dynamic_cast<I##Type*>(spDispatch.p);																				\
			if (!MappedName) throw "No object of the correct type is associated with handle '" + sz + "'";									\
		}																	\
	}

#define map_com_object_parameter(InputName, Type, MappedName)				\
	CComQIPtr<I##Type>					MappedName;							\
	{																		\
		CComPtr<IDispatch>				spDispatch;							\
		estring							sz;									\
		{																	\
			map_parameter(InputName, estring, MappedName);					\
			sz.assign(MappedName);											\
		}																	\
		sz.trim();															\
		if (_Module.GetSiriusApplication()->GetObject(sz, L#Type, (IDispatch**)&spDispatch)) CParameterMap::ThrowStringError(_Module.GetSiriusApplication().p);		\
		MappedName = spDispatch;																										\
		if (!MappedName) throw "No object of the correct type is associated with handle '" + sz + "'";									\
	}


/////////////////////////////////////////////////////////////////////////////
//	map_[com_]object_vector_parameter
//
//  These macros are used to map a vector of object handles (or keys) to a
//  std::vector of analytic or COM smart pointers. Usage restrictions and
//  interface requirements are as for map_com_to_analytic and map_analytic_to_com.
//
//	map_strikes_vector_parameter is a specialisation of map_object_vector_parameter
//	for the strikes case.
//
#define map_object_vector_parameter(InputName, Type, MappedName)			\
	std::vector<CAdapt<CComQIPtr<I##Type> > >	_s##MappedName;				\
	std::vector<MlEq##Type##Handle>				MappedName;					\
	{																		\
		CParameterMap	_pmInput;											\
		CComVariant		v;													\
		CParameterMap	_pm;												\
		if (_Module.GetSiriusApplication()->GetObjects(InputName, L#Type, &v)) CParameterMap::ThrowStringError(_Module.GetSiriusApplication().p);\
		if (_pmInput.SetValue(InputName)) ATLASSERT(false);					\
		if (_pm.SetValue(v)) ATLASSERT(false);								\
		if (_pm.IsBlank()){													\
			std::stringstream		ss;										\
			estring					sz(#InputName);							\
			ss << "No values have been entered into the array parameter'";	\
				if (sz.find('[') == std::string::npos){						\
				ss << sz;													\
			} else {														\
				ss << #MappedName;											\
			}																\
			ss << "'";														\
			throw ss.str();													\
		}																	\
		for (long nRow = 0; nRow < _pm.GetRows(); nRow++){					\
			for (long nCol = 0; nCol < _pm.GetCols(); nCol++){				\
				_pm.GetValue(nRow, nCol, &v);								\
				ATLASSERT(v.vt == VT_DISPATCH);								\
				CComPtr<IDispatch>	spDispatch = v.pdispVal;				\
				CComQIPtr<I##Type>	sp = spDispatch;						\
				if (!sp){													\
					std::string sz;											\
					_pmInput.GetValue(nRow, nCol, &sz);						\
					throw "No object of the correct type is associated with handle '" + sz + "'";		\
				}																\
				long				n;											\
				sp->get_AnalyticObject(&n);										\
				_s##MappedName.push_back(sp);									\
				MlEq##Type##Handle	h(reinterpret_cast<MlEq##Type##*>(n));		\
				MappedName.push_back(h);										\
			}																	\
		}																		\
	}

#define map_com_object_vector_parameter(InputName, Type, MappedName)		\
	std::vector<CAdapt<CComQIPtr<I##Type> > >		MappedName;				\
	{																		\
		CParameterMap	_pmInput;											\
		CComVariant		v;													\
		CParameterMap	_pm;												\
		if (_Module.GetSiriusApplication()->GetObjects(InputName, L#Type, &v)) CParameterMap::ThrowStringError(_Module.GetSiriusApplication().p);		\
		if (_pmInput.SetValue(InputName)) ATLASSERT(false);					\
		if (_pm.SetValue(v)) ATLASSERT(false);								\
		if (_pm.IsBlank()){													\
			std::stringstream		ss;										\
			estring					sz(#InputName);							\
			ss << "No values have been entered into the array parameter'";	\
				if (sz.find('[') == std::string::npos){						\
				ss << sz;													\
			} else {														\
				ss << #MappedName;											\
			}																\
			ss << "'";														\
			throw ss.str();													\
		}																	\
		for (long nRow = 0; nRow < _pm.GetRows(); nRow++){					\
			for (long nCol = 0; nCol < _pm.GetCols(); nCol++){				\
				_pm.GetValue(nRow, nCol, &v);								\
				ATLASSERT(v.vt == VT_DISPATCH);								\
				CComPtr<IDispatch>	spDispatch = v.pdispVal;				\
				CComQIPtr<I##Type>	sp = spDispatch;						\
				if (!sp){																			\
					std::string sz;																	\
					_pmInput.GetValue(nRow, nCol, &sz);												\
					throw "No object of the correct type is associated with handle '" + sz + "'";	\
				}																					\
				MappedName.push_back(sp);															\
			}																						\
		}																							\
	}

#define map_strikes_vector_parameter(InputName, MappedName)						\
	std::vector<std::vector<MlEqStrikeHandle> > 			MappedName;			\
	std::vector<CAdapt<CComQIPtr<IStrikes> > >		_s##MappedName;				\
	{																			\
		map_object_vector_parameter(InputName, Strikes, _##MappedName);			\
		MappedName.resize(_##MappedName.size());								\
		for (long n = 0; n < MappedName.size(); n++){							\
			MappedName[n] = _##MappedName[n]->m_hStrikes;						\
		}																		\
		_s##MappedName = _s_##MappedName;										\
	}


/////////////////////////////////////////////////////////////////////////////
//	[declare/implement]_member_[variable / string]
//
//	These macros are used to implement properties on serialisable objects
//	that form part of their key. The macros call the Rekey function in the
//	corresponding collection to update the collection Handle <-> Key tables.
//
// Note that declare_member_variable will work with analytic-level object
// properties (e.g. MlEqAsset).
//
#define declare_member_variable_rekey_ro(comtype, comproperty)	\
	STDMETHOD(get_##comproperty)(comtype* pVal);

#define declare_member_string_rekey(property)					\
	declare_member_variable_rekey(BSTR, property)

#define declare_member_variable_rekey(comtype, comproperty)		\
	declare_member_variable_rekey_ro(comtype, comproperty)		\
	STDMETHOD(put_##comproperty)(comtype newVal);

#define implement_member_variable_rekey_ro(singularname, collectionname, comtype, comproperty, analytictype, analyticproperty)	\
	STDMETHODIMP C##singularname::get_##comproperty(comtype* pVal)																\
	{																															\
		begin_function																											\
		*pVal = (comtype)m_h->Get##analyticproperty();																			\
		end_function																											\
	}

#define implement_member_string_rekey_ro(singularname, collectionname, comproperty, analyticproperty)							\
	STDMETHODIMP C##singularname::get_##comproperty(BSTR* pVal)																	\
	{																															\
		begin_function																											\
		return estring::GetBSTR(m_h->Get##analyticproperty(), pVal);															\
		end_function																											\
	}

#define implement_member_string_rekey(singularname, collectionname, comproperty, analyticproperty)								\
	implement_member_string_rekey_ro(singularname, collectionname, comproperty, analyticproperty)								\
	STDMETHODIMP C##singularname::put_##comproperty(BSTR newVal)																\
	{																															\
		begin_function																											\
		estring sz(newVal);																										\
		sz.trim(&sz);																											\
		if (sz == m_h->Get##analyticproperty()) return S_OK;																	\
		sz.trim(&sz);																											\
		std::string szOldKey = CComObjectCollectionSerialisableKey(CComPtr<I##singularname>(this));								\
		m_h->Put##analyticproperty(sz);																							\
		CComPtr<I##collectionname> spCollection;																				\
		if (g_pApplication->GetObjectManager().IIDToCollectionPtr(IID_I##singularname, spCollection)) throw "Unhandled exception in C"#singularname##"::put_"#comproperty;					\
		g_pApplication->GetObjectManager().Rekey(CComPtr<I##singularname>(dynamic_cast<I##singularname*>(this)), spCollection, szOldKey, NULL);		\
		end_function																											\
	}

#define implement_member_variable_rekey(singularname, collectionname, comtype, comproperty, analytictype, analyticproperty)		\
	implement_member_variable_rekey_ro(singularname, collectionname, comtype, comproperty, analytictype, analyticproperty)		\
	STDMETHODIMP C##singularname::put_##comproperty(comtype newVal)																\
	{																															\
		begin_function																											\
		if ((analytictype)newVal == m_h->Get##analyticproperty()) return S_OK;													\
		std::string szOldKey = CComObjectCollectionSerialisableKey(CComPtr<I##singularname>(this));								\
		m_h->Put##analyticproperty((analytictype)newVal);																		\
		CComPtr<I##collectionname> spCollection;																				\
		if (g_pApplication->GetObjectManager().IIDToCollectionPtr(IID_I##singularname, spCollection)) throw "Unhandled exception in C"#singularname##"::put_"#comproperty;		\
		g_pApplication->GetObjectManager().Rekey(CComPtr<I##singularname>(dynamic_cast<I##singularname*>(this)), spCollection, szOldKey, NULL);		\
		end_function																											\
	}


//////////////////////////////////////////////////////////////////////
//  propagate_error[_ex]
//
//  Looks at any stored ErrorInfo object and throws its string description.
//
//	The extended version will return the error text associated with an
//  HRESULT if the rich error text is not available.
//
#define propagate_error																\
	{																				\
		std::string szError;														\
		CParameterMap::ErrorHandler(_Module.GetSiriusApplication(), S_OK, &szError);\
		if (szError.size()) throw szError;											\
	}

#define propagate_error_ex(hr)														\
	{																				\
		std::string szError;														\
		CParameterMap::ErrorHandler(_Module.GetSiriusApplication(), hr, &szError);	\
		if (szError.size()) throw szError;											\
	}


/////////////////////////////////////////////////////////////////////////////
//	Old-style macros from this point.
//
//	Do not use since I intend to deprecate them. Use the lower case
//	equivalents since they enforce better disentanglement between analytics
//	and interface.
//
//	Use this to declare and associate an interface-level member variable 
//  with the COM interface. Note that you should use special case macros
//  for strings, variants and objects.
#define DECLARE_MEMBER_VARIABLE(type, member, property)		\
	public:													\
		STDMETHOD(get_##property)(type* pVal)				\
		{													\
			*pVal = member;									\
			return S_OK;									\
		}													\
		STDMETHOD(put_##property)(type newVal)				\
		{													\
			member = newVal;								\
			return S_OK;									\
		}													\
	protected:												\
		type member;

//	as for DECLARE_MEMBER_VARIABLE but for read only properties
#define DECLARE_MEMBER_VARIABLE_RO(type, member, property)	\
	public:													\
		STDMETHOD(get_##property)(type* pVal)				\
		{													\
			*pVal = member;									\
			return S_OK;									\
		}													\
	protected:												\
		type member;

//	this macro is the special case of DECLARE_MEMBER_VARIABLE for string properties
#define DECLARE_MEMBER_STRING(member, property)				\
	public:													\
		STDMETHOD(get_##property)(BSTR* pVal)				\
		{													\
			return member.GetBSTR(pVal);					\
		}													\
		STDMETHOD(put_##property)(BSTR newVal)				\
		{													\
			HRESULT hr;										\
			if (hr = member.Set(newVal)) return hr;			\
			member.trim();									\
			return S_OK;									\
		}													\
	protected:												\
		estring		member;

//	this macro is the special case of DECLARE_MEMBER_VARIABLE for a VARIANT_BOOL type property
#define DECLARE_MEMBER_BOOL(member, property)				\
	public:													\
		STDMETHOD(get_##property)(YesNoEnum* pVal)			\
		{													\
			*pVal = member ? Yes : No;						\
			return S_OK;									\
		}													\
		STDMETHOD(put_##property)(YesNoEnum newVal)			\
		{													\
			member = newVal ? true : false;					\
			return S_OK;									\
		}													\
	protected:												\
		bool		member;

//	this macro is the special case of DECLARE_MEMBER_VARIABLE_RO for a BSTR type property
#define DECLARE_MEMBER_STRING_RO(member, property)			\
	public:													\
		STDMETHOD(get_##property)(BSTR* pVal)				\
		{													\
			return member.GetBSTR(pVal);					\
		}													\
	protected:												\
		estring		member;

//	this macro is the special case of DECLARE_MEMBER_VARIABLE for a VARIANT type property
#define DECLARE_MEMBER_VARIANT(member, property)			\
	public:													\
		STDMETHOD(get_##property)(VARIANT* pVal)			\
		{													\
			return CComVariant(member).Detach(pVal);		\
		}													\
		STDMETHOD(put_##property)(VARIANT newVal)			\
		{													\
			return member.Copy(&newVal);					\
		}													\
	protected:												\
		CComVariant member;

//	this macro is used in product type implementation header files to declare a member object variable associated with a COM interface
#define DECLARE_MEMBER_OBJECT(Type, member, property)		\
	public:													\
		STDMETHOD(get_##property)(I##Type** pVal)			\
		{													\
			return _sp##member.CopyTo(pVal);				\
		}													\
		STDMETHOD(put_##property)(I##Type* newVal)			\
		{													\
			_sp##member = newVal;							\
			if (newVal){									\
				long n;										\
				_sp##member->get_AnalyticObject(&n);		\
				member = reinterpret_cast<MlEq##Type##*>(n);\
			} else {										\
				member = NULL;								\
			}												\
			return S_OK;									\
		}													\
	protected:												\
		CComPtr<I##Type>				_sp##member;		\
		MlEq##Type##Handle				member;

#define DECLARE_MEMBER_COM_OBJECT(type, member, property)	\
	public:													\
		STDMETHOD(get_##property)(type** pVal)				\
		{													\
			return member.CopyTo(pVal);						\
		}													\
		STDMETHOD(put_##property)(type* newVal)				\
		{													\
			member = newVal;								\
			return S_OK;									\
		}													\
	protected:												\
		CComPtr<type> member;

#define DECLARE_MEMBER_COM_OBJECT_RO(type, member, property)\
	public:													\
		STDMETHOD(get_##property)(type** pVal)				\
		{													\
			return member.CopyTo(pVal);						\
		}													\
	protected:												\
		CComPtr<type> member;


#endif
