//	objectmanager.h
//
//	Author :				David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifndef _OBJECTMANAGER_H
#define _OBJECTMANAGER_H

#include "resource.h"
#include "xmlstreamer.h"
#include "excelinterface.h"
#include "comobjectcollectionserialisablekey.h"
#include "guidex.h"

#undef max

//////////////////////////////////////////////////////////////////////
//	CSiriusSingularObject
//
//	Stores information about a sirius object like its CLISD and
//  corresponding collection.
//
class CSiriusSingularObject
{	
public:
	enum object_type {
		// these enumerations should be in sequence 1, 2, 4, 8, 16 etc.
		ProductType = 1,			// the object is a product type imported from a quant library
		MarketDataComponent = 2,	// the object is exposed to the application object via a market data object		
		StandardObject = 4			// any other object type						
	};
	
	CSiriusSingularObject(void);	
	CLSID								GetCLSID(void) const;
	bool								GetProductLike(void) const;
	IID									GetIID(void) const;
	IID									GetCollectionIID(void) const;
	CLSID								GetCollectionCLSID(void) const;
	const std::string&                  GetCollectionName(void) const;
	std::string							GetCollectionName_lc(void) const;
	GUID								GetLIBID(void) const;
	const std::string&					GetName(void) const;
	const std::string&					GetProperName(void) const;
	object_type							GetObjectType(void) const;
	CComPtr<IDispatch>					GetCollectionPtr(void) const;	
	template<class ICollection> HRESULT	Install(const GUID& libid, const CLSID& clsid, const IID& iid, object_type ot, bool bProductLike, const std::string& szName, CComPtr<ICollection>& spCollection, const CLSID& clsidCollection)
	{
		HRESULT hr = S_OK;		
		m_bProductLike = bProductLike;
		m_clsid = clsid;
		m_iid = iid;
		m_libid = libid;		
		m_ot = ot;		
		if (m_ot == StandardObject){
			// create the collection
			hr = spCollection.CoCreateInstance(clsidCollection);
			m_iidCollection = __uuidof(ICollection);
			m_clsidCollection = clsidCollection;
			ATLASSERT(clsidCollection != CLSID_NULL);
		} else if (m_ot == MarketDataComponent){
			// the collection will have already been created
			ATLASSERT(spCollection);
			m_iidCollection = __uuidof(ICollection);
			m_clsidCollection = clsidCollection;
			ATLASSERT(clsidCollection != CLSID_NULL);
		} else if (m_ot == ProductType){
			// a collection is never created in this case
			ATLASSERT(!spCollection);
			m_iidCollection = IID_NULL;	
			m_clsidCollection = CLSID_NULL;
		} else {
			ATLASSERT(false);			
		}
		if (m_spCollection = spCollection){
			CParameterMap::GetObjectName(m_spCollection, &m_szCollectionName);
		} else {
			m_szCollectionName.erase();
		}
		
		if (!szName.length() && !hr){
			hr = CParameterMap::CLSIDToObjectName(clsid, &m_szProperName);
		} else {
			m_szProperName = szName;
		}
		estring::StripWhiteSpace(&m_szProperName);
		m_szName = m_szProperName;
		estring::lc(&m_szName);		
		ATLTRACE("Installing object. Singular name is '%s', Proper case is '%s', Return code is 0x%x (zero denotes success).\n", m_szName.c_str(), m_szProperName.c_str(), hr);
		return hr;
	}	
	bool											IsValid(void) const;
												
protected:										
	bool											m_bProductLike;					// true if the object has a header row which we display in bold
	CLSID											m_clsid;						// CLSID of the singular object
	IID												m_iid;							// IID of the singular object
	GUID											m_libid;
	object_type										m_ot;
	std::string										m_szName;						// object singular name
	std::string										m_szProperName;					// object name in proper case	
	CComPtr<IDispatch>								m_spCollection;					// pointer to the collection interface (null for product types)
	IID												m_iidCollection;				// IID of the associated collection object (null for product types)
	CLSID											m_clsidCollection;				// CLSID of the assocaited collection object (null for product types)
	std::string										m_szCollectionName;				// name of the corresponding collection (null for product types)
};												


//////////////////////////////////////////////////////////////////////
//	ObjectManager
//
//	Functionality splits as follows:
//		1) General casts.
//		2) Global index table functions.
//		3) CLSID, IID, Name Mapping functions.
//		4) Object manager initialisation.
//		5) Object information functions.
//
class CObjectManager
{
// General casts - These allow you to return a vector of underlying object base class pointers. These macros
// will work if the underlying objects all inherit from that base class.
public:	
	#define GetObjectVectorStart(ParameterMap, SmartPointerVector, BaseClassPointerVector, BaseClassPointerType)			\
	{																														\
		typedef	BaseClassPointerType					handle_array_start_type;											\
		std::vector<BaseClassPointerType*>*				pBaseClassPointerVector = &BaseClassPointerVector;					\
		std::vector<CAdapt<CComQIPtr<IDispatch> > >*	pSmartPointerVector = &SmartPointerVector;							\
		const CParameterMap*					ppm = &ParameterMap;														\
		if (CObjectManager::GetObjectVector(ParameterMap, SmartPointerVector)) propagate_error;								\
		pBaseClassPointerVector->clear();																					\
		for (long nElement = 0; nElement < SmartPointerVector.size(); nElement++){											\
			IID iid;																										\
			if (CParameterMap::GetObjectIID(CComDispatchDriverEx((*pSmartPointerVector)[nElement].m_T), &iid)) ATLASSERT(false);\
			if (false){
	#define GetObjectVectorCast(Class)																	\
			} else if (iid == IID_I##Class){															\
				I##Class*				p_com;															\
				C##Class*				p_obj;															\
				p_com = dynamic_cast<I##Class*>((*pSmartPointerVector)[nElement].m_T.p);				\
				p_obj = (C##Class*)p_com;																\
				pBaseClassPointerVector->push_back(&*p_obj->m_h);
	#define GetObjectVectorEnd(Error, iid)																\
			} else {																					\
				std::string	sz;																			\
				ppm->GetValue(nElement, 0, &sz);														\
				CParameterMap::ThrowComErrorRS(Error, sz, iid);											\
			}																							\
		}																								\
	}
	static void Cast(const CParameterMap& pmIn, std::vector<CAdapt<CComQIPtr<IDispatch> > >& asp, std::vector<MlEqVolDataHandle>& apOut);

protected:
	template<class I>
	static HRESULT GetObjectVector(const CParameterMap& pmIn, std::vector<CAdapt<CComQIPtr<I> > >& asp)
	{		
		asp.clear();
		if (pmIn.GetCols() != 1) return CParameterMap::ReturnErrorR(IDS_COLUMNS_INVALID);
		for (long nRow = 0; nRow < pmIn.GetRows(); nRow++){
			CComPtr<IDispatch>			spDispatch;
			CComQIPtr<I>				sp;
			CComVariant					vElement;
			if (pmIn.GetValue(nRow, 0, &vElement)) return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
			if (vElement.vt == VT_DISPATCH){
				spDispatch = vElement.pdispVal;
			} else {
				// assume vElement is an index
				if (g_pApplication->GetObjectManager().GetObject(vElement, spDispatch)){
					return CParameterMap::ReturnErrorRV(IDS_NO_OBJECT_WITH_HANDLE, vElement);
				}
			}
			sp = spDispatch;
			asp.push_back(sp);
		}
		return S_OK;
	}

// Global index table.
protected:		
	struct index_and_count							// this structure is an index name with a reference counter
	{
		estring			szIndex;
		long			nCount;
	};
	typedef std::map<CAdapt<CComPtr<IDispatch> >, index_and_count>	ObjectToIndexType;
	ObjectToIndexType												m_mapObjectToIndex;
	
public:	
	template<class ISingular, class ICollection> void AddIndex(CComPtr<ISingular> spSingular, CComPtr<ICollection> spCollection, const std::string& szObjectName, estring szIndex)	
	{
		// Only ever call this function from the CComObjectCollection or its child classes
		ATLASSERT(spSingular);
		ATLASSERT(spCollection);
		if (!IsMaintainedCollection(spCollection)) return;
		if (szIndex.left("::").CompareNoCase(szObjectName)){
			// we need to prefix the index with the object name
			szIndex.assign(szObjectName + "::" + szIndex);
		}
		ObjectToIndexType::iterator		it;	
		if ((it = m_mapObjectToIndex.find(CComPtr<IDispatch>(spSingular.p))) == m_mapObjectToIndex.end()){
			// Object has not been added yet. Add it with a reference count of 1. szIndex 'wins' in the sense that we will always
			// use that name in the map even if this function is called with the same object but a different value of szIndex.=
			index_and_count iac;
			iac.nCount = 1;
			iac.szIndex.assign(szIndex);
			m_mapObjectToIndex[CComPtr<IDispatch>(spSingular.p)] = iac;		
		} else {
			// Object has already been added. All we do is increment the reference count. We don't change the index name (see above).			
			it->second.nCount++;
		}
	}
	bool GetIndex(CComPtr<IDispatch> spObject, std::string* pszIndex, bool bThrow) const;
	template <class ISingular> bool							HasIndex(CComPtr<ISingular> spSingular) const
	{
		return m_mapObjectToIndex.find(CComPtr<IDispatch>(spSingular.p)) != m_mapObjectToIndex.end();
	}
	void Reindex(CComPtr<IDispatch> spObject, const std::string& szOldIndex, const std::string& szNewIndex);	
	template<class ISingular, class ICollection> void RemoveIndex(CComPtr<ISingular> spSingular, CComPtr<ICollection> spCollection)	
	{
		// Only ever call this function from the CComObjectCollection or its child classes		
		ATLASSERT(spSingular);
		ATLASSERT(spCollection);
		if (!IsMaintainedCollection(spCollection)) return;		
		ObjectToIndexType::iterator								it;
		it = m_mapObjectToIndex.find(CComPtr<IDispatch>(spSingular.p));
		if (it == m_mapObjectToIndex.end()) return;
		if (it->second.nCount <= 1){
			// we can now delete the object
			m_mapObjectToIndex.erase(it);
		} else {
			it->second.nCount--;		
		}		
	}

// CLSID, IID, Name Mapping functions
public:	
	HRESULT											CLSIDToName(const CLSID& clsid, std::string* pszName) const;
	HRESULT											CLSIDToProductLike(const CLSID& clsid, bool* pbProductLike) const;
	HRESULT											CLSIDToProperName(const CLSID& clsid, std::string* pszName) const;
	HRESULT											CollectionCLSIDToCollectionName(const CLSID& clsid, std::string* pszName) const;
	HRESULT											CollectionNameToCollectionCLSID(std::string szName, CLSID* pclsid) const;
	HRESULT											CollectionNameToCollectionIID(std::string szName, IID* piid) const;
	HRESULT											CollectionNameToCollectionProgID(const std::string& szName, CComBSTR& sProgID) const;
	HRESULT											IIDToCLSID(const IID& iid, CLSID* pclsid) const;
	HRESULT											IIDToCollectionName(const IID& iid, std::string* pszName) const;
	template <class ICollection> HRESULT			IIDToCollectionPtr(const IID& iid, CComPtr<ICollection>& spCollection) const
	{
		HRESULT							hr;
		CComPtr<IDispatch>				sp;
		if (hr = IIDToCollectionPtr(iid, sp)) return hr;
		spCollection = dynamic_cast<ICollection*>(sp.p);
		return spCollection ? S_OK : E_FAIL;
	}
	HRESULT											IIDToLIBID(const IID& iid, GUID* plibid) const;
	HRESULT											IIDToName(const IID& iid, std::string* pszName) const;	
	HRESULT											NameToCLSID(std::string szName, CLSID* pclsid) const;
	HRESULT											NameToCollectionName(std::string szName, std::string* pszCollectionName) const;	
	HRESULT											NameToIID(const std::string& szName, IID* piid) const;	
	HRESULT											NameToProgID(std::string szName, CComBSTR& sProgID) const;
	std::string										NameToProperName(std::string szName) const;
	
protected:
	HRESULT											CollectionIIDToCollectionCLSID(const IID& iid, CLSID* pclsid) const;
	HRESULT											CollectionIIDToIID(const IID& iid, IID* piid) const;
	HRESULT											HandleToObjectIID(const IID& iidSingular, const CComVariant& vHandle, IID* piid) const;
	HRESULT											HandleToObjectIID(const IID& iidSingular, std::string szHandle, IID* piid) const;
	HRESULT											IIDToCollectionPtr(const IID& iid, CComPtr<IDispatch>& sp) const;
	HRESULT											NameToCollectionPtr(std::string szName, CComPtr<IDispatch>& sp) const;

// Object manager initialisation.
public:
	template<class ICollection> HRESULT				InstallMarketDataComponent(const GUID& libid, const CLSID& clsid, const IID& iid, CComPtr<ICollection> spCollection, const CLSID& clsidCollection, bool bProductLike)
	{
		HRESULT						hr;
		CSiriusSingularObject				sso;	
		if (hr = sso.Install(libid, clsid, iid, CSiriusSingularObject::MarketDataComponent, bProductLike, "", spCollection, clsidCollection)) return hr;
		return Install(sso);
	}
	HRESULT											InstallProductType(const GUID& libid, const CLSID& clsid, const IID& iid);
	HRESULT											InstallProductType(const std::string& szLIBID, const std::string& szCLSID, const std::string& szIID, const std::string& szName);	
	template<class ICollection> HRESULT				InstallStandardObject(const GUID& libid, const CLSID& clsid, const IID& iid, bool bProductLike, CComPtr<ICollection>& spCollection, const CLSID& clsidCollection)
	{		
		HRESULT hr;
		CSiriusSingularObject sso;		
		if (hr = sso.Install(libid, clsid, iid, CSiriusSingularObject::StandardObject, bProductLike, "", spCollection, clsidCollection)) return hr;
		return Install(sso);
	}		
protected:
	HRESULT											Install(const CSiriusSingularObject& sso);

// Object information functions.
public:
	HRESULT											CreateObject(const CLSID& clsid, CComPtr<IDispatch>& spObject) const;
	HRESULT											CreateObject(const std::string& szObject, CComPtr<IDispatch>& spObject) const;
	HRESULT											ImplementGetValue(CComPtr<IDispatch> spObject, VARIANT* pVal) const;
	HRESULT											ImplementGetValue(CComPtr<IDispatch> spObject, const estring& szProperty, VARIANT* pVal) const;
	HRESULT											GetCollectionIndex(CComPtr<IDispatch> spObjectContainingCollection, CComPtr<IDispatch> spIndexObject, const estring& szCollectionName, std::string* pszIndex) const;
	HRESULT											GetFileName(const IID& iid, const std::string& szLocation, std::string szName, DataSourceEnum ds, long nDate, std::string* pszName) const;
	template<class ISingular> HRESULT				GetObject(const CComVariant& vIn, CComPtr<ISingular>& spSingular) const
	{
		CComPtr<IDispatch>	spObject;
		HRESULT				hr;
		if (hr = GetObject(vIn, __uuidof(ISingular), spObject)) return hr;		
		spSingular = CComQIPtr<ISingular>(spObject);
		return spSingular ? S_OK : E_FAIL;
	}
	HRESULT											GetObject(const CComVariant& vIn, const IID& IID_ISingular, CComPtr<IDispatch>& spSingular) const;
	HRESULT											GetObject(const CComVariant& vIn, std::string& szObjectName, CComPtr<IDispatch>& spSingular) const;
	HRESULT											GetObjectListOfType(long nTypes, VARIANT* pVal) const;
	HRESULT											GetObjectProperties(INVOKEKIND inv, CComPtr<IDispatch> spObject, VARIANT* pVal, CParameterMap* ppm) const;
	HRESULT											GetObjectProperties(INVOKEKIND inv, const std::string& szObject, VARIANT* pVal, CParameterMap* ppm) const;
	bool											IsCLSIDSiriusCollection(const CLSID& clsid) const;
	// Returns true if the collection is one of those maintained in this class.
	template <class I> bool							IsMaintainedCollection(CComPtr<I> sp) const
	{
		return m_mapCollectionPtrToBool.find(CComPtr<IDispatch>(sp.p)) != m_mapCollectionPtrToBool.end();
	}
	//	Returns true if the object is savable to a Sirius database
	template<class ISingular> bool					IsSerialisable(CComPtr<ISingular> spObject) const
	{	
		// We test the object for the "Name" property. ToDo - we need a better way.
		try {
			DISPID									dispid;			
			if (!spObject->GetIDsOfNames(IID_NULL, &CComBSTR(L"Name"), 1, LOCALE_USER_DEFAULT, &dispid)){									
				return true;
			} else {
				return false;
			}
		} catch (...){
			ATLASSERT(false);
			return false;
		}						
	}
	bool											IsSiriusCollection(std::string sz) const;	
	bool											IsSiriusCollection(CComPtr<IDispatch> sp) const;
	bool											IsSiriusProduct(std::string sz) const;
	bool											IsSiriusSingular(std::string sz) const;
protected:	
	bool											IsIIDSiriusCollection(const IID& iid) const;
	bool											IsIIDSiriusSingular(const IID& iid) const;	
	bool											IsSiriusMarketDataObject(std::string sz) const;

// Operations on objects.
public:
	void											Clear(CSiriusSingularObject::object_type ot);
	HRESULT											ClearCollections(void) const;
	HRESULT											ImplementPutValue(const VARIANT& vIn, CComPtr<IDispatch> spObject);
	HRESULT											ImplementPutValue(std::vector<CParameterMap>& vpm, CComPtr<IDispatch> spObject);
	HRESULT											ImplementPutValue(CComPtr<IDispatch> spObject, const std::string& szProperty, const CComVariant& Value);
	// Adds an object to the appropriate maintained collection and returns an index
	// for future object identification.
	template<class ISingular> estring	InsertObject(CComPtr<ISingular> spSingular, bool bExternal)
	{						
		static std::map<std::string, long>		mapWatermark;
		std::string								szObjectName;		
		CComPtr<IDispatch>						spCollection;			// maintained collection corresponding to spSingular
		estring									szIndex;				// index returned
				
		// Return the existing handle if the object is already in a maintained collection.
		if (GetIndex(spSingular.p, &szIndex, false)){		
			return szIndex;			
		}
		if (CParameterMap::GetObjectName(spSingular, &szObjectName)) propagate_error;
		if (NameToCollectionPtr(szObjectName, spCollection)) propagate_error;		
	
		if (bExternal){
			// Insert an external (i.e. Excel-based) handle.			
			long						nCalculated;			
			szIndex = CExcelInterface::GetCellBasedHandle(szObjectName);
			
			// The cell based handle is of the form ObjectName::CellAddress. So we are variant in time,
			// we need to append a number denoting the number of times the handle has been recalculated.
			
			std::string					szMapKey = szObjectName + ">>" + szIndex;
			if (!(nCalculated = m_mapCalculation[szMapKey])){
				// Not yet calculated. Notice that we can assume that nCalculated is zero if the map does not contain szHandle.
				nCalculated = 1;
				// If the object is serialisable then we need to delete any item with the same key as spSingular.
				if (IsSerialisable(spSingular)){
					CComDispatchDriverEx(spCollection).Invoke1(L"Remove", &CComVariant(spSingular));
				}
			} else {
				// has been calculated - remove it (we check for any errors later when adding the object to the collection)
				nCalculated++;								
				CComDispatchDriverEx(spCollection).Invoke1(L"Remove", &szIndex.GetValue());				
			}

			if (CComDispatchDriverEx(spCollection).Invoke2(L"Add", &szIndex.GetValue(), &CComVariant(spSingular))) propagate_error;
			m_mapCalculation[szMapKey] = nCalculated;
			// append the calculation number to the handle
			szIndex += "::";
			szIndex += estring(nCalculated);						
			ATLTRACE("Inserting external handle '%s'.\n", szIndex.c_str());
			return szIndex;
		} else {
			// Insert an internal handle (if serialisable then the key (adorned with the object name) will do, else use ObjectName::Internal[Watermark]).
			if (IsSerialisable(spSingular)){				
				if (CComDispatchDriverEx(spCollection).Invoke2(L"Add", &CComVariant(), &CComVariant(spSingular))) propagate_error;
				szIndex = CComObjectCollectionSerialisableKey(spSingular.p).GetKeyAndObjectName(false, false);
			} else {
				if (mapWatermark.find(szObjectName) == mapWatermark.end()){
					// Not yet added. Therefore initialise it.
					mapWatermark[szObjectName] = 0;
				}	
				std::map<std::string, long>::iterator it = mapWatermark.find(szObjectName);
				if (++it->second == std::numeric_limits<long>::max()) throw "Fatal error! Internal handle limit (" + estring(it->second) + ") reached for object '" + szObjectName + "'. Sirius.dll must be restarted.";
				szIndex = szObjectName + "::Internal" + estring(it->second);
				if (CComDispatchDriverEx(spCollection).Invoke2(L"Add", &szIndex.GetValue(), &CComVariant(spSingular))) propagate_error;
			}
			ATLTRACE("Inserting internal handle '%s'.\n", szIndex.c_str());
			return szIndex;
		}		
	}
	void											Refresh(CComPtr<IDispatch> spObject, bool bRecursive) const;
	template<class ISingular, class ICollection> void	Rekey(CComPtr<ISingular> spSingular, CComPtr<ICollection> spCollection, const std::string& szOldKey, std::string* pszNewKey) const
	//	pszNewKey - (returned, nullable)
	{
		CComBSTR									sNewKey;
		if (spCollection->Rekey(spSingular, estring::GetBSTR(szOldKey), &sNewKey)) propagate_error;
		if (pszNewKey) pszNewKey->assign(estring(sNewKey));
	}
protected:
	HRESULT											ImplementPutValue(CComPtr<IDispatch> spObject, const std::string& szProperty, const CComVariant& Value, CComPtr<ITypeInfo> pti, FUNCDESC* pFuncDesc, MEMBERID memid);
	void											Refresh(std::set<CComPtr<IDispatch> >* pset, CComPtr<IDispatch> spObject, bool bRecursive, DataSourceEnum ds, long nDate) const;
	CComPtr<IDispatch>								RefreshGetNewObject(IID iid, const CComVariant& vObject, DataSourceEnum ds, long nDate) const;

// Data members.
protected:											
	std::map<std::string, long>						m_mapCalculation;						// calculation number map
	std::map<GuidEx, GuidEx>						m_mapIIDToCLSID;
	std::map<std::string, CSiriusSingularObject>	m_mapNameToObject;
	std::map<GuidEx, CComPtr<IDispatch> >			m_mapIIDToCollectionPtr;
	std::map<CAdapt<CComPtr<IDispatch> >, bool>		m_mapCollectionPtrToBool;				// we can make this map do something more useful if we need to
	std::map<GuidEx, GuidEx>						m_mapCollectionIIDToIID;
	std::map<GuidEx, GuidEx>						m_mapCollectionIIDToCollectionCLSID;
	std::map<GuidEx, GuidEx>						m_mapCollectionCLSIDToCollectionIID;
	std::map<GuidEx, std::string>					m_mapCollectionCLSIDToCollectionName;	// case insensitive name to CLSID		
	std::map<std::string, GuidEx>					m_mapCollectionNameToCollectionIID;		// case insensitive name to IID
	std::map<std::string, GuidEx>					m_mapCollectionNameToCollectionCLSID;	// case insensitive name to CLSID
	std::map<GuidEx, bool>							m_mapCLSIDToProductLike;
	std::map<GuidEx, std::string>					m_mapCLSIDToName;
	std::map<GuidEx, std::string>					m_mapCLSIDToProperName;
	std::map<GuidEx, std::string>					m_mapIIDToName;
	std::map<GuidEx, GuidEx>						m_mapIIDToLIBID;
};											
													
#endif
