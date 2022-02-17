//	objectmanager.cpp: implementation of the CObjectManager class.
//
//	Author :		   David Cuin
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "objectmanager.h"
#include "comobjectcollectionserialisablekey.h"
#include "siriusapplication.h"
#include "mleqdate.h"

////////////////////////////////////////////////////////////////////////////
//	CSiriusSingularObject implementation
//
CSiriusSingularObject::CSiriusSingularObject(void)
{
	m_bProductLike = false;
	m_clsid = m_clsidCollection = CLSID_NULL;
	m_iid = m_iidCollection = IID_NULL;
	m_libid = GUID_NULL;	
	m_ot = StandardObject;	
}
IID	CSiriusSingularObject::GetCollectionIID(void) const
{
	return m_iidCollection;
}
CLSID CSiriusSingularObject::GetCollectionCLSID(void) const
{
	return m_clsidCollection;
}
const std::string& CSiriusSingularObject::GetCollectionName(void) const
{
	return m_szCollectionName;
}
std::string CSiriusSingularObject::GetCollectionName_lc(void) const
{
	std::string szCollectionName_lc = m_szCollectionName;
	estring::lc(&szCollectionName_lc);
	return szCollectionName_lc;
}
CLSID CSiriusSingularObject::GetCLSID(void) const
{
	return m_clsid;
}
CComPtr<IDispatch> CSiriusSingularObject::GetCollectionPtr(void) const
{
	return m_spCollection;
}
bool CSiriusSingularObject::GetProductLike(void) const
{
	return m_bProductLike;
}
IID CSiriusSingularObject::GetIID(void) const
{
	return m_iid;
}
GUID CSiriusSingularObject::GetLIBID(void) const
{
	return m_libid;
}
CSiriusSingularObject::object_type	CSiriusSingularObject::GetObjectType(void) const
{
	return m_ot;
}
const std::string& CSiriusSingularObject::GetName(void) const
{
	return m_szName;
}
const std::string& CSiriusSingularObject::GetProperName(void) const
{
	return m_szProperName;
}
bool CSiriusSingularObject::IsValid(void) const
{
	return m_libid != GUID_NULL && m_clsid != CLSID_NULL && m_iid != IID_NULL && m_szName.length() != 0;
}


//////////////////////////////////////////////////////////////////////////////
//	CObjectManager implementation
//
bool CObjectManager::GetIndex(CComPtr<IDispatch> spObject, std::string* pszIndex, bool bThrow) const
{
	ObjectToIndexType::const_iterator it;
	it = m_mapObjectToIndex.find(spObject);
	if (it != m_mapObjectToIndex.end()){
		pszIndex->assign(it->second.szIndex);
		return true;
	}
	if (bThrow){
		std::string szName;
		CParameterMap::GetObjectName(spObject, &szName);
		throw "Could not find a handle for the object with name '" + szName + "'";
	}
	return false;
}
//
//	This should be called in response to a change in the key associated with a serialisable object. This
//  function replaces an entry (if any) in the handle table with a new handle name.
//
/*static*/ void CObjectManager::Reindex(CComPtr<IDispatch> spObject, const std::string& szOldIndex, const std::string& szNewIndex)
{
	ObjectToIndexType::iterator it;
	it = m_mapObjectToIndex.find(spObject);
	if (it == m_mapObjectToIndex.end()) return;
	
	// it->second.szIndex contains the old index along with the object name
	std::string szObjectName = it->second.szIndex.left("::");
	ATLASSERT(szObjectName.size());	
	if (it->second.szIndex.right("::") != szOldIndex) return;
	it->second.szIndex = szObjectName + "::" + szNewIndex;
}
//
//	clears the object manager records of all objects of a given type
//
void CObjectManager::Clear(CSiriusSingularObject::object_type ot)
{
	bool								bFinished;	
	do {
		bFinished = true;		
		for (std::map<std::string, CSiriusSingularObject>::iterator it = m_mapNameToObject.begin(); it != m_mapNameToObject.end(); ++it){
			if (it->second.GetObjectType() == ot){
				std::string szCollectionName_lc(it->second.GetCollectionName_lc());
				m_mapIIDToCLSID.erase(it->second.GetIID());
				m_mapIIDToCollectionPtr.erase(it->second.GetIID());
				m_mapCollectionPtrToBool.erase(it->second.GetCollectionPtr());
				m_mapCollectionIIDToIID.erase(it->second.GetCollectionIID());
				m_mapCollectionIIDToCollectionCLSID.erase(it->second.GetCollectionIID());
				m_mapCollectionCLSIDToCollectionIID.erase(it->second.GetCollectionCLSID());
				m_mapCollectionCLSIDToCollectionName.erase(it->second.GetCollectionCLSID());
				m_mapCollectionNameToCollectionIID.erase(szCollectionName_lc);
				m_mapCollectionNameToCollectionCLSID.erase(szCollectionName_lc);				
				m_mapCLSIDToProductLike.erase(it->second.GetCLSID());
				m_mapCLSIDToName.erase(it->second.GetCLSID());
				m_mapCLSIDToProperName.erase(it->second.GetCLSID());
				m_mapIIDToName.erase(it->second.GetIID());
				m_mapIIDToLIBID.erase(it->second.GetIID());
				m_mapNameToObject.erase(it);				
				bFinished = false;
				break;
			}
		}
	} while (!bFinished);
}
//
//	Delete all the handles in the object manager
//
HRESULT CObjectManager::ClearCollections(void) const
{
	std::map<std::string, CSiriusSingularObject>::const_iterator	it;	
	HRESULT													hr;
	for (it = m_mapNameToObject.begin(); it != m_mapNameToObject.end(); it++){
		CComPtr<IDispatch> sp = it->second.GetCollectionPtr();
		ATLTRACE("Clearing collection '%s'\n", it->first.c_str());
		if (sp && (hr = CComDispatchDriverEx(sp).Invoke0(L"Clear"))) return hr;
	}
	return S_OK;
}
//
//	Maps a singular object CLSID to its name.
//
HRESULT CObjectManager::CLSIDToName(const CLSID& clsid, std::string* pszName) const
{	
	std::map<GuidEx, std::string>::const_iterator it;
	if ((it = m_mapCLSIDToName.find(clsid)) == m_mapCLSIDToName.end()) return E_FAIL;
	pszName->assign(it->second);
	return S_OK;
}
//
//	Returns true if the object denoted by an input CLSID should be presented
//	in Excel in the same way as a product.
//
HRESULT CObjectManager::CLSIDToProductLike(const CLSID& clsid, bool* pbProductLike) const
{
	std::map<GuidEx, bool>::const_iterator it;	
	if ((it = m_mapCLSIDToProductLike.find(clsid)) == m_mapCLSIDToProductLike.end()) return E_FAIL;
	*pbProductLike = it->second;
	return S_OK;
}
//
//	Maps a singular object CLSID to its proper name.
//
HRESULT CObjectManager::CLSIDToProperName(const CLSID& clsid, std::string* pszName) const
{
	std::map<GuidEx, std::string>::const_iterator it;
	if ((it = m_mapCLSIDToProperName.find(clsid)) == m_mapCLSIDToProperName.end()) return E_FAIL;
	pszName->assign(it->second);
	return S_OK;
}
//
//	Maps an input collection CLSID to the corresponding collection name
//
HRESULT CObjectManager::CollectionCLSIDToCollectionName(const CLSID& clsid, std::string* pszName) const
{
	std::map<GuidEx, std::string>::const_iterator it;
	if ((it = m_mapCollectionCLSIDToCollectionName.find(clsid)) == m_mapCollectionCLSIDToCollectionName.end()) return E_FAIL;
	pszName->assign(it->second);
	return S_OK;
}
//
//	Maps an input IID of a Sirius collection to the corresponding CLSID.
//
HRESULT CObjectManager::CollectionIIDToCollectionCLSID(const IID& iid, CLSID* pclsid) const
{
	std::map<GuidEx, GuidEx>::const_iterator it;	
	if ((it = m_mapCollectionIIDToCollectionCLSID.find(iid)) == m_mapCollectionIIDToCollectionCLSID.end()){
		ATLASSERT(false);		// this should never fail
		return E_FAIL;
	}
	*pclsid = it->second;
	return S_OK;
}
//
//	maps an input IID of a Sirius collection to the corresponding singular IID
//
HRESULT CObjectManager::CollectionIIDToIID(const IID& iid, IID* piid) const
{
	std::map<GuidEx, GuidEx>::const_iterator it;
	if ((it = m_mapCollectionIIDToIID.find(iid)) == m_mapCollectionIIDToIID.end()) return E_FAIL;	// this can fail if iid is not a collection
	*piid = it->second;
	return S_OK;
}
//
//	Returns the collection CLSID associated with an input collection name.
//
HRESULT CObjectManager::CollectionNameToCollectionCLSID(std::string szName, CLSID* pclsid) const
{		
	std::map<std::string, GuidEx>::const_iterator it;
	estring::lc(&szName);
	if ((it = m_mapCollectionNameToCollectionCLSID.find(szName)) == m_mapCollectionNameToCollectionCLSID.end()){
		// the input name could be of the form Library.Object.Version
		std::vector<std::string>		vector;
		estring::Split(szName, ".", &vector);	
		if (vector.size() != 3) return E_FAIL;
		return CollectionNameToCollectionCLSID(vector[1], pclsid);
	}
	
	*pclsid = it->second;
	return S_OK;	
}
//
//	Returns the collection IID associated with an input collection name.
//
HRESULT CObjectManager::CollectionNameToCollectionIID(std::string szName, IID* piid) const
{			
	std::map<std::string, GuidEx>::const_iterator it;
	estring::lc(&szName);
	if ((it = m_mapCollectionNameToCollectionIID.find(szName)) == m_mapCollectionNameToCollectionIID.end()) return E_FAIL;	// this could fail if szName is not a collection
	*piid = it->second;
	return S_OK;	
}
//
//	Returns the collection ProgID associated with an input collection name.
//
HRESULT CObjectManager::CollectionNameToCollectionProgID(const std::string& szName, CComBSTR& sProgID) const
{
	CLSID								clsid;	
	HRESULT								hr;
	if (hr = CollectionNameToCollectionCLSID(szName, &clsid)) return hr;
	return CParameterMap::ProgIDFromCLSID(clsid, sProgID);
}
//
//	Creates an instance of a singular OR COLLECTION object.
//	We do NOT insert any created object into any collection.
//
HRESULT CObjectManager::CreateObject(const CLSID& clsid, CComPtr<IDispatch>& spObject) const
{		
	if (spObject.CoCreateInstance(clsid)){
		CComBSTR sName;
		CParameterMap::ProgIDFromCLSID(clsid, sName);		
		return CParameterMap::ReturnErrorRSR(IDS_CREATE_OBJECT, estring(sName), IDS_DLL_NOT_REGISTERED);
	}
	return S_OK;
}
HRESULT	CObjectManager::CreateObject(const std::string& szObject, CComPtr<IDispatch>& spObject) const
{
	CLSID clsid;
	if (NameToCLSID(szObject, &clsid) && CollectionNameToCollectionCLSID(szObject, &clsid)) return CParameterMap::ReturnErrorRSR(IDS_CREATE_OBJECT, szObject, IDS_DLL_NOT_REGISTERED);
	return CreateObject(clsid, spObject);
}
//
//	Returns the 'Name' component of the m_coll variable in the property
//	szCollectionName of the object associated with spCollection corresponding
//	to an input spIndexObject
//
HRESULT CObjectManager::GetCollectionIndex(CComPtr<IDispatch> spObjectContainingCollection, CComPtr<IDispatch> spIndexObject, const estring& szCollectionName, std::string* pszIndex) const
{
	HRESULT								hr;
	CComVariant							vCollectionObject;
	CComVariant							vIndex;

	if (!spObjectContainingCollection) return E_FAIL;
	if (CComDispatchDriverEx(spObjectContainingCollection).GetPropertyByName((CComBSTR)szCollectionName, &vCollectionObject)) return CParameterMap::ReturnErrorRS(IDS_PROPERTY_NOT_FOUND, szCollectionName);
	if (vCollectionObject.vt != VT_DISPATCH) return E_FAIL;
	if (hr = CComDispatchDriverEx(vCollectionObject.pdispVal).GetPropertyByName(L"Index", &CComVariant(spIndexObject.p), 1, &vIndex)) return hr;
	pszIndex->assign(estring(vIndex));
	return pszIndex->size() ? S_OK : E_FAIL;
}
//	Manufactures a file name from an input object, data source and date etc.
HRESULT CObjectManager::GetFileName(const IID& iid, const std::string& szLocation, std::string szName, DataSourceEnum ds, long nDate, std::string* pszName) const
{
	HRESULT								hr;
	std::string							szObjectName;					// name of the object that we are saving	
	std::stringstream					ssFileName;
	std::string							szFileRoot = _Module.GetFileSystemRoot();	

	// manufacture the file name
	if (hr = IIDToCollectionName(iid, &szObjectName)) return hr;
	ssFileName << szFileRoot << '\\';
	if (szLocation.size()) ssFileName << szLocation << '\\';
	if (nDate) ssFileName << MlEqDate(nDate).GetString() << '\\';	
	ssFileName << szObjectName << '\\';
	if (ds != NoDataSource) ssFileName << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << '\\';
	if (szName.size()){
		bool								bWildcardAtEnd = false;			// true if szName ends with '*'

		if (estring::right(&szName, 1) == "*"){
			bWildcardAtEnd = true;
			szName.erase(szName.end() - 1);
		}

		estring::ReplaceStrInStr(&szName, ">", "&gt;");
		estring::ReplaceStrInStr(&szName, "*", "&star;");
		estring::ReplaceStrInStr(&szName, ":", "&colon;");
		estring::ReplaceStrInStr(&szName, "?", "&question;");
		
		ssFileName << szName << (bWildcardAtEnd ? "*" : "") << ".xml";
	}
	pszName->assign(ssFileName.str());
	return S_OK;
}
//
//	Returns an object from an appropriate collection (a maintained collection if the input
//  identifier is an index).
//
HRESULT CObjectManager::GetObject(const CComVariant& vIn, const IID& IID_ISingular, CComPtr<IDispatch>& spSingular) const
{
	HRESULT								hr;
	CParameterMap						pmIn;
	std::string							szIn;
	IID									iid = IID_ISingular;
	CLSID								clsidSingular = CLSID_NULL;
	CComVariant							vObject;
					
	if (hr = pmIn.SetValue(vIn)) return hr;		// we must do this process the Excel range case

	// Case 1: vIn is a handle
	if (!pmIn.GetString(&szIn)){
		if (iid == IID_IDispatch) HandleToObjectIID(IID_ISingular, szIn, &iid);
		if (iid != IID_IDispatch){
			CComPtr<IDispatch>			spCollection;
			std::string					sz = szIn;
			if (!IIDToCLSID(iid, &clsidSingular) && !IIDToCollectionPtr(iid, spCollection)){															
				CParameterMap::RemoveCalculationNumberFromHandle(clsidSingular, &sz);
				if (!CComDispatchDriverEx(spCollection.p).GetPropertyByName(L"Item", &CComVariant(sz.c_str()), 1, &vObject)){
					ATLASSERT(vObject.vt == VT_DISPATCH);
					spSingular = vObject.pdispVal;
					return S_OK;
				}
				// Another possibility is that the user has entered a handle of the form [ObjectName]::[ObjectInstanceName] (e.g. Asset::.FTSE).
				// We allow for this as it is an extremely powerful notation.
				std::vector<std::string>	vector;
				estring::Split(sz, "::", &vector);
				if (vector.size() == 2){
					IID		iidHandle = IID_NULL;
					NameToIID(vector[0], &iidHandle);
					if (iidHandle == iid){
						if (!CComDispatchDriverEx(spCollection.p).GetPropertyByName(L"Item", &CComVariant(vector[1].c_str()), 1, &vObject)){
							ATLASSERT(vObject.vt == VT_DISPATCH);
							spSingular = vObject.pdispVal;
							return S_OK;
						}
					}
				}
			}
		}
	}

	// Case 2: vIn is a dispatch pointer to a Sirius object	
	if (!pmIn.GetObject(&vObject)){	
		if (!CParameterMap::GetObjectIID(vObject.pdispVal, &iid) && IsIIDSiriusSingular(iid)){
			ATLASSERT(vObject.vt == VT_DISPATCH);
			spSingular = vObject.pdispVal;
			return S_OK;
		}		
	}		

	// Case 3: vIn, and therefore pmIn, are the data for the Sirius object IID_ISingular
	if (IID_ISingular != IID_IDispatch){
		CComPtr<IDispatch>		spObject;
		CComVariant				v;			// value of pmIn
		if (clsidSingular == CLSID_NULL) IIDToCLSID(IID_ISingular, &clsidSingular);
		if (clsidSingular != CLSID_NULL && !spObject.CoCreateInstance(clsidSingular) && !pmIn.GetValue(&v)){			
			if (!CComDispatchDriverEx(spObject).PutPropertyByName(L"Value", &v)){
				spSingular = spObject;
				return S_OK;				
			}
		}
	}

	// We were unsuccessful if this point is reached
	return CParameterMap::ReturnErrorRV(IDS_NO_OBJECT_WITH_HANDLE, vIn);
}
HRESULT CObjectManager::GetObject(const CComVariant& vIn, std::string& szObjectName, CComPtr<IDispatch>& spSingular) const
{
	HRESULT								hr;
	IID									iid;	

	if (!szObjectName.size()) return GetObject(vIn, spSingular);
	if (hr = NameToIID(szObjectName, &iid)) return hr;
	return GetObject(vIn, iid, spSingular);	
}
//
//	return an array of installed object names and their GUIDs for a certain
//	object type
//
HRESULT CObjectManager::GetObjectListOfType(long nTypes, VARIANT* pVal) const
//	nTypes - total of the various object::type enumerations specifying which objects we display
{
	SAFEARRAY*							psa;
	long								nIndex[2];						// element in the variant array
	SAFEARRAYBOUND						sab[2];							// safe array size specification vector
	CComVariant							v;
	CComBSTR							b;
	unsigned long						cElements = 0;
	std::map<std::string, CSiriusSingularObject>::const_iterator itr;
	
	// get the number of elements	
	for (itr = m_mapNameToObject.begin(); itr != m_mapNameToObject.end(); ++itr){
		if (!(itr->second.GetObjectType() & nTypes)) continue;
		cElements++;
	}	
	
	// create the safe array
	sab[0].lLbound = 1;	
	sab[0].cElements = cElements;
	sab[1].lLbound = 1;
	sab[1].cElements = 2;
	psa = ::SafeArrayCreate(VT_BSTR, 2, sab);
	if (!psa) return E_FAIL;			// ToDo - better error
			
	// populate the safe array		
	nIndex[0] = 0;
	for (itr = m_mapNameToObject.begin(); itr != m_mapNameToObject.end(); ++itr){
		if (!(itr->second.GetObjectType() & nTypes)) continue;
		nIndex[0]++;
		nIndex[1] = 1;
		b = estring(itr->second.GetProperName());
		::SafeArrayPutElement(psa, nIndex, b);		
		nIndex[1] = 2;		
		b = estring(itr->second.GetCLSID());
		::SafeArrayPutElement(psa, nIndex, b);		
	}
		
	// ToDo - newer versions of ATL have a CComVariant constructor / assignment operator that takes a SAFEARRAY pointer
	v.vt = VT_ARRAY | VT_BSTR;
	v.parray = psa;
	return v.Detach(pVal);
}
//
//	Returns a list of parameter names that the product type implementation
//  object expects. We do not return any indexed properties.
//
HRESULT CObjectManager::GetObjectProperties(INVOKEKIND inv, CComPtr<IDispatch> spObject, VARIANT* pVal, CParameterMap* ppm) const
//	inv - specifies which type(s) we include
//	spObject - name of object whose properties are required
//	pVal - 2D variant array returned (nullable)
//	ppm - parameter map form of pVal (returned, nullable)
{
	HRESULT								hr;
	CComPtr<ITypeInfo>					pti;
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc;
	MEMBERID							memid;	
	CParameterMap						pmRet;
	long								nRow = 0;
	
	if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;
	if (pti->GetTypeAttr(&pTypeAttr)) return hr;	
	if (hr = pmRet.SetSize(pTypeAttr->cFuncs, 1)) return hr;			// this is the maximum possible number of rows
	for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){
		if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){
			if (pFuncDesc->invkind & inv){				
				if (pFuncDesc->invkind == INVOKE_PROPERTYGET && pFuncDesc->cParams != 0 || pFuncDesc->invkind == INVOKE_PROPERTYPUT && pFuncDesc->cParams != 1) continue;	// we omit parameterised properties
				if (pFuncDesc->wFuncFlags & FUNCFLAG_FHIDDEN) continue;
				memid = pFuncDesc->memid;
				// get the method name
				CComBSTR sProperty;
				if (!(pti->GetDocumentation(memid, &sProperty, NULL, NULL, NULL))){
					pmRet.SetValue(nRow++, 0, sProperty);
				}
			}
			pti->ReleaseFuncDesc(pFuncDesc);
		}	
	}	
	pti->ReleaseTypeAttr(pTypeAttr);
	pmRet.SetRows(nRow);
	if (ppm && (hr = ppm->Attach(&pmRet))) return hr;
	if (pVal && (hr = pmRet.GetValue(pVal))) return hr;
	return S_OK;
}
HRESULT CObjectManager::GetObjectProperties(INVOKEKIND inv, const std::string& szObject, VARIANT* pVal, CParameterMap* ppm) const
//	pVal - 2D variant array returned (nullable)
//	ppm - parameter map form of pVal (returned, nullable)
{
	HRESULT								hr;	
	CComPtr<IDispatch>					spObject;						// product type model created
	
	if (!szObject.size()) return CParameterMap::ReturnErrorRR(IDS_DATA_ELEMENT_NOT_FOUND, IDS_HEADER_PRODUCT_TYPE, IID_ISiriusApplication);
	if (hr = CreateObject(szObject, spObject)) return hr;
	return GetObjectProperties(inv, spObject, pVal, ppm);
}
//
//	Returns the IID of an object associated with an input handle.
//
HRESULT CObjectManager::HandleToObjectIID(const IID& iidSingular, const CComVariant& vHandle, IID* piid) const
{
	if (vHandle.vt == VT_DISPATCH){
		return CParameterMap::GetObjectIID(CComPtr<IDispatch>(vHandle.pdispVal), piid);
	}	
	return HandleToObjectIID(iidSingular, estring(vHandle), piid);
}

HRESULT CObjectManager::HandleToObjectIID(const IID& iidSingular, std::string szHandle, IID* piid) const
{
	std::vector<std::string>	vector;
	CLSID						clsidSingular = CLSID_NULL;
	HRESULT						hr;

	if (iidSingular != IID_IDispatch){
		if (hr = IIDToCLSID(iidSingular, &clsidSingular)) return hr;	
	}
	CParameterMap::RemoveCalculationNumberFromHandle(clsidSingular, &szHandle);
	// from here, the handle must be of the form a::b
	estring::Split(szHandle, "::", &vector);
	if (vector.size() != 2) return E_FAIL;	
	return NameToIID(vector[0], piid);
}
//
//	Maps a singular object IID to its corresponding CLSID.
//
HRESULT CObjectManager::IIDToCLSID(const IID& iid, CLSID* pclsid) const
{			
	std::map<GuidEx, GuidEx>::const_iterator it;	
	if ((it = m_mapIIDToCLSID.find(iid)) == m_mapIIDToCLSID.end()) return E_FAIL;
	*pclsid = it->second;
	return S_OK;
}
//
//	Maps a singular object IID to the name of its corresponding collection.
//
HRESULT CObjectManager::IIDToCollectionName(const IID& iid, std::string* pszName) const
{
	HRESULT								hr;
	CComPtr<IDispatch>					sp;

	if (hr = IIDToCollectionPtr(iid, sp)) return hr;
	return CParameterMap::GetObjectName(sp, pszName);
}
//
//	Maps a singular object IID to the pointer associated with the corresponding
//	collection.
//
HRESULT	CObjectManager::IIDToCollectionPtr(const IID& iid, CComPtr<IDispatch>& sp) const
{
	std::map<GuidEx, CComPtr<IDispatch> >::const_iterator it;
	if ((it = m_mapIIDToCollectionPtr.find(iid)) == m_mapIIDToCollectionPtr.end()){
		// this should never fail
		ATLASSERT(false);
		return E_FAIL;
	}
	sp = it->second;
	return S_OK;
}
//
//	Maps a singular object IID to the GUID of the library in which the
//	object is implemented.
//
HRESULT CObjectManager::IIDToLIBID(const IID& iid, GUID* plibid) const
{
	ATLASSERT(false);			// ToDo - check
	std::map<GuidEx, GuidEx>::const_iterator it;
	if ((it = m_mapIIDToLIBID.find(iid)) == m_mapIIDToLIBID.end()) return E_FAIL;	
	*plibid = it->second;
	return S_OK;
}
//
//	Maps a singular object IID to its name.
//
HRESULT CObjectManager::IIDToName(const IID& iid, std::string* pszName) const
{					
	ATLASSERT(false);		// ToDo - check
	std::map<GuidEx, std::string>::const_iterator it;
	if ((it = m_mapIIDToName.find(iid)) == m_mapIIDToName.end()) return E_FAIL;
	pszName->assign(it->second);
	return S_OK;
}
//
//	Populates a variant with the available properties of an object (other than
//	the Value property). The variant always effectively contains only one
//	parameter map.
//
HRESULT CObjectManager::ImplementGetValue(CComPtr<IDispatch> spObject, VARIANT* pVal) const
{
	HRESULT								hr;
	CComPtr<ITypeInfo>					pti;
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc;
	MEMBERID							memid;	
	CParameterMap						pmRet;
	long								nRow = 0;
	
	if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;
	if (pti->GetTypeAttr(&pTypeAttr)) return hr;	
	if (hr = pmRet.SetSize(16, 2)) return hr;
	for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){
		if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){	
			if (pFuncDesc->invkind == DISPATCH_PROPERTYPUT && !(pFuncDesc->wFuncFlags & FUNCFLAG_FHIDDEN)){
				// Having DISPATCH_PROPERTYPUT here instead of DISPATCH_PROPERTYGET is not a mistake.
				// We want to maintain the symmetry between ImplementGetValue and ImplementPutValue.
				memid = pFuncDesc->memid;
				// get the method name
				CComBSTR sProperty;
				if (!(pti->GetDocumentation(memid, &sProperty, NULL, NULL, NULL))){
					if (CStringResource(IDS_HEADER_VALUE).CompareNoCase(sProperty) && pFuncDesc->cParams == 1){
						// We omit the Value property and all properties with parameters.
						CComVariant v;
						CComDispatchDriverEx(spObject).GetPropertyByName(sProperty, &v);
						if (nRow >= pmRet.GetRows()){							
							// add more rows
							pmRet.SetRows(16 + nRow * 2);	// the 16 is there in case nRow is zero (it could be if pm.GetRows() below evaluates to zero)
						}						
						pmRet.SetValue(nRow, 0, sProperty);
						VARTYPE vt = pFuncDesc->lprgelemdescParam->tdesc.vt;												
						if (vt == VT_USERDEFINED && v.vt == VT_I4){
							HREFTYPE				refType = pFuncDesc->lprgelemdescParam->tdesc.hreftype;
							CComPtr<ITypeInfo>		spTypeInfoRef;
							HRESULT					hrLocal = pti->GetRefTypeInfo(refType, &spTypeInfoRef);
							CComVariant				vPropertyValue;
							if (!hrLocal){								
								std::string		szValue;
								if (!(hrLocal = CEnumMap::GetString(spTypeInfoRef, LIBID_Sirius, v.lVal, &szValue))){
									pmRet.SetValue(nRow++, 1, szValue);
								}								
							}
							if (hrLocal) pmRet.SetValue(nRow++, 1, v);							
						} else if (v.vt == VT_DISPATCH){
							if (!v.pdispVal){
								pmRet.SetValue(nRow++, 1, v);
							} else if (IsSiriusCollection(CComPtr<IDispatch>(v.pdispVal))){
								// If v.pdispVal is a collection then we add in a set of object handles and notionals.
								CComVariant vCollection;
								CParameterMap pmCollection;
								if (CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Value", &vCollection) || pmCollection.SetValue(vCollection)){
									pmRet.SetValue(nRow++, 1, v);
								} else {
									// add pmCollection to pmRet
									pmRet.RemoveRow(nRow);
									pmRet.RemoveBlanksAtEnd();
									pmRet.AddToEnd(pmCollection);	
									nRow += pmCollection.GetRows();
								}
							} else {																																
								pmRet.SetValue(nRow++, 1, v);
							}
						} else {
							pmRet.SetValue(nRow++, 1, v);
						}
					}
				}
			}
			pti->ReleaseFuncDesc(pFuncDesc);
		}	
	}	
	pti->ReleaseTypeAttr(pTypeAttr);
	pmRet.RemoveBlanksAtEnd();
	return pmRet.GetValue(pVal);
}
//
//	Returns a value for a property exposed by an object. We convert enumerations
//	to strings and objects to handle names (if possible). We also decompose
//	collections into their constituent handles and notionals.
//
HRESULT CObjectManager::ImplementGetValue(CComPtr<IDispatch> spObject, const estring& szProperty, VARIANT* pVal) const
{		
	CComVariant							v;
	HRESULT								hr;	
		
	if (hr = CComDispatchDriverEx(spObject).GetPropertyByName(szProperty.GetBSTR(), &v)) return hr;
	if (v.vt == VT_DISPATCH){
		// first try mapping to an object index
		estring szIndex;		
		if (GetIndex(v.pdispVal, &szIndex, false)){
			return szIndex.GetValue(pVal);
		}		
		// now try the collection case
		if (IsSiriusCollection(CComPtr<IDispatch>(v.pdispVal))){
			// return a set of object handles and notionals.			
			return CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Value", pVal);			
		}		
	} else if (v.vt == VT_I4){
		// we may have an enumeration case that we need to decompose into its corresponding string
		CComPtr<ITypeInfo>				spti;		
		FUNCDESC*						pFuncDesc = NULL;		
		TYPEATTR*						pTypeAttr = NULL;					
		if (!spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &spti)){
			spti->GetTypeAttr(&pTypeAttr);
			for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){
				if (!(spti->GetFuncDesc(nIndex, &pFuncDesc))){
					if (pFuncDesc->invkind == DISPATCH_PROPERTYGET){
						CComBSTR		sProperty;
						if (!(spti->GetDocumentation(pFuncDesc->memid, &sProperty, NULL, NULL, NULL))){
							if (!szProperty.CompareNoCase(sProperty)){
								if (pFuncDesc->elemdescFunc.tdesc.vt == VT_USERDEFINED){										
									CComPtr<ITypeInfo>				spTypeInfoRef;										
									std::map<long, std::string>		map;										
									if (!spti->GetRefTypeInfo(pFuncDesc->elemdescFunc.tdesc.hreftype, &spTypeInfoRef)){										
										estring					szValue;
										if (!CEnumMap::GetString(spTypeInfoRef, LIBID_Sirius, v.lVal, &szValue)){
											szValue.GetValue(&v);												
										}										
									}
								}
							}
						}
					}
					spti->ReleaseFuncDesc(pFuncDesc);
				}
			}								
			spti->ReleaseTypeAttr(pTypeAttr);
		}
	}
	return v.Detach(pVal);
}
//
//	These nice little functions assume that an object is entirely representable by
//  a set of put property functions.
//
HRESULT CObjectManager::ImplementPutValue(CComPtr<IDispatch> spObject, const std::string& szProperty, const CComVariant& vValue, CComPtr<ITypeInfo> pti, FUNCDESC* pFuncDesc, MEMBERID memid)
{
	VARTYPE								vt = pFuncDesc->lprgelemdescParam->tdesc.vt;
	HRESULT								hr = S_OK;
	CComDispatchDriverEx				ddObject(spObject);

	if (vt == VT_PTR){
		// Object is expected. We need to determine the object needed, and
		// see if the user has passed in a handle.
		ELEMDESC* ped = (ELEMDESC*)pFuncDesc->lprgelemdescParam->tdesc.lptdesc;
		if (ped->tdesc.vt == VT_USERDEFINED){
			CComPtr<ITypeInfo>	pti_ref;									
			if (!(hr = pti->GetRefTypeInfo(ped->tdesc.hreftype, &pti_ref))){
				TYPEATTR*	pTypeAttrRef;
				if (!(hr = pti_ref->GetTypeAttr(&pTypeAttrRef))){
					GUID guid = pTypeAttrRef->guid;
					// get / create the appropriate object,
					// wrap it in a variant and and put the property.
					CComPtr<IDispatch> sp;
					if (vValue.vt == VT_EMPTY || vValue.vt == VT_BSTR && !::SysStringLen(vValue.bstrVal)){
						// set to null pointer																								
						hr = ddObject.PutProperty(memid, &CComVariant(CComPtr<IDispatch>()));
					} else if (!(hr = GetObject(vValue, guid, sp))){
						// If sp is a basket asset and vValue is a string that does not contain an explicit pay currency then we error.
						// We do this since there is ambigutiy in the assumed pay currency of a basket asset which does not explicitly
						// define the pay currency.												
						if (guid == IID_IAsset && sp && vValue.vt == VT_BSTR){							
							estring szKey(vValue.bstrVal);
							bool bSpreadsheetHandle;
							CParameterMap::RemoveCalculationNumberFromHandle(CLSID_Asset, &szKey, &bSpreadsheetHandle);
							if (!bSpreadsheetHandle){
								// i.e. szKey is an asset key.
								if (szKey.find(">") == std::string::npos){									
									CComPtr<IAsset> spAsset = dynamic_cast<IAsset*>(sp.p);
									AssetTypeEnum at = InvalidAsset;
									spAsset->get_AssetType(&at);
									if (at == QuantoBasket || at == CompositeBasket){
										throw "You need to supply an explicit pay currency for '" + szKey + "' since it is a basket asset";
									}
								}
							}
						}
						if (hr = ddObject.PutProperty(memid, &CComVariant(sp))){
							std::string szObjectName;
							CParameterMap::GetObjectName(spObject, &szObjectName);
							estring szError;
							szError.SetError(hr);
							CParameterMap::ReturnErrorS("Error '" + szError + "' encountered when assigning a value to the property '" + szProperty + "' of object '" + szObjectName + "'");
						}

					}
					pti_ref->ReleaseTypeAttr(pTypeAttrRef);
				}									
			}
			return hr;
		} else {
			// error - ToDo - better description
			return E_FAIL;
		}		
	} else if (vt & VT_ARRAY){								
		// not supported
		ATLASSERT(false);
		return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);		
	} else if (vt == VT_USERDEFINED){
		// user defined type with a string input																
		HREFTYPE						refType = pFuncDesc->lprgelemdescParam->tdesc.hreftype;
		CComPtr<ITypeInfo>				spTypeInfoRef;		
		if (hr = pti->GetRefTypeInfo(refType, &spTypeInfoRef)) return hr;			
		std::map<std::string, long>		map;
		std::string						sz;
		long							n = 0;
		CParameterMap					pmElement;
		if (!pmElement.SetValue(vValue) && pmElement.GetValue(&n)){
			// assume the string form of the enumerator is given
			pmElement.GetString(&sz);
			if (!sz.size()){
				// here we are assuming that the null string always corresponds to the zero enumerator value
				n = 0;
			} else if (hr = CEnumMap::GetEnum(spTypeInfoRef, LIBID_Sirius, sz, &n)){
				return CParameterMap::ReturnErrorSRS(sz, IDS_MSG_INVALID_ENUMERATOR, szProperty);
			}
		}
		return ddObject.PutProperty(memid, &CComVariant(n));
	} else {
		// straightforward value is expected
		if (!(hr = ddObject.PutProperty(memid, const_cast<CComVariant*>(&vValue)))){
			return S_OK;
		} else {
			std::string						szValue;
			CParameterMap					pmElement;			
			pmElement.SetValue(vValue);
			pmElement.GetString(&szValue);
			if (szValue.size()){				
				return CParameterMap::ReturnErrorSRS(szValue, IDS_MSG_INVALID_PARAMETER, szProperty);
			} else {
				return CParameterMap::ReturnErrorRS(IDS_INVALID_VALUE_FOR_PARAMETER, szProperty);
			}
		}
	}	
}
HRESULT CObjectManager::ImplementPutValue(const VARIANT& vIn, CComPtr<IDispatch> spObject)
{	
	HRESULT								hr;
	std::vector<CParameterMap>			vpm;
	
	if (hr = CParameterMap::ArrayToVector(vIn, NULL, &vpm)) return hr;
	hr = ImplementPutValue(vpm, spObject);
	return hr;
}
HRESULT CObjectManager::ImplementPutValue(std::vector<CParameterMap>& vpm, CComPtr<IDispatch> spObject)
{		
	HRESULT								hr;	
	long								nRowUsed;						// row of the parameter map		
	std::vector<bool>					vectorRowsUsed;					// list of parameter map rows used			
	IID									iid;							// this is the IID of spObject (we use this only for error handling)
	CComPtr<ITypeInfo>					pti;
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc;
	MEMBERID							memid;
	int									nRowsUsed = 0;					// number of rows used in the parameter maps
	CComDispatchDriverEx				ddObject(spObject);
	CParameterMap						pmIndex;						// collection handle indexes (populated if bHasCollectionIndexExtension is true).
	bool								bHasCollectionIndexExtension = false;
				
	if (hr = CParameterMap::GetObjectIID(spObject, &iid)) return hr;
	if (vpm.size() == 1){
		if (vpm[0].GetCols() != 2 && vpm[0].GetCols() != 3) return CParameterMap::ReturnErrorR(IDS_COLUMNS_INVALID, iid);
		CParameterMap pmCol;
		vpm[0].GetColumn(1, &pmCol, false);		// i.e. we do not remove blanks from the column - this supports blank values at the top and bottom of the original matrix.
		vpm[0].RemoveColumn(1);
		vpm.push_back(pmCol);
		if (!vpm[0].IsVector()){
			// attempt to extract pmIndex
			if (vpm[0].GetCols() == 2){
				vpm[0].GetColumn(1, &pmIndex, false);	// i.e. we do not remove blanks from the column - we only expect values for collection properties
				vpm[0].RemoveColumn(1);
				bHasCollectionIndexExtension = true;
			}
		}
	}
	
	if (vpm.size() != 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, iid);
	if (!vpm[0].IsVector()) return CParameterMap::ReturnErrorR(IDS_VECTOR_NAMES_INVALID);
	if (!vpm[1].IsVector()) return CParameterMap::ReturnErrorR(IDS_VECTOR_VALUES_INVALID);
	if (vpm[0].GetCols() != 1 && vpm[0].Transpose()) ATLASSERT(false);
	if (vpm[1].GetCols() != 1 && vpm[1].Transpose()) ATLASSERT(false);
	if (vpm[0].GetRows() < vpm[1].GetRows()) return CParameterMap::ReturnErrorR(IDS_VECTOR_VALUES_INVALID);
	vpm[1].SetRows(vpm[0].GetRows());
	vectorRowsUsed.resize(vpm[0].GetRows(), false);	
	
	// For each writable property exposed by the object, we need to attempt to set its value from the input.
	// We check that every parameter in the input has been used.
	if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;
	if (!(pti->GetTypeAttr(&pTypeAttr))){
		for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){
			if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){				
				if (pFuncDesc->invkind == DISPATCH_PROPERTYPUT && !(pFuncDesc->wFuncFlags & FUNCFLAG_FHIDDEN)){
					// only process if it's a non-hidden put property function
					memid = pFuncDesc->memid;
					// get the method name
					CComBSTR sProperty;
					if (!(pti->GetDocumentation(memid, &sProperty, NULL, NULL, NULL))){
						// get the parameter map row with this method name
						CComVariant vElement;
						estring		szProperty(sProperty);						
						if (!vpm[0].FindCell(szProperty, true, CParameterMap::LeftColumnOnly, &nRowUsed, NULL) && !vpm[1].GetValue(nRowUsed, 0, &vElement)){
							nRowsUsed++;
							vectorRowsUsed[nRowUsed] = true;							
							hr = ImplementPutValue(spObject, szProperty, vElement, pti, pFuncDesc, memid);
						} else {
							// szProperty could not be found in the parameter map. However, szProperty may denote a collection.
							// In which case, there may be elements in the parameter map that we can set to the the collection.
							if (pFuncDesc->lprgelemdescParam->tdesc.vt == VT_PTR){
								ELEMDESC* ped = (ELEMDESC*)pFuncDesc->lprgelemdescParam->tdesc.lptdesc;
								if (ped->tdesc.vt == VT_USERDEFINED){
									CComPtr<ITypeInfo>	pti_ref;
									if (!(hr = pti->GetRefTypeInfo(ped->tdesc.hreftype, &pti_ref))){
										TYPEATTR*	pTypeAttrRef;
										if (!(hr = pti_ref->GetTypeAttr(&pTypeAttrRef))){																					
											IID			iidSingular;											
											if (!CollectionIIDToIID(pTypeAttrRef->guid, &iidSingular)){
												// szProperty is a collection. The operation must be successful from here.
												CLSID		clsidCollection;												
												if (!(hr = CollectionIIDToCollectionCLSID(pTypeAttrRef->guid, &clsidCollection))){													
													CComPtr<IDispatch> spCollection;
													if (!(hr = spCollection.CoCreateInstance(clsidCollection))){
														// add any elements in the parameter map that are objects of type iidSingular should be added to spCollection
														for (long nRow = 0, nIndex = 0; nRow < vectorRowsUsed.size(); ++nRow){
															if (vectorRowsUsed[nRow]) continue;
															CComVariant v;
															if (!vpm[0].GetValue(nRow, 0, &v)){
																// We only proceed from here if v maps to an iid equal to iidSingular.
																IID iidThis;
																if (!HandleToObjectIID(iidSingular, v, &iidThis) && iidThis == iidSingular){
																	CComPtr<IDispatch> spObject;
																	if (!GetObject(v, iidSingular, spObject)){
																		double fNotional;
																		if (vpm[1].GetValue(nRow, 0, &fNotional)){
																			CComVariant vElement;
																			vpm[1].GetValue(nRow, 0, &vElement);
																			CParameterMap::ReturnErrorRV(IDS_INVALID_PARAMETER, vElement);
																			hr = E_FAIL;
																			break;
																		}
																		// Add the object to our collection.
																		estring szHandle;
																		if (!IsSerialisable(spObject)){
																			if (bHasCollectionIndexExtension){																				
																				pmIndex.GetValue(nRow, 0, &szHandle);
																				if (szHandle.islong()) szHandle.clear();	// Don't use in case it clashes with an nIndex value.
																			}
																			if (!szHandle.size()){
																				szHandle = ++nIndex;
																			}
																		}
																		if (hr = CComDispatchDriver(spCollection.p).Invoke2(DISPID_ADD, &szHandle.GetValue(), &CComVariant(spObject))) break;
																		CComVariant	params[] = {fNotional, spObject};
																		if (hr = CComDispatchDriverEx(spCollection.p).PutPropertyByName(L"Notional", params, 2)) break;
																		vectorRowsUsed[nRow] = true;
																		nRowsUsed++;
																	}
																}
															}
														}
														// put spCollection if all is well
														if (!hr) hr = ddObject.PutProperty(memid, &CComVariant(spCollection));
													}
												}
											}											
											pti_ref->ReleaseTypeAttr(pTypeAttrRef);
										}
									}
								}
							}
						}
					}
				}
				pti->ReleaseFuncDesc(pFuncDesc);	
			}
			if (hr) break;
		}
		pti->ReleaseTypeAttr(pTypeAttr);
	}
	if (hr) return hr;

	// return an error if there are unused lines in the input parameter maps
	if (vpm[0].GetRows() != nRowsUsed){
		// Find at least one of the unused parameters. Note that we allow blank parameter names if their values are also blank.
		std::string sz;
		for (long nRow = 0; nRow < vectorRowsUsed.size(); nRow++){
			if (!vectorRowsUsed[nRow]){
				// this row will do				
				vpm[0].GetValue(nRow, 0, &sz);				
				if (!sz.size()){
					// We error if the corresponding value is not blank
					if (!vpm[1].IsBlank(nRow, 0)){
						CComVariant v;
						vpm[1].GetValue(nRow, 0, &v);
						return CParameterMap::ReturnErrorRV(IDS_VALUE_NOT_USED, v);
					}
				} else {
					return CParameterMap::ReturnErrorRS(IDS_PARAMETER_NOT_USED, sz, iid);
				}
			}
		}		
	}	
	
	// all is well if this point is reached
	return S_OK;
}
HRESULT CObjectManager::ImplementPutValue(CComPtr<IDispatch> spObject, const std::string& szProperty, const CComVariant& Value)
{
	HRESULT								hr;
	CComPtr<ITypeInfo>					pti;	
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc = NULL;
	MEMBERID							memid;	
	bool								bFound = false;
	
	// If the property type is an enumerator then we need to map any enumerator string to the corresponding enumerator value.
	if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;
	if (pti->GetTypeAttr(&pTypeAttr)) return hr;
	for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs && !bFound; nIndex++){
		if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){
			if (pFuncDesc->invkind != INVOKE_PROPERTYPUT && pFuncDesc->cParams != 1) continue;	// we omit parameterised properties
			memid = pFuncDesc->memid;
			// get the method name
			CComBSTR sProperty;
			if (!(pti->GetDocumentation(memid, &sProperty, NULL, NULL, NULL))){
				if (!estring(sProperty).CompareNoCase(szProperty)){
					// found the function
					bFound = true;										
					hr = ImplementPutValue(spObject, szProperty, CComVariant(Value), pti, pFuncDesc, memid);
				}
			}
			pti->ReleaseFuncDesc(pFuncDesc);
		}	
	}	
	pti->ReleaseTypeAttr(pTypeAttr);
	if (!bFound) return CParameterMap::ReturnErrorRS(IDS_PROPERTY_NOT_FOUND, szProperty);	
	return hr;
}
//
//	Installs a Sirius object (i.e. one which exposes a handle or a product
//  type) into the object manager maps.
//
/*protected*/ HRESULT CObjectManager::Install(const CSiriusSingularObject& sso)
{
	m_mapNameToObject[sso.GetName()] = sso;				
	m_mapIIDToCLSID[sso.GetIID()] = GuidEx(sso.GetCLSID());
	m_mapIIDToCollectionPtr[sso.GetIID()] = sso.GetCollectionPtr();
	m_mapCollectionPtrToBool[sso.GetCollectionPtr()] = true;
	m_mapCLSIDToProductLike[sso.GetCLSID()] = sso.GetProductLike();	
	m_mapCLSIDToName[sso.GetCLSID()] = sso.GetName();
	m_mapCLSIDToProperName[sso.GetCLSID()] = sso.GetProperName();	
	m_mapIIDToName[sso.GetIID()] = sso.GetName();
	m_mapIIDToLIBID[sso.GetIID()] = sso.GetLIBID();				
	if (sso.GetCollectionIID() != IID_NULL){
		m_mapCollectionIIDToIID[sso.GetCollectionIID()] = sso.GetIID();
		m_mapCollectionIIDToCollectionCLSID[sso.GetCollectionIID()] = sso.GetCollectionCLSID();
		m_mapCollectionCLSIDToCollectionIID[sso.GetCollectionCLSID()] = sso.GetCollectionIID();
		ATLASSERT(sso.GetCollectionName().size());		
		m_mapCollectionNameToCollectionIID[sso.GetCollectionName_lc()] = sso.GetCollectionIID();
		m_mapCollectionNameToCollectionCLSID[sso.GetCollectionName_lc()] = sso.GetCollectionCLSID();
		m_mapCollectionCLSIDToCollectionName[sso.GetCollectionCLSID()] = sso.GetCollectionName();
	}
	return S_OK;
}

//
//	Installs a product type (e.g. Cliquet, Napoleon, Call) into the object manager.
//	
HRESULT CObjectManager::InstallProductType(const GUID& libid, const CLSID& clsid, const IID& iid)
{
	HRESULT						hr;	
	CSiriusSingularObject				sso;
	CComPtr<IDispatch>			spDummy;

	if (hr = sso.Install(libid, clsid, iid, CSiriusSingularObject::ProductType, true, "", spDummy, CLSID_NULL)) return hr;
	return Install(sso);
}
HRESULT CObjectManager::InstallProductType(const std::string& szLIBID, const std::string& szCLSID, const std::string& szIID, const std::string& szName)
{
	HRESULT						hr;
	GUID						libid;
	CLSID						clsid;
	IID							iid;	
	OLECHAR*					s;
	CSiriusSingularObject				sso;
	CComPtr<IDispatch>			spDummy;
	
	libid = GUID_NULL;
	clsid = CLSID_NULL;
	iid = IID_NULL;
	s = estring::AnsiToUnicode(szLIBID.c_str());
	IIDFromString(s, &libid);	
	delete s;	
	s = estring::AnsiToUnicode(szCLSID.c_str());
	CLSIDFromString(s, &clsid);
	delete s;
	s = estring::AnsiToUnicode(szIID.c_str());
	IIDFromString(s, &iid);
	delete s;		
	if (hr = sso.Install(libid, clsid, iid, CSiriusSingularObject::ProductType, true, szName, spDummy, CLSID_NULL)) return hr;
	return Install(sso);		
}
//
//	See IsSiriusCollection
//
bool CObjectManager::IsCLSIDSiriusCollection(const CLSID& clsid) const
{	
	std::map<GuidEx, GuidEx>::const_iterator it = m_mapCollectionCLSIDToCollectionIID.find(clsid);
	if (it != m_mapCollectionCLSIDToCollectionIID.end()) return true;
	return false;
}
//
//	See IsIIDSiriusColletion
//
bool CObjectManager::IsIIDSiriusCollection(const IID& iid) const
{
	std::map<GuidEx, GuidEx>::const_iterator it = m_mapCollectionIIDToCollectionCLSID.find(iid);
	return it != m_mapCollectionIIDToCollectionCLSID.end();
}

//
//	Returns true if the interface IID is the default interface for a 
//  sirius non-colletion object.
//
bool CObjectManager::IsIIDSiriusSingular(const IID& iid) const
{
	std::map<GuidEx, GuidEx>::const_iterator it = m_mapIIDToCLSID.find(iid);
	return it != m_mapIIDToCLSID.end();
}
//
//	Returns true if the input object pointer / object name is a Sirius
//	collection. NOTE THAT THIS DOESN'T MEAN THAT THE COLLECTION IS MAINTAINED.
//
bool CObjectManager::IsSiriusCollection(std::string sz) const
{
	estring::lc(&sz);	
	std::map<std::string, GuidEx>::const_iterator it = m_mapCollectionNameToCollectionIID.find(sz);
	if (it != m_mapCollectionNameToCollectionIID.end()) return true;
	// the input name could be of the form Library.Object.Version
	std::vector<std::string>		vector;
	estring::Split(sz, ".", &vector);	
	if (vector.size() != 3) return false;
	return IsSiriusCollection(vector[1]);
}
bool CObjectManager::IsSiriusCollection(CComPtr<IDispatch> sp) const
{
	IID									iid;
	if (CParameterMap::GetObjectIID(sp, &iid)) return false;
	return IsIIDSiriusCollection(iid);
}
//
//	Returns true if the input object is a Sirius market data object.
//
bool CObjectManager::IsSiriusMarketDataObject(std::string sz) const
{
	estring::lc(&sz);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(sz)) == m_mapNameToObject.end()) return false;
	return (it->second.GetObjectType() & CSiriusSingularObject::MarketDataComponent) ? true : false;
}
//
//	Returns true if the input string is the name of a product that is
//  accessible to the current user.
//
bool CObjectManager::IsSiriusProduct(std::string sz) const
{
	estring::lc(&sz);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(sz)) == m_mapNameToObject.end()) return false;
	return it->second.GetObjectType() & CSiriusSingularObject::ProductType ? true : false;
}
//
//	Returns true if the input object is a Sirius singular object.
//
bool CObjectManager::IsSiriusSingular(std::string sz) const
{
	estring::lc(&sz);
	std::map<std::string, CSiriusSingularObject>::const_iterator it = m_mapNameToObject.find(sz);
	return it != m_mapNameToObject.end();
}
//
//	Maps a singular object name to its corresponding CLSID.
//
HRESULT CObjectManager::NameToCLSID(std::string szName, CLSID* pclsid) const
{	
	estring::lc(&szName);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(szName)) == m_mapNameToObject.end()){
		// the input name could be of the form Library.Object.Version
		std::vector<std::string>		vector;
		estring::Split(szName, ".", &vector);	
		if (vector.size() != 3) return E_FAIL;
		return NameToCLSID(vector[1], pclsid);
	}
	
	*pclsid = it->second.GetCLSID();
	return S_OK;
}
//
//	Maps a singular object name to the name of its corresponding collection
//
HRESULT CObjectManager::NameToCollectionName(std::string szName, std::string* pszCollectionName) const
{
	estring::lc(&szName);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(szName)) == m_mapNameToObject.end()) return E_FAIL;	
	if (!it->second.GetCollectionName().size()) return E_FAIL;
	pszCollectionName->assign(it->second.GetCollectionName());
	return S_OK;
}
//
//	Maps a singular object name to the object pointer of the corresponding
//	collection.
//
HRESULT CObjectManager::NameToCollectionPtr(std::string szName, CComPtr<IDispatch>& sp) const
{
	estring::lc(&szName);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(szName)) == m_mapNameToObject.end()) return E_FAIL;
	sp = it->second.GetCollectionPtr();
	return S_OK;
}
//
//	Maps a singular object name to its corresponding IID.
//
HRESULT CObjectManager::NameToIID(const std::string& szName, IID* piid) const
{
	estring szName_lc(szName);
	estring::lc(&szName_lc);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(szName_lc)) == m_mapNameToObject.end()) return CParameterMap::ReturnErrorS("Invalid object name '" + szName + "'");
	*piid = it->second.GetIID();
	return S_OK;
}
//
//	Maps a singular object name to its corresponding ProgID.
//
HRESULT CObjectManager::NameToProgID(std::string szName, CComBSTR& sProgID) const
{		
	CLSID								clsid;	
	HRESULT								hr;
	if (hr = NameToCLSID(szName, &clsid)) return hr;
	return CParameterMap::ProgIDFromCLSID(clsid, sProgID);
}
//
//	Returns the proper case form of an input object name
//
std::string CObjectManager::NameToProperName(std::string szName) const
{		
	estring::lc(&szName);
	std::map<std::string, CSiriusSingularObject>::const_iterator it;
	if ((it = m_mapNameToObject.find(szName)) == m_mapNameToObject.end()) return szName;
	return it->second.GetProperName();
}
//
//	Examines all the puttable properties in an input object and if they
//  are Sirius serialisable objects, then they are replaced with objects
//  from the appropriate maintained collection (if available).
//
void CObjectManager::Refresh(CComPtr<IDispatch> spObject, bool bRecursive) const
{			
	long								nDate = 0L;
	DataSourceEnum						ds = NoDataSource;
	std::set<CComPtr<IDispatch> >       set;							// This is used as a recursion blocker.
		
	if (!IsSerialisable(spObject)){
		throw "Input object is not serialisable";
	} else {
		CComVariant						vDataSource, vDate;
		if (!CComDispatchDriverEx(spObject).GetPropertyByName(L"DataSource", &vDataSource)){
			ds = CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, vDataSource, NoDataSource);
		}
		if (!CComDispatchDriverEx(spObject).GetPropertyByName(L"Date", &vDate) && !vDate.ChangeType(VT_I4)){
			nDate = vDate.lVal;
		}
		if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
		if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	}
		
	Refresh(&set, spObject, bRecursive, ds, nDate);
}
/*protected*/ void CObjectManager::Refresh(std::set<CComPtr<IDispatch> >* pset, CComPtr<IDispatch> spObject, bool bRecursive, DataSourceEnum ds, long nDate) const
{
	CComPtr<ITypeInfo>					pti;
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc;
	std::string							szError;
		
	if (pset->find(spObject) != pset->end()){
		// We have already considered this object.
		// ToDo - on March 3, 2005 I decided that this was sufficient to block infinite recursion. 
		// Let's keep an eye on this for a year or so.
		return;				
	} else {
		pset->insert(spObject);
	}
			
	pin_date(ds, nDate);
	
	if (spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) throw "Unhandled exception in CObjectManager::Refresh";
	if (!(pti->GetTypeAttr(&pTypeAttr))){		
		for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){
			if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){
				if (pFuncDesc->invkind == DISPATCH_PROPERTYPUT && !(pFuncDesc->wFuncFlags & FUNCFLAG_FHIDDEN) && pFuncDesc->cParams == 1){	// We omit parameterised properties.
					if (pFuncDesc->lprgelemdescParam->tdesc.vt == VT_PTR){
						// This point is reached only if are considering a non-hidden non-parameterised put property 
						// function whose property type is an object.
						ELEMDESC* ped = (ELEMDESC*)pFuncDesc->lprgelemdescParam->tdesc.lptdesc;
						if (ped->tdesc.vt == VT_USERDEFINED){
							CComPtr<ITypeInfo>	pti_ref;									
							if (!pti->GetRefTypeInfo(ped->tdesc.hreftype, &pti_ref)){
								TYPEATTR* pTypeAttrRef;
								if (!pti_ref->GetTypeAttr(&pTypeAttrRef)){																		
									MEMBERID				memid = pFuncDesc->memid;
									CComVariant				vOldObject;
#									ifdef _DEBUG
										{	
											CComBSTR sProperty;
											std::string szObjectName;
											if (IsSerialisable(spObject)){
												szObjectName = CComObjectCollectionSerialisableKey(spObject.p);
											} else {
												CParameterMap::GetObjectName(spObject, &szObjectName);
											}
											pti->GetDocumentation(memid, &sProperty, NULL, NULL, NULL);
											ATLTRACE("CObjectManager::Refresh: Object name is '%s'. Property under consideration is '%s'.\n", szObjectName.c_str(), estring(sProperty).c_str());
										}
#									endif									
									CComDispatchDriverEx(spObject).GetProperty(memid, &vOldObject);
									if (vOldObject.vt == VT_DISPATCH && vOldObject.pdispVal){
										// Only consider replacing an object if there is already one defined.
										GUID				iid = pTypeAttrRef->guid;
										if (IsIIDSiriusSingular(iid)){
											CComPtr<IDispatch> spNewObject;											
											try {
												spNewObject = RefreshGetNewObject(iid, vOldObject, ds, nDate);
											} catch (const std::string& sz){
												szError.assign(sz);
											}
											if (!szError.size() && bRecursive){
												Refresh(pset, spNewObject, bRecursive, ds, nDate);
											}
											if (!szError.size() && spNewObject.p != vOldObject.pdispVal){
												if (CComDispatchDriverEx(spObject).PutProperty(memid, &CComVariant(spNewObject))) CParameterMap::ErrorHandler(_Module.GetSiriusApplication(), S_OK, &szError);
											}											
										} else if (IsIIDSiriusCollection(iid)){
											// Refresh every item in the collection.
											// This is nasty due to the fact that we can't use late binding (or we get large
											// switch statements) and we can't call Remove and Add on the collection to switch
											// objects since that would mess up the container map.
											IID										  iidSingular = IID_NULL;
											std::vector<CAdapt<CComPtr<IDispatch> >	> aspNewObjects;
											std::vector<CComVariant>			      avNotionals;
											CComDispatchDriverEx				      ddCollection(vOldObject.pdispVal);
											CComVariant							      vItems;
											
											CollectionIIDToIID(iid, &iidSingular);
											ddCollection.GetPropertyByName(L"Count", &vItems);
											if (vItems.vt == VT_I4){
												for (CComVariant vItem(1L); vItem.lVal <= vItems.lVal && !szError.size(); vItem.lVal++){
													CComVariant				vOldSingularObject, vNotional(1.0);
													CComPtr<IDispatch>		spNewSingularObject;
													ddCollection.GetPropertyByName(L"Item", &vItem, 1, &vOldSingularObject);
													ddCollection.GetPropertyByName(L"Notional", &vItem, 1, &vNotional);
													try {
														spNewSingularObject = RefreshGetNewObject(iidSingular, vOldSingularObject, ds, nDate);
													} catch (const std::string& sz){
														szError.assign(sz);
													}
													if (!szError.size()){
														aspNewObjects.push_back(spNewSingularObject);
														avNotionals.push_back(vNotional);
														if (bRecursive){
															Refresh(pset, spNewSingularObject, bRecursive, ds, nDate);
														}
													}														
												}
												if (!szError.size()){
													// We can now clear ddCollection and add the aspNewObjects set to it.
													ddCollection.Invoke0(L"Clear");																													
													for (long n = 0; n < aspNewObjects.size() && !szError.size(); n++){
														if (ddCollection.Invoke2(L"Add", &CComVariant(), &CComVariant(aspNewObjects[n].m_T))){
															CParameterMap::ErrorHandler(_Module.GetSiriusApplication(), S_OK, &szError);
														} else {
															CComVariant av[2] = {avNotionals[n], aspNewObjects[n].m_T};	
															if (ddCollection.PutPropertyByName(L"Notional", av, 2)){
																CParameterMap::ErrorHandler(_Module.GetSiriusApplication(), S_OK, &szError);
															}
														}
													}														
												}												
											}
										}
									} else {
										// We don't do anything for the null object case - I certainly never intend to
										// support reputting a null.										
									}
									pti_ref->ReleaseTypeAttr(pTypeAttrRef);
								}
							}
						}
					}
				}
				pti->ReleaseFuncDesc(pFuncDesc);	
			}
			if (szError.size()) break;				
		}
		pti->ReleaseTypeAttr(pTypeAttr);
	}
	if (szError.size()) throw szError;
		
	// Now set the date and data source on the object spObject to ds and nDate.
	// Intuitively, this seems unnecessary but some objects (e.g. assets) might 
	// perform processing when a date is put - even if it is unchanged.
#	ifdef _DEBUG
	{
		// The date and datasource on the object should not change as a result
		// of the PutPropertyByName functions.
		CComVariant vDate_Current, vDS_Current;
		CComDispatchDriverEx(spObject).GetPropertyByName(L"Date", &vDate_Current);
		CComDispatchDriverEx(spObject).GetPropertyByName(L"DataSource", &vDS_Current);
		if (vDate_Current.ChangeType(VT_I4) || vDS_Current.ChangeType(VT_I4)) ATLASSERT(false);
		ATLASSERT(vDate_Current.lVal == nDate && vDS_Current.lVal == ds);
	}
#	endif

	CComDispatchDriverEx(spObject).PutPropertyByName(L"Date", &CComVariant(nDate));
	CComDispatchDriverEx(spObject).PutPropertyByName(L"DataSource", &CComVariant(ds));
}

CComPtr<IDispatch> CObjectManager::RefreshGetNewObject(IID iid, const CComVariant& vObject, DataSourceEnum ds, long nDate) const
{
	CComVariant							vName;
	CComDispatchDriverEx(vObject.pdispVal).GetPropertyByName(L"Name", &vName);
	if (vName.ChangeType(VT_BSTR)) throw "Invalid name property value on an object encountered";
	
	CComObjectCollectionSerialisableKey key(vName.bstrVal, nDate, ds);											
	CComPtr<IDispatch>					spNewObject;
	if (GetObject(key, iid, spNewObject)) propagate_error;
	return spNewObject;
}
