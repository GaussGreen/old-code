//	ComObjectCollectionSerialisable.h : This implements the value property for a collection.
//
//  Author :							David Cuin
//
//	This class extends the collection concept to use the {Name, Date, DataSource}
//	as an alternative means of representing an object. I call this the Key. The
//	'proper' container map (m_coll) is still Handle -> Object as before
//
//  I have introduced two other maps Key -> Handle and Handle -> Key.
//	The user can therfore specify either the Handle or the Key to identify an object.
//  This greatly facilitates referencing these objects in VBA etc.
//
//	As far as this class is concerned, the Index can be regarded as the collective
//	name for a Handle or Key. (Remember that Index can also refer to a number denoting
//  an arbitrary position in a collection).

/////////////////////////////////////////////////////////////////////////////

#ifndef _COMOBJECTCOLLECTIONSERIALISABLE_H
#define _COMOBJECTCOLLECTIONSERIALISABLE_H

#include "comobjectcollectionserialisablekey.h"
#include "siriuscomobjectcollection.h"
#include "MlEqDate.h"


template <class ISingular, class CCollection, class ICollection, const CLSID* pclsidCollection, const IID* piidCollection, const GUID* plibidCollection, const CLSID* pclsidSingular, const UINT idr>
	class ATL_NO_VTABLE CComObjectCollectionSerialisable:
	public CSiriusComObjectCollection<ISingular, CCollection, ICollection, pclsidCollection, piidCollection, plibidCollection, pclsidSingular, idr>
{
protected:	
	typedef std::map<CAdapt<CComBSTR>, CAdapt<CComBSTR> >		HandleToKeyType;
	typedef std::map<CAdapt<CComBSTR>, CAdapt<CComBSTR> >		KeyToHandleType;
	HandleToKeyType												m_HandleToKey;
	KeyToHandleType												m_KeyToHandle;							// this is used for quick deletion of m_HandleToKey elements
	
	// Virtual Functions. Override these in base classes for specific functionality.
	virtual void CheckAddHandle(std::string* pszHandle)													{return;}							// Override to change the input handle value HandleIn.
	virtual HRESULT GetClone(const std::string& szKey, CComPtr<IDispatch>& spObject)					{ATLASSERT(false)/*You should override this if IsCloneable is overridden*/; return E_FAIL;}	// Override this if we can derive an object from an existing one - this function is only called in the maintained collection case.
	virtual bool IsClonable(void) const																	{return false;}						// If an object is clonable then return true here.	
	virtual void KeyAdded(const std::string& szKey)														{return;}							// Called when a key has been inserted into the maintained collection.
	virtual void KeyCleared(const std::string& szKey)													{return;}							// Called when a key has been cleared from the maintained collection.	
	virtual void KeyRekeyed(const std::string& szKeyOld, const std::string& szKeyNew)					{return;}							// Called when a key has been renamed in the maintained collection
	virtual void KeysAllCleared(void)																	{return;}							// Called when all the keys have been cleared from the maintained collection.
	virtual HRESULT LoadGDAObject(const std::string& szName, long nDate, CComPtr<ISingular>& spObject)	{return E_FAIL;}					// Override this if an object can be obtained from a GDA instance.
	virtual bool ShouldAddToMaintainedCollection(CComPtr<IDispatch> spObject)						    {return true;}						// Override this if we need to prevent certain objects from being added to the maintained collections (e.g. cloned assets). We cannot support this for the general collection or Flatten() breaks.
			
	// Override this to provide specific functionality for adding a key, handle pair.
	// Return S_FALSE if the special add is unnecessary, S_OK if the add was processed
	// or a fail value if the special add was considered and failed.	
	virtual HRESULT SpecialAdd(CComPtr<ISingular> spSingular, const CComBSTR& sKey)						{return S_FALSE;}
	
public:	
	//	Add a new item into the collection.
	STDMETHOD(Add)(VARIANT HandleIn, ISingular* Value)
	{				
		begin_function
		HRESULT							hr;
		CComPtr<ISingular>				spValue(Value);
		CComDispatchDriverEx			ddValue(Value);				
		CComBSTR						sHandle;
		estring							szHandle;
		CComBSTR						sKey;
		estring							szKey;
		std::string						szName;							// name of the object being added
		long							nDate;							// date of the object being added
		DataSourceEnum					ds;								// data source of the object being added
		bool							bIsSpreadsheetHandle;
		bool							bIsValidHandle;
		bool							bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICollection>(this));
		
		// Get the handles, sHandle or szHandle.
		{
			CComVariant						Handle(HandleIn);			// this scoping is deliberate
			if (Handle.vt == VT_ERROR && Handle.scode == DISP_E_PARAMNOTFOUND){
				Handle.Attach(&CComVariant());
			}
			if (Handle.vt == VT_EMPTY){
				Handle.ChangeType(VT_BSTR);
			}
			if (hr = MapIndex(&Handle)) return hr;
			if (hr = VCUE::GenericCopy<BSTR, VARIANT>::copy(&sHandle, &Handle)) return hr;
			szHandle = estring(sHandle);
			szHandle.trim();
			CheckAddHandle(&szHandle);
									
			if (!szHandle.size()){
				// A blank handle is OK - we just stamp a GUID (as the user never sees the handles!)
				szHandle = GetNewGuid();
				bIsSpreadsheetHandle = false;
				bIsValidHandle = true;				
			} else {
				CParameterMap::RemoveCalculationNumberFromHandle(*pclsidSingular, &szHandle, &bIsSpreadsheetHandle);
				bIsValidHandle = bIsSpreadsheetHandle;
			}
			sHandle = szHandle;
		}

		// Get the key (by referencing the properties associated with the object)
		sKey = GetKey(spValue, &szName, &nDate, &ds);
		szKey = estring(sKey);
				
		if (bIsMaintainedCollection && !ShouldAddToMaintainedCollection(Value) && !bIsSpreadsheetHandle){
			// Note that we emit spreadsheet handles to the maintained collection.
			// The reason for this is quite deep: If bIsSpreadsheetHandle is set to true, then
			// that means that a handle associated with this object is being sent to the
			// Excel worksheet. The user is therefore free to recalculate the function returning
			// that handle, which therefore means that the current object in the maintained
			// collection will be replaced.
			return CParameterMap::ReturnErrorS("Object '" + szKey + "' cannot be added to the maintained collection");
		}

		// Unless szHandle, when respelt and unabridged, equals szKey (whereupon
		// we substitute a guid), the only valid form of szHandle is a spreadsheet-
		// style handle.
		if (!bIsValidHandle){
			std::string				szRespeltKey;
			long					nRespeltDate;
			DataSourceEnum			dsRespelt;
			std::string				szRespeltName;
			CComObjectCollectionSerialisableKey::Respell(m_szSingularName, szHandle, &szRespeltKey, NULL, &szRespeltName, &nRespeltDate, &dsRespelt, NULL, NULL);
			if (szRespeltKey == szKey){
				szHandle = GetNewGuid();
				sHandle = szHandle;
			} else if (dsRespelt == ds && !nDate && nRespeltDate == MlEqDate::GetCurrentDate()){
				// This can arise if the object has a zero date. (The implied behaviour of this is always to compare
				// against the current system date.)
				szHandle = GetNewGuid();
				sHandle = szHandle;
			} else {
				throw "The handle '" + estring(HandleIn) + "' is invalid";
			}
		}
		
		// SpecialAdd specific processing				
		if ((hr = SpecialAdd(Value, sKey)) != S_FALSE) return hr;
		
		// Normal item addition
		// If an item is already defined against this key then remove it.
		if (m_KeyToHandle.find(sKey) != m_KeyToHandle.end()){
			if (hr = Remove(CComVariant(sKey))) return hr;
		}
		// But if the handle is already in use then we error.
		if (m_coll.find(sHandle) != m_coll.end()){
			return CParameterMap::ReturnErrorRS(IDS_HANDLE_ALREADY_IN_USE, sHandle, __uuidof(ICollection));
		}

		// Add the item to the maps.				
		g_pApplication->GetObjectManager().AddIndex(CComPtr<ISingular>(Value), CComPtr<ICollection>(this), m_szSingularName, szKey);
#		ifdef _DEBUG				
		if (bIsMaintainedCollection){				
			std::string szKey = CComObjectCollectionSerialisableKey(Value);
			ATLTRACE("Adding (%s) %s to the maintained collection.\n", m_szSingularName.c_str(), szKey.c_str());
		}
#		endif
		m_coll[sHandle] = Value;									// Handle -> Object
		m_HandleToKey[sHandle] = sKey;								// Handle -> Key
		m_KeyToHandle[sKey] = sHandle;								// Key -> Handle
		m_ObjectToNotional[CComPtr<IDispatch>(Value)] = 1.0;		// Object -> Notional
		m_ObjectToLocalIndex[CComPtr<IDispatch>(Value)] = szKey;	// Object -> LocalIndex
		if (bIsMaintainedCollection) KeyAdded(szKey);
		
		end_function
	}
	
	//	Empty the collection
	STDMETHOD(Clear)()
	{
		begin_function
		bool bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICollection>(this));
		m_HandleToKey.clear();
		m_KeyToHandle.clear();						
		CComObjectCollection<ISingular, CCollection, ICollection, pclsidCollection, piidCollection, plibidCollection, idr>::Clear();
		if (bIsMaintainedCollection) KeysAllCleared();
		end_function
	}

	// Retrieve an item from the collection
	STDMETHOD(get_Item)(VARIANT IndexIn, ISingular** pVal)
	//	Index - can be either a number (1-based collection item counter), a Handle or a Key (in general, partially spelt)
	{		
		begin_function				
		HRESULT							hr;
		CComVariant						Index(IndexIn);
		estring							szIndex;
		bool							bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICollection>(this));
		
		if (!pVal) return E_POINTER;										
		if (hr = MapIndex(&Index)) return hr;
		if (Index.vt != VT_BSTR){
			// If the index isn't a string, but can be converted to a long value, we use ATL's (default) implementation
			if (hr = Index.ChangeType(VT_I4)) return hr;
			CComVariant vItem;
			if (Index.lVal <= 0) return E_INVALIDARG;		// CollectionType::get_Item does not catch this as an error - it just returns the first item
			if (hr = CollectionType::get_Item(Index.lVal, &vItem)) return hr;
			if (vItem.vt != VT_DISPATCH) return E_FAIL;
			*pVal = (ISingular*)vItem.pdispVal;
			(*pVal)->AddRef();
			return S_OK;						
		} else {
			// Index is a string
			szIndex = Index.bstrVal;
			CParameterMap::RemoveCalculationNumberFromHandle(*pclsidSingular, &szIndex);
		}
		
		// Respell szIndex and get any name, date and data source associated with it (we may use them later!)
		long								nDate;
		DataSourceEnum						ds;
		std::string							szName;
		estring								szKey;
		estring								szShortKey;
		
		// We need to record whether or not the default date and data source stacks have been used to
		// obtain a date and data source. This is so any subsequent load calls default the date and
		// data source according to other predefined rules (e.g. Taurus data acquisition policy).
		
		bool								bDateDefaulted = false;
		bool								bDataSourceDefaulted = false;

		CComObjectCollectionSerialisableKey::Respell(m_szSingularName, szIndex, &szKey, &szShortKey, &szName, &nDate, &ds, &bDateDefaulted, &bDataSourceDefaulted);
		
		ContainerType::const_iterator		itHandleToObject = m_coll.end();
		KeyToHandleType::const_iterator		itKeyToHandle;
				
		// We try the full key first (important as they change frequently).
		if (itHandleToObject == m_coll.end()){
			if ((itKeyToHandle = m_KeyToHandle.find(szKey.GetBSTR())) != m_KeyToHandle.end()){
				// Key found.
				itHandleToObject = m_coll.find(itKeyToHandle->second);
			}
		}
		
		// If the date has been defaulted then we try to find an item with a dateless key.
		if (itHandleToObject == m_coll.end() && bDateDefaulted){			
			if ((itKeyToHandle = m_KeyToHandle.find(szShortKey.GetBSTR())) != m_KeyToHandle.end()){
				// Key found
				itHandleToObject = m_coll.find(itKeyToHandle->second);
			}
		}

		// Try the handle.
		if (itHandleToObject == m_coll.end()){
			itHandleToObject = m_coll.find(szIndex.GetBSTR());
		}
		
		if (bIsMaintainedCollection && itHandleToObject == m_coll.end()){
			// Neither an object corresponding to a key or handle has been found.
			// If "this" is the maintained collection we can try to get the object
			// from an external source as follows.
			// (If this collection is not maintained then we can't do this since 
			// acquiring an object always inserts it into a maintained collection 
			// therefore breaking symmetry).
									
			CComPtr<IDispatch>			spObject;
			// Try obtaining a GDA object which has been created in process with a GDA[...] function call.			
			if (!spObject && ds == Last && !bDataSourceDefaulted){
				CComPtr<ISingular>		sp;
				LoadGDAObject(szName, nDate, sp);
				spObject = sp;
				// Notice here that we don't add this object to the maintained collection. This is the
				// expected behaviour. The user can now change the GDA object as they wish and Sirius
				// will pick up this object whenever it's required.
			}
			// Try to load the object from an external source (which can obviously include the Sirius database).
			if (!spObject){				
				// We have a slight complexity here with clonable objects:
				// 1) The first call to _Module.Load should not consider creating a clone.
				// 2) We then make a call to GetClone (if the object is clonable).
				// 3) We then make a further call to _Module.Load; this time with cloning enabled.
				// These three steps best mirror the situation of only accessing an object from the 
				// database if we really need to.
				bool bClonable = IsClonable();
				try {
					// We don't propagate defaulted data sources - we delegate data source selection
					// to the relevant load function.
					_Module.Load(*pclsidSingular, szName, bDataSourceDefaulted ? "" : CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds), nDate, false, spObject);
				} catch (const std::string& szError){					
					// If the load fails (for example, we may not have database access enabled), then
					// we may be able to derive spObject item from an existing one. (A good example of
					// this is with assets with differing composite and pay currencies to an existing
					// asset.)
					if (bClonable){
						GetClone(szKey, spObject);
						if (!spObject){
							_Module.Load(*pclsidSingular, szName, bDataSourceDefaulted ? "" : CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds), nDate, true, spObject);
						}						
					} else {
						throw szError;
					}										
				}
			}			
			if (spObject){				
				// return the object					
				return spObject.CopyTo((IDispatch**)pVal);
			}
		}

		if (itHandleToObject == m_coll.end()){
			// error from here
			estring szHandle(Index);
			if (szHandle.find("::") == szHandle.npos){
				szHandle.assign(m_szSingularName + "::" + szHandle);
			}
			return CParameterMap::ReturnErrorRS(IDS_NO_OBJECT_WITH_HANDLE, szHandle, __uuidof(ICollection));
		}
		
		// Item found if this point is reached.
		CComVariant v(itHandleToObject->second);
		if (v.vt != VT_DISPATCH) return E_FAIL;
		*pVal = (ISingular*)v.pdispVal;		
		(*pVal)->AddRef();
		end_function
	}

	STDMETHOD(get_Xml)(/*[out, retval]*/ BSTR* pVal)
	{
		HRESULT hr;
		xmlstreamer							ssXML;
			
		// We want to use the weak form of CXmlStreamer::GetXML since we don't want to adjust the data sources.		
		hr = CXmlStreamer::GetXML(CComVariant(this), ssXML);
		return estring(ssXML).GetBSTR(pVal);
	}

	bool IsInCollection(const CComObjectCollectionSerialisableKey& key) const
	{
		return m_KeyToHandle.find((CComBSTR)key) != m_KeyToHandle.end();
	}
	
	//	Resynchronise the key associated with an object with the data in the object.
	//  Any assertions in this function are usually indicative of memory leaks in other classes! YOU MUST FIX THESE!
	STDMETHOD(Rekey)(ISingular* Object, BSTR OldKey, BSTR* NewKey)
	//	NewKey - (returned, nullable)
	{
		// We cannot simply remove and re-add the item as doing so would
		// break any current iterations through m_coll via IEnumVariant.
		// The trick here is to keep m_coll intact.
		if (!g_pApplication->GetObjectManager().HasIndex(CComPtr<ISingular>(Object))){
			return S_OK;
		}
		
		HRESULT								hr;
		CComBSTR							sOldKey(OldKey);
		CComBSTR							sNewKey;
		KeyToHandleType::iterator			itKeyToHandle;
		HandleToKeyType::iterator			itHandleToKey;
		ContainerType::const_iterator		itHandleToObject;
		ObjectToLocalIndexType::iterator	itObjectToLocalIndex;
		bool								bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICollection>(this));

		if ((itKeyToHandle = m_KeyToHandle.find(sOldKey)) == m_KeyToHandle.end()){
			return E_FAIL;
		}
		if ((itHandleToObject = m_coll.find(itKeyToHandle->second)) == m_coll.end()){
			ATLASSERT(false);
			return E_FAIL;
		}
		if ((sNewKey = GetKey(dynamic_cast<ISingular*>(itHandleToObject->second.pdispVal))) == sOldKey) return S_OK;
		// If sNewKey already points to an object then we need to delete the old one.
		if (m_KeyToHandle.find(sNewKey) != m_KeyToHandle.end()){
			if (hr = Remove(CComVariant(sNewKey))) return hr;
		}
		if ((itHandleToKey = m_HandleToKey.find(itKeyToHandle->second)) == m_HandleToKey.end()){
			ATLASSERT(false);
			return E_FAIL;
		}
		if ((itObjectToLocalIndex = m_ObjectToLocalIndex.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal))) == m_ObjectToLocalIndex.end()){
			ATLASSERT(false);
			return E_FAIL;
		}
		ATLASSERT(itObjectToLocalIndex->second == estring(sOldKey));

		// Adjust the maps
		m_KeyToHandle.erase(itKeyToHandle);
		m_KeyToHandle[sNewKey] = itHandleToKey->first;
		itHandleToKey->second = sNewKey;
		itObjectToLocalIndex->second = estring(sNewKey);
		g_pApplication->GetObjectManager().Reindex(itHandleToObject->second.pdispVal, estring(sOldKey), estring(sNewKey));
		ATLASSERT(m_coll.size() == m_HandleToKey.size() && m_coll.size() == m_KeyToHandle.size() && m_coll.size() == m_ObjectToNotional.size() && m_coll.size() == m_ObjectToLocalIndex.size());	// Failure here is often indicative of memory leaks in other classes! YOU MUST FIX THIS!
		if (bIsMaintainedCollection) KeyRekeyed(estring(sOldKey), estring(sNewKey));
		return S_OK;
	}
	
	//	Removes an item from the collection
	//  Any assertions in this function are usually indicative of memory leaks in other classes! YOU MUST FIX THESE!
	STDMETHOD(Remove)(VARIANT IndexIn)
	{
		begin_function
		HRESULT							hr;
		CComVariant						Index(IndexIn);
		estring							szIndex;
		estring							szKey;
		estring							szShortKey;
		bool							bDateDefaulted = false;
		bool							bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICollection>(this));

		if (hr = MapIndex(&Index)) return hr;		
		if (Index.vt == VT_DISPATCH){
			// It's quite nice to support this. What we do is get the correctly spelt
			// key associated with this object.
			CComPtr<ISingular>			spSingular(dynamic_cast<ISingular*>(Index.pdispVal));
			std::string					szObjectName;

			CParameterMap::GetObjectName(spSingular, &szObjectName);
			if (szObjectName != m_szSingularName) throw "Cannot remove objects of type '" + szObjectName + "' from a collection of objects of type '" + m_szSingularName + "'";
			Index = CComVariant((CComBSTR)CComObjectCollectionSerialisableKey(spSingular));
		}

		if (Index.vt != VT_BSTR){
			// If the index isn't a string, but can be converted to a long value, we use the long implementation of Remove
			if (hr = Index.ChangeType(VT_I4)) return hr;
			return Remove(Index.lVal);
		} else {
			// Index is a string
			szIndex = Index.bstrVal;
			CParameterMap::RemoveCalculationNumberFromHandle(*pclsidSingular, &szIndex);
		}

		CComObjectCollectionSerialisableKey::Respell(m_szSingularName, szIndex, &szKey, &szShortKey, NULL, NULL, NULL, &bDateDefaulted, NULL);

		// Get itHandleToObject, itHandleToKey and itKeyToHandle.
		ContainerType::iterator					itHandleToObject;
		HandleToKeyType::iterator				itHandleToKey;
		KeyToHandleType::iterator				itKeyToHandle;
		
		if ((itKeyToHandle = m_KeyToHandle.find(szKey.GetBSTR())) != m_KeyToHandle.end()){
			// Key found.	
			if ((itHandleToKey = m_HandleToKey.find(itKeyToHandle->second)) == m_HandleToKey.end()){
				ATLASSERT(false);
				return E_FAIL;
			}
			if ((itHandleToObject = m_coll.find(itKeyToHandle->second)) == m_coll.end()){
				ATLASSERT(false);
				return E_FAIL;
			}
		} else if (bDateDefaulted && (itKeyToHandle = m_KeyToHandle.find(szKey.GetBSTR())) != m_KeyToHandle.end()){
			// Short key found
			ATLASSERT(false);	// ToDo - test and then incorporate this statement into the first if with an ||.
			if ((itHandleToKey = m_HandleToKey.find(itKeyToHandle->second)) == m_HandleToKey.end()){
				ATLASSERT(false);
				return E_FAIL;
			}
			if ((itHandleToObject = m_coll.find(itKeyToHandle->second)) == m_coll.end()){
				ATLASSERT(false);
				return E_FAIL;
			}
		} else if ((itHandleToObject = m_coll.find(szIndex.GetBSTR())) != m_coll.end()){
			// Handle found
			if ((itHandleToKey = m_HandleToKey.find(itHandleToObject->first)) == m_HandleToKey.end()){
				ATLASSERT(false);
				return E_FAIL;
			}
			if ((itKeyToHandle = m_KeyToHandle.find(itHandleToKey->second)) == m_KeyToHandle.end()){
				ATLASSERT(false);
				return E_FAIL;
			}
		} else {
			throw "No object associated with handle '" + estring(Index) + "'";
		}
		
		// And get the other maps
		ObjectToNotionalType::iterator			itObjectToNotional;
		ObjectToLocalIndexType::iterator		itObjectToLocalIndex;
		
		itObjectToNotional = m_ObjectToNotional.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToNotional == m_ObjectToNotional.end()){
			ATLASSERT(false);
			return E_FAIL;				
		}
		itObjectToLocalIndex = m_ObjectToLocalIndex.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToLocalIndex == m_ObjectToLocalIndex.end()){
			ATLASSERT(false);
			return E_FAIL;
		}
		
		// Now remove the item from the maps. Note that we could have just erased on the
		// handles and keys but it's not clear what erase does if the item isn't present.
		// Erasing on iterators is safer and allows us in previous code to inform the
		// client of failure.		
		g_pApplication->GetObjectManager().RemoveIndex(CComPtr<ISingular>(dynamic_cast<ISingular*>(itHandleToObject->second.pdispVal)), CComPtr<ICollection>(this));
		if (bIsMaintainedCollection) KeyCleared(estring(itKeyToHandle->first.m_T));
		m_coll.erase(itHandleToObject);															// this removes Name -> Object
		m_HandleToKey.erase(CComBSTR(itKeyToHandle->second));									// this removes Handle -> Name
		m_KeyToHandle.erase(itKeyToHandle);														// this removes Name -> Handle
		m_ObjectToNotional.erase(itObjectToNotional);											// this removes Object -> Notional
		m_ObjectToLocalIndex.erase(itObjectToLocalIndex);										// this removes Object -> LocalIndex
		ATLASSERT(m_coll.size() == m_HandleToKey.size() && m_coll.size() == m_KeyToHandle.size() && m_coll.size() == m_ObjectToNotional.size() && m_coll.size() == m_ObjectToLocalIndex.size());				
		end_function
	}

	// Remove an item from the collection index by number.
	STDMETHOD(Remove)(size_t Index)
	{
		begin_function
		ContainerType::iterator					itHandleToObject;
		HandleToKeyType::iterator				itHandleToKey;
		HandleToKeyType::iterator				itKeyToHandle;
		ObjectToNotionalType::iterator			itObjectToNotional;
		ObjectToLocalIndexType::iterator		itObjectToLocalIndex;
		bool									bIsMaintainedCollection = g_pApplication->GetObjectManager().IsMaintainedCollection(CComPtr<ICollection>(this));
				
		if (Index <= 0 || Index > m_coll.size()) return E_INVALIDARG;
		itHandleToObject = m_coll.begin();
		std::advance(itHandleToObject, Index - 1);
		itHandleToKey = m_HandleToKey.find(CComBSTR(itHandleToObject->first));
		if (itHandleToKey == m_HandleToKey.end()) return E_FAIL;
		itKeyToHandle = m_KeyToHandle.find(itHandleToKey->second);
		if (itKeyToHandle == m_KeyToHandle.end()) return E_FAIL;		
		itObjectToNotional = m_ObjectToNotional.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToNotional == m_ObjectToNotional.end()) return E_FAIL;
		itObjectToLocalIndex = m_ObjectToLocalIndex.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToLocalIndex == m_ObjectToLocalIndex.end()) return E_FAIL;				
		
		// Now remove the item from the maps.		
		g_pApplication->GetObjectManager().RemoveIndex(CComPtr<ISingular>(dynamic_cast<ISingular*>(itHandleToObject->second.pdispVal)), CComPtr<ICollection>(this));
		if (bIsMaintainedCollection) KeyCleared(estring(itKeyToHandle->first.m_T));
		m_coll.erase(itHandleToObject);															// this removes Handle -> Object
		m_HandleToKey.erase(itHandleToKey);														// this removes Handle -> Key
		m_KeyToHandle.erase(itKeyToHandle);														// this removes Key -> Handle
		m_ObjectToNotional.erase(itObjectToNotional);											// this removes Object -> Notional
		m_ObjectToLocalIndex.erase(itObjectToLocalIndex);										// this removes Object -> LocalIndex
		ATLASSERT(m_coll.size() == m_HandleToKey.size() && m_coll.size() == m_KeyToHandle.size() && m_coll.size() == m_ObjectToNotional.size() && m_coll.size() == m_ObjectToLocalIndex.size());	// Failure here is often indicative of memory leaks in other classes! YOU MUST FIX THIS!		
		end_function
	}

protected:
	//  Returns the key value associated with an object
	CComBSTR GetKey(CComPtr<ISingular> spSingular, std::string* pszNameOut = NULL, long* pnDateOut = NULL, DataSourceEnum* pdsOut = NULL)
	// pszNameOut - object's name (returned, nullable)
	{
		CComBSTR						sName;
		DATE							date;
		DataSourceEnum					ds;
		
		spSingular->get_Name(&sName);
		spSingular->get_Date(&date);
		spSingular->get_DataSource(&ds);
		if (pszNameOut) pszNameOut->assign(estring(sName));
		if (pnDateOut) *pnDateOut = date;
		if (pdsOut) *pdsOut = ds;
		return CComObjectCollectionSerialisableKey(estring(sName), date, ds);
	}
	
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal)
	{
		begin_function
		HRESULT								hr;	
		CComBSTR							sOut;							// return string
		for (ContainerType::iterator it = m_coll.begin(); it != m_coll.end(); it++){
			CComVariant	vObject(it->second);
			CComBSTR	s;
			if (vObject.vt != VT_DISPATCH) return !S_OK;		
			if (hr = ((ISingular*)vObject.pdispVal)->Save(&s)) return hr;
			if (sOut.Length()) sOut += L", ";
			sOut += s;
		}
		if (pVal){
			return sOut.CopyTo(pVal);	
		} else {
			return S_OK;
		}
		end_function
	}	
};


#endif
