//	comobjectcollection.h : This is a template object collection class. It requires:
//
//							ISingular - the type of interface that you're calling.
//							            i.e. IObjects collects IObject so ISingular would be IObject since IObject is the singular interface
//
//							CCollection  - the actual collection class (i.e. CObjects)
//
//							ICollection	 - the interface of the collection class (i.e. IObjects)
//
//							pclsidCollection, piidCollection - pointers to the class ID, the interface ID for the collection object
//
//							plibid - pointer to the library ID of the singular AND collection
//
//							pclisdSingular - pointer to the class ID for the singular object
//
//							idr	- resource identifier for the class (usually of form IDR_Collection)
//
//	Author :				David Cuin
//
//	Notes :					ATL's collection implementation only allows numeric indexing (due to container
//							independence) so we have to provide our own implementation to allow string indexing
//						    into the map.
//
//							The bass class variable, m_coll, holds Handle -> Object.
//
//							We also support an Object -> Notional map. This provides a method of associating a quantity
//							with each element in the collection.
//
//							
/////////////////////////////////////////////////////////////////////////////


#ifndef _COMOBJECTCOLLECTION_H
#define _COMOBJECTCOLLECTION_H

#include <guidex.h>
#include <VCUE_Collection.h>
#include <VCUE_Copy.h>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#undef max


//////////////////////////////////////////////////////////////////////////////
//	CComObjectCollection
//
//	main class
//
typedef std::map<CAdapt<CComBSTR>, CComVariant>				ContainerType;
typedef VARIANT												ExposedType;
typedef IEnumVARIANT										EnumeratorInterface;				
typedef VCUE::MapCopy<ContainerType, ExposedType>			CopyType;
typedef std::map<CAdapt<CComPtr<IDispatch> >, double>		ObjectToNotionalType;
typedef std::map<CAdapt<CComPtr<IDispatch> >, estring>		ObjectToLocalIndexType;				// used by get_Index
#define collection_max_size_increase						256									// Maximum amount by which we can increase a collection size.

typedef CComEnumOnSTL<EnumeratorInterface, &__uuidof(EnumeratorInterface), ExposedType, CopyType, ContainerType > EnumeratorType;

template <class ISingular, class CCollection, class ICollection, const CLSID* pclsidCollection, const IID* piidCollection, const GUID* plibid, const UINT idr>
	class ATL_NO_VTABLE CComObjectCollection : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCollection, pclsidCollection>,
	public ISupportErrorInfo,
	public IDispatchImpl<VCUE::ICollectionOnSTLCopyImpl< ICollection, ContainerType, ExposedType, CopyType, EnumeratorType >,  piidCollection, plibid>
{		
protected:
	// Virtual Functions. Override these in base classes for specific functionality.		
	virtual HRESULT MapIndex(CComVariant* pIndex)						{return S_OK;}			// Override this to enforce any Index spelling conventions.
	
	// These functions allow the association of object indices with a global index
	virtual void AddGlobal(const CComVariant& Index, ISingular* Value)	{return;}
	virtual void ClearGlobal(void)										{return;}
	virtual void RemoveGlobal(ISingular* Value)							{return;}

	// This function allows further respelling of an input index; e.g. removal of a calculation adornment.
	virtual ContainerType::const_iterator ImplementGetItem(const CComVariant& Index) const
	{		
		CComVariant						vHandle;
		ContainerType::const_iterator   itHandleToObject;
		
		if (vHandle.ChangeType(VT_BSTR, &Index)) throw "No object associated with handle";
		
		estring							szHandle(vHandle.bstrVal);

		if ((itHandleToObject = m_coll.find(szHandle.GetBSTR())) == m_coll.end()){
			throw "No object associated with handle '" + szHandle + "'";
		}
		return itHandleToObject;		
	}

public:	
	CComObjectCollection(){}
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	DECLARE_REGISTRY_RESOURCEID(idr)
	BEGIN_COM_MAP(CCollection)
		COM_INTERFACE_ENTRY(ICollection)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
					
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid)
	{
		static const IID* arr[] = { piidCollection };
		for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
			if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
		}
		return S_FALSE;
	}
		
	HRESULT FinalRelease(void)
	{
		return Clear();					// to remove Object -> Index (which is, of course, a global map)
	}

	//	Adds a list of input objects to this collection, avoiding duplicates
	void MergeCollection(CComPtr<ICollection> spCollection)
	{									
		CCollection*						pCollection = dynamic_cast<CCollection*>(spCollection.p);
		
		for (ContainerType::const_iterator it = pCollection->m_coll.begin(); it != pCollection->m_coll.end(); it++){
			// check that it->second.pdispVal is not in the parent
			bool bFound = false;
			for (ContainerType::const_iterator itThis = m_coll.begin(); itThis != m_coll.end() && !bFound; itThis++){
				if (itThis->second.pdispVal == it->second.pdispVal){
					// Already in collection - increment the notional.
					double fNotional_This, fNotional_Other;
					get_Notional(itThis->second, &fNotional_This);
					get_Notional(it->second, &fNotional_Other);
					fNotional_This += fNotional_Other;
					put_Notional(itThis->second, fNotional_This);					
					// We need to check other elements in itThis and add their notional contributions.
					// This is why I have this curious bFound logic; which appears ostensibly that I've forgotten about
					// break and continue statements in for loops!
					bFound = true;
				}
			}
			if (!bFound && Add(CComVariant()/*i.e. use the generated key*/, (ISingular*)it->second.pdispVal)) propagate_error;
		}		
	}
	
	void MergeSingular(CComPtr<ISingular> spSingular)
	{
		if (!spSingular) return;
		for (ContainerType::const_iterator itThis = m_coll.begin(); itThis != m_coll.end(); itThis++){
			if (itThis->second.pdispVal == spSingular.p){
				// Already in collection - increment the notional by 1.				
				double fNotional;
				if (get_Notional(itThis->second, &fNotional)) propagate_error;
				fNotional += 1.0;
				if (put_Notional(itThis->second, fNotional)) propagate_error;				
			}
		}
		// we need to add this item		
		if (Add(CComVariant()/*i.e. use the generated key*/, spSingular)) propagate_error;
	}
		

// The protection of all member functions is deliberate. We only want them to be usable if
// the idl file exposes them. The parent class then publicises these functions.
protected:							
	ObjectToNotionalType				m_ObjectToNotional;				// this allows us to associate a notional amount with each element added
	ObjectToLocalIndexType				m_ObjectToLocalIndex;			// this is used to implement get_Index	
	GuidEx								m_guid;							// used to allocate 'random' handles
		
	typedef VCUE::ICollectionOnSTLCopyImpl< ICollection, ContainerType, ExposedType, CopyType, EnumeratorType > CollectionType;	
	
	//	Add a new item into the collection.
	STDMETHOD(Add)(VARIANT IndexIn, ISingular* Value)
	{		
		CComDispatchDriverEx			ddValue(Value);
		HRESULT							hr;				
		CComBSTR						sHandle;
		CComVariant						Index(IndexIn);		
						
		if (Index.vt != VT_EMPTY){
			if (hr = MapIndex(&Index)) return hr;
			if (hr = VCUE::GenericCopy<BSTR, VARIANT>::copy(&sHandle, &Index)) return hr;
		}				
		{
			estring szHandle(sHandle);
			szHandle.trim();
			if (!szHandle.size()){
				szHandle = GetNewGuid();
			}
			sHandle = szHandle;
		}

		// Fail if an item is already defined against sHandle.
		if (m_coll.find(sHandle) != m_coll.end()){
			return CParameterMap::ReturnErrorSS("Cannot create duplicate name", estring(sHandle), __uuidof(ICollection));
		}
				
		// Add the item to the maps				
		AddGlobal(Index, Value);		
		m_coll[sHandle] = Value;															// Handle -> Object
		m_ObjectToNotional[CComPtr<IDispatch>(Value)] = 1.0;								// Object -> Notional
		m_ObjectToLocalIndex[CComPtr<IDispatch>(Value)] = estring(sHandle);					// Object -> LocalIndex
		return S_OK;
	}

	//	Empty the collection
	STDMETHOD(Clear)()
	{				
		begin_function
		ClearGlobal();		
		m_ObjectToNotional.clear();
		m_ObjectToLocalIndex.clear();
		m_coll.clear();
		end_function
	}
	
	STDMETHOD(get_Index)(ISingular* Item, BSTR* Index)
	{
		ObjectToLocalIndexType::const_iterator	it;
		if ((it = m_ObjectToLocalIndex.find(CComPtr<IDispatch>(Item))) == m_ObjectToLocalIndex.end()) return CParameterMap::ReturnErrorS("Could not find object", __uuidof(ICollection));
		return it->second.GetBSTR(Index);
	}

	// Retrieve an item from the collection
	STDMETHOD(get_Item)(VARIANT IndexIn, ISingular** pVal)
	//	Index - can be either a number (1-based collection item counter) or a Handle (we strip of the calculation counter).
	{		
		begin_function		
		HRESULT							hr;		
		ContainerType::const_iterator	itHandleToObject;
		CComVariant						Index(IndexIn);				
		
		if (pVal == NULL) return E_POINTER;
		if (hr = MapIndex(&Index)) return hr;
		if (Index.vt != VT_BSTR){
			// If the index isn't a string, but can be converted to a long value, we use ATL's (default) implementation.
			CComVariant vHandle;			
			if (!(hr = vHandle.ChangeType(VT_I4, &Index))){
				CComVariant vItem;
				if (vHandle.lVal <= 0) return E_INVALIDARG;		// CollectionType::get_Item does not catch this as an error - it just returns the first item
				if (hr = CollectionType::get_Item(vHandle.lVal, &vItem)) return hr;
				if (vItem.vt != VT_DISPATCH) return E_FAIL;				
				CComPtr<ISingular> spSingular(dynamic_cast<ISingular*>(vItem.pdispVal));
				return spSingular.CopyTo(pVal);				
			}
		}				
		itHandleToObject = ImplementGetItem(Index);

		// This point is reached only if an item has been found.
		CComVariant vItem(itHandleToObject->second);		
		if (vItem.vt != VT_DISPATCH) return E_FAIL;		
		return CComPtr<ISingular>(dynamic_cast<ISingular*>(vItem.pdispVal)).CopyTo(pVal);
		end_function
	}

	std::string GetNewGuid(void)
	{		
		return ++m_guid;
	}

	STDMETHOD(get_Notional)(VARIANT IndexIn, double* pf)
	{
		CComPtr<ISingular>						sp;
		HRESULT									hr;
		ObjectToNotionalType::const_iterator	it;
		CComVariant								Index(IndexIn);			// Sometimes it is better to copy a variant (especially if one might modify it). This function, and others in this class are good examples.
		
		if (hr = MapIndex(&Index)) return hr;
		if (Index.vt == VT_DISPATCH){
			// an object has been passed in
			sp = (ISingular*)Index.pdispVal;
		} else {
			if (hr = get_Item(Index, &sp)) return hr;
		}
		if ((it = m_ObjectToNotional.find(CComPtr<IDispatch>(sp.p))) == m_ObjectToNotional.end()) return CParameterMap::ReturnErrorSV("Could not find object", Index, __uuidof(ICollection));
		*pf = it->second;
		return S_OK;			
	}

	STDMETHOD(put_Notional)(VARIANT IndexIn, double f)
	{
		CComPtr<ISingular>					sp;
		HRESULT								hr;
		ObjectToNotionalType::iterator		it;
		CComVariant							Index(IndexIn);

		if (hr = MapIndex(&Index)) return hr;
		if (Index.vt == VT_DISPATCH){
			// an object has been passed in
			sp = (ISingular*)Index.pdispVal;
		} else {
			if (hr = get_Item(Index, &sp)) return hr;
		}
		if ((it = m_ObjectToNotional.find(CComPtr<IDispatch>(sp.p))) == m_ObjectToNotional.end()) return CParameterMap::ReturnErrorSV("Could not find object", Index, __uuidof(ICollection));			
		it->second = f;	
		return S_OK;	
	}
		
	//	Removes an item from the collection
	STDMETHOD(Remove)(VARIANT IndexIn)
	{
		HRESULT									hr;
		CComVariant								vHandle;
		ContainerType::iterator					itHandleToObject;		
		ObjectToNotionalType::iterator			itObjectToNotional;
		ObjectToLocalIndexType::iterator		itObjectToLocalIndex;
		CComVariant								Index(IndexIn);
						
		if (Index.vt != VT_BSTR){
			// If the index isn't a string, but can be converted to a long value, we use the long implementation of Remove
			if (!vHandle.ChangeType(VT_I4, &Index)){
				return Remove(vHandle.lVal);
			}
		} else {
			if (hr = MapIndex(&Index)) return hr;
		}
		
		if (hr = vHandle.ChangeType(VT_BSTR, &Index)) return hr;		
		itHandleToObject = m_coll.find(CComBSTR(vHandle.bstrVal));						
		if (itHandleToObject == m_coll.end()) return E_FAIL;
		itObjectToNotional = m_ObjectToNotional.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToNotional == m_ObjectToNotional.end()) return E_FAIL;				
		itObjectToLocalIndex = m_ObjectToLocalIndex.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToLocalIndex == m_ObjectToLocalIndex.end()) return E_FAIL;						
		RemoveGlobal(dynamic_cast<ISingular*>(itHandleToObject->second.pdispVal));
		m_coll.erase(itHandleToObject);															// this removes Handle -> Object
		m_ObjectToNotional.erase(itObjectToNotional);											// this removes Object -> Notional
		m_ObjectToLocalIndex.erase(itObjectToLocalIndex);										// this removes Object -> LocalIndex
		ATLASSERT(m_coll.size() == m_ObjectToNotional.size() && m_coll.size() == m_ObjectToLocalIndex.size());	// Failure here is often indicative of memory leaks in other classes! YOU MUST FIX THIS!
		return S_OK;
	}

	// Remove an item from the collection index by number.
	STDMETHOD(Remove)(size_t Index)
	{
		ContainerType::iterator					itHandleToObject;		
		ObjectToNotionalType::iterator			itObjectToNotional;
		ObjectToLocalIndexType::iterator		itObjectToLocalIndex;
				
		if (Index <= 0 || Index > m_coll.size()) return E_INVALIDARG;
		itHandleToObject = m_coll.begin();
		std::advance(itHandleToObject, Index - 1);
		itObjectToNotional = m_ObjectToNotional.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToNotional == m_ObjectToNotional.end()) return E_FAIL;
		itObjectToLocalIndex = m_ObjectToLocalIndex.find(CComPtr<IDispatch>(itHandleToObject->second.pdispVal));
		if (itObjectToLocalIndex == m_ObjectToLocalIndex.end()) return E_FAIL;						
		RemoveGlobal(dynamic_cast<ISingular*>(itHandleToObject->second.pdispVal));
		m_coll.erase(itHandleToObject);															// this removes Handle -> Object
		m_ObjectToNotional.erase(itObjectToNotional);											// this removes Object -> Notional
		m_ObjectToLocalIndex.erase(itObjectToLocalIndex);										// this removes Object -> LocalIndex
		ATLASSERT(m_coll.size() == m_ObjectToNotional.size() && m_coll.size() == m_ObjectToLocalIndex.size());	// Failure here is often indicative of memory leaks in other classes! YOU MUST FIX THIS!
		return S_OK;
	}	
};


#endif