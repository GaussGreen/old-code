//	SiriusComObjectCollection.h : This is a template object collection class designed specifically for
//							      Sirius objects that are properies of SiriusApplication.
//
//								  The class contains maintained collection logic.
//
//
//	Author :					  David Cuin
//
/////////////////////////////////////////////////////////////////////////////


#ifndef _SIRIUSCOMOBJECTCOLLECTION_H
#define _SIRIUSCOMOBJECTCOLLECTION_H

#include <comobjectcollection.h>
#include "guidex.h"
#include "siriusapplication.h"

template <class ISingular, class CCollection, class ICollection, const CLSID* pclsidCollection, const IID* piidCollection, const GUID* plibidCollection, const CLSID* pclsidSingular, const UINT idr>
	class ATL_NO_VTABLE CSiriusComObjectCollection:
	public CComObjectCollection<ISingular, CCollection, ICollection, pclsidCollection, piidCollection, plibidCollection, idr>
{
public:
	HRESULT FinalConstruct(void)
	{
		// This is used to set up m_szSingularName
		HRESULT							hr;
		CComBSTR						sProgID;
		std::vector<std::string>		vector;
		estring							szSingularProgID;
														
		if (hr = CParameterMap::ProgIDFromCLSID(*pclsidSingular, sProgID)) return hr;
		szSingularProgID.Set(sProgID);		
		// m_szSingularName occurs between the two periods in szSingularProgID
		szSingularProgID.Split(".", &vector);
		if (vector.size() != 3) return E_FAIL;
		m_szSingularName.assign(vector[1]);		
		return S_OK;
	}

protected:
	ContainerType::const_iterator ImplementGetItem(const CComVariant& Index) const
	{		
		CComVariant							vHandle;
		ContainerType::const_iterator       itHandleToObject;
		
		if (vHandle.ChangeType(VT_BSTR, &Index)) throw "No object associated with handle";
		
		estring							szHandle(vHandle.bstrVal);
		CParameterMap::RemoveCalculationNumberFromHandle(*pclsidSingular, &szHandle);

		if ((itHandleToObject = m_coll.find(szHandle.GetBSTR())) == m_coll.end()){
			// error
			estring		szHandle(Index);
			if (szHandle.find("::") == szHandle.npos){
				szHandle.assign(m_szSingularName + "::" + szHandle);
			}
			throw "No object associated with handle '" + szHandle + "'";
		}
		return itHandleToObject;		
	}

protected:
	estring								m_szSingularName;				// name of the singular object associated with the collection	

protected:	
// The protection of all member functions is deliberate. We only want them to be usable if
// sirius.idl exposes them. The parent class then publicises these functions.

	void AddGlobal(const CComVariant& Index, ISingular* Value)
	{
		g_pApplication->GetObjectManager().AddIndex(CComPtr<ISingular>(Value), CComPtr<ICollection>(this), m_szSingularName, estring(Index));				
	}

	void ClearGlobal(void)
	{				
		// Clear the relevant objects from the Object -> Index map (only if this object is maintained)		
		for (ContainerType::const_iterator it = m_coll.begin(); it != m_coll.end(); it++){						
			g_pApplication->GetObjectManager().RemoveIndex(CComPtr<ISingular>(dynamic_cast<ISingular*>(it->second.pdispVal)), CComPtr<ICollection>(this));
		}		
	}

	void RemoveGlobal(ISingular* Value)
	{		
		g_pApplication->GetObjectManager().RemoveIndex(CComPtr<ISingular>(Value), CComPtr<ICollection>(this));		
	}

protected:
	STDMETHOD(get_Size)(/*[out, retval]*/ long* pn)
	{
		begin_function
		*pn = m_coll.size();			
		end_function
	}

	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal)
	{
		// Return a parameter map with columns Name, Object, Notional
		begin_function
		CParameterMap							pm(m_coll.size(), 2);
		long									nRow = 0;
		ObjectToNotionalType::const_iterator	itObjectToNotional;

		for (ContainerType::const_iterator it = m_coll.begin(); it != m_coll.end(); it++){
			pm.SetValue(nRow, 0, it->second);
			double		fNotional = 1.0;			
			if ((itObjectToNotional = m_ObjectToNotional.find(CComPtr<IDispatch>(it->second.pdispVal))) != m_ObjectToNotional.end()){
				fNotional = itObjectToNotional->second;
			}			
			pm.SetValue(nRow, 1, fNotional);
			nRow++;
		}						
		return pm.GetValue(pVal);		
		end_function
	}

	STDMETHOD(put_Size)(/*[in]*/ long nNewSize)
	{		
		begin_function
		if (nNewSize == m_coll.size()){
			// do nothing			
		} else if (nNewSize < 0){
			throw estring(nNewSize) + " is an invalid value for the Size property";
		} else if (!nNewSize){
			return Clear();
		} else if (nNewSize < m_coll.size()){
			while (nNewSize < m_coll.size()){
				if (Remove(CComVariant((long)m_coll.size()))) propagate_error(_Module.GetSiriusApplication());
			}		
		} else if (nNewSize > m_coll.size() + collection_max_size_increase){
			throw "The increment in size (up to " + estring(nNewSize - m_coll.size()) + " elements) is too big. The limit is " + estring(collection_max_size_increase);
		} else {
			// Increase from m_coll.size() to nNewSize
			while (nNewSize > m_coll.size()){
				CComPtr<ISingular> spSingular;
				spSingular.CoCreateInstance(*pclsidSingular);
				if (Add(CComVariant(), spSingular)) propagate_error;
			}
		}
		end_function
	}
	
	
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal)
	{
		// This is and has to be symmetric with get_Value.
		// ToDo - avoid setting member variables until the end.
		CParameterMap					pm;		
		HRESULT							hr;

		pm.SetValue(newVal);
		if (pm.GetCols() != 2) throw "The number of columns is invalid in the value variant of a collection of objects of type '" + m_szSingularName + "'";
		if (hr = Clear()) return hr;
		for (long nRow = 0; nRow < pm.GetRows(); nRow++){
			CComVariant					vHandle;
			double						fNotional;
			CComPtr<IDispatch>			spObject;
						
			if (hr = pm.GetValue(nRow, 0, &vHandle)) return hr;			
			if (hr = g_pApplication->GetObjectManager().GetObject(vHandle, spObject)) return hr;
			if (hr = pm.GetValue(nRow, 1, &fNotional)) return hr;									
			if (hr = Add(CComVariant(), (ISingular*)spObject.p)) return hr;
			if (hr = put_Notional(CComVariant(spObject.p), fNotional)) return hr;
		}		
		return S_OK;
	}
};


#endif