//	comobjectcollectionfunctions.h
//
//	Implements functions common to more than one, but not all, collection
//  classes. This is preferable to multiple inheritance for disparate collections.
//
//	Author :			David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _COMOBJECTCOLLECTIONFUNCTIONS_H
#define _COMOBJECTCOLLECTIONFUNCTIONS_H

typedef std::map<CAdapt<CComBSTR>, CComVariant>				ContainerType;
typedef std::map<CAdapt<CComPtr<IDispatch> >, double>		ObjectToNotionalType;
#include "xmlstreamer.h"


template <class ISingular>
class CComObjectCollectionFunctions
{
protected:
	const ContainerType*							m_pcoll;
	const ObjectToNotionalType*						m_pObjectToNotional;
	const std::string*								m_pszError;						// used by implement load collection
	CComQIPtr<IEvaluatable>							m_spEvaluatable;				// used by GetZeroCurves

public:	
	CComObjectCollectionFunctions(void) : m_pszError(NULL), m_pcoll(NULL), m_pObjectToNotional(NULL) {}
	explicit CComObjectCollectionFunctions(const ContainerType* pcoll) : m_pszError(NULL), m_pcoll(pcoll), m_pObjectToNotional(NULL) {}
	explicit CComObjectCollectionFunctions(const ContainerType* pcoll, const ObjectToNotionalType* pObjectToNotional) : m_pszError(NULL), m_pcoll(pcoll), m_pObjectToNotional(pObjectToNotional) {}
	explicit CComObjectCollectionFunctions(const std::string* pszError) : m_pszError(pszError), m_pcoll(NULL), m_pObjectToNotional(NULL) {}
	explicit CComObjectCollectionFunctions(CComQIPtr<IEvaluatable> spEvaluatable) : m_pszError(NULL), m_pcoll(NULL), m_pObjectToNotional(NULL), m_spEvaluatable(spEvaluatable){}

	//	Evaluates all the products in a collection (even if the collection does not
	//  have an explicit products collection). The singular object must implement
	//	an evaluate function.
	HRESULT Evaluate(const BSTR& Calculate, IResult** pVal)
	{
		CComPtr<IResult>						spResultOut;
		HRESULT									hr;
		ObjectToNotionalType::const_iterator	itObjectToNotional;

		ATLASSERT(m_pObjectToNotional);
		if (hr = spResultOut.CoCreateInstance(CLSID_Result)) return hr;
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComQIPtr<ISingular>			spSingular(it->second.pdispVal);
			CComPtr<IResult>				spResult;
			CComVariant						vResult;
			CComDispatchDriverEx			ddSingular(spSingular);
			
			if (hr = ddSingular.Invoke1(L"Evaluate", &CComVariant(Calculate), &vResult)) return hr;
			if (vResult.vt != VT_DISPATCH) return E_FAIL;
			if (!(spResult = dynamic_cast<IResult*>(vResult.pdispVal))) return E_FAIL;

			if ((itObjectToNotional = m_pObjectToNotional->find(CComPtr<IDispatch>(it->second.pdispVal))) != m_pObjectToNotional->end()){
				spResult->Multiply(itObjectToNotional->second);
			}
			if (hr = spResultOut->Add(spResult)) return hr;
		}
		return spResultOut.CopyTo(pVal);
	}

	//	Returns a date that is the (consistent) date for each underlying. An exception is
	//  thrown if the dates for the singular objects differ.
	long GetDate(void)
	{
		long							nDateRet = 0;
		
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			DATE						date = 0.0;			
			CComQIPtr<ISingular>		spSingular(it->second.pdispVal);
			if (spSingular->get_Date(&date)) propagate_error;
			long						nDate = date;			
			if (!nDateRet){
				nDateRet = nDate;
			} else if (nDateRet != nDate){
				throw "Inconsistent dates " + MlEqDate(nDateRet).GetString() + " and " + MlEqDate(nDate).GetString() + " encountered when accessing the date property for a collection object";
			}

		}
		
		return nDateRet;
	}
	
	//  Returns a collection of FX assets on which the object depends. The singular object
	//  must implement a get_FxUnderlyings function
	HRESULT GetFxUnderlyings(IAssets** pVal)
	{
		HRESULT								hr;
		CComPtr<IAssets>					spAssetsReturn;
		CAssets*							pAssetsReturn;

		if (hr = spAssetsReturn.CoCreateInstance(CLSID_Assets)) return hr;
		pAssetsReturn = dynamic_cast<CAssets*>(spAssetsReturn.p);
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){			
			CComQIPtr<ISingular>			spSingular(it->second.pdispVal);
			CComPtr<IAssets>				spAssets;		
			ATLASSERT(spSingular);
			if (hr = spSingular->GetFxUnderlyings(&spAssets)) return hr;
			pAssetsReturn->MergeCollection(spAssets);
		}		
		return spAssetsReturn.CopyTo(pVal);
	}
	
	//	Returns a collection of assets on which the object depends. The singular object
	//  must implement a GetUnderlyings function.
	HRESULT GetUnderlyings(IAssets** pVal)
	{
		HRESULT								hr;
		CComPtr<IAssets>					spAssetsReturn;
		CAssets*							pAssetsReturn;

		if (hr = spAssetsReturn.CoCreateInstance(CLSID_Assets)) return hr;
		pAssetsReturn = dynamic_cast<CAssets*>(spAssetsReturn.p);
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){			
			CComQIPtr<ISingular>			spSingular(it->second.pdispVal);
			CComPtr<IAssets>				spAssets;		
			ATLASSERT(spSingular);
			if (hr = spSingular->GetUnderlyings(&spAssets)) return hr;
			pAssetsReturn->MergeCollection(spAssets);
		}		
		return spAssetsReturn.CopyTo(pVal);
	}
	
	//  Returns the volatility structures on which the collection of assets on which the
	//  object depends. The singular object must implement a GetUnderlyings function.
	HRESULT GetVolatilityStructures(IVolatilityStructures** pVal)
	{
		HRESULT								hr;
		CComPtr<IAssets>					spAssets;
		CComPtr<IVolatilityStructures>		spVolatilityStructures;

		if (hr = GetUnderlyings(&spAssets)) return hr;
		if (hr = spAssets->GetVolatilityStructures(&spVolatilityStructures)) return hr;
		return spVolatilityStructures.CopyTo(pVal);
	}

	//  Returns the zero curves on which the collection of FX underlyings on which the
	//  object depends. The m_spEvaluatable pointer must be valid.
	HRESULT GetZeroCurves(IZeroCurves** pVal)
	{
		HRESULT								hr;		
		CComPtr<IAssets>					spAssets;
		CComPtr<IZeroCurves>				spZeroCurvesReturn;
		CZeroCurves*						pZeroCurvesReturn;
		long								nCount;

		if (hr = spZeroCurvesReturn.CoCreateInstance(CLSID_ZeroCurves)) return hr;
		pZeroCurvesReturn = dynamic_cast<CZeroCurves*>(spZeroCurvesReturn.p);
		if (hr = m_spEvaluatable->GetFxUnderlyings(&spAssets)) return hr;				
		spAssets->get_Count(&nCount);
		for (CComVariant vItem = 1L; vItem.lVal <= nCount; vItem.lVal++){
			CComPtr<IAsset>					spAsset;
			CComPtr<IZeroCurve>				spZeroCurve;
			spAssets->get_Item(vItem, &spAsset);
			spAsset->get_CurrencyCurve(&spZeroCurve);
			pZeroCurvesReturn->MergeSingular(spZeroCurve);
		}
		return spZeroCurvesReturn.CopyTo(pVal);
	}
		
	template <class ICollection>
	HRESULT ImplementLoadCollection_NTFS(CComPtr<ICollection> spCollection, UINT nConstraint, const estring& szConstraintValue, const std::string& szFirstLocation, const std::string& szSecondLocation, DataSourceEnum ds, long nDate, HRESULT (*load)(const std::string&, DataSourceEnum, long, CComPtr<ISingular>&)) const
	//	nConstraint - resource identifier of a property name string
	//  load - this is the address of the load function for the corresponding singular object (e.g. CPosition::Load or CDeal::Load)
	{
		HRESULT								hr;		
		std::map<std::string, bool>			mapObjects;						// objects to load
		std::string							szDirectory;					// directory continaing the collection
		estring								szConstraint;
		
		// Get the union of the files in szFirstLocation and szSecondLocation
		if (szFirstLocation.size()){
			if (hr = g_pApplication->GetObjectManager().GetFileName(__uuidof(ISingular), szFirstLocation, "", ds, nDate, &szDirectory)) return hr;
			if (hr = GetObjectFilesInDirectory(szDirectory, &mapObjects)) return hr;
		}
		if (szSecondLocation.size()){
			if (hr = g_pApplication->GetObjectManager().GetFileName(__uuidof(ISingular), szSecondLocation, "", ds, nDate, &szDirectory)) return hr;
			if (hr = GetObjectFilesInDirectory(szDirectory, &mapObjects)) return hr;
		}
	
		szConstraint.LoadString(nConstraint);	
		for (std::map<std::string, bool>::const_iterator it = mapObjects.begin(); it != mapObjects.end(); it++){
			// Note that the methodology here assumes that criteria for
			// the selection of szFirstLocation and szSecondLocation is
			// identical in this function to that pointed to by load().
			CComPtr<ISingular>			spSingular;
			if (hr = load(it->first, ds, nDate, spSingular)) return hr;

			// check the constraint is satisfied
			bool						bConstraintOK = true;
			if (szConstraintValue.size()){
				CComDispatchDriverEx	dd(spSingular);
				CComVariant				vConstraintValue;
				if (hr = dd.GetPropertyByName(szConstraint.GetBSTR(), &vConstraintValue)) bConstraintOK = false;
				if (szConstraintValue != vConstraintValue) bConstraintOK = false;
			}				
			// add the singular object to the collection if the constraint is satisfied
			if (bConstraintOK){								
				if (hr = CComDispatchDriverEx(spCollection).Invoke2(DISPID_ADD, &CComVariant(), &CComVariant(spSingular))) return hr;
				// We need an explicit insertion because spCollection, in general, is not
				// the maintained colletion and an arbitrary load function is not
				// guaranteed to insert either.
				g_pApplication->GetObjectManager().InsertObject(spSingular, false);
			}	
		}
		return S_OK;
	}	

	//	Enables the loading of a collection from a SQL Server database
	template<class ICollection>
	HRESULT ImplementLoadCollection_SQL(const std::stringstream& ssQuery, DataSourceEnum ds, CComPtr<ICollection>& spCollection)
	// ssQuery - Must return ID and Value as the zeroth and first parameters.
	//           We allow it to return duplicate IDs. The first one is the
	//           only one that is processed.
	{
		_ConnectionPtr									spConnection;
		_RecordsetPtr									prs;							// ADO recordset				
		HRESULT											hr;
		std::map<std::string, bool>						mapAdded;						// true if we have already considered an ID		
		std::map<std::string, bool>::const_iterator		it;

		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw *m_pszError;
			while (!prs->adoEOF){
				CComVariant			vID;				
				if ((vID = prs->GetFields()->GetItem(0L)->GetValue()).vt != VT_BSTR) return E_FAIL;
				estring				sz(vID);
				if (mapAdded.find(sz) == mapAdded.end()){
					// not yet added
					CComVariant		vField;
					mapAdded[sz] = true;
					if ((vField = prs->GetFields()->GetItem(1L)->GetValue()).vt != VT_BSTR) return E_FAIL;					
					CComPtr<ISingular>					spSingular;					
					if (hr = CXmlStreamer::GetObject((char*)_bstr_t(vField.bstrVal), ds, spSingular)) return hr;
					if (hr = spCollection->Add(CComVariant(), spSingular)) return hr;
					// We need an explicit insertion because spCollection is not, in general, the maintained colletion
					g_pApplication->GetObjectManager().InsertObject(spSingular, false);	
				}
				prs->MoveNext();
			}
		} catch (_com_error& e){
			throw estring(e);
		}
		return S_OK;
	}
	
	static void ImplementLoad_NTFS(const IID& iid, const std::string& szFirstLocation, const std::string& szSecondLocation, const std::string& szName, DataSourceEnum ds, long nDate, CComPtr<ISingular>& spObject)
	//  szSecondLocation - if non-blank then we try this location
	{
		std::string							szFileName;						// root directory of the Sirius file system database				
		std::string							szError;

		// try szFirstLocation
		if (!g_pApplication->GetObjectManager().GetFileName(iid, szFirstLocation, szName, ds, nDate, &szFileName)){
			try {
				ImplementLoad_NTFS(szFileName, ds, spObject);
			} catch (const std::string& sz){
				szError.assign(sz);
			}
		}
		// try szSecondLocation
		if (szSecondLocation.size() && szSecondLocation != szFirstLocation){
			if (g_pApplication->GetObjectManager().GetFileName(iid, szSecondLocation, szName, ds, nDate, &szFileName)) propagate_error;
			ImplementLoad_NTFS(szFileName, ds, spObject);
		} else {
			if (szError.size()) throw szError;
		}
	}

	static HRESULT ImplementSave_NTFS(const xmlstreamer& ssXML, const IID& iid, const std::string& szLocation, const std::string& szName, DataSourceEnum ds, long nDate, BSTR* pVal)
	//	pVal - returned, nullable
	{
		HRESULT								hr;
		std::string							szFileName;

		if (hr = g_pApplication->GetObjectManager().GetFileName(iid, szLocation, szName, ds, nDate, &szFileName)) return hr;
		hr = ImplementSave_NTFS(ssXML, szFileName, pVal);
		return hr;
	}

	// Puts the data source property to a given value for all objects in an input collection.
	template<class ICollection> void PutDataSource(CComPtr<ICollection> spCollection, DataSourceEnum ds)
	{
		long nCount = 0L;
		spCollection->get_Count(&nCount);
		for (CComVariant v = 1L; v.lVal <= nCount; v.lVal++){
			CComPtr<ISingular> spSingular;
			spCollection->get_Item(v, &spSingular);
			spSingular->put_DataSource(ds);
		}
	}
	
	//	Sets the DataSource property for all objects in the collection.
	void PutDataSource(DataSourceEnum ds)
	{
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComQIPtr<ISingular>		spSingular(it->second.pdispVal);
			if (spSingular->PutDataSource(ds)) propagate_error;
		}		
	}

	//	Sets the Date property for all objects in the collection.
	void PutDate(long nDate)
	{
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComQIPtr<ISingular>		spSingular(it->second.pdispVal);
			if (spSingular->PutDate(nDate)) propagate_error;
		}		
	}

	//  Calls the Refresh method for every object in the collection.
	void Refresh(/*[in, defaultvalue(0L)]*/ VARIANT_BOOL Recursive)
	{
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComVariant	vObject(it->second);
			ATLASSERT(vObject.vt == VT_DISPATCH);
			CComPtr<ISingular> spSingular = dynamic_cast<ISingular*>(vObject.pdispVal);
			if (spSingular->Refresh(Recursive)){
				propagate_error(spSingular);
			}
		}
	}

	//	Calls the Reset method for every object in the collection.
	void Reset(void)
	{
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComVariant	vObject(it->second);
			ATLASSERT(vObject.vt == VT_DISPATCH);
			CComPtr<ISingular> spSingular = dynamic_cast<ISingular*>(vObject.pdispVal);
			if (spSingular->Reset()){
				propagate_error(spSingular);
			}
		}
	}

	//	Calls the Shift method for every object in the collection.
	void Shift(double Amount)
	{
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComVariant	vObject(it->second);
			ATLASSERT(vObject.vt == VT_DISPATCH);
			CComPtr<ISingular> spSingular = dynamic_cast<ISingular*>(vObject.pdispVal);
			if (spSingular->Shift(Amount)){
				propagate_error(spSingular);
			}
		}
	}
	
	//	Calls the Shift method for every object in the collection.
	void Stick(void)
	{
		for (ContainerType::const_iterator it = m_pcoll->begin(); it != m_pcoll->end(); it++){
			CComVariant	vObject(it->second);
			ATLASSERT(vObject.vt == VT_DISPATCH);
			CComPtr<ISingular> spSingular = dynamic_cast<ISingular*>(vObject.pdispVal);
			if (spSingular->Stick()){
				propagate_error(spSingular);
			}
		}
	}

	// Emulates sp_user_update_pl for an object class
	template <class ICollection>
	static HRESULT UpdatePL_NTFS(long nDate, HRESULT (*load)(const std::string&, DataSourceEnum, long, CComPtr<ICollection>&))
	{
		CComPtr<ICollection>	spCollection;
		HRESULT					hr;
		long					nItems = 0;
		if (!load("", Last, nDate, spCollection)){
			spCollection->get_Count(&nItems);
			for (long nItem = 1; nItem <= nItems; nItem++){
				CComPtr<ISingular>	spSingular;
				if (hr = spCollection->get_Item(CComVariant(nItem), &spSingular)) return hr;
				if (hr = spSingular->put_DataSource(PL)) return hr;
				if (hr = spSingular->Save(NULL)) return hr;
			}
		}
		return S_OK;
	}

protected:
	//	Adds all the files that define Sirius objects in an input directory to a string-bool map.
	static HRESULT GetObjectFilesInDirectory(const std::string& szDirectory, std::map<std::string, bool>* pmap)
	{
		HRESULT								hr;
		char								szCurrentDirectory[MAX_PATH];	// active directory when this function is called
		WIN32_FIND_DATA						fd;
		HANDLE								h = NULL;						// return of FindFirstFile

		if (!::GetCurrentDirectory(sizeof(szCurrentDirectory), szCurrentDirectory)) return E_FAIL;
		if (!::SetCurrentDirectory(szDirectory.c_str())) throw "Invalid path '" + szDirectory + "'";
		hr = ((h = ::FindFirstFile("*", &fd)) == INVALID_HANDLE_VALUE) ? E_FAIL : S_OK;
		while (!hr){
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)){
				estring						szFileName(fd.cFileName);
				szFileName.resize(__max(szFileName.length() - 4, 0));		// remove extension
				(*pmap)[szFileName] = true;
			}
			if (!::FindNextFile(h, &fd)) break;
		}
		if (h != INVALID_HANDLE_VALUE && (!::FindClose(h))) ATLASSERT(false);
		if (!::SetCurrentDirectory(szCurrentDirectory)) ATLASSERT(false);
		return hr;
	}
	
	static void ImplementLoad_NTFS(estring szFileName, DataSourceEnum ds, CComPtr<ISingular>& spObject)
	// szFileName - could contain a wildcard in the file name
	{
		std::fstream						file;
		xmlstreamer							ssXML;

		if (szFileName.find('*') != szFileName.npos){
			char							szCurrentDirectory[MAX_PATH];		// active directory when this function is called			
			WIN32_FIND_DATA					fd;
			HANDLE							h = NULL;							// return of FindFirstFile
			std::string						szDirectory;
			
			std::string::size_type st = szFileName.rfind('\\');
			if (st != szFileName.npos){
				szDirectory = szFileName.left(st);
				szFileName = szFileName.mid(st + 1);
			}

			if (!::GetCurrentDirectory(sizeof(szCurrentDirectory), szCurrentDirectory)) throw "Unhandled exception in CComObjectCollectionFunctions::ImplementLoad_NTFS";
			if (!::SetCurrentDirectory(szDirectory.c_str())) throw "Invalid path '" + szDirectory + "'";
			h = ::FindFirstFile(szFileName.data(), &fd);
			if (!::SetCurrentDirectory(szCurrentDirectory)) ATLASSERT(false);
			if (h == INVALID_HANDLE_VALUE) throw "Invalid or missing file '" + szFileName + "'";				
			szFileName.assign(szDirectory + '\\' + fd.cFileName);
		}
		
		file.open(szFileName.c_str(), std::ios::in);
		if (file.fail()) throw "Invalid or missing file '" + szFileName + "'";
		while (!file.fail() && !file.eof()){
			std::string szLine;
			std::getline(file, szLine);
			ssXML << szLine;
		}	
		if (CXmlStreamer::GetObject((char*)ssXML, ds, spObject)) propagate_error;
	}

	//	Implmemets saving of an object from the local file system.
	static HRESULT ImplementSave_NTFS(const xmlstreamer& ssXML, const std::string& szFileName, BSTR* pVal)
	//	pVal - returned, nullable
	{
		std::fstream						file;
		estring								szDirectory;
		std::string::size_type				nPos;

		if ((nPos = szFileName.rfind('\\')) != std::string::npos){
			szDirectory = szFileName.substr(0, nPos);
			if (::GetFileAttributes(szDirectory.c_str()) == 0xFFFFFFFF){
				// The directory does not exist - we need to create it.
				// Unfortunately the ::CreateDirectory function does not
				// support creating a whole path in one go. Therefore we
				// have to create each sub-directory in turn.
				std::vector<std::string>	asz;
				szDirectory.Split("\\", &asz);
				if (asz.size() > 1){
					szDirectory = asz[0] + "\\";
					for (long n = 1; n < asz.size(); n++){
						szDirectory += asz[n];
						szDirectory += "\\";
						::CreateDirectory(szDirectory.c_str(), NULL);
					}
				}
			}
		}

		file.open(szFileName.c_str(), std::ios::out | std::ios::trunc);
		if (file.fail()) throw "Invalid file name '" + szFileName + "'";
		file << (char*)ssXML;
		file.close();
		if (pVal){
			return CComBSTR(szFileName.c_str()).CopyTo(pVal);
		} else {
			return S_OK;
		}
	}
};


#endif