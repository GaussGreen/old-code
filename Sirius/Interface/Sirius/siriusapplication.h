//	siriusapplication.h : Declaration of CSiriusApplication
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SIRIUSAPPLICATION_H_
#define __SIRIUSAPPLICATION_H_

class CAssets;
#include "resource.h"
#include "objectmanager.h"
#include "SiriusCP.h"

class ATL_NO_VTABLE CSiriusApplication : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSiriusApplication, &CLSID_SiriusApplication>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISiriusApplication, &IID_ISiriusApplication, &LIBID_Sirius>,
	public CProxy_ISiriusApplicationEvents<CSiriusApplication>,
	public IConnectionPointContainerImpl<CSiriusApplication>
{	
protected:
	// collections - all collection objects that are exposed to Excel via handles must be declared using the DECLARE_COLLECTION macro	
	DECLARE_COLLECTION(ARPropVolDataCollection)
	DECLARE_COLLECTION(ARPropFitVolDataCollection)
	DECLARE_COLLECTION(Arrays)
	DECLARE_COLLECTION(Dates)
	DECLARE_COLLECTION(DateSchedules)
	DECLARE_COLLECTION(Deals)
	DECLARE_COLLECTION(HermiteVolDataCollection)
	DECLARE_COLLECTION(HullAndWhiteCollection)
	DECLARE_COLLECTION(HullVolDataCollection)	
	DECLARE_COLLECTION(Interpolators)
	DECLARE_COLLECTION(JumpWingVolDataCollection)
	DECLARE_COLLECTION(JumpWingFitVolDataCollection)
	DECLARE_COLLECTION(Matrices)
	DECLARE_COLLECTION(MultiplierVolDataCollection)
	DECLARE_COLLECTION(MonteCarloCollection)
	DECLARE_COLLECTION(ParameterLists)
	DECLARE_COLLECTION(Parameters)
	DECLARE_COLLECTION(Products)
	DECLARE_COLLECTION(Positions)
	DECLARE_COLLECTION(RamFitVolDataCollection)
	DECLARE_COLLECTION(RamVolDataCollection)
	DECLARE_COLLECTION(Results)
	DECLARE_COLLECTION(Scenarios)
	DECLARE_COLLECTION(StrikesCollection)
	DECLARE_COLLECTION(SVLVolDataCollection)
	DECLARE_COLLECTION(SVLFitVolDataCollection)	
	DECLARE_COLLECTION(VolDataCollection)
	DECLARE_COLLECTION(PdeCollection)
		
public:			
	CSiriusApplication(){}

	DECLARE_CLASSFACTORY_SAFEPROC_SINGLETON(CSiriusApplication)		// uses custom singleton implementation	
	DECLARE_REGISTRY_RESOURCEID(IDR_SIRIUSAPPLICATION)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSiriusApplication)
		COM_INTERFACE_ENTRY(ISiriusApplication)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY_IMPL(IConnectionPointContainer)
	END_COM_MAP()
	
	// Event maps
	BEGIN_CONNECTION_POINT_MAP(CSiriusApplication)
		CONNECTION_POINT_ENTRY(__DIID__ISiriusApplicationEvents)
	END_CONNECTION_POINT_MAP()
	
	STDMETHOD(Clear)(void);
	STDMETHOD(CreateObject)(/*[in]*/ BSTR ObjectName, /*[in, optional]*/ VARIANT Arg1, /*[in, optional]*/ VARIANT Arg2, /*[in, optional]*/ VARIANT Arg3, /*[in, optional]*/ VARIANT Arg4, /*[in, optional]*/ VARIANT Arg5, /*[in, optional]*/ VARIANT Arg6, /*[in, optional]*/ VARIANT Arg7, /*[in, optional]*/ VARIANT Arg8, /*[in, optional]*/ VARIANT Arg9, /*[in, optional]*/ VARIANT Arg10, /*[in, optional]*/ VARIANT Arg11, /*[in, optional]*/ VARIANT Arg12, /*[in, optional]*/ VARIANT Arg13, /*[in, optional]*/ VARIANT Arg14, /*[in, optional]*/ VARIANT Arg15, /*[in, optional]*/ VARIANT Arg16, /*[in, optional]*/ VARIANT Arg17, /*[in, optional]*/ VARIANT Arg18, /*[in, optional]*/ VARIANT Arg19, /*[in, optional]*/ VARIANT Arg20, /*[in, optional]*/ VARIANT Arg21, /*[in, optional]*/ VARIANT Arg22, /*[in, optional]*/ VARIANT Arg23, /*[in, optional]*/ VARIANT Arg24, /*[in, optional]*/ VARIANT Arg25, /*[in, optional]*/ VARIANT Arg26, /*[in, optional]*/ VARIANT Arg27, /*[in, optional]*/ VARIANT Arg28, /*[in, optional]*/ VARIANT Arg29, /*[in, optional]*/ VARIANT Arg30, /*[in, optional]*/ VARIANT Arg31, /*[in, optional]*/ VARIANT Arg32, /*[out, retval]*/ IDispatch** pVal);
	STDMETHOD(get_MarketData)(IMarketData** pVal);
	STDMETHOD(GetObject)(/*[in]*/ VARIANT Handle, /*[in]*/ BSTR ObjectName, /*[out, retval]*/ IDispatch** pVal);
	STDMETHOD(GetObjects)(/*[in]*/ VARIANT Handles, /*[in]*/ BSTR out, /*[ObjectName, retval]*/ VARIANT* pObjectsArray);
	STDMETHOD(GetObjectTypes)(/*[out, retval]*/ VARIANT* pVal);			
	STDMETHOD(GetProductTypes)(/*[out, retval]*/ VARIANT* pVal);	
	STDMETHOD(InsertObject)(/*[in]*/ IDispatch* Object, /*[in, defaultvalue(TRUE)]*/ VARIANT_BOOL ExternalHandle, /*[out, retval]*/ BSTR* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);	
	STDMETHOD(Load)(/*[in]*/ BSTR ObjectName, /*[in]*/ BSTR Identifier, /*[in, optional]*/ VARIANT DataSourceOpt, /*[in, optional]*/ VARIANT DateOpt, /*[out, retval]*/ IDispatch** pVal);
	STDMETHOD(put_MarketData)(IMarketData* newVal);
	STDMETHOD(Save)(/*[in]*/ BSTR Handle);	
	STDMETHOD(StartErrors)(void);
	STDMETHOD(AddError)(/*[in]*/ BSTR Error);
	STDMETHOD(StopErrors)(void);
	STDMETHOD(get_Errors)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_FirstError)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_LastError)(/*[out, retval]*/ BSTR* pVal);	
	STDMETHOD(get_Setup)(/*[out, retval]*/ ISetup** pVal);		
	STDMETHOD(get_Locations)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_FunctionsDisabled)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(put_FunctionsDisabled)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(GetVersion)(/*[in]*/ BSTR ModuleName, /*[out, retval]*/ BSTR* pVersion);
	STDMETHOD(EnumToString)(/*[in]*/ BSTR EnumName, /*[in]*/ long Value, /*[out, retval]*/ BSTR* psName);
	STDMETHOD(StringToEnum)(/*[in]*/ BSTR EnumName, /*[in]*/ BSTR Name, /*[out, retval]*/ long* pnValue);
	STDMETHOD(Enums)(/*[in]*/ BSTR EnumName, /*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(UpdatePL)(/*[in]*/ DATE CopyDateOpt);
	STDMETHOD(CreateObjectFromXml)(/*[in]*/ BSTR Xml, /*[out, retval]*/ IDispatch** pVal);
	
	HRESULT								FinalConstruct(void);	
	void								FinalRelease(void);			
	HRESULT								GetVersion(const std::string& szModuleName, std::string* psz) const;
	static HRESULT						GetAssets(CComPtr<IAssets>& sp);
	static CAssets*						GetAssets(void);
	double								GetCorrelation(const std::string& szAsset1, const std::string& szAsset2, const std::string& szDataSourceOpt, long nDate) const;
	static CComPtr<ICorrelationMatrix>	GetCorrelationMatrix(DataSourceEnum ds, long nDate, bool bCreateIfNecessary);
	static MlEqCorrelationMatrixHandle	GetCorrelationMatrixHandle(DataSourceEnum ds, long nDate, bool bCreateIfNecessary);
	HRESULT								GetCreateParameterNames(CComDispatchDriverEx& ddObject, const std::vector<CParameterMap>& vpm, std::vector<std::string>* pvector, std::map<std::string, long>* pmap, std::map<std::string, long>* pmapCase) const;
	bool								GetFunctionsDisabled(void) const;
	_ConnectionPtr						GetLocations(FileSystemEnum fs, const CSiriusComModule::connection_parameters* pcp, std::vector<estring>* pasz, std::map<estring, estring>* pmap, bool bReturnLowerCase, bool bThrow) const;
	CObjectManager&						GetObjectManager(void);
	HRESULT								GetOneParameter(CComDispatchDriverEx& ddObject, std::string szParameter, CComVariant* pv);
	HRESULT								GetProperties(CComPtr<IDispatch> spObject, VARIANT* pResult);
	void								GetLastError(VARIANT* pResult);
	static HRESULT						GetZeroCurves(CComPtr<IZeroCurves>& sp);
	static void							InsertCorrelation(const std::string& sz1, const std::string& sz2, DataSourceEnum ds, long nDate, bool bThrow);
	HRESULT								LoadProductTypes(void);
	void								PutFunctionsDisabled(bool b);
				
protected:	
	HRESULT								CreateObject(const std::string& szObjectName, IDispatch** pVal, int nArgs, /*VARIANT*/...);	
	HRESULT								InstallStandardObjects(void);

	CObjectManager						m_ObjectManager;				// object manager
	CComPtr<IMarketData>				m_spMarketData;					// market data handle
	std::vector<estring>				m_vectorErrors;					// list of errors that have occurred
	bool								m_bErrorLogging;				// true if we are logging any errors that arise
	bool								m_bFunctionsDisabled;			// true if we disable all Sirius functions
};													

#endif
