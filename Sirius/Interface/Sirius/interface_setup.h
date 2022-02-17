//	interface_setup.h : Declaration of the CSetup
//
//				        This class is effectively the COM interface to 
//					    the setup dialog.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SETUP_H_
#define __SETUP_H_

#include "resource.h"

class ATL_NO_VTABLE CSetup : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSetup, &CLSID_Setup>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISetup, &IID_ISetup, &LIBID_Sirius>
{
public:
	CSetup(){}	
	HRESULT								FinalConstruct();
	HRESULT								FinalRelease();

	DECLARE_CLASSFACTORY_SINGLETON(CSetup)
	DECLARE_REGISTRY_RESOURCEID(IDR_SETUP)
	DECLARE_NOT_AGGREGATABLE(CSetup)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSetup)
		COM_INTERFACE_ENTRY(ISetup)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

public:	
	STDMETHOD(get_DisplaySiriusStatus)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(put_DisplaySiriusStatus)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(get_UseYesterdaysClose)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(put_UseYesterdaysClose)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(get_EnablePublishing)(/*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(put_EnablePublishing)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(get_FileSystem)(/*[out, retval]*/ FileSystemEnum* pVal);
	STDMETHOD(put_FileSystem)(/*[in]*/ FileSystemEnum newVal);
	STDMETHOD(get_ProductTypesFile)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_ProductTypesFile)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_SiriusDataSourceName)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_SiriusDataSourceName)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_FileSystemRoot)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_FileSystemRoot)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_ConnectionMode)(/*[out, retval]*/ DatabaseModeEnum* pVal);
	STDMETHOD(put_ConnectionMode)(/*[in]*/ DatabaseModeEnum newVal);
	STDMETHOD(get_UseBwlForSpotSchedule)(/*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(put_UseBwlForSpotSchedule)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(ShowSetupDialog)();
	STDMETHOD(Bounce)(void);
	STDMETHOD(get_Location)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_Location)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_LocationFirstTry)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_LocationFirstTry)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_LocationFirstTryEnabled)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(put_LocationFirstTryEnabled)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(get_DefaultDataSource)(/*[out, retval]*/ DataSourceEnum* pVal);
	STDMETHOD(put_DefaultDataSource)(/*[in]*/ DataSourceEnum newVal);
	STDMETHOD(get_SiriusDatabaseUserName)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_SiriusDatabaseUserName)(/*[in]*/ BSTR newVal);
	STDMETHOD(put_SiriusDatabasePassword)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_BwlDataSourceName)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_BwlDataSourceName)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_BwlDatabaseUserName)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_BwlDatabaseUserName)(/*[in]*/ BSTR newVal);
	STDMETHOD(put_BwlDatabasePassword)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_BwlInterpolationRule)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_BwlInterpolationRule)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_BwlSpotScheduleSource)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(put_BwlSpotScheduleSource)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_RevertOnDestruction)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(put_RevertOnDestruction)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(SaveToRegistry)(void);
	STDMETHOD(LoadFromRegistry)(void);

protected:
	bool								m_bRevertOnExit;
};

#endif
