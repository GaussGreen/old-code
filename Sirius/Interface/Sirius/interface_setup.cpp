//	interface_setup.cpp : Implementation of CSetup
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Sirius.h"
#include "interface_setup.h"
#include "dlgsetup.h"

STDMETHODIMP CSetup::Bounce(void)
{
	begin_function
	_Module.ResetAllConnections();
	end_function
}

HRESULT CSetup::FinalConstruct()
{
	 m_bRevertOnExit = false;
	 return S_OK;
}

HRESULT CSetup::FinalRelease()
{
	if (m_bRevertOnExit){
		_Module.Revert();
	}
	return S_OK;
}

STDMETHODIMP CSetup::get_BwlDatabaseUserName(BSTR *pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetBwlUserName(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_BwlDataSourceName(BSTR *pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetBwlSQLServerDSN(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_BwlInterpolationRule(BSTR *pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetBwlInterpolationRule(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_BwlSpotScheduleSource(BSTR* pVal)
{	
	begin_function
	return estring::GetBSTR(_Module.GetBwlSpotScheduleSource(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_ConnectionMode(DatabaseModeEnum* pVal)
{
	begin_function
	*pVal = _Module.GetDatabaseMode();
	end_function
}

STDMETHODIMP CSetup::get_DefaultDataSource(DataSourceEnum* pVal)
{
	begin_function
	*pVal = _Module.GetDefaultDataSource();
	end_function
}

STDMETHODIMP CSetup::get_DisplaySiriusStatus(VARIANT_BOOL *pVal)
{
	begin_function
	*pVal = _Module.GetDisplaySiriusStatus();
	end_function
}

STDMETHODIMP CSetup::get_EnablePublishing(VARIANT_BOOL* pVal)
{	
	begin_function
	*pVal = _Module.GetEnablePublishing();	
	end_function
}

STDMETHODIMP CSetup::get_FileSystem(FileSystemEnum* pVal)
{
	begin_function
	*pVal = _Module.GetFileSystem();
	end_function
}

STDMETHODIMP CSetup::get_FileSystemRoot(BSTR* pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetFileSystemRoot(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_Location(BSTR* pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetLocation(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_LocationFirstTry(BSTR* pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetLocationFirstTry(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_LocationFirstTryEnabled(VARIANT_BOOL *pVal)
{
	begin_function
	*pVal = _Module.GetLocationFirstTryEnabled();
	end_function
}

STDMETHODIMP CSetup::get_ProductTypesFile(BSTR* pVal)
{	
	begin_function
	return estring::GetBSTR(_Module.GetProductTypesFileName(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_RevertOnDestruction(VARIANT_BOOL *pVal)
{
	begin_function
	*pVal = m_bRevertOnExit ? VARIANT_TRUE : VARIANT_FALSE;
	end_function
}

STDMETHODIMP CSetup::get_SiriusDatabaseUserName(BSTR* pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetSiriusUserName(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_SiriusDataSourceName(BSTR* pVal)
{
	begin_function
	return estring::GetBSTR(_Module.GetSiriusSQLServerDSN(), pVal);
	end_function
}

STDMETHODIMP CSetup::get_UseBwlForSpotSchedule(VARIANT_BOOL* pVal)
{	
	begin_function
	*pVal = _Module.GetUseBwlForSpotSchedule();	
	end_function
}

STDMETHODIMP CSetup::get_UseYesterdaysClose(VARIANT_BOOL *pVal)
{
	begin_function
	*pVal = _Module.GetUseYesterdaysClose();
	end_function
}

STDMETHODIMP CSetup::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISetup };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CSetup::LoadFromRegistry(void)
{
	begin_function
	_Module.Revert();
	end_function
}

STDMETHODIMP CSetup::put_BwlDatabaseUserName(BSTR newVal)
{
	begin_function
	_Module.PutBwlUserName(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_BwlDatabasePassword(BSTR newVal)
{
	begin_function
	_Module.PutBwlPassword(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_BwlDataSourceName(BSTR newVal)
{
	begin_function
	_Module.PutBwlSQLServerDSN(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_BwlInterpolationRule(BSTR newVal)
{
	begin_function
	_Module.PutBwlInterpolationRule(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_BwlSpotScheduleSource(BSTR newVal)
{
	begin_function
	_Module.PutBwlSpotScheduleSource(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_ConnectionMode(DatabaseModeEnum newVal)
{	
	begin_function
	_Module.PutDatabaseMode(newVal);
	end_function
}

STDMETHODIMP CSetup::put_DisplaySiriusStatus(VARIANT_BOOL newVal)
{
	begin_function
	_Module.PutDisplaySiriusStatus(newVal ? true : false);
	end_function
}

STDMETHODIMP CSetup::put_DefaultDataSource(DataSourceEnum newVal)
{
	begin_function
	_Module.PutDefaultDataSource(newVal);
	end_function
}

STDMETHODIMP CSetup::put_EnablePublishing(VARIANT_BOOL newVal)
{
	begin_function
	_Module.PutEnablePublishing(newVal ? true : false);
	end_function
}

STDMETHODIMP CSetup::put_FileSystem(FileSystemEnum newVal)
{
	begin_function
	_Module.PutFileSystem(newVal);
	end_function
}

STDMETHODIMP CSetup::put_FileSystemRoot(BSTR newVal)
{
	begin_function
	_Module.PutFileSystemRoot(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_Location(BSTR newVal)
{
	begin_function
	_Module.PutLocation(&estring(newVal), NULL, NULL);
	end_function
}

STDMETHODIMP CSetup::put_LocationFirstTry(BSTR newVal)
{
	begin_function
	_Module.PutLocationFirstTry(&estring(newVal), NULL, NULL);
	end_function
}

STDMETHODIMP CSetup::put_LocationFirstTryEnabled(VARIANT_BOOL newVal)
{
	begin_function
	_Module.PutLocationFirstTryEnabled(newVal ? true : false);
	end_function
}

STDMETHODIMP CSetup::put_ProductTypesFile(BSTR newVal)
{
	begin_function
	_Module.PutProductTypesFileName(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_RevertOnDestruction(VARIANT_BOOL newVal)
{
	begin_function
	m_bRevertOnExit = newVal ? true : false;
	end_function
}

STDMETHODIMP CSetup::put_SiriusDataSourceName(BSTR newVal)
{
	begin_function
	_Module.PutSiriusSQLServerDSN(estring(newVal));	
	end_function
}

STDMETHODIMP CSetup::put_SiriusDatabasePassword(BSTR newVal)
{
	begin_function
	_Module.PutSiriusPassword(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_SiriusDatabaseUserName(BSTR newVal)
{
	begin_function
	_Module.PutSiriusUserName(estring(newVal));
	end_function
}

STDMETHODIMP CSetup::put_UseBwlForSpotSchedule(VARIANT_BOOL newVal)
{
	begin_function
	_Module.PutUseBwlForSpotSchedule(newVal ? true : false);
	end_function
}

STDMETHODIMP CSetup::put_UseYesterdaysClose(VARIANT_BOOL newVal)
{
	begin_function
	_Module.PutUseYesterdaysClose(newVal ? true : false);
	end_function
}

STDMETHODIMP CSetup::SaveToRegistry(void)
{
	begin_function
	_Module.MembersToRegistry();
	end_function
}

STDMETHODIMP CSetup::ShowSetupDialog()
{
	begin_function
	CDlgSetup().DoModal();
	end_function
}
