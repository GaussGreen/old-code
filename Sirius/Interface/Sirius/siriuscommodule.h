// siriuscommodule.h: interface for the CSiriusComModule class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _SIRIUSCOMMODULE_H
#define _SIRIUSCOMMODULE_H

#pragma once

class										CComDispatchDriverEx;
class										estring;

static const IID							IID_IBwlEngine				= { 0x875db756, 0x6cc6, 0x11d6, { 0xb9, 0xbc, 0x0, 0xb0, 0xd0, 0xab, 0xf2, 0xc7 } };
static const CLSID							CLSID_BwlEngine				= { 0x875db757, 0x6cc6, 0x11d6, { 0xb9, 0xbc, 0x0, 0xb0, 0xd0, 0xab, 0xf2, 0xc7 } };
static const CLSID							CLSID_SiriusDataManager     = { 0x3fc33402,	0x792e,	0x4afb,	{ 0x8d, 0xf4, 0x2b,0x68, 0x3d, 0x33, 0xe0, 0x7f } };
static const CLSID							CLSID_Taurus				= { 0xac180cbd, 0x1b19, 0x43fb, { 0x87, 0x6b, 0x8a,0xe4, 0xab, 0x8b, 0x1c, 0xc4 } };

typedef IDispatch							IBwlEngine;
typedef IDispatch							ISiriusDataManager;
typedef IDispatch							ITaurus;

class CSiriusApplication;
class CSiriusComModule : public CComModule
{
public:
	struct connection_parameters {
		std::string							m_szDSN;
		std::string							m_szUserName;
		std::string							m_szPassword;
	};
protected:
	struct state {
		DatabaseModeEnum					m_dbm;
		DataSourceEnum						m_dsDefault;
		std::string							m_szDataSourceDefault;			// In attempt to prevent string copies, I maintain the string form of m_dsDefault
		FileSystemEnum						m_fs;
		std::string							m_szBwlInterpolationRule;
		std::string							m_szBwlPassword;
		std::string							m_szBwlSpotScheduleSource;
		std::string							m_szBwlSQLServerDSN;
		std::string							m_szBwlUserName;		
		bool								m_bDisplayProductOnly;
		bool								m_bDisplaySiriusStatus;
		bool								m_bEnablePublishing;		
		std::string							m_szFileSystemRoot;
		std::string							m_szLocation;
		std::string							m_szLocationFirstTry;
		bool								m_bLocationFirstTryEnabled;
		std::string							m_szProductTypesFileName;		
		std::string							m_szSiriusPassword;
		std::string							m_szSiriusSQLServerDSN;
		std::string							m_szSiriusUserName;
		bool								m_bUseBwlForSpotSchedule;
		bool								m_bUseYesterdaysClose;				
	} m_state, m_state_init;

	enum BwlUsable {
		Unknown = 0,
		Usable = 1,
		Unusable = 2
	};	

public:
	std::string								GetBwlInterpolationRule(void) const;
	std::string								GetBwlPassword(void) const;	
	std::string								GetBwlSpotScheduleSource(void) const;
	std::string								GetBwlSQLServerDSN(void) const;
	std::string								GetBwlUserName(void) const;			
	DatabaseModeEnum						GetDatabaseMode(void) const;
	DataSourceEnum							GetDefaultDataSource(void) const;
	const std::string&						GetDefaultDataSourceStr(void) const;
	bool									GetDisplayProductOnly(void) const;
	bool									GetDisplaySiriusStatus(void) const;
	bool									GetEnablePublishing(void) const;		
	FileSystemEnum							GetFileSystem(void) const;
	std::string								GetFileSystemRoot(void) const;
	std::string								GetLocationFirstTry(void) const;
	bool									GetLocationFirstTryEnabled(void) const;
	std::string								GetLocation(void) const;
	std::string								GetProductTypesFileName(void) const;
	std::string								GetSiriusPassword(void) const;
	std::string								GetSiriusSQLServerDSN(void) const;		
	std::string								GetSiriusUserName(void) const;
	bool									GetUseBwlForSpotSchedule(void) const;
	static bool								GetUseYesterdaysClose(void);					// we take the address of this function in places - hence it is static
			
	void									PutBwlInterpolationRule(const std::string& sz);
	void									PutBwlPassword(const std::string& sz);	
	void									PutBwlSpotScheduleSource(const std::string& sz);
	void									PutBwlSQLServerDSN(const std::string& sz);
	void									PutBwlUserName(const std::string& sz);	
	void									PutDatabaseMode(DatabaseModeEnum dbm);
	void									PutDefaultDataSource(DataSourceEnum dsDefault);
	void									PutDisplayProductOnly(bool b);
	void									PutDisplaySiriusStatus(bool b);
	void									PutEnablePublishing(bool b);	
	void									PutFileSystem(FileSystemEnum fs);
	void									PutFileSystemRoot(const std::string& sz);
	void									PutLocation(std::string* psz, const connection_parameters* pcp, FileSystemEnum* pfs);
	void									PutLocationFirstTry(std::string* psz, const connection_parameters* pcp, FileSystemEnum* pfs);
	void									PutLocationFirstTryEnabled(bool b);
	void									PutProductTypesFileName(const std::string& sz);
	void									PutSiriusSQLServerDSN(const std::string& sz);
	void									PutSiriusPassword(const std::string& sz);
	void									PutSiriusUserName(const std::string& sz);
	void									PutUseBwlForSpotSchedule(bool b);
	void									PutUseYesterdaysClose(bool b);

public:
	CSiriusComModule(void);
	virtual									~CSiriusComModule(void){}
	HRESULT									AddError(const estring& sz) const;
	HRESULT									DeclareExcel(const VARIANT& ExcelApplicationInstance);
	HRESULT									DisplayExplorer(void);
	HRESULT									GetConnection(_ConnectionPtr& spConnection, const connection_parameters* pcp = NULL) const;
	bool									GetEnableLoad(void) const;
	HRESULT									GetExcel(CComDispatchDriverEx& ddExcel) const;
	HRESULT									GetMarketDataRoot(std::string* psz) const;
	std::string								GetName(void) const;
	bool									GetNTAuthentication(const std::string& szDSN) const;
	CComVariant								GetObjectSignature(const std::string& szObjectName, const std::string& szIdentifier, const std::string& szDataSource, long nDate) const;
	CSiriusApplication*						_GetSiriusApplication(void) const;
	CComPtr<ISiriusApplication>				GetSiriusApplication(void) const;
	HRESULT									GetSiriusCaption(CComVariant* pv) const;
	static bool								IsReutersRunning(void);
	HRESULT									Load(const std::string& szObjectName, const std::string& szIdentifier, const std::string& szDataSource, DATE date, bool bAllowCloning, CComPtr<IDispatch>& spObject);
	HRESULT									Load(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, bool bAllowCloning, CComPtr<IDispatch>& spObject);
	double									LoadBwlSpot(const std::string& szCode, const std::string& szCodeType, const std::string& szDatabase, const std::string& szDataSource, const std::string& szInterpolation);
	HRESULT									LoadFromBwl(CLSID clsid, const std::string& szIdentifier, CComPtr<IDispatch>& spObject);
	void									MembersToRegistry(void);
	void									ResetAllConnections(void);
	void									Revert(void);	
	void									Term(void);

protected:	
	void									CheckBwlConnection(void);
	static void								Crypt(TCHAR* inp, DWORD inplen, const TCHAR* key = "", DWORD keylen = 0);
	void									GetBwlSpotScheduleSource(std::vector<std::string>* pasz) const;
	HRESULT									GetExcel(CComPtr<IDispatch>& spExcel) const;
	HRESULT									GetSiriusCaption(std::string* psz) const;
	HRESULT									_Load(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, CComPtr<IDispatch>& spObject);
	HRESULT									_LoadAsset(const std::string& szIdentifier, const std::string& szDataSource, DATE date, bool bAllowCloning, CComPtr<IDispatch>& spAsset);
	HRESULT									_LoadDirect(CLSID clsid, const std::string& szIdentifier, DataSourceEnum ds, DATE date, CComPtr<IDispatch>& spObject) const;
	HRESULT									_LoadFromServer(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, CComPtr<IDispatch>& spObject);
	HRESULT									_LoadFromTaurus(CLSID clsid, const std::string& szIdentifier, const std::string& szDataSource, DATE date, CComPtr<IDispatch>& spObject);
	static HRESULT							MemberToRegistry(HKEY hKey, const std::string& szKey, std::string szValue, bool bSilent, bool bEncrypt);
	static HRESULT							RegistryToMember(HKEY hKey, const std::string& szName, std::string* psz, bool bEncrypted);
	void									RegistryToMembers(state* pstate, bool bRefresh) const;
	void									ResetBwlDatabaseConnection(void);	
	void									ResetSiriusDatabaseConnection(void);
	HRESULT									TaurusAttach(bool bAlwaysAttach);	
	std::string								ValidateLocation(const std::string& sz, const connection_parameters* pcp, FileSystemEnum* pfs) const;
	
public:
	static const char						s_szODBCDriverName[];
	HFONT									m_hFont;
	HIMAGELIST								m_hImageList;		
	
protected:
	static const char						s_sz_registry_database_mode[];
	static const char						s_sz_registry_path[];
	static const char						s_sz_registry_file_system[];
	static const char						s_sz_registry_sirius_sql_server[];
	static const char						s_sz_registry_file_root[];
	static const char						s_sz_registry_product_types[];
	static const char						s_sz_registry_sirius_user_name[];
	static const char						s_sz_registry_sirius_password[];
	static const char						s_sz_registry_product_only[];
	static const char						s_sz_registry_enable_publishing[];
	static const char						s_sz_registry_use_bwl_for_spot_schedule[];
	static const char						s_sz_registry_config_file_path[];
	static const char						s_sz_registry_use_yesterday_close[];
	static const char						s_sz_registry_display_sirius_status[];
	static const char						s_sz_registry_location[];
	static const char						s_sz_registry_location_first_try[];
	static const char						s_sz_registry_location_first_try_enabled[];
	static const char						s_sz_default_data_source[];
	static const char						s_sz_registry_bwl_sql_server[];
	static const char						s_sz_registry_bwl_user_name[];
	static const char						s_sz_registry_bwl_password[];
	static const char						s_sz_registry_bwl_interpolation[];
	static const char						s_sz_registry_bwl_spot_schedule_source[];

	BwlUsable								m_BwlUsable;							// true if we can use bwl database directly
	
	CComPtr<ISiriusDataManager>				m_spSiriusServer;						// ToDo - remove this and all _LoadFromServer in this file and the corresponding cpp.
	CComPtr<ITaurus>						m_spTaurus;

	CComPtr<IBwlEngine>						m_spBwlEngine;							// direct bwl connection		
	CComPtr<IDispatch>						m_spExcel;								// interface for Excel.Application (required for functions interacting with Excel)
	CComPtr<IDispatch>						m_spExplorer;							// Sirius Explorer pointer

private:
	mutable CComPtr<ISiriusApplication>		m_spApplication;
	mutable _ConnectionPtr					m_spConnection;							// Sirius database connection
};

#endif
	