//	dlgsetup.h : Declaration of the CDlgSetup
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DLGSETUP_H_
#define __DLGSETUP_H_

#include "resource.h"

class CDlgSetup : public CAxDialogImpl<CDlgSetup>
{
public:
	BEGIN_MSG_MAP(CDlgSetup)
		MESSAGE_HANDLER(WM_INITDIALOG, OnInitDialog)
		COMMAND_ID_HANDLER(IDOK, OnOK)
		COMMAND_ID_HANDLER(IDCANCEL, OnCancel)
		COMMAND_HANDLER(IDC_RADIO_LOCAL, BN_CLICKED, OnClickedRadioLocal)
		COMMAND_HANDLER(IDC_RADIO_SQL, BN_CLICKED, OnClickedRadioSQL)		
		COMMAND_HANDLER(IDC_BROWSE, BN_CLICKED, OnClickedBrowse)
		COMMAND_HANDLER(IDC_COMBO_DSN, CBN_SELCHANGE, OnSelchangeComboDSN)
		COMMAND_HANDLER(IDC_BROWSE_ROOT, BN_CLICKED, OnClickedBrowseRoot)
		COMMAND_HANDLER(IDC_BOUNCE_DATABASES, BN_CLICKED, OnClickedBounceDatabases)
		COMMAND_HANDLER(IDC_SAVE, BN_CLICKED, OnClickedSave)
		COMMAND_HANDLER(IDC_REVERT, BN_CLICKED, OnClickedRevert)
		NOTIFY_HANDLER(IDC_TAB, TCN_SELCHANGE, OnSelChangeTab)
		COMMAND_HANDLER(IDC_COMBO_DSN_BWL, CBN_SELCHANGE, OnSelchangeComboDSNBwl)
		COMMAND_HANDLER(IDC_FIRST_TRY, BN_CLICKED, OnClickedLocationFirstTry)
		COMMAND_HANDLER(IDC_COMBO_LOCATION, CBN_DROPDOWN, OnDropDownLocation)
		COMMAND_HANDLER(IDC_COMBO_LOCATION_TRY, CBN_DROPDOWN, OnDropDownLocationTry)
	END_MSG_MAP()

	CDlgSetup(){}
	~CDlgSetup(){}	
	enum { IDD = IDD_DLGSETUP };
	
	LRESULT								OnClickedRevert(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnClickedSave(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);	
	LRESULT								OnChangeEditSave(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnClickedBrowse(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnClickedBrowseRoot(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnClickedRadioLocal(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnClickedRadioSQL(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT								OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnSelchangeComboDSN(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnSelchangeComboDSNBwl(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnSelChangeTab(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);
	LRESULT								OnClickedBounceDatabases(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);	
	LRESULT								OnClickedLocationFirstTry(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnDropDownLocation(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnDropDownLocationTry(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
		
protected:	
	enum { Sirius = 0, SiriusDatabase = 1, BwlDatabase = 2 };
	static int							s_nCurSel;						// currently selected tab
	
	FileSystemEnum						m_fs;
	estring								m_szDSN;						// data source name selected
	estring								m_szFileRoot;					// file root path selected
	estring								m_szProductTypesFile;			// product type file name
	estring								m_szUserName;					// login ID (for SQL server authentication)
	estring								m_szPassword;					// user password (for SQL server authentication)
	estring								m_szLocation;					// location (used for market data access)
	estring								m_szLocationFirstTry;			// first location to try when loading
	bool								m_bLocationFirstTryEnabled;
	bool								m_bDisplaySiriusStatus;
	bool								m_bEnablePublishing;
	bool								m_bUseBwlForSpotSchedule;
	bool								m_bUseYesterdaysClose;
	bool								m_bFunctionsDisabled;
	DatabaseModeEnum					m_dbm;
	DataSourceEnum						m_dsDefault;

	estring								m_szBwlDSN;
	estring								m_szBwlUserName;
	estring								m_szBwlPassword;
	estring								m_szBwlInterpolation;
	estring								m_szBwlSpotScheduleSource;
	
		
	void								DialogToMembers(void);
	BOOL								EndDialog(int nRetCode);
	void								GetConnectionParameters(CSiriusComModule::connection_parameters* pcp, FileSystemEnum* pfs) const;
	UINT								GetDlgItemText(int nID, std::string* psz) const;
	void								GetUniversalPath(TCHAR* szPath, std::string* pszPathOut) const;
	DWORD								GetUniversalPath(TCHAR* szPath, DWORD* pdwBufferSize, std::string* pszPathOut) const;
	void								MembersToDialog(void);	
	void								MembersToModule(void);
	void								ModuleToMembers(void);
	void								PopulateDSNCombo(int nControlID, const std::string& szPrefix) const;
	void								PopulateLocationCombos(void);
	void								ShowControl(int nID, bool bShow);
};

#endif
