//	dlgsetup.cpp : Implementation of CDlgSetup
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "dlgsetup.h"
#include "siriusapplication.h"

/*static*/ int								CDlgSetup::s_nCurSel = CDlgSetup::Sirius;

//	update member variables with dialog contents
void CDlgSetup::DialogToMembers(void)
{			
	estring								szDatabaseMode;
	estring								szDefaultDataSource;
	
	m_fs = IsDlgButtonChecked(IDC_RADIO_SQL) ? fsSQLServer : fsNTFS;
	GetDlgItemText(IDC_COMBO_DSN, &m_szDSN);
	GetDlgItemText(IDC_EDIT_ROOT, &m_szFileRoot);
	GetDlgItemText(IDC_EDIT_PRODUCT_TYPES, &m_szProductTypesFile);
	GetDlgItemText(IDC_EDIT_USER, &m_szUserName);
	GetDlgItemText(IDC_EDIT_PASSWORD, &m_szPassword);
	GetDlgItemText(IDC_COMBO_LOCATION, &m_szLocation);
	GetDlgItemText(IDC_COMBO_LOCATION_TRY, &m_szLocationFirstTry);
	GetDlgItemText(IDC_COMBO_DSN_BWL, &m_szBwlDSN);
	GetDlgItemText(IDC_EDIT_USER_BWL, &m_szBwlUserName);
	GetDlgItemText(IDC_EDIT_PASSWORD_BWL, &m_szBwlPassword);
	GetDlgItemText(IDC_EDIT_INTERPOLATION_BWL, &m_szBwlInterpolation);
	GetDlgItemText(IDC_EDIT_SPOT_SCHEDULE_SOURCE_BWL, &m_szBwlSpotScheduleSource);
	m_bLocationFirstTryEnabled = IsDlgButtonChecked(IDC_FIRST_TRY) ? true : false;
	m_bDisplaySiriusStatus = IsDlgButtonChecked(IDC_DISPAY_SIRIUS_STATUS) ? true : false;
	m_bEnablePublishing = IsDlgButtonChecked(IDC_CHECK_PUBLISHING) ? true : false;
	m_bUseYesterdaysClose = IsDlgButtonChecked(IDC_YESTERDAY_CLOSE) ? true : false;
	m_bFunctionsDisabled = IsDlgButtonChecked(IDC_DISPLAY_FUNCTIONS) ? true : false;
	m_bUseBwlForSpotSchedule = IsDlgButtonChecked(IDC_USE_BWL) ? true : false;
	GetDlgItemText(IDC_COMBO_DATABASE_MODE, &szDatabaseMode);
	CEnumMap::GetEnum("DatabaseModeEnum", LIBID_Sirius, szDatabaseMode, &m_dbm);
	GetDlgItemText(IDC_COMBO_DATASOURCE, &szDefaultDataSource);
	CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szDefaultDataSource, &m_dsDefault);
}

BOOL CDlgSetup::EndDialog(int nRetCode)
{
	return CAxDialogImpl<CDlgSetup>::EndDialog(nRetCode);
}

void CDlgSetup::GetConnectionParameters(CSiriusComModule::connection_parameters* pcp, FileSystemEnum* pfs) const
{
	GetDlgItemText(IDC_COMBO_DSN, &pcp->m_szDSN);
	GetDlgItemText(IDC_EDIT_USER, &pcp->m_szUserName);
	GetDlgItemText(IDC_EDIT_PASSWORD, &pcp->m_szPassword);
	*pfs = IsDlgButtonChecked(IDC_RADIO_SQL) ? fsSQLServer : fsNTFS;
}

//	Wraps CAxDialogImpl<CDlgSetup>::GetDlgItemText
UINT CDlgSetup::GetDlgItemText(int nID, std::string* psz) const
{
	UINT								nRet;
	char								szValue[MAX_PATH];
		
	nRet = CAxDialogImpl<CDlgSetup>::GetDlgItemText(nID, szValue, sizeof(szValue));
	psz->assign(szValue);
	estring::trim(psz);
	return nRet;
}

void CDlgSetup::GetUniversalPath(TCHAR* szPath, std::string* pszPathOut) const
{
	DWORD								dwBufferSize;
	DWORD								dwRet;
			
	dwBufferSize = MAX_PATH;			// probably sufficient
	if ((dwRet = GetUniversalPath(szPath, &dwBufferSize, pszPathOut)) == ERROR_MORE_DATA){
		// Try again with the new buffer size
		if (GetUniversalPath(szPath, &dwBufferSize, pszPathOut) != NO_ERROR) ATLASSERT(false);
	} else if (dwRet == ERROR_NOT_CONNECTED){
		// This network connection does not exist
		pszPathOut->assign(szPath);
	} else if (dwRet != NO_ERROR){
		ATLASSERT(false);
		pszPathOut->assign(szPath);
	}
	estring::trim(pszPathOut);
}
DWORD CDlgSetup::GetUniversalPath(TCHAR* szPath, DWORD* pdwBufferSize, std::string* pszPathOut) const
{				
	TCHAR*								pszUniversalPath = new TCHAR[*pdwBufferSize];
	UNIVERSAL_NAME_INFO*				puni;			// universal name structure
	DWORD								dwRet;
	
	puni = (UNIVERSAL_NAME_INFO*)pszUniversalPath;	
	if ((dwRet = WNetGetUniversalName(szPath, UNIVERSAL_NAME_INFO_LEVEL, pszUniversalPath, pdwBufferSize)) == NO_ERROR) pszPathOut->assign(puni->lpUniversalName);
	delete pszUniversalPath;
	return dwRet;
}

//	update the dialog contents with the member variables
void CDlgSetup::MembersToDialog(void)
{	
	CheckDlgButton(m_fs == fsSQLServer ? IDC_RADIO_SQL : IDC_RADIO_LOCAL, BST_CHECKED);	
	CheckDlgButton(m_fs != fsSQLServer ? IDC_RADIO_SQL : IDC_RADIO_LOCAL, BST_UNCHECKED);	
	CheckDlgButton(IDC_CHECK_PUBLISHING, m_bEnablePublishing ? BST_CHECKED : BST_UNCHECKED);	// do this first before any text control Change messages can be sent
	CheckDlgButton(IDC_DISPAY_SIRIUS_STATUS, m_bDisplaySiriusStatus ? BST_CHECKED : BST_UNCHECKED);
	CheckDlgButton(IDC_USE_BWL, m_bUseBwlForSpotSchedule ? BST_CHECKED : BST_UNCHECKED);
	CheckDlgButton(IDC_YESTERDAY_CLOSE, m_bUseYesterdaysClose ? BST_CHECKED : BST_UNCHECKED);
	CheckDlgButton(IDC_DISPLAY_FUNCTIONS, m_bFunctionsDisabled ? BST_CHECKED : BST_UNCHECKED);
	CheckDlgButton(IDC_FIRST_TRY, m_bLocationFirstTryEnabled ? BST_CHECKED : BST_UNCHECKED);
		
	if (SendMessage(GetDlgItem(IDC_COMBO_DATABASE_MODE), CB_SELECTSTRING, 0, (LPARAM)CEnumMap::GetString("DatabaseModeEnum", LIBID_Sirius, m_dbm).c_str()) == CB_ERR){
		// select the first item
		SendMessage(GetDlgItem(IDC_COMBO_DATABASE_MODE), CB_SETCURSEL, 0, 0);
	}
	if (SendMessage(GetDlgItem(IDC_COMBO_DATASOURCE), CB_SELECTSTRING, 0, (LPARAM)CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_dsDefault).c_str()) == CB_ERR){
		// select the first item
		SendMessage(GetDlgItem(IDC_COMBO_DATASOURCE), CB_SETCURSEL, 0, 0);
	}
	SetDlgItemText(IDC_EDIT_PRODUCT_TYPES, m_szProductTypesFile.c_str());
	SetDlgItemText(IDC_EDIT_ROOT, m_szFileRoot.c_str());
	if (SendMessage(GetDlgItem(IDC_COMBO_DSN), CB_SELECTSTRING, 0, (LPARAM)m_szDSN.c_str()) == CB_ERR){
		// select the first item
		SendMessage(GetDlgItem(IDC_COMBO_DSN), CB_SETCURSEL, 0, 0);
	}
	if (SendMessage(GetDlgItem(IDC_COMBO_DSN_BWL), CB_SELECTSTRING, 0, (LPARAM)m_szBwlDSN.c_str()) == CB_ERR){
		// select the first item
		SendMessage(GetDlgItem(IDC_COMBO_DSN_BWL), CB_SETCURSEL, 0, 0);
	}
	SetDlgItemText(IDC_EDIT_USER, m_szUserName.c_str());
	SetDlgItemText(IDC_EDIT_PASSWORD, m_szPassword.c_str());	
	SetDlgItemText(IDC_COMBO_LOCATION, m_szLocation.c_str());
	SetDlgItemText(IDC_COMBO_LOCATION_TRY, m_szLocationFirstTry.c_str());
	SetDlgItemText(IDC_COMBO_DSN_BWL, m_szBwlDSN.c_str());
	SetDlgItemText(IDC_EDIT_USER_BWL, m_szBwlUserName.c_str());
	SetDlgItemText(IDC_EDIT_PASSWORD_BWL, m_szBwlPassword.c_str());
	SetDlgItemText(IDC_EDIT_INTERPOLATION_BWL, m_szBwlInterpolation.c_str());
	SetDlgItemText(IDC_EDIT_SPOT_SCHEDULE_SOURCE_BWL, m_szBwlSpotScheduleSource.c_str());

	// Set the tab pages
	s_nCurSel = ::SendMessage(GetDlgItem(IDC_TAB), TCM_GETCURSEL, 0L, 0L);
	
	// Sirius Database
	ShowControl(IDC_RADIO_SQL, s_nCurSel == SiriusDatabase || s_nCurSel == BwlDatabase);
	ShowControl(IDC_RADIO_LOCAL, s_nCurSel == SiriusDatabase || s_nCurSel == BwlDatabase);
	::EnableWindow(GetDlgItem(IDC_RADIO_SQL), s_nCurSel == SiriusDatabase);
	::EnableWindow(GetDlgItem(IDC_RADIO_LOCAL), s_nCurSel == SiriusDatabase);
	ShowControl(IDC_EDIT_ROOT, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_TXT_DATASOURCE, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_TXT_LOCATION, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_TXT_LOCATION_TRY, s_nCurSel == SiriusDatabase);
	::EnableWindow(GetDlgItem(IDC_TXT_LOCATION_TRY), m_bLocationFirstTryEnabled == true);
	ShowControl(IDC_COMBO_LOCATION_TRY, s_nCurSel == SiriusDatabase);
	::EnableWindow(GetDlgItem(IDC_COMBO_LOCATION_TRY), m_bLocationFirstTryEnabled == true);
	ShowControl(IDC_COMBO_DATASOURCE, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_FIRST_TRY, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_COMBO_LOCATION, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_BROWSE_ROOT, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_TXT_SAVE, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_COMBO_DSN, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_TXT_USER, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_EDIT_USER, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_TXT_PASSWORD, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_EDIT_PASSWORD, s_nCurSel == SiriusDatabase);
	ShowControl(IDC_ICON_32, s_nCurSel == SiriusDatabase || s_nCurSel == BwlDatabase);
	if (s_nCurSel == SiriusDatabase){
		// hide / show controls according to the file system type (NTFS or SQL Server)
		if (IsDlgButtonChecked(IDC_RADIO_SQL)){
			// display SQL DSN
			SetDlgItemText(IDC_TXT_SAVE, _T("&Data source name:"));
			::ShowWindow(GetDlgItem(IDC_COMBO_DSN), SW_SHOW);
			::ShowWindow(GetDlgItem(IDC_EDIT_ROOT), SW_HIDE);
			::ShowWindow(GetDlgItem(IDC_BROWSE_ROOT), SW_HIDE);
			bool bNT = _Module.GetNTAuthentication(m_szDSN);
			::EnableWindow(GetDlgItem(IDC_EDIT_USER), !bNT);
			::EnableWindow(GetDlgItem(IDC_EDIT_PASSWORD), !bNT);
			// disable changes to 'Use Bwl if local file system'
			::EnableWindow(GetDlgItem(IDC_USE_BWL), false);
		} else {
			// display local file system
			SetDlgItemText(IDC_TXT_SAVE, "&File system root:");
			::ShowWindow(GetDlgItem(IDC_COMBO_DSN), SW_HIDE);
			::ShowWindow(GetDlgItem(IDC_EDIT_ROOT), SW_SHOW);
			::ShowWindow(GetDlgItem(IDC_BROWSE_ROOT), SW_SHOW);
			::EnableWindow(GetDlgItem(IDC_EDIT_USER), false);
			::EnableWindow(GetDlgItem(IDC_EDIT_PASSWORD), false);
			// enable changes to 'Use Bwl if local file system'
			::EnableWindow(GetDlgItem(IDC_USE_BWL), true);
		}		
	}

	// Bwl Database
	ShowControl(IDC_TXT_DSN_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_COMBO_DSN_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_TXT_USER_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_EDIT_USER_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_TXT_PASSWORD_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_EDIT_PASSWORD_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_TXT_INTERPOLATION_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_EDIT_INTERPOLATION_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_TXT_SPOT_SCHEDULE_SOURCE_BWL, s_nCurSel == BwlDatabase);
	ShowControl(IDC_EDIT_SPOT_SCHEDULE_SOURCE_BWL, s_nCurSel == BwlDatabase);
	if (s_nCurSel == BwlDatabase){
		bool bNT = _Module.GetNTAuthentication(m_szBwlDSN);
		::EnableWindow(GetDlgItem(IDC_EDIT_USER_BWL), !bNT);
		::EnableWindow(GetDlgItem(IDC_EDIT_PASSWORD_BWL), !bNT);
	}

	// Sirius
	ShowControl(IDC_ICON_SETUP, s_nCurSel == Sirius);
	ShowControl(IDC_ICON_SETUP, s_nCurSel == Sirius);
	ShowControl(IDC_COMBO_DATABASE_MODE, s_nCurSel == Sirius);
	ShowControl(IDC_TXT_DATABASE, s_nCurSel == Sirius);
	ShowControl(IDC_EDIT_PRODUCT_TYPES, s_nCurSel == Sirius);
	ShowControl(IDC_TXT_PRODUCT_TYPES, s_nCurSel == Sirius);
	ShowControl(IDC_BROWSE, s_nCurSel == Sirius);
	ShowControl(IDC_DISPAY_SIRIUS_STATUS, s_nCurSel == Sirius);
	ShowControl(IDC_DISPLAY_FUNCTIONS, s_nCurSel == Sirius);
	ShowControl(IDC_USE_BWL, s_nCurSel == Sirius);
	ShowControl(IDC_YESTERDAY_CLOSE, s_nCurSel == Sirius);
	ShowControl(IDC_CHECK_PUBLISHING, s_nCurSel == Sirius);
}

void CDlgSetup::MembersToModule(void)
{	
	// Process the locations first since they require extra validation.	
	if (_Module.GetFileSystem() != m_fs ||
		_Module.GetSiriusSQLServerDSN() != m_szDSN ||
		_Module.GetLocation() != m_szLocation ||
		_Module.GetLocationFirstTry() != m_szLocationFirstTry){		
		CSiriusComModule::connection_parameters cp;
		FileSystemEnum fs;
		GetConnectionParameters(&cp, &fs);
		_Module.PutLocation(&m_szLocation, &cp, &fs);
		_Module.PutLocationFirstTry(&m_szLocationFirstTry, &cp, &fs);
		// Write m_szLocation and m_szLocationFirstTry back to the dialog since location validation can change these values
		SetDlgItemText(IDC_COMBO_LOCATION, m_szLocation.c_str());
		SetDlgItemText(IDC_COMBO_LOCATION_TRY, m_szLocationFirstTry.c_str());
	}

	if (_Module.GetFileSystem() != m_fs) _Module.PutFileSystem(m_fs);
	if (_Module.GetSiriusSQLServerDSN() != m_szDSN) _Module.PutSiriusSQLServerDSN(m_szDSN);	
	if (_Module.GetDatabaseMode() != m_dbm) _Module.PutDatabaseMode(m_dbm);
	if (_Module.GetDefaultDataSource() != m_dsDefault) _Module.PutDefaultDataSource(m_dsDefault);
	if (_Module.GetDisplaySiriusStatus() != m_bDisplaySiriusStatus) _Module.PutDisplaySiriusStatus(m_bDisplaySiriusStatus);	
	if (_Module.GetLocationFirstTryEnabled() != m_bLocationFirstTryEnabled) _Module.PutLocationFirstTryEnabled(m_bLocationFirstTryEnabled);	
	if (_Module.GetFileSystemRoot() != m_szFileRoot) _Module.PutFileSystemRoot(m_szFileRoot);
	if (_Module.GetProductTypesFileName() != m_szProductTypesFile) _Module.PutProductTypesFileName(m_szProductTypesFile);
	if (_Module.GetSiriusUserName() != m_szUserName) _Module.PutSiriusUserName(m_szUserName);
	if (_Module.GetSiriusPassword() != m_szPassword) _Module.PutSiriusPassword(m_szPassword);
	if (_Module.GetEnablePublishing() != m_bEnablePublishing) _Module.PutEnablePublishing(m_bEnablePublishing);
	if (_Module.GetUseYesterdaysClose() != m_bUseYesterdaysClose) _Module.PutUseYesterdaysClose(m_bUseYesterdaysClose);
	if (_Module.GetUseBwlForSpotSchedule() != m_bUseBwlForSpotSchedule) _Module.PutUseBwlForSpotSchedule(m_bUseBwlForSpotSchedule);
	if (_Module.GetBwlSQLServerDSN() != m_szBwlDSN) _Module.PutBwlSQLServerDSN(m_szBwlDSN);
	if (_Module.GetBwlUserName() != m_szBwlUserName) _Module.PutBwlUserName(m_szBwlUserName);
	if (_Module.GetBwlPassword() != m_szBwlPassword) _Module.PutBwlPassword(m_szBwlPassword);
	if (_Module.GetBwlInterpolationRule() != m_szBwlInterpolation) _Module.PutBwlInterpolationRule(m_szBwlInterpolation);
	if (_Module.GetBwlSpotScheduleSource() != m_szBwlSpotScheduleSource) _Module.PutBwlSpotScheduleSource(m_szBwlSpotScheduleSource);
	if (g_pApplication->GetFunctionsDisabled() != m_bFunctionsDisabled) g_pApplication->PutFunctionsDisabled(m_bFunctionsDisabled);	
}

void CDlgSetup::ModuleToMembers(void)
{
	m_fs = _Module.GetFileSystem();
	m_szProductTypesFile = _Module.GetProductTypesFileName();
	m_szFileRoot = _Module.GetFileSystemRoot();
	m_szDSN = _Module.GetSiriusSQLServerDSN();
	m_szUserName = _Module.GetSiriusUserName();
	m_szPassword = _Module.GetSiriusPassword();
	m_bDisplaySiriusStatus = _Module.GetDisplaySiriusStatus();
	m_bEnablePublishing = _Module.GetEnablePublishing();
	m_bUseYesterdaysClose = _Module.GetUseYesterdaysClose();
	m_bFunctionsDisabled = g_pApplication->GetFunctionsDisabled();
	m_szLocation = _Module.GetLocation();	
	m_szLocationFirstTry = _Module.GetLocationFirstTry();
	m_bLocationFirstTryEnabled = _Module.GetLocationFirstTryEnabled();
	m_bUseBwlForSpotSchedule = _Module.GetUseBwlForSpotSchedule();
	m_dbm = _Module.GetDatabaseMode();
	m_dsDefault = _Module.GetDefaultDataSource();
	m_szBwlDSN = _Module.GetBwlSQLServerDSN();
	m_szBwlUserName = _Module.GetBwlUserName();
	m_szBwlPassword = _Module.GetBwlPassword();
	m_szBwlInterpolation = _Module.GetBwlInterpolationRule();
	m_szBwlSpotScheduleSource = _Module.GetBwlSpotScheduleSource();
}

LRESULT CDlgSetup::OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{		
	EndDialog(wID);
	return 0;
}

//	Reestablish all database connections:
LRESULT CDlgSetup::OnClickedBounceDatabases(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{	
	_Module.ResetAllConnections();
	return 0;
}

//	browse for product type files
LRESULT CDlgSetup::OnClickedBrowse(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	OPENFILENAME						ofn;	
	char								szFilter[]="XML Files (*.xml)\0*.xml\0All Files (*.*)\0*.*\0";
	char								szFile[MAX_PATH];

	DialogToMembers();
	memset(&ofn, NULL, sizeof(OPENFILENAME));
	memset(szFile, NULL, sizeof(szFile));
	strcpy(szFile, m_szProductTypesFile.c_str());
	ofn.lStructSize = sizeof(OPENFILENAME);
	ofn.hwndOwner = NULL;	
	ofn.lpstrFilter = szFilter;
	ofn.nFilterIndex = 1;
	ofn.lpstrFile = szFile;
	ofn.nMaxFile = sizeof(szFile);
	ofn.Flags = OFN_HIDEREADONLY | OFN_FILEMUSTEXIST;
	ofn.lpstrTitle = "Browse Product Types Files";
	if (GetOpenFileName(&ofn) == IDOK){		
		GetUniversalPath(ofn.lpstrFile, &m_szProductTypesFile);			
		MembersToDialog();
	}		
	return 0;
}

LRESULT CDlgSetup::OnClickedBrowseRoot(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	BROWSEINFO							bi;								// browse information structrue
	TCHAR								szDisplayName[MAX_PATH];
	LPITEMIDLIST						pidl;
	TCHAR								szPath[MAX_PATH];
	CComPtr<IMalloc>					spMalloc;
									
	// setup and show the browse for folder
	DialogToMembers();
	SHGetMalloc(&spMalloc);
	ATLASSERT(spMalloc);	
	memset(&bi, 0, sizeof(BROWSEINFO));
	bi.hwndOwner = ::GetActiveWindow();
	bi.pszDisplayName = szDisplayName;
	bi.lpszTitle = TEXT("&File system root:");

	// show browse for folder
	if (!(pidl = SHBrowseForFolder(&bi))) return 0;		// user cancelled
	
	// extract the folder path chosen and get the universal path
	SHGetPathFromIDList(pidl, szPath);
	spMalloc->Free(pidl);	
	GetUniversalPath(szPath, &m_szFileRoot);
	
	// update dialog and return	
	MembersToDialog();
	return 0;
}

LRESULT CDlgSetup::OnClickedRadioLocal(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{			
	DialogToMembers();
	MembersToDialog();	
	return 0;
}

LRESULT CDlgSetup::OnClickedLocationFirstTry(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	DialogToMembers();
	MembersToDialog();
	return 0;
}
	
LRESULT CDlgSetup::OnClickedRadioSQL(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	DialogToMembers();
	MembersToDialog();	
	return 0;
}

LRESULT CDlgSetup::OnClickedRevert(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	_Module.Revert();
	ModuleToMembers();
	MembersToDialog();
	return 0;
}

LRESULT CDlgSetup::OnClickedSave(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{		
	try {
		DialogToMembers();
		MembersToModule();
		_Module.MembersToRegistry();
	} catch (_com_error& e){
		CParameterMap::DisplayError(estring(e), MB_ICONSTOP);	
	} catch (const std::string& sz){	
		CParameterMap::DisplayError(sz, MB_ICONSTOP);
	} catch (...){
		ATLASSERT(false);
	}	
	return 0;
}

LRESULT CDlgSetup::OnDropDownLocation(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	PopulateLocationCombos();
	return 0;
}

LRESULT CDlgSetup::OnDropDownLocationTry(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	PopulateLocationCombos();
	return 0;
}

LRESULT CDlgSetup::OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	TCITEM								tci;
		
	// set up the tab control
	::memset(&tci, NULL, sizeof(tci));
	tci.mask = TCIF_TEXT;
	tci.pszText = "Sirius Database";
	CWindow(GetDlgItem(IDC_TAB)).SendMessage(TCM_INSERTITEM, SiriusDatabase, (LPARAM)&tci);
	tci.pszText = "Bwl Database";
	CWindow(GetDlgItem(IDC_TAB)).SendMessage(TCM_INSERTITEM, BwlDatabase, (LPARAM)&tci);
	tci.pszText = "Sirius Setup";
	CWindow(GetDlgItem(IDC_TAB)).SendMessage(TCM_INSERTITEM, Sirius, (LPARAM)&tci);
	::SendMessage(GetDlgItem(IDC_TAB), TCM_SETCURSEL, s_nCurSel, 0L);

	// load member variables		
	ModuleToMembers();

	// limit edit / combo  control text
	CWindow(GetDlgItem(IDC_EDIT_PRODUCT_TYPES)).SendMessage(EM_LIMITTEXT, MAX_PATH - 1, 0);
	CWindow(GetDlgItem(IDC_EDIT_ROOT)).SendMessage(EM_LIMITTEXT, MAX_PATH - 1, 0);
	CWindow(GetDlgItem(IDC_EDIT_USER)).SendMessage(EM_LIMITTEXT, 32, 0);
	CWindow(GetDlgItem(IDC_EDIT_USER_BWL)).SendMessage(EM_LIMITTEXT, 32, 0);
	CWindow(GetDlgItem(IDC_EDIT_PASSWORD)).SendMessage(EM_LIMITTEXT, 32, 0);
	CWindow(GetDlgItem(IDC_EDIT_PASSWORD_BWL)).SendMessage(EM_LIMITTEXT, 32, 0);
	CWindow(GetDlgItem(IDC_COMBO_LOCATION)).SendMessage(CB_LIMITTEXT, 16, 0);
	CWindow(GetDlgItem(IDC_EDIT_INTERPOLATION_BWL)).SendMessage(EM_LIMITTEXT, 32, 0);
	CWindow(GetDlgItem(IDC_EDIT_SPOT_SCHEDULE_SOURCE_BWL)).SendMessage(EM_LIMITTEXT, 32, 0);

	// populate the database mode combo
	std::vector<std::string>			aszDatabaseModes;
	CEnumMap::GetEnumList("DatabaseModeEnum", LIBID_Sirius, &aszDatabaseModes);
	for (std::vector<std::string>::const_iterator it = aszDatabaseModes.begin(); it != aszDatabaseModes.end(); it++){
		SendMessage(GetDlgItem(IDC_COMBO_DATABASE_MODE), CB_ADDSTRING, 0, (LPARAM)(it->c_str()));
	}

	// populate the DSN combos
	PopulateDSNCombo(IDC_COMBO_DSN, "sirius_");
	PopulateDSNCombo(IDC_COMBO_DSN_BWL, "bwl_");

	// populate the default data source combo
	std::vector<std::string>			aszDataSource;
	CEnumMap::GetEnumList("DataSourceEnum", LIBID_Sirius, &aszDataSource);
	for (std::vector<std::string>::const_iterator it = aszDataSource.begin(); it != aszDataSource.end(); it++){
		SendMessage(GetDlgItem(IDC_COMBO_DATASOURCE), CB_ADDSTRING, 0, (LPARAM)(it->c_str()));
	}

	// initialise the dialog controls
	MembersToDialog();
	
	// return 1 so the operating system sets the appropriate control focus
	return 1;
}

LRESULT CDlgSetup::OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{					
	DialogToMembers();
	try {
		MembersToModule();
		EndDialog(wID);
	} catch (_com_error& e){
		CParameterMap::DisplayError(estring(e), MB_ICONSTOP);
	} catch (const std::string& sz){	
		CParameterMap::DisplayError(sz, MB_ICONSTOP);	
	} catch (...){
		ATLASSERT(false);
	}
	return 0;
}

LRESULT CDlgSetup::OnSelchangeComboDSN(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{	
	DialogToMembers();	
	MembersToDialog();	
	return 0;
}

LRESULT CDlgSetup::OnSelchangeComboDSNBwl(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	DialogToMembers();
	MembersToDialog();
	return 0;
}

LRESULT CDlgSetup::OnSelChangeTab(int idCtrl, LPNMHDR pnmh, BOOL& bHandled)
{
	DialogToMembers();
	MembersToDialog();
	return 0;
}

void CDlgSetup::PopulateDSNCombo(int nControlID, const std::string& szPrefix) const
//	nControlID - e.g. IDC_COMBO_DSN
//	szPrefix - e.g. "sirius_", "bwl_"
{
	SQLHANDLE				hEnv = NULL;							// SQL environment handle
	
	do {		
		SQLRETURN			rc ;									// function return value
		SQLCHAR				szSource[SQL_MAX_DSN_LENGTH + 1];		// data source name			
		SQLCHAR				szDesc[255];							// data source description
		SQLSMALLINT			nBuff, nDesc;							// string buffer lengths
		SQLUSMALLINT		nDir = SQL_FETCH_FIRST_SYSTEM;			// source / description page to return

		// allocate BWL_DRIVER_NAME, 0, sizeof(SQLHANDLE));
		if ((rc = SQLAllocHandle(SQL_HANDLE_ENV, SQL_NULL_HANDLE, &hEnv)) != SQL_SUCCESS) break;
		if ((rc = SQLSetEnvAttr(hEnv, SQL_ATTR_ODBC_VERSION, (void*)SQL_OV_ODBC3, SQL_IS_INTEGER)) != SQL_SUCCESS) break;

		// get the sources
		while ((rc = SQLDataSources(hEnv, nDir, szSource, SQL_MAX_DSN_LENGTH + 1, &nBuff, szDesc, 255, &nDesc) != SQL_NO_DATA_FOUND)){
			if (!strcmp((char*)szDesc, _Module.s_szODBCDriverName) && !strncmp((char*)szSource, szPrefix.c_str(), szPrefix.size())){
				SendMessage(GetDlgItem(nControlID), CB_ADDSTRING, 0, (LPARAM)szSource);
			}			
			nDir = SQL_FETCH_NEXT;
		}
		
	} while (false);
	if (hEnv) SQLFreeHandle(SQL_HANDLE_ENV, hEnv);
}

// populate the location and try location combos
void CDlgSetup::PopulateLocationCombos(void)
{	
	CSiriusComModule::connection_parameters		cp;
	std::vector<estring>						aszLocations;
	FileSystemEnum								fs;
	estring										szLocation;
	estring										szLocationFirstTry;
	
	GetDlgItemText(IDC_COMBO_LOCATION, &szLocation);
	GetDlgItemText(IDC_COMBO_LOCATION_TRY, &szLocationFirstTry);		
	GetConnectionParameters(&cp, &fs);
	
	SendMessage(GetDlgItem(IDC_COMBO_LOCATION), CB_RESETCONTENT, 0, 0);
	SendMessage(GetDlgItem(IDC_COMBO_LOCATION_TRY), CB_RESETCONTENT, 0, 0);

	g_pApplication->GetLocations(fs, &cp, &aszLocations, NULL, false, false);
	for (std::vector<estring>::const_iterator it = aszLocations.begin(); it != aszLocations.end(); it++){
		SendMessage(GetDlgItem(IDC_COMBO_LOCATION), CB_ADDSTRING, 0, (LPARAM)(it->c_str()));
		SendMessage(GetDlgItem(IDC_COMBO_LOCATION_TRY), CB_ADDSTRING, 0, (LPARAM)(it->c_str()));
	}
	
	SetDlgItemText(IDC_COMBO_LOCATION, szLocation.c_str());
	SetDlgItemText(IDC_COMBO_LOCATION_TRY, szLocationFirstTry.c_str());
}

void CDlgSetup::ShowControl(int nID, bool bShow)
{
	::EnableWindow(GetDlgItem(nID), bShow ? TRUE : FALSE);
	::ShowWindow(GetDlgItem(nID), bShow ? SW_SHOWNA : SW_HIDE);	
}
