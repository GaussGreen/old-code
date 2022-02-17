//	dlgmoduleload.cpp : Implementation of CDlgModuleLoad
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "dlgmoduleload.h"
#include "excelinterface.h"
#include "lmcons.h"
#include "siriusapplication.h"

/*static*/ WINDOWPLACEMENT				CDlgModuleLoad::s_wp = {sizeof(WINDOWPLACEMENT), 0, SW_SHOWNORMAL, {-1, -1}, {-1, -1}, {0, 0, 0, 0}};


//	Close button is pressed
LRESULT CDlgModuleLoad::OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	EndDialog(wID);
	return 0;
}

//	Delete the currently selected module from the database
LRESULT CDlgModuleLoad::OnClickedDelete(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	std::stringstream					ssQuery;						// SQL query to execute		
	_ConnectionPtr						spConnection;
	
	// m_nModuleSelected holds the module key to delete
	if (!m_nModuleSelected){
		::MessageBeep(MB_ICONEXCLAMATION);		
		CParameterMap::DisplayError(CStringResource(IDS_NO_MODULE_SELECTED).str() + ".", MB_ICONEXCLAMATION);
		return 0;
	}

	// prompt
	if (::MessageBox(::GetActiveWindow(), CStringResource(IDS_DELETE_CONFIRM).str().c_str(), _T("Delete Module"), MB_YESNOCANCEL | MB_ICONQUESTION) != IDYES) return 0;
		
	ssQuery << "sp_user_delete_vbmodule " << m_nModuleSelected;
	try {		
		_Module.GetConnection(spConnection);
		spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);				
	} catch (_com_error& e){
		CParameterMap::DisplayError(estring(e), MB_ICONEXCLAMATION);
		return 0;
	}

	// successful if this point is reached - delete the corresponding tree control item		
	SendMessage(GetDlgItem(IDC_TREE_MODULE), TVM_DELETEITEM, 0, (LPARAM)m_hModuleNodeSelected);	
	return 0;
}

//	A user node is being expanded - pull in the module details associated with that user
LRESULT CDlgModuleLoad::OnItemExpandingTreeModule(int idCtrl, LPNMHDR pnmh, BOOL& bHandled)
{
	LPNMTREEVIEW						pnmtv = (LPNMTREEVIEW)pnmh;
	TCHAR								szUser[UNLEN + 1];
	_ConnectionPtr						spConnection;
	_RecordsetPtr						prs;
	CComVariant							vKey;
	CComVariant							vModuleName;
	CComVariant							vModuleType;
	std::string							szModuleName;
	TVINSERTSTRUCT						tvis;
	std::stringstream					ssQuery;						// SQL query to execute	
	
	if (TreeView_GetChild(GetDlgItem(IDC_TREE_MODULE), pnmtv->itemNew.hItem)){
		// don't refresh if the item already has children
		return 0;
	}

	pnmtv->itemNew.mask |= TVIF_TEXT;
	pnmtv->itemNew.pszText = szUser;
	pnmtv->itemNew.cchTextMax = UNLEN + 1;
	TreeView_GetItem(GetDlgItem(IDC_TREE_MODULE), &pnmtv->itemNew);		// szUser now contains the tree text	

	memset(&tvis, NULL, sizeof(tvis));
	tvis.hParent = pnmtv->itemNew.hItem;					
	tvis.hInsertAfter = TVI_SORT;
	tvis.item.mask = TVIF_TEXT | TVIF_IMAGE | TVIF_SELECTEDIMAGE | TVIF_PARAM | TVIF_CHILDREN;
	//tvis.item.iImage = tvis.item.iSelectedImage = 0;

	// read in the module names associated with the user
	ssQuery << "sp_user_get_vbmodule_names '" << szUser << "'";
	try {		
		_Module.GetConnection(spConnection);
		prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		while (!prs->adoEOF){
			vKey = prs->GetFields()->GetItem(0L)->GetValue();
			vModuleType = prs->GetFields()->GetItem(2L)->GetValue();
			ATLASSERT(vKey.vt == VT_I4);			
			ATLASSERT(vModuleType.vt == VT_I4);
			CParameterMap pmModuleName;
			pmModuleName.SetValue((CComVariant)prs->GetFields()->GetItem(1L)->GetValue());
			pmModuleName.GetString(&szModuleName);
			tvis.item.pszText = (LPTSTR)szModuleName.c_str();
			tvis.item.lParam = vKey.lVal;
			tvis.item.iImage = tvis.item.iSelectedImage = MapImageModuleType(NULL, &vModuleType.lVal);
			SendMessage(GetDlgItem(IDC_TREE_MODULE), TVM_INSERTITEM, 0, (LPARAM)&tvis);
			prs->MoveNext();
		}	
	} catch (_com_error& e){
		CParameterMap::DisplayError(estring(e), MB_ICONSTOP);
		EndDialog(IDCANCEL);
		return 0;
	}

	// Remove the expandable icon if we fail to get any children
	if (!TreeView_GetChild(GetDlgItem(IDC_TREE_MODULE), pnmtv->itemNew.hItem)){
		pnmtv->itemNew.cChildren = 0;
		::SendMessage(GetDlgItem(IDC_TREE_MODULE), TVM_SETITEM, 0, (LPARAM)&pnmtv->itemNew);
	}		
	return 0;
}

//	Dialog Startup
LRESULT CDlgModuleLoad::OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	_ConnectionPtr						spConnection;
	_RecordsetPtr						prs;
	CComVariant							vUser;
	std::string							szUser;
	TVINSERTSTRUCT						tvis;
	
	if (_Module.GetExcel(m_ddExcel)){
		ATLASSERT(false);
		EndDialog(NULL);
		return 0;
	}	
	m_nModuleSelected = 0;
	m_hModuleNodeSelected = NULL;
	memset(&tvis, NULL, sizeof(tvis));
	tvis.hParent = TVI_ROOT;					
	tvis.hInsertAfter = TVI_SORT;
	tvis.item.mask = TVIF_TEXT | TVIF_IMAGE | TVIF_SELECTEDIMAGE | TVIF_PARAM | TVIF_CHILDREN;
	tvis.item.iImage = tvis.item.iSelectedImage = 4;	
	tvis.item.cChildren = 1;
	SendMessage(GetDlgItem(IDC_EDIT_PREVIEW), WM_SETFONT, (WPARAM)_Module.m_hFont, MAKELPARAM(FALSE, 0));
	TreeView_SetImageList(GetDlgItem(IDC_TREE_MODULE), _Module.m_hImageList, TVSIL_NORMAL);

	// Set the module tree from the database (we only read in the top level items).
	try {
		_Module.GetConnection(spConnection);
		prs = spConnection->Execute(L"sp_user_get_vbmodule_users", NULL, -1);
		while (!prs->adoEOF){
			CParameterMap pmUser;
			pmUser.SetValue((CComVariant)prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue());			
			pmUser.GetString(&szUser);
			tvis.item.pszText = (LPTSTR)szUser.c_str();			
			SendMessage(GetDlgItem(IDC_TREE_MODULE), TVM_INSERTITEM, 0, (LPARAM)&tvis);
			prs->MoveNext();
		}	
	} catch (_com_error& e){
		CParameterMap::DisplayError(estring(e), MB_ICONSTOP);		
		EndDialog(IDCANCEL);
		return 0;				
	}
	if (!SendMessage(GetDlgItem(IDC_TREE_MODULE), TVM_GETCOUNT, 0, 0)){		
		CParameterMap::DisplayError(CStringResource(IDS_NO_MODULES_FOUND).str() + ".", MB_ICONINFORMATION);		
		EndDialog(IDCANCEL);
	}

	// Set the download combo box.
	CComDispatchDriverEx				ddWorkbooks, ddWorkbook;
	CComVariant							WorkbooksCount;
	CComVariant							Name;
	CParameterMap						pmName;
	std::string							szName, szActiveName;
	COMBOBOXEXITEM						cbei;

	memset(&cbei, NULL, sizeof(cbei));
	cbei.mask = CBEIF_IMAGE | CBEIF_SELECTEDIMAGE | CBEIF_TEXT;
	cbei.iImage = cbei.iSelectedImage = 5;
	SendMessage(GetDlgItem(IDC_COMBOBOXEX), CBEM_SETIMAGELIST, 0, (LPARAM)_Module.m_hImageList);	
	m_ddExcel.GetPropertyByName(L"ActiveWorkbook", ddWorkbook);
	if (ddWorkbook) ddWorkbook.GetPropertyByName(L"Name", &Name);
	pmName.SetValue(Name);
	pmName.GetString(&szActiveName);
	if (!m_ddExcel.GetPropertyByName(L"Workbooks", ddWorkbooks)){
		if (!ddWorkbooks.GetPropertyByName(L"Count", &WorkbooksCount)){
			ATLASSERT(WorkbooksCount.vt == VT_I4);
			for (long nWorkbook = 1; nWorkbook <= WorkbooksCount.lVal; nWorkbook++){
				if (ddWorkbooks.GetPropertyByName(L"Item", &CComVariant(nWorkbook), 1, ddWorkbook)) continue;
				if (ddWorkbook.GetPropertyByName(L"Name", &Name)) continue;
				pmName.SetValue(Name);
				pmName.GetString(&szName);
				cbei.iItem = nWorkbook - 1;
				cbei.pszText = (LPTSTR)szName.c_str();
				SendMessage(GetDlgItem(IDC_COMBOBOXEX), CBEM_INSERTITEM, 0, (LPARAM)&cbei);
				if (szName == szActiveName){
					// select the active workbook
					::SendMessage(GetDlgItem(IDC_COMBOBOXEX), CB_SETCURSEL, nWorkbook - 1, 0);
				}
			}		
		}
	}

	// Error if there are no workbooks.
	if (!::SendMessage((HWND)::SendMessage(GetDlgItem(IDC_COMBOBOXEX), CBEM_GETCOMBOCONTROL, 0, 0), CB_GETCOUNT, 0, 0)){
		CParameterMap::DisplayError(CStringResource(IDS_NO_WORKBOOKS_FOUND).str() + ".", MB_ICONINFORMATION);
		EndDialog(IDCANCEL);
	}
		
	// Get initial window placements.
	m_wp_init.length = m_wp_tree_init.length = m_wp_preview_init.length = m_wp_edit_init.length = m_wp_combo_init.length = m_wp_downloadto_init.length = m_wp_delete_init.length = m_wp_load_init.length = m_wp_close_init.length = sizeof(WINDOWPLACEMENT);
	GetWindowPlacement(&m_wp_init);
	::GetWindowPlacement(GetDlgItem(IDC_TREE_MODULE), &m_wp_tree_init);
	::GetWindowPlacement(GetDlgItem(IDC_PREVIEW), &m_wp_preview_init);
	::GetWindowPlacement(GetDlgItem(IDC_EDIT_PREVIEW), &m_wp_edit_init);
	::GetWindowPlacement(GetDlgItem(IDC_COMBOBOXEX), &m_wp_combo_init);
	::GetWindowPlacement(GetDlgItem(IDC_DOWNLOAD_TO), &m_wp_downloadto_init);
	::GetWindowPlacement(GetDlgItem(IDC_DELETE), &m_wp_delete_init);
	::GetWindowPlacement(GetDlgItem(1), &m_wp_load_init);
	::GetWindowPlacement(GetDlgItem(2), &m_wp_close_init);

	// Resize to memory
	if (s_wp.rcNormalPosition.right){
		SetWindowPlacement(&s_wp);
	}
	return 1;
}

//	Import the module denoted by m_nModuleSelected
LRESULT CDlgModuleLoad::OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{			
	_ConnectionPtr						spConnection;
	_RecordsetPtr						prs;
	std::stringstream					ssQuery;						// SQL query to execute	
	CComVariant							vName;
	CComVariant							vText;
	CComVariant							vModuleType;
	CComDispatchDriverEx				ddWorkbook;						// workbook to set	
		
	if (!m_nModuleSelected){
		::MessageBeep(MB_ICONEXCLAMATION);
		CParameterMap::DisplayError(CStringResource(IDS_NO_MODULE_SELECTED).str() + ".", MB_ICONEXCLAMATION);
		return 0;
	}

	// get the selected workbook
	long nIndex;
	if ((nIndex = SendMessage(GetDlgItem(IDC_COMBOBOXEX), CB_GETCURSEL, 0, 0)) != CB_ERR){
		long nTextLen = SendMessage(GetDlgItem(IDC_COMBOBOXEX), CB_GETLBTEXTLEN, nIndex, 0) + 1;		
		LPSTR lpszBuffer = new TCHAR[nTextLen];
		SendMessage(GetDlgItem(IDC_COMBOBOXEX), CB_GETLBTEXT, nIndex, (LPARAM)lpszBuffer);
		m_ddExcel.GetPropertyByName(L"Workbooks", &CComVariant(lpszBuffer), 1, ddWorkbook);
		delete lpszBuffer;
	}
	if (!ddWorkbook){
		CParameterMap::DisplayError(CStringResource(IDS_NO_WORKBOOK_SELECTED).str() + ".", MB_ICONSTOP);
		EndDialog(IDCANCEL);
		return 0;
	}

	// get the module name and the text from the database
	ssQuery << "sp_user_get_vbmodule_code " << m_nModuleSelected;
	try {
		_Module.GetConnection(spConnection);
		prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		vName = prs->GetFields()->GetItem(1L)->GetValue();
		vText = prs->GetFields()->GetItem(2L)->GetValue();
		vModuleType = prs->GetFields()->GetItem(3L)->GetValue();
	} catch (_com_error& e){
		CParameterMap::DisplayError(estring(e), MB_ICONSTOP);
		EndDialog(IDCANCEL);
	}
			
	// add the module to the active workbook
	CComDispatchDriverEx ddVBProject, ddVBComponents, ddVBComponent;	
	if (!ddWorkbook.GetPropertyByName(L"VBProject", ddVBProject)){
		if (!ddVBProject.GetPropertyByName(L"VBComponents", ddVBComponents)){
			if (!ddVBComponents.Invoke1(L"Add", &vModuleType, ddVBComponent)){
				// Set the name property. The base string is vName. Try setting this as the 
				// name property of ddVBComponent. If we fail then we suffix a number. 
				// We give up after a few tries								
				for (long nSuffix = 0; nSuffix < 65536; nSuffix++){
					estring szName(vName);
					if (nSuffix) szName += estring(nSuffix);
					if (!ddVBComponent.PutPropertyByName(L"Name", &CComVariant(szName))){
						CComDispatchDriverEx ddCodeModule;
						if (!ddVBComponent.GetPropertyByName(L"CodeModule", ddCodeModule)){
							if (!ddCodeModule.Invoke1(L"AddFromString", &vText)){
								if (vModuleType.lVal == CExcelInterface::vbext_ct_ClassModule){
									// Set the PublicNotCreatable (=2) value if we can.
									// We're not worried about errors here.
									CComDispatchDriverEx ddProperties;
									if (!ddVBComponent.GetPropertyByName(L"Properties", ddProperties)){										
										CComDispatchDriverEx ddProperty;																				
										// "Item" is, unconventionally, a method for the interface on ddProperties.
										if (!ddProperties.Invoke1(L"Item", &CComVariant(L"Instancing"), ddProperty)){										
											ddProperty.PutPropertyByName(L"Value", &CComVariant(2L));
										}
									}
								}
								
								// Grey the load button (which gives the user indication that the load has worked) and move the focus.
								SetFocus();		// we need this as the module import occasionally drops the focus							
								GotoDlgCtrl(GetDlgItem(IDCANCEL));
								::EnableWindow(GetDlgItem(IDOK), FALSE);
								return 0;
							}
						}
						break;
					}
				}																
			}
		}
	}

	// fail if this point is reached   	
	CParameterMap::DisplayError(MB_ICONSTOP);		
	EndDialog(IDCANCEL);
	return 0;
}

// The module tree control selection has changed - preview any new module
LRESULT CDlgModuleLoad::OnSelChangedTreeModule(int idCtrl, LPNMHDR pnmh, BOOL& bHandled)
{
	::EnableWindow(GetDlgItem(IDOK), TRUE);	
	LPNMTREEVIEW						pnmtv = (LPNMTREEVIEW)pnmh;	
	int									nKey = pnmtv->itemNew.lParam;
	_ConnectionPtr						spConnection;
	_RecordsetPtr						prs;
	std::stringstream					ssQuery;						// SQL query to execute	
	CComVariant							vName;
	CComVariant							vText;
	
	if (nKey){
		m_nModuleSelected = nKey;	
		m_hModuleNodeSelected = pnmtv->itemNew.hItem;
	}
	if (!m_nModuleSelected) return 0;
	ssQuery << "sp_user_get_vbmodule_code " << m_nModuleSelected;
	try {
		_Module.GetConnection(spConnection);
		prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		vName = prs->GetFields()->GetItem(0L)->GetValue();
		vText = prs->GetFields()->GetItem(2L)->GetValue();
		SetDlgItemText(IDC_EDIT_PREVIEW, (LPCTSTR)_bstr_t(vText.bstrVal));		
	} catch (_com_error& e){
		// only report this problem if nKey is non-zero (or we break the delete module functionality)
		if (nKey){
			CParameterMap::DisplayError(estring(e), MB_ICONSTOP);
			EndDialog(IDCANCEL);
		}
	}
		
	// enable delete button if the current user matches vName
	estring szName = CParameterMap::GetUser();
	
	if (!szName.CompareNoCaseAndSpace(estring(vName))){
		::EnableWindow(GetDlgItem(IDC_DELETE), TRUE);
	} else {
		if (::GetFocus() == GetDlgItem(IDC_DELETE)) NextDlgCtrl();
		::EnableWindow(GetDlgItem(IDC_DELETE), FALSE);
	}
	return 0;
}

LRESULT CDlgModuleLoad::OnSize(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{			
	WINDOWPLACEMENT wpuse = {sizeof(WINDOWPLACEMENT), 0, SW_SHOWNORMAL, {-1, -1}, {-1, -1}, {0, 0, 0, 0}};	
	WINDOWPLACEMENT wp;

	wp.length = sizeof(WINDOWPLACEMENT);
	GetWindowPlacement(&wp);

	s_wp = wp;

	long dy = (wp.rcNormalPosition.bottom - wp.rcNormalPosition.top) - (m_wp_init.rcNormalPosition.bottom - m_wp_init.rcNormalPosition.top);	
	long dx = (wp.rcNormalPosition.right - wp.rcNormalPosition.left) - (m_wp_init.rcNormalPosition.right - m_wp_init.rcNormalPosition.left);
	
	// Initial amounts are t_i, e_i, c_i, cx_i.
	// t - width of tree
	// e - width of edit
	// cx - width of dialog
	// c - total width of non-control space
	long t_i = m_wp_tree_init.rcNormalPosition.right - m_wp_tree_init.rcNormalPosition.left;
	long e_i = m_wp_edit_init.rcNormalPosition.right - m_wp_edit_init.rcNormalPosition.left;
	long cx_i = m_wp_init.rcNormalPosition.right - m_wp_init.rcNormalPosition.left;
	long c_i = cx_i - t_i - e_i;	
	long c = c_i;
	long cx = wp.rcNormalPosition.right - wp.rcNormalPosition.left;
	long t = std::min(((double)cx - (double)c) / ((double)cx_i - (double)c_i) * (double)t_i, 250.0);
	long e = cx - c - t;

	// Module tree
	wpuse.rcNormalPosition = m_wp_tree_init.rcNormalPosition;
	wpuse.rcNormalPosition.bottom += dy;		
	wpuse.rcNormalPosition.right = wpuse.rcNormalPosition.left + t;	
	CWindow(GetDlgItem(IDC_TREE_MODULE)).SetWindowPlacement(&wpuse);

	// Preview window static
	wpuse.rcNormalPosition = m_wp_preview_init.rcNormalPosition;
	::OffsetRect(&wpuse.rcNormalPosition, t - t_i, 0);
	CWindow(GetDlgItem(IDC_PREVIEW)).SetWindowPlacement(&wpuse);

	// Preview window
	wpuse.rcNormalPosition = m_wp_edit_init.rcNormalPosition;
	wpuse.rcNormalPosition.bottom += dy;	
	wpuse.rcNormalPosition.left += (t - t_i);
	wpuse.rcNormalPosition.right += (cx - cx_i);		
	CWindow(GetDlgItem(IDC_EDIT_PREVIEW)).SetWindowPlacement(&wpuse);

	// Download to combo
	wpuse.rcNormalPosition = m_wp_combo_init.rcNormalPosition;
	::OffsetRect(&wpuse.rcNormalPosition, t - t_i, dy);
	wpuse.rcNormalPosition.right += ((cx - cx_i) - (t - t_i));	
	CWindow(GetDlgItem(IDC_COMBOBOXEX)).SetWindowPlacement(&wpuse);

	// Download to static
	wpuse.rcNormalPosition = m_wp_downloadto_init.rcNormalPosition;
	::OffsetRect(&wpuse.rcNormalPosition, t - t_i, dy);
	CWindow(GetDlgItem(IDC_DOWNLOAD_TO)).SetWindowPlacement(&wpuse);

	// Delete button
	wpuse.rcNormalPosition = m_wp_delete_init.rcNormalPosition;
	::OffsetRect(&wpuse.rcNormalPosition, dx, dy);
	CWindow(GetDlgItem(IDC_DELETE)).SetWindowPlacement(&wpuse);

	// Load button
	wpuse.rcNormalPosition = m_wp_load_init.rcNormalPosition;
	::OffsetRect(&wpuse.rcNormalPosition, dx, dy);
	CWindow(GetDlgItem(1)).SetWindowPlacement(&wpuse);

	// Save button
	wpuse.rcNormalPosition = m_wp_close_init.rcNormalPosition;
	::OffsetRect(&wpuse.rcNormalPosition, dx, dy);
	CWindow(GetDlgItem(2)).SetWindowPlacement(&wpuse);

	RedrawWindow();	
	return 0;
}

LRESULT CDlgModuleLoad::OnSizing(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	LPRECT prect = reinterpret_cast<LPRECT>(lParam);					
	// Minimum size
	if (prect->right - prect->left < 270){
		prect->right = prect->left + 270;
	}
	if (prect->bottom - prect->top < 190){
		prect->bottom = prect->top + 190;
	}	
	
	// hello there
	return 1L;
}










//	Asserts and returns -1 if failed.
long CDlgModuleLoad::MapImageModuleType(int* pnImage, long* pnModuleType)
//	pnImage - if not null then we return a module type corresponding to the image value
//  pnModuleType - if not null then we return an image corresponding to the module type value
{
	if (pnImage && pnModuleType){
		ATLASSERT(false);
		return -1;
	}	
	if (pnImage){
		if (*pnImage == 1){
			return CExcelInterface::vbext_ct_StdModule;
		} else if (*pnImage == 7){
			return CExcelInterface::vbext_ct_ClassModule;
		} else {		
			ATLASSERT(false);
			return -1;
		}
	}

	if (pnModuleType){
		if (*pnModuleType == CExcelInterface::vbext_ct_StdModule){
			return 1;
		} else if (*pnModuleType == CExcelInterface::vbext_ct_ClassModule){
			return 7;
		} else {
			ATLASSERT(false);
			return -1;
		}
	}

	ATLASSERT(false);
	return -1;
}
