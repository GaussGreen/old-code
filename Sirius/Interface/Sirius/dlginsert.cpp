//	dlginsert.cpp : Implementation of CDlgInsert
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "dlginsert.h"
#include "siriusapplication.h"
#include "excelinterface.h"

//	if m_szAt is blank then create any workbooks and worksheets as necessary
HRESULT CDlgInsert::CheckInsertAt(void)
{
	CComVariant vWorkbooks, vCount, vActiveWorkbook, vWorksheets, vWorksheet, vRange;
		
	if (m_szAt.size()) return S_OK;		
	// create a new workbook if necessary	
	if (m_ddExcel.GetPropertyByName(L"Workbooks", &vWorkbooks)) ATLASSERT(false);
	if (CComDispatchDriverEx(vWorkbooks.pdispVal).GetPropertyByName(L"Count", &vCount)) ATLASSERT(false);		
	ATLASSERT(vCount.vt == VT_I4);
	if (!vCount.lVal){
        CComDispatchDriverEx(vWorkbooks.pdispVal).Invoke0(L"Add");
	}
    // always create a new worksheet
	if (m_ddExcel.GetPropertyByName(L"ActiveWorkbook", &vActiveWorkbook)) ATLASSERT(false);
	if (CComDispatchDriverEx(vActiveWorkbook.pdispVal).GetPropertyByName(L"Worksheets", &vWorksheets)) ATLASSERT(false);
	CComDispatchDriverEx(vWorksheets.pdispVal).Invoke0(L"Add", &vWorksheet);	
	if (m_ddExcel.GetPropertyByName(L"Range", &CComVariant(L"A1"), 1, &vRange)) ATLASSERT(false);	
	return CExcelInterface::GetAddress(vRange.pdispVal, &m_szAt);		
}

//	Populate the list control with the product or object list as appropriate
HRESULT CDlgInsert::LoadList(void)
{	
	CComVariant							vList;
	HRESULT								hr;	
	std::string							szItem;							// item to insert into the list control
	
	switch (m_dt){
	case insert_object:
		if (hr = g_pApplication->GetObjectTypes(&vList)) return hr;
		break;
	case insert_product:
		if (hr = g_pApplication->GetProductTypes(&vList)) return hr;
		break;
	default:
		ATLASSERT(false);
	}
	if (hr = m_pmList.SetValue(vList)) return hr;
	if (hr = m_pmList.RemoveColumn(1)) return hr;
	SendMessage(GetDlgItem(IDC_LIST), LB_RESETCONTENT, 0, 0);	
	for (long nRow = m_pmList.GetRows() - 1; nRow >=0; nRow--){
		m_pmList.GetValue(nRow, 0, &szItem);
		SendMessage(GetDlgItem(IDC_LIST), LB_INSERTSTRING, 0, (LPARAM)szItem.c_str());
	}

	// select the first item
	SendMessage(GetDlgItem(IDC_LIST), LB_SELECTSTRING, -1, (LPARAM)szItem.c_str());
	return S_OK;
}

LRESULT CDlgInsert::OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	EndDialog(wID);
	return 0;
}

//	Reload the product types
LRESULT CDlgInsert::OnClickedRefresh(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	g_pApplication->LoadProductTypes();
	LoadList();
	return 0;
}

//	Double clicking the list control has the same effect as clicking the OK button
LRESULT CDlgInsert::OnDblClkList(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	return OnOK(wNotifyCode, wID, hWndCtl, bHandled);	
}

//	Called just before the dialog is displayed
LRESULT CDlgInsert::OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{	
	bool								bDisplayProductOnly = FALSE;		
			
	if (_Module.GetExcel(m_ddExcel)){
		ATLASSERT(false);
		EndDialog(NULL);
		return 0;
	}	
	
	bDisplayProductOnly = _Module.GetDisplayProductOnly();
	switch (lParam){
	case insert_object:
		SetWindowText(_T("Insert Object"));
		::ShowWindow(GetDlgItem(IDC_PRODUCT_ONLY), SW_HIDE);		
		::EnableWindow(GetDlgItem(IDC_REFRESH), false);
		m_dt = insert_object;
		break;
	case insert_product:
		SetWindowText(_T("Insert Product"));
		::ShowWindow(GetDlgItem(IDC_PRODUCT_ONLY), SW_SHOW);			
		SendMessage(GetDlgItem(IDC_PRODUCT_ONLY), BM_SETCHECK, bDisplayProductOnly ? BST_CHECKED : BST_UNCHECKED, 0);
		::EnableWindow(GetDlgItem(IDC_REFRESH), true);
		m_dt = insert_product;
		break;	
	default:
		ATLASSERT(false);
		EndDialog(NULL);
		return 0;
	}
	
	// get m_szAt	
	CComVariant vActiveCell, vWorksheet, vName, vSelection, vAddress;	
	m_szAt.erase();
	if (m_ddExcel.GetPropertyByName(L"ActiveCell", &vActiveCell)) ATLASSERT(false);
	if (vActiveCell.pdispVal){
		if (CComDispatchDriverEx(vActiveCell.pdispVal).GetPropertyByName(L"Worksheet", &vWorksheet)) ATLASSERT(false);
		if (CComDispatchDriverEx(vWorksheet.pdispVal).GetPropertyByName(L"Name", &vName)) ATLASSERT(false);
		if (m_ddExcel.GetPropertyByName(L"Selection", &vSelection)) ATLASSERT(false);
		if (CComDispatchDriverEx(vSelection.pdispVal).GetPropertyByName(L"AddressLocal", &vAddress)) ATLASSERT(false);
		m_szAt = "'" + estring(vName.bstrVal) + "'!" + estring(vAddress.bstrVal);
	}

	// load the list control data
	LoadList();
	return 1;	// i.e. let the system set the focus
}

LRESULT CDlgInsert::OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled)
{
	estring								szItem;							// list control item
	int									nItem;							// index of the selected item
			
	if ((nItem = SendMessage(GetDlgItem(IDC_LIST), LB_GETCURSEL, 0, 0)) == LB_ERR){
		::MessageBeep(MB_ICONEXCLAMATION);
		return 0;
	}
	if (m_pmList.GetValue(nItem, 0, &szItem)){
		::MessageBeep(MB_ICONEXCLAMATION);
		return 0;
	}	
	    
	switch (m_dt){
	case insert_object:
		if (szItem == IDS_HEADER_PRODUCT){
			// transform the form to the insert product case			
			BOOL bDummy;
			OnInitDialog(0, 0, insert_product, bDummy);
			return 0;
		} else {
			CheckInsertAt();
			if (CExcelInterface::InsertObject(szItem, m_szAt)){
				CParameterMap::DisplayError(MB_ICONEXCLAMATION);
				return 0;
			}
		}
		break;
	case insert_product:				
		_Module.PutDisplayProductOnly(IsDlgButtonChecked(IDC_PRODUCT_ONLY) ? true : false);
		CheckInsertAt();
		if (CExcelInterface::InsertProduct(szItem, m_szAt)){
			CParameterMap::DisplayError(MB_ICONEXCLAMATION);
			return 0;
		}
		break;
	default:
		ATLASSERT(false);
	}    
	EndDialog(wID);	
	return 0;
}