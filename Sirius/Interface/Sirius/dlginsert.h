//	dlginsert.h : Declaration of the CDlgInsert
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DLGINSERT_H_
#define __DLGINSERT_H_

#include "resource.h"

class CDlgInsert :  public CAxDialogImpl<CDlgInsert>
{
public:
	enum dialog_type {
		insert_object = 1,
		insert_product = 2
	};

	BEGIN_MSG_MAP(CDlgInsert)
		MESSAGE_HANDLER(WM_INITDIALOG, OnInitDialog)
		COMMAND_ID_HANDLER(IDOK, OnOK)
		COMMAND_ID_HANDLER(IDCANCEL, OnCancel)
		COMMAND_HANDLER(IDC_REFRESH, BN_CLICKED, OnClickedRefresh)
		COMMAND_HANDLER(IDC_LIST, LBN_DBLCLK, OnDblClkList)
	END_MSG_MAP()
	
	CDlgInsert(){}
	~CDlgInsert(){}
	enum { IDD = IDD_DLGINSERT };

protected:
	HRESULT								CheckInsertAt(void);	
	HRESULT								LoadList(void);
	LRESULT								OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnClickedRefresh(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);	
	LRESULT								OnDblClkList(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT								OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT								OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	
	dialog_type							m_dt;
	std::string							m_szAt;							// excel range address at which the item is inserted
	CComDispatchDriverEx				m_ddExcel;
	CParameterMap						m_pmList;						// items currently in the list box
};

#endif
