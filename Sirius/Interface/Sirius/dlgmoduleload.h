//	dlgmoduleload.h : Declaration of the CDlgModuleLoad
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DLGMODULELOAD_H_
#define __DLGMODULELOAD_H_

#include "resource.h"

class CDlgModuleLoad : public CAxDialogImpl<CDlgModuleLoad>
{
public:
	BEGIN_MSG_MAP(CDlgModuleLoad)
		MESSAGE_HANDLER(WM_INITDIALOG, OnInitDialog)
		MESSAGE_HANDLER(WM_SIZE, OnSize)
		MESSAGE_HANDLER(WM_SIZING, OnSizing)
		COMMAND_ID_HANDLER(IDOK, OnOK)
		COMMAND_ID_HANDLER(IDCANCEL, OnCancel)
		NOTIFY_HANDLER(IDC_TREE_MODULE, TVN_ITEMEXPANDING, OnItemExpandingTreeModule)
		NOTIFY_HANDLER(IDC_TREE_MODULE, TVN_SELCHANGED, OnSelChangedTreeModule)
		COMMAND_HANDLER(IDC_DELETE, BN_CLICKED, OnClickedDelete)
	END_MSG_MAP()	
	CDlgModuleLoad(){}
	~CDlgModuleLoad(){}
	enum { IDD = IDD_DLGMODULELOAD };

protected:	
	LRESULT OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);	
	LRESULT OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);	
	LRESULT OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT OnItemExpandingTreeModule(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);
	LRESULT OnSelChangedTreeModule(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);
	LRESULT OnClickedDelete(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT OnSize(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnSizing(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	long MapImageModuleType(int* pnImage, long* pnModuleType);

	CComDispatchDriverEx				m_ddExcel;	
	int									m_nModuleSelected;	
	HTREEITEM							m_hModuleNodeSelected;
	WINDOWPLACEMENT						m_wp_init, m_wp_tree_init, m_wp_preview_init, m_wp_edit_init, m_wp_combo_init, m_wp_downloadto_init, m_wp_delete_init, m_wp_load_init, m_wp_close_init;
	static WINDOWPLACEMENT				s_wp;							// to remember the size
};

#endif
