//	dlgmodulesaveload.h : Implements the save and load module dialogs.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DLGMODULESAVELOAD_H_
#define __DLGMODULESAVELOAD_H_

#include "resource.h"

class CDlgModuleSaveLoad : public CAxDialogImpl<CDlgModuleSaveLoad>
{
public:
	BEGIN_MSG_MAP(CDlgModuleSaveLoad)
		MESSAGE_HANDLER(WM_INITDIALOG, OnInitDialog)
		MESSAGE_HANDLER(WM_SIZE, OnSize)
		MESSAGE_HANDLER(WM_SIZING, OnSizing)
		COMMAND_ID_HANDLER(IDOK, OnOK)
		COMMAND_ID_HANDLER(IDCANCEL, OnCancel)
		NOTIFY_HANDLER(IDC_TREE_MODULE, TVN_ITEMEXPANDING, OnItemExpandingTreeModule)
		NOTIFY_HANDLER(IDC_TREE_MODULE, TVN_SELCHANGED, OnSelChangedTree)
		COMMAND_HANDLER(IDC_DELETE, BN_CLICKED, OnClickedDelete)
	END_MSG_MAP()	
	CDlgModuleSaveLoad(){}
	~CDlgModuleSaveLoad(){}
	enum { IDD = IDD_DLGMODULESAVELOAD };

	LRESULT OnCancel(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	virtual LRESULT OnClickedDelete(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	virtual LRESULT OnItemExpandingTreeModule(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);
	virtual LRESULT OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	virtual LRESULT OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled) = 0;
	virtual LRESULT OnSelChangedTree(int idCtrl, LPNMHDR pnmh, BOOL& bHandled) = 0;
	virtual LRESULT OnSize(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnSizing(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);				

protected:
	long MapImageModuleType(int* pnImage, long* pnModuleType);

	CComDispatchDriverEx				m_ddExcel;		
	WINDOWPLACEMENT						m_wp_init, m_wp_tree_init, m_wp_preview_init, m_wp_edit_init, m_wp_saveload_init, m_wp_close_init;
	static WINDOWPLACEMENT				s_wp;							// to remember the size
};

class CDlgModuleSave : public CDlgModuleSaveLoad
{
public:
	LRESULT OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);	
	LRESULT OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT OnSelChangedTree(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);

protected:
	struct module_details {
		int								m_nModuleSelected;				// (one-based) element of m_vector corresponding to the module to save
		estring							m_szModuleText;
		int								m_nModuleType;					// e.g. vbext_ct_StdModule, vbext_ct_ClassModule
	};

	std::vector<CAdapt<CComPtr<IDispatch> >	>	m_vector;				// list of code module pointers	
	module_details								m_md;	
};

class CDlgModuleLoad : public CDlgModuleSaveLoad
{
public:
	LRESULT OnClickedDelete(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT OnInitDialog(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnItemExpandingTreeModule(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);
	LRESULT OnOK(WORD wNotifyCode, WORD wID, HWND hWndCtl, BOOL& bHandled);
	LRESULT OnSelChangedTree(int idCtrl, LPNMHDR pnmh, BOOL& bHandled);
	LRESULT OnSize(UINT uMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);

protected:
	int									m_nModuleSelected;	
	HTREEITEM							m_hModuleNodeSelected;	
	WINDOWPLACEMENT						m_wp_combo_init, m_wp_downloadto_init, m_wp_delete_init;
};

#endif
