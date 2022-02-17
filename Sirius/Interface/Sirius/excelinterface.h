//	excelinterface.h: interface for the CExcelInterface class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _EXCELINTERFACE_H
#define _EXCELINTERFACE_H

#pragma once

class CExcelInterface  
{
protected:
	class open_workbook : public CComDispatchDriverEx
	{
	public:
		open_workbook(IDispatch* lp) : CComDispatchDriverEx(lp){}				
		virtual ~open_workbook()
		{			
			this->PutPropertyByName(L"Saved", &CComVariant(VARIANT_TRUE));			
			CComVariant avClose[] = {VARIANT_MISSING, VARIANT_MISSING, VARIANT_FALSE};			
			this->InvokeN(L"Close", avClose, 3);
		}
	};
		
	class range_offset
	{
	public:		
		range_offset(void) : m_bValid(true)
		{
			m_nRowOffset = m_nColOffset = m_nRows = m_nCols = m_nRowLimit = m_nColLimit = 0;			
		}
		range_offset(int nRowOffset, int nColOffset, int nRows, int nCols, int nRowLimit, int nColLimit) : m_bValid(true)
		{			
			m_nRowOffset = nRowOffset;
			m_nColOffset = nColOffset;
			m_nRows = nRows;
			m_nCols = nCols;
			m_nRowLimit = nRowLimit;
			m_nColLimit = nColLimit;
		}
		range_offset(int nRowOffset, int nColOffset, int nRowLimit, int nColLimit) : m_bValid(true)
		{			
			m_nRowOffset = nRowOffset;
			m_nColOffset = nColOffset;
			m_nRows = 1;
			m_nCols = 1;
			m_nRowLimit = nRowLimit;
			m_nColLimit = nColLimit;
		}
		range_offset(bool bValid) : m_bValid(bValid)
		{
			m_nRowOffset = m_nColOffset = m_nRows = m_nCols = m_nRowLimit = m_nColLimit = 0;
		}
		HRESULT							MapRange(const CComDispatchDriverEx& ddRangeIn, CComDispatchDriverEx& ddRangeOut) const;
		HRESULT							GetFormula(const CComDispatchDriverEx& ddRange, std::string* pszOut) const;
		static HRESULT					GetCreateFormula(const CComDispatchDriverEx& ddRange, std::vector<range_offset>& vtrRangeOffset, const CLSID& CLSID_Singular, std::string* pszOut);
		HRESULT							GetCreateFormula(const CComDispatchDriverEx& ddRange, const CLSID& CLSID_Singular, std::string* pszOut) const;
		int								GetRows(void) const;
		int								GetColumns(void) const;
		HRESULT							GetShortAddress(const CComDispatchDriverEx& ddRange, std::string* pszOut) const;		
		void							RemoveTopRow(void);
		range_offset					TopLeft(void) const;
		HRESULT							Union(const range_offset& ro);
	
	protected:						
		bool							m_bValid;						// provides a mechanism to invalidate a range offset object
		int								m_nRowOffset;
		int								m_nColOffset;
		int								m_nRows;
		int								m_nCols;
		int								m_nRowLimit;
		int								m_nColLimit;		
	};
	
public:
	static const long					vbext_ct_StdModule;
	static const long					vbext_ct_ClassModule;
	
	CExcelInterface();
	virtual								~CExcelInterface();

	static void							Flatten(CComDispatchDriverEx& ddObject, VARIANT* pResult);	
	static void							Flatten(const CComVariant& vIn, VARIANT* pResult);
	template<class ISingular> static void Flatten(CComPtr<ISingular> sp, VARIANT* pResult)
	{
		Flatten(sp.p, pResult);
	}	
	template<class T> static void		Flatten(const std::vector<T>& v, VARIANT* pResult)
	{
		CParameterMap pm;
		pm.SetValue(v);		
		pm.GetValue(pResult);
	}	
	static HRESULT						GetAddress(const CComDispatchDriverEx& ddRange, std::string* pszOut);
	static HRESULT						GetAddressOfRange(const CComDispatchDriverEx& ddRange, std::string* pszOut);
	static HRESULT						GetCallerCellValue(CComVariant* pResult);
	static std::string					GetCellBasedHandle(const std::string& szPrefix);
	static HRESULT						GetWorkbookName(std::string* pszName);
	static HWND							GetWindow(void);
	static HRESULT						InitialiseExcel(const VARIANT& ExcelApplicationInstance);
	static HRESULT						InsertObject(const std::string& szObject, const std::string& szAt);
	static HRESULT						InsertProduct(const std::string& szProductType, const std::string& szAt);
	static HRESULT						VbaKey(const std::string& szFileName, std::string* pszKey);
	static HRESULT						UpdateSiriusStatusCaption();
	static HRESULT						Resize(const VARIANT& v, VARIANT* pv);
	
protected:
	static const char					s_szSiriusMenuTitle[];
	static const long					xlAutomatic;
	static const long					xlContinuous;		
	static const long					xlDiagonalDown;
	static const long					xlDiagonalUp;
	static const long					xlEdgeBottom;
	static const long					xlEdgeLeft;
	static const long					xlEdgeRight;
	static const long					xlEdgeTop;	
	static const long					xlInsideHorizontal;
	static const long					xlInsideVertical;
	static const long					xlMedium;
	static const long					xlNone;
	static const long					xlThin;
	static const long					msoControlPopup;
	static const long					msoControlButton;
	static const long					xlValidateList;
	static const long					xlValidAlertInformation;
	static const long					xlBetween;
	static const long					xlValues;
	static const long					xlPart;
	static const long					xlByRows;
	static const long					xlNext;
	static const long					xlShared;
	static const long					xlWBATWorksheet;
	static const long					xlFormulas;
	static const long					xlFormats;
	static const long					xlUp;
	static const long					xlRight;
	static const long					xlCenter;
	static const long					xlLeft;
	static const long					xlBottom;
	static const long					xlThick;
	static const long					xlDouble;
	
	static HRESULT						Autofit(const CComDispatchDriverEx& ddRange, const range_offset& ro);
	static HRESULT						CleanUpMenu(void);
	static HRESULT						DrawBorder(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nColourIndex);	
	static HRESULT						DrawLines(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nEdge, long nLineStyle, long nWeight = xlThin, long nColourIndex = xlAutomatic);
	static HRESULT						DrawLinesHelper(const CComDispatchDriverEx& ddRange, long nEdge, long nLineStyle, long nWeight = xlThin, long nColourIndex = xlAutomatic);		
	static void							FormatLikeObject(CComDispatchDriverEx& ddRange, long nRows, long nCols, long nRangeRows, long nRangeCols);
	static void							FormatLikeProduct(CComDispatchDriverEx& ddRange, long nRows, long nCols, long nRangeRows, long nRangeCols);
	static HRESULT						GetBoundary(CComPtr<IDispatch> spFormulaCell, CComPtr<IDispatch>& spBoundary);
	static HRESULT						GetLineStyle(CComPtr<IDispatch> spRange, long nRowOffset, long nColOffset, long nBorder, long* pnLineStyle);	
	static HRESULT						GetWritableRange(const std::string& szAddress, CComDispatchDriverEx& ddRangeOut, long* pnRangeRows, long* pnRangeCols);	
	static HRESULT						InitialiseExcelDrawMenuOption(/*const*/ CComDispatchDriverEx/*&*/ddControls, long nType, BOOL BeginGroup, long FaceID, LPCOLESTR Caption, LPCOLESTR OnAction, bool bEnabled, CAdapt<CComPtr<IDispatch> >* pout);
	static HRESULT						InsertObject(const CLSID& clsid, CComDispatchDriverEx& ddRange, long nRangeRows, const long nRangeCols, long* pnRowOffset, long nColOffset, long nParentRowOffset = -1, long nParentColOffset = -1);	
	static HRESULT						_Resize(const VARIANT& v, VARIANT* pv);
	static HRESULT						SetAlignment(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nHorizontalAlignment, long nVerticalAlignment);
	static HRESULT						SetBackColour(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nColour);
	static HRESULT						SetFont(const CComDispatchDriverEx& ddRange, const range_offset& ro, CComVariant vFontName, long nSize, long nColour, bool bBold);
	static void							VbaKeyDecryptKey(const char *encrypt_pass, char *s);
	static LPSTR						VbaKeyFindEncryptKey(const char *filename);
	static int							VbaKeyMemStr(const char *buf, const char *s, size_t size);		
	static HRESULT						Write(const CComDispatchDriverEx& ddRange, const range_offset& ro, const CComVariant& vIn);	
	static HRESULT						WriteComment(const CComDispatchDriverEx& ddRange, const range_offset& ro, const CComVariant& vIn);
	static HRESULT						WriteValidation(const CComDispatchDriverEx& ddRange, const range_offset& ro, const CComBSTR& sParameterName, const std::string& szList, const std::string& szInitialValue);
};

#endif
