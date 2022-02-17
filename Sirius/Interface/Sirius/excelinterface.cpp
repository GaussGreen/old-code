//	excelinterface.cpp: implementation of the CExcelInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "excelinterface.h"
#include "siriusapplication.h"
#include "objectmanager.h"

/*static*/ const char					CExcelInterface::s_szSiriusMenuTitle[] = "&Sirius";
/*static*/ const long					CExcelInterface::xlContinuous = 1;
/*static*/ const long					CExcelInterface::xlMedium = -4138;
/*static*/ const long					CExcelInterface::xlAutomatic = -4105;
/*static*/ const long					CExcelInterface::xlEdgeRight = 10;
/*static*/ const long					CExcelInterface::xlEdgeBottom = 9;
/*static*/ const long					CExcelInterface::xlEdgeLeft = 7;
/*static*/ const long					CExcelInterface::xlEdgeTop = 8;
/*static*/ const long					CExcelInterface::xlNone = -4142;
/*static*/ const long					CExcelInterface::xlDiagonalDown = 5;
/*static*/ const long					CExcelInterface::xlDiagonalUp = 6;
/*static*/ const long					CExcelInterface::xlInsideVertical = 11;
/*static*/ const long					CExcelInterface::xlInsideHorizontal = 12;
/*static*/ const long					CExcelInterface::xlThin = 2;
/*static*/ const long					CExcelInterface::msoControlPopup = 10;
/*static*/ const long					CExcelInterface::msoControlButton = 1;
/*static*/ const long					CExcelInterface::vbext_ct_StdModule = 1;
/*static*/ const long					CExcelInterface::vbext_ct_ClassModule = 2;
/*static*/ const long					CExcelInterface::xlValidateList = 3;
/*static*/ const long					CExcelInterface::xlValidAlertInformation = 3;
/*static*/ const long					CExcelInterface::xlBetween = 1;
/*static*/ const long					CExcelInterface::xlValues = -4163;
/*static*/ const long					CExcelInterface::xlPart = 2;
/*static*/ const long					CExcelInterface::xlByRows = 1;
/*static*/ const long					CExcelInterface::xlNext = 1;
/*static*/ const long					CExcelInterface::xlShared = 2;
/*static*/ const long					CExcelInterface::xlWBATWorksheet = -4167;
/*static*/ const long					CExcelInterface::xlFormulas = -4123;
/*static*/ const long					CExcelInterface::xlFormats = -4122;
/*static*/ const long					CExcelInterface::xlUp = -4162;
/*static*/ const long					CExcelInterface::xlRight = -4152;
/*static*/ const long					CExcelInterface::xlCenter = -4108;
/*static*/ const long					CExcelInterface::xlLeft = -4131;
/*static*/ const long					CExcelInterface::xlBottom = -4107;
/*static*/ const long					CExcelInterface::xlThick = 4;
/*static*/ const long					CExcelInterface::xlDouble = -4119;


#ifdef BLOCK_SIZE
#error macro 'BLOCK_SIZE' already defined
#endif

#define BLOCK_SIZE 32768


////////////////////////////////////////////////////////////////////////////
//	range_offset implementation
//
//	transforms an input range object to an output range object according to
//	the data in this object
//
HRESULT CExcelInterface::range_offset::MapRange(const CComDispatchDriverEx& ddRangeIn, CComDispatchDriverEx& ddRangeOut) const
{	
	CComVariant							v1;
	CComVariant							v2;
	HRESULT								hr;
	CComVariant							params[2];
		
	if (m_nRowOffset < 0 || m_nRowOffset >= m_nRowLimit || m_nColOffset < 0 || m_nColOffset >= m_nColLimit) return E_FAIL;
			
	// implement "Range.Offset"
	params[0].vt = VT_I4;
	params[0].lVal = m_nColOffset;
	params[1].vt = VT_I4;
	params[1].lVal = m_nRowOffset;
	if (hr = ddRangeIn.GetPropertyByName(L"Offset", params, 2, &v1)) return hr;
	if (v1.vt != VT_DISPATCH) return E_FAIL;
	
	// implement "Range.Resize"
	params[0].vt = VT_I4;
	params[0].lVal = __max(0, __min(m_nCols, m_nColLimit - m_nColOffset));
	params[1].vt = VT_I4;
	params[1].lVal = __max(0, __min(m_nRows, m_nRowLimit - m_nRowOffset));
	if (hr = CComDispatchDriverEx(v1.pdispVal).GetPropertyByName(L"Resize", params, 2, &v2)) return hr;
	if (v2.vt != VT_DISPATCH) return E_FAIL;
	
	// set pddRangeOut	
	ddRangeOut = v2.pdispVal;	
	return S_OK;
}
//
//	Returns the short address associated with the data in this object and an input origin range.
//	This function is used for developing formulae that take a range as an argument.
//
HRESULT CExcelInterface::range_offset::GetShortAddress(const CComDispatchDriverEx& ddRange, std::string* pszOut) const
{
	if (!m_bValid){
		pszOut->clear();
	} else {	
		CComDispatchDriverEx				ddRangeMap;	
		HRESULT								hr;
		CComVariant							v;
		CComVariant							params[5];						// parameter list of the Address property on the excel range object
				
		if (hr = MapRange(ddRange, ddRangeMap)) return hr;	
		params[0].vt = VT_EMPTY;								// relative to			
		params[1].vt = VT_BOOL;	params[1].boolVal = false;		// external reference
		params[2].vt = VT_I2; params[2].iVal = 1;				// reference style
		params[3].vt = VT_BOOL;	params[3].boolVal = false;		// column absolute
		params[4].vt = VT_BOOL;	params[4].boolVal = false;		// row absolute
		if (hr = ddRangeMap.GetPropertyByName(L"Address", params, 5, &v)) return hr;
		if (v.vt != VT_BSTR) return E_FAIL;
		pszOut->assign((LPCTSTR)_bstr_t(v.bstrVal));
	}
	return S_OK;
}
//
//	Returns a string that represents a formula for the address of the range
//	implied by pRangeIn and the data in this object
//
HRESULT CExcelInterface::range_offset::GetFormula(const CComDispatchDriverEx& ddRange, std::string* pszOut) const
{	
	HRESULT								hr;
	std::string							szAddress;
	
	if (hr = GetShortAddress(ddRange, &szAddress)) szAddress = "#REF!";	
	pszOut->assign("=" + szAddress);
	return hr;
}
//
//	Returns a string that represents the create handle formula for a range. I.e. it returns MLCreate(ObjectName, ...)
//
HRESULT CExcelInterface::range_offset::GetCreateFormula(const CComDispatchDriverEx& ddRange, const CLSID& CLSID_Singular, std::string* pszOut) const
{			
	HRESULT								hrRet = S_OK;
	HRESULT								hr;
	std::string							szAddress;
	std::string							szObject;						// object name associated with the input CLSID

	if (hr = GetShortAddress(ddRange, &szAddress)){
		hrRet = hr;
		szAddress = "#REF!";
	}
	if (hr = CParameterMap::CLSIDToObjectName(CLSID_Singular, &szObject)){
		hrRet = hr;
		szObject = "?";
	}
	pszOut->assign("=MLCreate(\"" + szObject + "\", " + szAddress + ")");
	return hrRet;
}
//
//	Returns a string that represents the create handle formula for a multiargument handle creation function.
//
/*static*/ HRESULT CExcelInterface::range_offset::GetCreateFormula(const CComDispatchDriverEx& ddRange, std::vector<range_offset>& vtrRangeOffset, const CLSID& CLSID_Singular, std::string* pszOut)
//	ddRange - defines the origin against which all range_offsets are calculated
//	vtrRangeOffset - s.t. each element along with ddRange defines the area of a handle parameter
//	pszOut - create function (returned)
{
	HRESULT								hrRet = S_OK;
	HRESULT								hr;
	std::string							szObject;						// object name associated with the input CLSID
	std::stringstream					ss;

	ss << "=MLCreate";
	if (hr = CParameterMap::CLSIDToObjectName(CLSID_Singular, &szObject)){
		hrRet = hr;
		szObject = "?";
	}
	ss << szObject << "(";	
	for (std::vector<range_offset>::iterator it = vtrRangeOffset.begin(); it < vtrRangeOffset.end(); it++){
		std::string szAddress;
		if (hr = it->GetShortAddress(ddRange, &szAddress)){
			hrRet = hr;
			szAddress = "#REF!";
		}
		if (it != vtrRangeOffset.begin()) ss << ",";		
		ss << szAddress;
	}
	ss << ")";
	pszOut->assign(ss.str());	
	return hrRet;
}
int CExcelInterface::range_offset::GetRows(void) const
{
	return m_nRows;
}
int CExcelInterface::range_offset::GetColumns(void) const
{
	return m_nCols;
}
//	Removes the top row from the object
void CExcelInterface::range_offset::RemoveTopRow(void)
{
	if (!m_nRows) return;
	m_nRowOffset++;
	m_nRows--;
}
//	Returns a range_offset that represents the top left of an input offset
CExcelInterface::range_offset CExcelInterface::range_offset::TopLeft(void) const
{
	range_offset						ro;
	
	ro.m_nRowOffset = m_nRowOffset;
	ro.m_nColOffset = m_nColOffset;
	ro.m_nRows = 1;
	ro.m_nCols = 1;
	ro.m_nRowLimit = m_nRowLimit;
	ro.m_nColLimit = m_nColLimit;
	return ro;
}
//
//	Sets the current object as the union area of the current object and an input object
//
HRESULT CExcelInterface::range_offset::Union(const range_offset& ro)
{
	if (m_nRowLimit != ro.m_nRowLimit || m_nColLimit != m_nColLimit) return E_FAIL;
	// left boundary
	if (m_nColOffset > ro.m_nColOffset){
		ATLASSERT(false);	// ToDo - test
		m_nCols += (m_nColOffset - ro.m_nColOffset);
		m_nColOffset = ro.m_nColOffset;
	}
	// top boundary
	if (m_nRowOffset > ro.m_nRowOffset){
		ATLASSERT(false);	// ToDo - test
		m_nRows += (m_nRowOffset - ro.m_nRowOffset);
		m_nRowOffset = ro.m_nRowOffset;
	}
	// right boundary
	if (m_nColOffset + m_nCols < ro.m_nColOffset + ro.m_nCols){
		m_nCols = ro.m_nColOffset + ro.m_nCols - m_nColOffset;
	}		
	// bottom boundary
	if (m_nRowOffset + m_nRows < ro.m_nRowOffset + ro.m_nRows){
		m_nRows = ro.m_nRowOffset + ro.m_nRows - m_nRowOffset;
	}
	return S_OK;
}


////////////////////////////////////////////////////////////////////////////
//	CExcelInterface implementation
//
//	Construction / Destruction
//
CExcelInterface::CExcelInterface()
{
	// Don't ever create one of these objects - just use the statics
	ATLASSERT(false);
}
CExcelInterface::~CExcelInterface()
{
	// Don't ever create one of these objects - just use the statics
	ATLASSERT(false);
}
//
//	autofits the columns in an input Excel range constrained by a given range offset
//
/*static*/ HRESULT CExcelInterface::Autofit(const CComDispatchDriverEx& ddRange, const range_offset& ro)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;						// range to write to	
	CComVariant							vEntireColumn;
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;			
	if (hr = ddRangeUse.GetPropertyByName(L"EntireColumn", &vEntireColumn)) return hr;				
	CComDispatchDriverEx ddEntireColumn(vEntireColumn.pdispVal);	
	return ddEntireColumn.Invoke0(L"AutoFit");
}
//
//	Removes the 'Sirius' menu item
//
/*static*/ HRESULT CExcelInterface::CleanUpMenu(void)
{
	HRESULT								hr;
	long								nOld, nNew;	
	CComDispatchDriverEx				ddExcel;	
	CComVariant							vCommandBar, vControls, vCount, vControl;
		
	if (hr = _Module.GetExcel(ddExcel)) return E_FAIL;	
	if (ddExcel.GetPropertyByName(L"CommandBars", &CComVariant(L"Worksheet Menu Bar"), 1, &vCommandBar)) return E_FAIL;
	do {
		nOld = nNew = 0;
		if (CComDispatchDriverEx(vCommandBar.pdispVal).GetPropertyByName(L"Controls", &vControls)) return E_FAIL;
		if (CComDispatchDriverEx(vControls.pdispVal).GetPropertyByName(L"Count", &vCount)) return E_FAIL;
		if (vCount.vt != VT_I4) return E_FAIL;
		nOld = vCount.lVal;		
		if (!CComDispatchDriverEx(vCommandBar.pdispVal).GetPropertyByName(L"Controls", &CComVariant(s_szSiriusMenuTitle), 1, &vControl)){
			CComDispatchDriverEx(vControl.pdispVal).Invoke0(L"Delete");
		}
		if (CComDispatchDriverEx(vControls.pdispVal).GetPropertyByName(L"Count", &vCount)) return E_FAIL;
		if (vCount.vt != VT_I4) return E_FAIL;
		nNew = vCount.lVal;					
	} while (nOld != nNew);

	return S_OK;
}
//
//	draws a thick border around an input Excel range constrained by a given
//	range offset
//
/*static*/ HRESULT CExcelInterface::DrawBorder(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nColourIndex)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;
	CComVariant							params[1];
	CComVariant							vBorders;						// "Borders" property return value
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;
	DrawLinesHelper(ddRangeUse, xlDiagonalDown, xlNone);
	DrawLinesHelper(ddRangeUse, xlDiagonalUp, xlNone);	
	DrawLinesHelper(ddRangeUse, xlEdgeLeft, xlContinuous, xlMedium, nColourIndex);
	DrawLinesHelper(ddRangeUse, xlEdgeTop, xlContinuous, xlMedium, nColourIndex);
	DrawLinesHelper(ddRangeUse, xlEdgeBottom, xlContinuous, xlMedium, nColourIndex);
	DrawLinesHelper(ddRangeUse, xlEdgeRight, xlContinuous, xlMedium, nColourIndex);
	DrawLinesHelper(ddRangeUse, xlInsideVertical, xlNone);
	DrawLinesHelper(ddRangeUse, xlInsideHorizontal, xlNone);
	return S_OK;
}
//
//	Draws line(s) in an input range
//
/*static*/ HRESULT CExcelInterface::DrawLines(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nEdge, long nLineStyle, long nWeight/* = xlThin*/, long nColourIndex/* = xlAutomatic*/)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;
	CComVariant							params[1];
	CComVariant							vBorders;						// "Borders" property return value
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;
	DrawLinesHelper(ddRangeUse, nEdge, nLineStyle, nWeight, nColourIndex);
	return S_OK;
}
//
//	draws a border line in an input range
//
/*static*/ HRESULT CExcelInterface::DrawLinesHelper(const CComDispatchDriverEx& ddRange, long nEdge, long nLineStyle, long nWeight/* = xlThin*/, long nColourIndex/* = xlAutomatic*/)
{
	HRESULT								hr;
	CComVariant							vParameter(nEdge, VT_I4);
	CComVariant							vBorders;						// "Borders" property return value
		
	if (hr = ddRange.GetPropertyByName(L"Borders", &vParameter, 1, &vBorders)) return hr;
	if (vBorders.vt != VT_DISPATCH) return E_FAIL;
	CComDispatchDriverEx ddBorders(vBorders.pdispVal);
	if (hr = ddBorders.PutPropertyByName(L"LineStyle", &CComVariant(nLineStyle, VT_I4))) return hr;
	if (nLineStyle != xlNone){		
		if (hr = ddBorders.PutPropertyByName(L"Weight", &CComVariant(nWeight, VT_I4))) return hr;
		if (hr = ddBorders.PutPropertyByName(L"ColorIndex", &CComVariant(nColourIndex, VT_I4))) return hr;
	}
	return S_OK;
}
//
//  The Flatten functions are intended to be used to return results back to Excel.
//
void CExcelInterface::Flatten(CComDispatchDriverEx& ddObject, VARIANT* pResult)
{
	//  Here we flatten an object into the parameter list of the appropriate MLCreate[...] function.
	std::vector<CParameterMap>			vpm;
	std::vector<CComVariant>			vv;
	CComVariant							v;

	if (ddObject.GetPropertyByName(L"Value", &v)) propagate_error;	
	if (CParameterMap::ArrayToVector(v, &vv, &vpm, false/*i.e. don't decompose the 1 vector element case*/)) propagate_error;
	if (vv.size() == 1){
		CExcelInterface::Flatten(vv[0], pResult);
		return;
	}
	
	std::vector<std::string>			aszParameterNames;
	CParameterMap						pmResult;
	CParameterMap						pmNames;
	long								nParameter = 0;

	pmNames.SetAutoGrow(true);
	g_pApplication->GetCreateParameterNames(ddObject, vpm, &aszParameterNames, NULL, NULL);
	for (std::vector<CComVariant>::iterator it = vv.begin(); it != vv.end(); it++, nParameter++){
		CComVariant vOut;
		Flatten(*it, &vOut);
		
		CParameterMap pm;
		pm.SetValue(vOut);
		if (nParameter >= aszParameterNames.size()){
			std::stringstream ss;
			ss << "Parameter " << nParameter + 1;
			pmNames.SetValue(pmResult.GetRows(), 0, ss);
		} else {
			pmNames.SetValue(pmResult.GetRows(), 0, aszParameterNames[nParameter]);			
		}		
		pmResult.AddToEnd(pm);
	}
	if (pmResult.IsBlank()) throw "No data found";	
	pmNames.AddToRHS(pmResult);	
	pmResult.Attach(&pmNames);	
	pmResult.ReplaceBlanksWithSpaces();	// Don't just apply to pmNames. (e.g. Extra zeros will be added in the correlation matrix GetValue).
	if (pmResult.GetValue(pResult)) propagate_error;
}
/*static*/ void CExcelInterface::Flatten(const CComVariant& vIn, VARIANT* pResult)
{
	CParameterMap						pmResult;
		
	if (vIn.vt == VT_DISPATCH){		
		if (g_pApplication->GetObjectManager().IsSiriusCollection(vIn.pdispVal)){
			// Return an array of handles.			
			CComDispatchDriverEx		ddCollection(vIn.pdispVal);			
			CComVariant					vCount;
			
			if (ddCollection.GetPropertyByName(L"Count", &vCount)) propagate_error;
			if (vCount.vt != VT_I4 || pmResult.SetSize(vCount.lVal, 1)) throw "Unhandled exception in CExcelInterface::Flatten";
			if (!vCount.lVal) throw "No data found";
			for (long nItem = 1; nItem <= vCount.lVal; nItem++){
				CComVariant vObject;							
				if (ddCollection.GetPropertyByName(L"Item", &CComVariant(nItem), 1, &vObject)) propagate_error;
				estring szHandle = g_pApplication->GetObjectManager().InsertObject(CComPtr<IDispatch>(vObject.pdispVal), false);
				if (pmResult.SetValue(nItem - 1, 0, szHandle)) propagate_error;
			}
			if (pmResult.GetValue(pResult)) propagate_error;
		} else {						
			// Return one handle
			g_pApplication->GetObjectManager().InsertObject(CComPtr<IDispatch>(vIn.pdispVal), true).GetValue(pResult);
		}		
		return;
	}

	pmResult.SetValue(vIn);
	pmResult.Substitute("\t", "    ", true, false, -1, -1);			// Excel does not display tabs at the start of strings. Therefore replace with spaces.

	// If pmResult is a scalar string then split it up into "\r\n".	
	if (pmResult.IsScalar()){
		estring sz;
		if (!pmResult.GetValue(&sz)){
			std::vector<std::string> vector;				
			sz.Split("\r\n", &vector);
			if (vector.size() > 1){	
				pmResult.SetValue(vector);
				if (pmResult.GetValue(pResult)) propagate_error;
				return;
			}
		}
	}

	for (long nRow = 0; nRow < pmResult.GetRows(); nRow++){
		for (long nCol = 0; nCol < pmResult.GetCols(); nCol++){
			if (pmResult.GetElementPtr(nRow, nCol)->vt == VT_DISPATCH){
				if (pmResult.GetElementPtr(nRow, nCol)->pdispVal){
					try {
						bool bExternal = pmResult.IsScalar() ? true : false;
						estring sz = g_pApplication->GetObjectManager().InsertObject(CComPtr<IDispatch>(pmResult.GetElementPtr(nRow, nCol)->pdispVal), bExternal);
						pmResult.GetElementPtr(nRow, nCol)->Attach(&sz.GetValue());
					} catch (...){
						// Try a key (the object might have been declared as non-insertable into the maintained collection)						
						if (g_pApplication->GetObjectManager().IsSerialisable(CComPtr<IDispatch>(pmResult.GetElementPtr(nRow, nCol)->pdispVal))){
							pmResult.GetElementPtr(nRow, nCol)->Attach(&CComObjectCollectionSerialisableKey(pmResult.GetElementPtr(nRow, nCol)->pdispVal).GetKeyAndObjectName(false, false).GetValue());
						}
					}
				} else {
					// object is null
					CComVariant("").Detach(pmResult.GetElementPtr(nRow, nCol));
				}
			} else if (pmResult.GetElementPtr(nRow, nCol)->vt == VT_EMPTY){
				// Set to a blank string (since, by default, Excel represents empty variants as zeros).
				CComVariant("").Detach(pmResult.GetElementPtr(nRow, nCol));
			}
		}
	}

	if (pmResult.GetValue(pResult)) propagate_error;
}
//
//	places Sirius-type formatting around an input object
//
/*static*/ void CExcelInterface::FormatLikeObject(CComDispatchDriverEx& ddRange, long nRows, long nCols, long nRangeRows, long nRangeCols)
{
	DrawBorder(ddRange, range_offset(1, 0, nRows, nCols, nRangeRows, nRangeCols), 5);
	DrawLines(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), xlEdgeBottom, xlDouble, xlThick, 0);	
	SetBackColour(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), 5);
	SetBackColour(ddRange, range_offset(2, 0, nRows - 1, nCols, nRangeRows, nRangeCols), 2);
	SetFont(ddRange, range_offset(0, 0, nRows + 2, nCols, nRangeRows, nRangeCols), L"Tahoma", 8, 0, false);
	SetFont(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), L"Tahoma", 8, 6, true);	
	Autofit(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols));
}
//
//	places Sirius-type formatting around an input object
//
/*static*/ void CExcelInterface::FormatLikeProduct(CComDispatchDriverEx& ddRange, long nRows, long nCols, long nRangeRows, long nRangeCols)
{
	DrawBorder(ddRange, range_offset(1, 0, nRows, nCols, nRangeRows, nRangeCols), 5);
	DrawLines(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), xlEdgeBottom, xlDouble, xlThick, 0);	
	DrawLines(ddRange, range_offset(2, 1, nRows - 1, nCols - 1, nRangeRows, nRangeCols), xlInsideHorizontal, xlContinuous, xlThin, 0);		
	DrawLines(ddRange, range_offset(2, 1, nRows - 1, nCols - 1, nRangeRows, nRangeCols), xlEdgeLeft, xlContinuous, xlMedium, 0);		
	SetBackColour(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), 5);
	SetBackColour(ddRange, range_offset(2, 0, nRows - 1, 1, nRangeRows, nRangeCols), 37);
	SetBackColour(ddRange, range_offset(2, 1, nRows - 1, nCols - 1, nRangeRows, nRangeCols), 2);
	SetFont(ddRange, range_offset(0, 0, nRows + 2, nCols, nRangeRows, nRangeCols), L"Tahoma", 8, 0, false);
	SetAlignment(ddRange, range_offset(0, 1, nRows + 1, nCols - 1, nRangeRows, nRangeCols), xlRight, xlBottom);
	SetFont(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), L"Tahoma", 8, 6, true);
	SetAlignment(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols), xlLeft, xlCenter);
	Autofit(ddRange, range_offset(1, 0, 1, nCols, nRangeRows, nRangeCols));
}
//
//	returns the long address of an input range
//
/*static*/ HRESULT CExcelInterface::GetAddress(const CComDispatchDriverEx& ddRange, std::string* pszOut)
{	
	HRESULT								hr;
	CComVariant							v;
	CComVariant							params[5];						// parameter list of the Address property on the excel range object
			
	params[0].vt = VT_EMPTY;								// relative to			
	params[1].vt = VT_BOOL;	params[1].boolVal = true;		// external reference
	params[2].vt = VT_I2; params[2].iVal = 1;				// reference style
	params[3].vt = VT_BOOL;	params[3].boolVal = true;		// column absolute
	params[4].vt = VT_BOOL;	params[4].boolVal = true;		// row absolute
	if (hr = ddRange.GetPropertyByName(L"Address", params, 5, &v)) return hr;
	if (v.vt != VT_BSTR) return E_FAIL;
	pszOut->assign((LPCTSTR)_bstr_t(v.bstrVal));
	return S_OK;
}
//
//	returns the long address of an input range
//
/*static*/ HRESULT CExcelInterface::GetAddressOfRange(const CComDispatchDriverEx& ddRange, std::string* pszOut)
{	
	HRESULT								hr;
	CComVariant							v;
	CComVariant							params[5];						// parameter list of the Address property on the excel range object
			
	params[0].vt = VT_EMPTY;								// relative to			
	params[1].vt = VT_BOOL;	params[1].boolVal = true;		// external reference
	params[2].vt = VT_I2; params[2].iVal = 1;				// reference style
	params[3].vt = VT_BOOL;	params[3].boolVal = true;		// column absolute
	params[4].vt = VT_BOOL;	params[4].boolVal = true;		// row absolute
	if (hr = ddRange.GetPropertyByName(L"Address", params, 5, &v)) return hr;
	if (v.vt != VT_BSTR) return E_FAIL;
	pszOut->assign((LPCTSTR)_bstr_t(v.bstrVal));
	return S_OK;
}
//
//	Takes an input range containing a single cell. We then examine any cells
//  bounded by a rectangular border underneath this cell. We return the smallest
//	rectangular union of this range along with the input range.
//
/*static*/ HRESULT CExcelInterface::GetBoundary(CComPtr<IDispatch> spFormulaCell, CComPtr<IDispatch>& spBoundary)
{
	HRESULT								hr;
	long								nLineStyle;
	long								nColumn = 1;
	long								nRow = 1;

	if (hr = GetLineStyle(spFormulaCell, 1, 0, xlEdgeTop, &nLineStyle)) return hr;
	if (nLineStyle == xlNone) return S_OK;	// spBoundary will be null

	long nLineStyleTop, nLineStyleRight;
	do {		
		if (hr = GetLineStyle(spFormulaCell, 1, nColumn, xlEdgeTop, &nLineStyleTop)) return hr;
		if (hr = GetLineStyle(spFormulaCell, 1, nColumn++, xlEdgeRight, &nLineStyleRight)) return hr;		
	} while (nLineStyleTop == nLineStyle && nLineStyleRight != nLineStyle);

	long nLineStyleLeft, nLineStyleBottom;
	do {		
		if (hr = GetLineStyle(spFormulaCell, nRow, 0, xlEdgeLeft, &nLineStyleLeft)) return hr;
		if (hr = GetLineStyle(spFormulaCell, nRow++, 0, xlEdgeBottom, &nLineStyleBottom)) return hr;		
	} while (nLineStyleLeft == nLineStyle && nLineStyleBottom != nLineStyle);

	CComVariant							vResize;
	CComVariant							avResize[] = {nColumn, nRow};
	if (hr = CComDispatchDriverEx(spFormulaCell).GetPropertyByName(L"Resize", avResize, 2, &vResize)) return hr;
	spBoundary = vResize.pdispVal;
	return S_OK;
}
//
//	Returns the formula in the cell used ultimately to call this function.
//
/*static*/ HRESULT CExcelInterface::GetCallerCellValue(CComVariant* pResult)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddExcel;
	CComDispatchDriverEx				ddCaller;
	CComVariant							v;

	if (hr = _Module.GetExcel(ddExcel)) return hr;
	if (hr = ddExcel.GetPropertyByName(L"Caller", ddCaller)) return hr;
	if (hr = ddCaller.GetPropertyByName(L"Value", &v)) return hr;
	return v.Detach(pResult);
}
//
//	Returns a handle of the form szPrefix::RangeAddress.
//	This is only valid for an Excel process.
//
/*static*/ std::string CExcelInterface::GetCellBasedHandle(const std::string& szPrefix)
{
	CComVariant							v;
	CComVariant							params[5];						// parameter list of the Address property on the excel range object	
	CComDispatchDriverEx				ddExcel;
	
	if (_Module.GetExcel(ddExcel) || ddExcel.GetPropertyByName(L"Caller", &v) || v.vt != VT_DISPATCH) throw CStringResource(IDS_VALID_ONLY_IN_EXCEL);
	params[0].vt = VT_EMPTY;								// relative to			
	params[1].vt = VT_BOOL;	params[1].boolVal = true;		// external reference
	params[2].vt = VT_I2; params[2].iVal = 1;				// reference style	
	params[3].vt = VT_BOOL;	params[3].boolVal = true;		// column absolute
	params[4].vt = VT_BOOL;	params[4].boolVal = true;		// row absolute
	if (CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Address", params, 5, &v) || v.vt != VT_BSTR) throw "Unhandled exception in CExcelInterface::GetCellBasedHandle";	
	return estring(szPrefix + "::" + estring(v));	
}
//
//	Returns the line style of a border.
//
/*static*/ HRESULT CExcelInterface::GetLineStyle(CComPtr<IDispatch> spRange, long nRowOffset, long nColOffset, long nBorder, long* pnLineStyle)
//	spRange - input cell
//	nRowOffset, nColOffset - such that the border linestyle is tested for spRange offset by these values
//	nBorder - cell border to examine (e.g. xlEdgeTop)
//	pnLineStyle - value returned
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRange(spRange);
	CComVariant							vOffset;
	CComVariant							vBorder;
	CComVariant							vLineStyle;
	CComVariant							avOffset[] = {nColOffset, nRowOffset};

	if (hr = ddRange.GetPropertyByName(L"Offset", avOffset, 2, &vOffset)) return hr;
	if (hr = CComDispatchDriverEx(vOffset.pdispVal).GetPropertyByName(L"Borders", &CComVariant(nBorder), 1, &vBorder)) return hr;
	if (hr = CComDispatchDriverEx(vBorder.pdispVal).GetPropertyByName(L"LineStyle", &vLineStyle)) return hr;
	ATLASSERT(vLineStyle.vt == VT_I4);
	*pnLineStyle = vLineStyle.lVal;
	return S_OK;
}
//
//	Returns the handle of the Excel window (zero if the dll is not attached to Excel)
//
/*static*/ HWND CExcelInterface::GetWindow(void)
{
	CComDispatchDriverEx				ddExcel;	
	CComVariant							v;
		
	if (_Module.GetExcel(ddExcel)) return NULL;
	if (ddExcel.GetPropertyByName(L"HWND", &v)) return NULL;
	if (v.ChangeType(VT_I4)) return NULL;
	return reinterpret_cast<HWND>(v.lVal);
}
//
//	Returns the name of the active workbook
//
/*static*/ HRESULT CExcelInterface::GetWorkbookName(std::string* pszName)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddExcel;
	CComDispatchDriverEx				ddWorkbook;
	CComVariant							vFullName;
	
	if (hr = _Module.GetExcel(ddExcel)) return hr;
	if (hr = ddExcel.GetPropertyByName(L"ActiveWorkbook", ddWorkbook)) return hr;
	if (hr = ddWorkbook.GetPropertyByName(L"FullName", &vFullName)) return hr;
	pszName->assign(estring(vFullName));
	return S_OK;
}
//
//	Takes an input variant that represents the address of a range in the
//	excel session attached to this process. We return a range object
//	corresponding to this address along with the number of rows and columns
//	in this range.
//
//	This is useful in determining the boundaries of any output to a worksheet.
//
/*static*/ HRESULT CExcelInterface::GetWritableRange(const std::string& szAddress, CComDispatchDriverEx& ddRangeOut, long* pnRangeRows, long* pnRangeCols)
//	szAddress - must resolve to a valid excel workbook address
//	ddRangeOut - range corresponding to the output address (returned)
//	ppnRangeRows - number of rows in ddRangeOut (returned)
//	ppnRangeCols - number of columns in ddRangeOut (returned)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddExcel;	
	CComVariant							v;

	if (hr = _Module.GetExcel(ddExcel)) return hr;			
	if (hr = ddExcel.GetPropertyByName(L"Range", &CComVariant(szAddress.c_str()), 1, &v)) return hr;
	if (v.vt != VT_DISPATCH) return E_FAIL;
	ddRangeOut = v.pdispVal;		

	// get nRangeRows and nRangeCols
	if (hr = ddRangeOut.GetPropertyByName(L"Rows", &v)) return hr;		
	if (hr = CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Count", &v)) return hr;				
	*pnRangeRows = v.lVal;

	if (hr = ddRangeOut.GetPropertyByName(L"Columns", &v)) return hr;
	if (hr = CComDispatchDriverEx(v.pdispVal).GetPropertyByName(L"Count", &v)) return hr;
	*pnRangeCols = v.lVal;

	// if the range is only one cell then set nRangeRows and nRangeCols to the area of the worksheet with the range being the top left corner
	if (*pnRangeRows == 1 && *pnRangeCols == 1){
		if (hr = CParameterMap::ExcelRangeToRowCol(ddRangeOut, pnRangeRows, pnRangeCols, NULL, NULL)) return hr;
		*pnRangeRows = 65536 - *pnRangeRows;
		*pnRangeCols = 256 - *pnRangeCols;
	}	
	return S_OK;
}
//
//	Draws the Excel Sirius menu and stores an input pointer to the current Excel process interface
//
/*static*/ HRESULT CExcelInterface::InitialiseExcel(const VARIANT& ExcelApplicationInstance)
{
	HRESULT								hr = S_OK;
	CComDispatchDriverEx				ddExcel;
	CComVariant							vCommandBar, vControls, vControl, vWindow, vIndex;
	CComBSTR							bstrVersionStatus;
	CComVariant							varTitle;

		
	if (!(hr = _Module.DeclareExcel(ExcelApplicationInstance))){
		if (!(hr = _Module.GetExcel(ddExcel))){			
			UpdateSiriusStatusCaption();
			if (!(hr = CleanUpMenu())){
				if (!(hr = ddExcel.GetPropertyByName(L"CommandBars", &CComVariant(L"Worksheet Menu Bar"), 1, &vCommandBar))){
					if (!(hr = CComDispatchDriverEx(vCommandBar.pdispVal).GetPropertyByName(L"Controls", &vControls))){						
						if (!(hr = CComDispatchDriverEx(vCommandBar.pdispVal).GetPropertyByName(L"Controls", &CComVariant(L"Window"), 1, &vWindow))){
							if (!(hr = CComDispatchDriverEx(vWindow.pdispVal).GetPropertyByName(L"Index", &vIndex))){
								CComVariant av[] = {true/*temporary*/, vIndex/*before*/, 0L/*parameter*/, 1L/*ID*/, msoControlPopup/*type*/};
								if (!(hr = CComDispatchDriverEx(vControls.pdispVal).InvokeN(L"Add", av, 5, &vControl))){									
									if (!(hr = CComDispatchDriverEx(vControl.pdispVal).PutPropertyByName(L"Caption", &CComVariant(s_szSiriusMenuTitle)))){
										if (!(hr = CComDispatchDriverEx(vControl.pdispVal).GetPropertyByName(L"Controls", &vControls))){
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, FALSE, 25, L"&Explore...", L"MLDisplayExplorer", true, NULL);
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, FALSE, 1786, L"&Clear Handles", L"MLClearHandles", true, NULL);
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, TRUE, 473, L"Insert &Object...", L"MLInsertObject", true, NULL);
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, FALSE, 0, L"Insert &Product...", L"MLInsertProduct", true, NULL);
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, TRUE, 23, L"&Load VBA Module...", L"MLLoadVBAModule", true, NULL);
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, FALSE, 3, L"&Save VBA Module...", L"MLSaveVBAModule", true, NULL);
											
											// Retrieve market data options
											{
												CAdapt<CComPtr<IDispatch> >	menu;
												InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlPopup, FALSE, 0, L"&Retrieve Market Data", L"", true, &menu);
												InitialiseExcelDrawMenuOption(menu.m_T.p, msoControlButton, FALSE, 0, L"&From Handle", L"MLRetrieveFromHandle", false, NULL);
												InitialiseExcelDrawMenuOption(menu.m_T.p, msoControlButton, TRUE, 0, L"&Asset", L"MLRetrieveAsset", false, NULL);
												InitialiseExcelDrawMenuOption(menu.m_T.p, msoControlButton, FALSE, 0, L"&Dividend schedule", L"MLRetrieveDividendSchedule", false, NULL);
												InitialiseExcelDrawMenuOption(menu.m_T.p, msoControlButton, FALSE, 0, L"&Volatility structure", L"MLRetrieveVolatilityStructure", false, NULL);
											}
																						
											InitialiseExcelDrawMenuOption(vControls.pdispVal, msoControlButton, TRUE, 0, L"&Setup...", L"MLSetup", true, NULL);
										}
									}								
								}
							}
						}
					}
				}
			}
		}
	}

	if (hr){
		ATLASSERT(false);
		CParameterMap::DisplayError(CStringResource(IDS_UNHANDLED_EXCEPTION).str(), MB_ICONSTOP);
		return E_FAIL;
	}
	return S_OK;
}
/*static*/ HRESULT CExcelInterface::InitialiseExcelDrawMenuOption(/*const*/ CComDispatchDriverEx/*&*/ ddControls, long nType, BOOL BeginGroup, long FaceID, LPCOLESTR Caption, LPCOLESTR OnAction, bool bEnabled, CAdapt<CComPtr<IDispatch> >* pout)
{
	CComVariant							vControl;
	HRESULT								hr;
	CComDispatchDriverEx				ddControl;

	if (hr = ddControls.Invoke1(L"Add", &CComVariant(nType), &vControl)) return hr;
	ddControl = vControl.pdispVal;
	if (hr = ddControl.PutPropertyByName(L"BeginGroup", &CComVariant(BeginGroup))) return hr;
	if (nType == msoControlButton && (hr = ddControl.PutPropertyByName(L"FaceID", &CComVariant(FaceID)))) return hr;
	if (hr = ddControl.PutPropertyByName(L"Caption", &CComVariant(Caption))) return hr;		
	if (!bEnabled && (hr = ddControl.PutPropertyByName(L"Enabled", &CComVariant(0L)))) return hr;
	if (nType == msoControlButton && (hr = ddControl.PutPropertyByName(L"OnAction", &CComVariant(OnAction)))) return hr;
	if (pout){
		CComVariant						vControls;
		if (hr = ddControl.GetPropertyByName(L"Controls", &vControls)) return hr;		
		pout->m_T = vControls.pdispVal;
	}
	return S_OK;
}
//
//	this is called exclusively by Insert[...] functions
//
/*static*/ HRESULT CExcelInterface::InsertObject(const CLSID& clsid, CComDispatchDriverEx& ddRange, long nRangeRows, long nRangeCols, long* pnRowOffset, const long nColOffset, long nParentRowOffset/* = -1*/, long nParentColOffset/* = -1 */)
//	clsid - CLSID of the object to create
//	ddRange - contains the Excel range to which we can write
//	nRangeRows, nRangeCols - number of rows and columns from the top left of ddRange to which we can write
//	pnRowOffset (in, out) - row of ddRange where we write the object's data
//	nColOffset - column of ddRange where we write the object's data
//	nParentRowOffset, nParentColOffset - if the object has a parent then this is the region in ddRange that we can set handle data in the parent
{
	HRESULT								hr;	
	CComPtr<IDispatch>					spObject;						// interface to the object associated with ObjectName	
	std::string							szFormula;						// excel range formula
	CComVariant							vValue;							// default "Value" for the object
	bool								bProductLike;					// true if the format is like a product type
	long*								pnSize = NULL;
	short int							nDimensions;	
	std::string							szProperName;
	
	if (hr = g_pApplication->GetObjectManager().CreateObject(clsid, spObject)) return hr;
	if (hr = g_pApplication->GetObjectManager().CLSIDToProperName(clsid, &szProperName)) return hr;

	// If we have information about a parent object, then try the name property first.
	if (nParentRowOffset >= 0 && nParentColOffset >= 0){
		if (!CComDispatchDriverEx(spObject).GetPropertyByName(L"Name", &vValue)){
			Write(ddRange, range_offset(nParentRowOffset, nParentColOffset, nRangeRows, nRangeCols), vValue);
			return S_OK;
		}
	}

	// Use the Value property in all other cases. Note that if we have information about a parent object
	// then we attach a handle to the space written for the object.
	if (hr = g_pApplication->GetObjectManager().CLSIDToProductLike(clsid, &bProductLike)) return hr;
	if (hr = CComDispatchDriverEx(spObject).GetPropertyByName(L"Value", &vValue)) return hr;
	// examine the size and create an appropriate area on the spreadsheet	
	nDimensions = CParameterMap::GetArrayDimensions(vValue, &pnSize);
	if (nDimensions == 2 || nDimensions == 0){
		// normal single argument case
		long nRows = pnSize ? pnSize[0] + 1 : 8;
		long nCols = pnSize ? pnSize[1] : 2;
		// write this object's default data
		Write(ddRange, range_offset(*pnRowOffset + 1, nColOffset, 1, 1, nRangeRows, nRangeCols), szProperName.c_str());
		if (nCols > 1) Write(ddRange, range_offset(*pnRowOffset + 1, nColOffset + 1, 1, 1, nRangeRows, nRangeCols), estring(20, ' '));	// for a reasonable column width
		Write(ddRange, range_offset(*pnRowOffset + 2, nColOffset, nRows - 1, nCols, nRangeRows, nRangeCols), vValue);

		// If possible, abstract any enumerator-type properties and present the values as cell validation.
		if (bProductLike){
			// Display any enum lists
			CParameterMap		pmValue;
			CComPtr<ITypeInfo>	spti;
			TYPEATTR*			pTypeAttr;
			FUNCDESC*			pFuncDesc;
			long				nRow;			
			
			pmValue.SetValue(vValue);
			if (!CComDispatchDriverEx(spObject)->GetTypeInfo(0, LOCALE_USER_DEFAULT, &spti)){
				if (!(spti->GetTypeAttr(&pTypeAttr))){
					for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){						
						if (!(spti->GetFuncDesc(nIndex, &pFuncDesc))){
							MEMBERID memid = pFuncDesc->memid;
							if (pFuncDesc->invkind == DISPATCH_PROPERTYGET){
								// Get the method name.
								CComBSTR s;
								if (!(spti->GetDocumentation(memid, &s, NULL, NULL, NULL))){
									// If this method is in the parameter list then we can write the validation.
									if (!pmValue.FindCell(estring(s), true, CParameterMap::LeftColumnOnly, &nRow, NULL)){
										// method found - now check enumerator
										VARTYPE vt = pFuncDesc->elemdescFunc.tdesc.vt;
										if (vt == VT_USERDEFINED){
											// user defined type
											HREFTYPE				refType = pFuncDesc->elemdescFunc.tdesc.hreftype;
											CComPtr<ITypeInfo>		spTypeInfoRef;
											HRESULT					hrLocal = spti->GetRefTypeInfo(refType, &spTypeInfoRef);
											if (!hrLocal){
												std::string			szList;					// comma-separated list of possible enumerator values
												if (!(CEnumMap::GetEnumList(spTypeInfoRef, LIBID_Sirius, &szList))){
													WriteValidation(ddRange, range_offset(*pnRowOffset + 2 + nRow, nColOffset + 1, 1, 1, nRangeRows, nRangeCols), s, szList, "");
												}
											}
										}
									}
								}						
							}
							spti->ReleaseFuncDesc(pFuncDesc);
						}
					}
					spti->ReleaseTypeAttr(pTypeAttr);
				}
			}
		}
		
		// set a formula at the top of the object's data (we do this after the default data has been written to support auto-calculation)
		range_offset(*pnRowOffset + 2, nColOffset, nRows - 1, nCols, nRangeRows, nRangeCols).GetCreateFormula(ddRange, clsid, &szFormula);		
		Write(ddRange, range_offset(*pnRowOffset , nColOffset, nRangeRows, nRangeCols), szFormula.c_str());			
		
		// format the object
		CComDispatchDriverEx ddRangeFormat;
		range_offset(*pnRowOffset, nColOffset, nRows, nCols, nRangeRows, nRangeCols).MapRange(ddRange, ddRangeFormat);				
		if (bProductLike){			
			if (ddRangeFormat) FormatLikeProduct(ddRangeFormat, nRows, nCols, nRangeRows, nRangeCols);
		} else {
			if (ddRangeFormat) FormatLikeObject(ddRangeFormat, nRows, nCols, nRangeRows, nRangeCols);
		}
								
		// if we have a parent object then add the handle creation formula
		if (nParentRowOffset >= 0 && nParentColOffset >= 0){			
			range_offset(*pnRowOffset, nColOffset, nRangeRows, nRangeCols).GetFormula(ddRange, &szFormula);
			Write(ddRange, range_offset(nParentRowOffset, nParentColOffset, nRangeRows, nRangeCols), szFormula.c_str());
		}
		(*pnRowOffset) += nRows;
	} else if (nDimensions == 1){
		// multiple argument case
		SAFEARRAY*					psa;
		long						nl, nh;
		long						nIndex;
		CComVariant					vElement;
		long						nRowOffset_start = *pnRowOffset;
		std::vector<range_offset>	vtrRangeOffset;
		range_offset				roUnion(nRowOffset_start, nColOffset, nRangeRows, nRangeCols);
					
		if (vValue.vt & VT_BYREF){
			psa = *(vValue.pparray);
		} else {
			psa = vValue.parray;
		}
		if (!hr) hr = ::SafeArrayGetLBound(psa, 1, &nl);
		if (!hr) hr = ::SafeArrayGetUBound(psa, 1, &nh);
								
		Write(ddRange, range_offset(*pnRowOffset + 1, nColOffset, 1, 1, nRangeRows, nRangeCols), szProperName.c_str());
		for (nIndex = nl; nIndex <= nh && !hr; nIndex++){
			if (hr = CParameterMap::SafeArrayElementToVariant(psa, &nIndex, &vElement)) break;
			short int				nDimensionsElement;
			long*					pnSizeElement = NULL;
			nDimensionsElement = CParameterMap::GetArrayDimensions(vElement, &pnSizeElement);
			if (nDimensionsElement == 2 || nDimensionsElement == 0){
				long nRows = pnSizeElement ? pnSizeElement[0] : 1;
				long nCols = pnSizeElement ? pnSizeElement[1] : 1;
				// write the data
				range_offset ro(*pnRowOffset + 2, nColOffset, nRows, nCols, nRangeRows, nRangeCols);
				vtrRangeOffset.push_back(ro);
				if (hr = roUnion.Union(ro)) break;					
				Write(ddRange, ro, vElement);					
				SetFont(ddRange, range_offset(*pnRowOffset + 2, nColOffset, 1, nCols, nRangeRows, nRangeCols), L"Tahoma", 8, 0, false);
				SetBackColour(ddRange, ro, 2);
				(*pnRowOffset) += nRows;
			} else {
				hr = E_FAIL;
			}
			delete pnSizeElement;				
		}
		// set a formula at the top of the object's data (we do this after the default data has been written to support auto-calculation)
		range_offset::GetCreateFormula(ddRange, vtrRangeOffset, clsid, &szFormula);
		Write(ddRange, range_offset(nRowOffset_start, nColOffset, nRangeRows, nRangeCols), szFormula.c_str());
		
		// autofit the area
		Autofit(ddRange, roUnion);
		{
			CComDispatchDriverEx ddRangeWrite;
			roUnion.MapRange(ddRange, ddRangeWrite);
			FormatLikeObject(ddRangeWrite, roUnion.GetRows() - 1, roUnion.GetColumns(), nRangeRows, nRangeCols);
		}
		if (!hr && nParentRowOffset >= 0 && nParentColOffset >= 0){
			// if we have a parent object then add the handle creation formula				
			range_offset(nRowOffset_start, nColOffset, nRangeRows, nRangeCols).GetFormula(ddRange, &szFormula);
			Write(ddRange, range_offset(nParentRowOffset, nParentColOffset, nRangeRows, nRangeCols), szFormula.c_str());
		}
	} else {
		hr = E_FAIL;			
	}
	delete pnSize;
	return hr;
}
//
//	Callable only from Excel. Create a new object of type ObjectName
//	and write it to the Excel spreadsheet at location At
//
/*static*/ HRESULT CExcelInterface::InsertObject(const std::string& szObject, const std::string& szAt)
{
	HRESULT								hr;
	CLSID								clsid;							// CLSID associated with szObject
	CComDispatchDriverEx				ddRange;						// range object associated with the top-left corner of At
	long								nRangeRows;						// number of rows that we can write in the spreadsheet
	long								nRangeCols;						// number of columns that we can write in the spreadsheet
	long								nRowOffset = 0;					// number of rows written to the object
	CParameterMap						pm;
	
	if (hr = g_pApplication->GetObjectManager().NameToCLSID(szObject, &clsid)) return hr;
	if (hr = GetWritableRange(szAt, ddRange, &nRangeRows, &nRangeCols)) return hr;
	return InsertObject(clsid, ddRange, nRangeRows, nRangeCols, &nRowOffset, 0, -1, -1);
}
//
//	Callable only from Excel. Create a new product object of type ProductType
//	and write it to the Excel spreadsheet at location At
//
/*static*/ HRESULT CExcelInterface::InsertProduct(const std::string& szProductType, const std::string& szAt)
{	
	HRESULT								hr;
	CComDispatchDriverEx				ddRange;						// range object associated with the top-left corner of At
	long								nRangeRows;						// number of rows that we can write in the spreadsheet
	long								nRangeCols;						// number of columns that we can write in the spreadsheet
	CComPtr<IDispatch>					spProductType;					// interface to the product type object associated with ProductType
	long								nProductRows = 0;				// number of rows written in the product object
	long								nObjectRows = 0;				// number of rows written in the column of the worksheet associated with objects referenced by the principal product object		
	bool								bDisplayProductOnly = _Module.GetDisplayProductOnly();	// true if we only display the product object, false if we display any child objects
	
	if (hr = GetWritableRange(szAt, ddRange, &nRangeRows, &nRangeCols)) return hr;	
	if (hr = g_pApplication->GetObjectManager().CreateObject(szProductType, spProductType)) return hr;
	if (hr = Write(ddRange, range_offset(nProductRows + 1, 0, nRangeRows, nRangeCols), g_pApplication->GetObjectManager().NameToProperName(szProductType).c_str())) return hr;
	if (hr = Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), "")) return hr;
	
	// Iterate through all the set properties.
	CComPtr<ITypeInfo>					pti;	
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc;
	MEMBERID							memid;		
	std::string							szFormula;						// formula to insert into an Excel spreadsheet cell
		
	nProductRows++;		
	if (!(hr = CComDispatchDriverEx(spProductType)->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti))){
		if (!(pti->GetTypeAttr(&pTypeAttr))){
			for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs; nIndex++){
				if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){
					memid = pFuncDesc->memid;
					if (pFuncDesc->invkind == DISPATCH_PROPERTYGET){
						// get the method name
						CComBSTR bstr;
						CComBSTR bstrComment;
						if (!(pti->GetDocumentation(memid, &bstr, &bstrComment, NULL, NULL))){
							if (!Write(ddRange, range_offset(nProductRows + 1, 0, nRangeRows, nRangeCols), CComVariant(bstr))){
								WriteComment(ddRange, range_offset(nProductRows + 1, 0, nRangeRows, nRangeCols), CComVariant(bstrComment));
								// now pull in the starting property values
								VARTYPE vt = pFuncDesc->elemdescFunc.tdesc.vt;
								if (vt == VT_PTR){
									if (!bDisplayProductOnly){
										// Object is expected. Get its GUID.
										ELEMDESC* ped = (ELEMDESC*)pFuncDesc->elemdescFunc.tdesc.lptdesc;
										HRESULT hrLocal = (ped->tdesc.vt == VT_USERDEFINED ? S_OK : E_FAIL);
										if (!hrLocal){
											CComPtr<ITypeInfo>		pti_ref;										
											if (!(hrLocal = pti->GetRefTypeInfo(ped->tdesc.hreftype, &pti_ref))){
												TYPEATTR*	pTypeAttrRef;
												if (!(hrLocal = pti_ref->GetTypeAttr(&pTypeAttrRef))){
													// Create an instance of this object
													CLSID clsid;
													if (!(hrLocal = g_pApplication->GetObjectManager().IIDToCLSID(pTypeAttrRef->guid, &clsid))){
														hrLocal = InsertObject(clsid, ddRange, nRangeRows, nRangeCols, &nObjectRows, 3, nProductRows + 1, 1);
														nObjectRows += 2;	// To allow for the object handle and one blank row.
													}
													pti_ref->ReleaseTypeAttr(pTypeAttrRef);
												}										
											}
										}
										if (hrLocal){
											Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), "#unknown object");
										}
									}
								} else if (vt & VT_ARRAY){
									// ToDo - array variant
									Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), "#unknown array");
								} else if (vt == VT_USERDEFINED){
									// user defined type
									HREFTYPE				refType = pFuncDesc->elemdescFunc.tdesc.hreftype;
									CComPtr<ITypeInfo>		spTypeInfoRef;
									HRESULT					hrLocal = pti->GetRefTypeInfo(refType, &spTypeInfoRef);
									CComVariant				vPropertyValue;

									CComDispatchDriverEx(spProductType).GetProperty(memid, &vPropertyValue);
									if (!hrLocal){										
										std::string			szInitialValue;
										std::string			szList;					// comma-separated list of possible enumerator values
										ATLASSERT(vPropertyValue.vt == VT_I4);
										if (!(hrLocal = CEnumMap::GetString(spTypeInfoRef, LIBID_Sirius, vPropertyValue.lVal, &szInitialValue))){
											if (!(hrLocal = CEnumMap::GetEnumList(spTypeInfoRef, LIBID_Sirius, &szList))){
												WriteValidation(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), bstr, szList, szInitialValue);
											}
										}
									}
									if (hrLocal){
										Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), "#unknown user defined type");
									}
								} else if (vt == VT_I1){
									// we take this to be the Boolean type									
									CComVariant		vPropertyValue;
									if (!CComDispatchDriverEx(spProductType).GetProperty(memid, &vPropertyValue)){
										Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), vPropertyValue.boolVal ? L"True" : L"False");
									}
								} else {
									// normal variant
									CComVariant		vPropertyValue;									
									if (!CComDispatchDriverEx(spProductType).GetProperty(memid, &vPropertyValue)){
										Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), vPropertyValue);
									}
								}
								nProductRows++;
							}							
						}
					}
					pti->ReleaseFuncDesc(pFuncDesc);
				}
				if (hr) break;
			}
			pti->ReleaseTypeAttr(pTypeAttr);
		}		
	}
	
	// Set a formula at the top of the object's data (we do this after the default data has been written to support auto-calculation)	
	std::vector<range_offset>	vtrRangeOffset;
	vtrRangeOffset.push_back(range_offset(1, 0, 1, 1, nRangeRows, nRangeCols));					// product type	
	vtrRangeOffset.push_back(range_offset(false));												// optional parameter
	vtrRangeOffset.push_back(range_offset(false));												// optional parameter	
	vtrRangeOffset.push_back(range_offset(2, 0, nProductRows - 1, 1, nRangeRows, nRangeCols));	// parameters
	vtrRangeOffset.push_back(range_offset(2, 1, nProductRows - 1, 1, nRangeRows, nRangeCols));	// values	
	range_offset::GetCreateFormula(ddRange, vtrRangeOffset, CLSID_Product, &szFormula);
	Write(ddRange, range_offset(0, 0, nRangeRows, nRangeCols), szFormula.c_str());	
		
	// Set the MLEvaluate formula
	range_offset(0, 0, nRangeRows, nRangeCols).GetShortAddress(ddRange, &szFormula);	
	Write(ddRange, range_offset(nProductRows + 1, 0, nRangeRows, nRangeCols), "Result:");
	Write(ddRange, range_offset(nProductRows + 1, 1, nRangeRows, nRangeCols), std::string("=MLEvaluate(" + szFormula + ")").c_str());
	
	// Format the area		
	FormatLikeProduct(ddRange, nProductRows, 2, nRangeRows, nRangeCols);	
	return hr;
}


////////////////////////////////////////////////////////////////////////////
//	Resize
//
//	Emulates the GdaResize() function
//
/*static*/ HRESULT CExcelInterface::Resize(const VARIANT& v, VARIANT* pv)
{
	if (_Resize(v, pv)){
		// failure
		CComVariant vOut(v);
		return vOut.Detach(pv);
	}
	return S_OK;
}
HRESULT CExcelInterface::_Resize(const VARIANT& v, VARIANT* pv)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddExcel;
	CComDispatchDriverEx				ddCaller;
	CComDispatchDriverEx				ddRows, ddColumns;
	CComVariant							vRows, vColumns;
	CParameterMap						pm;
	
	if (hr = _Module.GetExcel(ddExcel)) return hr;
	if (hr = ddExcel.GetPropertyByName(L"Caller", ddCaller)) return hr;
	if (hr = ddCaller.GetPropertyByName(L"Rows", ddRows)) return hr;
	if (hr = ddCaller.GetPropertyByName(L"Columns", ddColumns)) return hr;
	if (hr = ddRows.GetPropertyByName(L"Count", &vRows)) return hr;
	if (hr = ddColumns.GetPropertyByName(L"Count", &vColumns)) return hr;
	if (hr = pm.SetValue(v)) return hr;
	if (hr = vRows.ChangeType(VT_I4)) return hr;
	if (hr = vColumns.ChangeType(VT_I4)) return hr;		
	if (hr = pm.Resize(vRows.lVal, vColumns.lVal, "")) return hr;
	if (hr = pm.SetColumns(vRows.lVal)) return hr;
	return pm.GetValue(pv);
}
//
//	sets the horizontal and vertical alignment properties of an Excel range constrained by a given range offset
//
/*static*/ HRESULT CExcelInterface::SetAlignment(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nHorizontalAlignment, long nVerticalAlignment)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;						// range to write to
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;	
	if (hr = ddRangeUse.PutPropertyByName(L"HorizontalAlignment", &CComVariant((long)nHorizontalAlignment))) return hr;
	if (hr = ddRangeUse.PutPropertyByName(L"VerticalAlignment", &CComVariant((long)nVerticalAlignment))) return hr;
	return S_OK;
}
//
//	set the background colour of an Excel range constrained by a given range offset
//
/*static*/ HRESULT CExcelInterface::SetBackColour(const CComDispatchDriverEx& ddRange, const range_offset& ro, long nColour)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;						// range to write to	
	CComVariant							vInterior;
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;		
	if (hr = ddRangeUse.GetPropertyByName(L"Interior", &vInterior)) return hr;
	CComDispatchDriverEx ddInterior(vInterior.pdispVal);	
	return ddInterior.PutPropertyByName(L"ColorIndex", &CComVariant((long)nColour, VT_I4));
}
//
//	set the background colour of an Excel range constrained by a given range offset
//	
/*static*/ HRESULT CExcelInterface::SetFont(const CComDispatchDriverEx& ddRange, const range_offset& ro, CComVariant vFontName, long nSize, long nColour, bool bBold)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;						// range to adjust	
	CComVariant							vFont;
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;		
	if (hr = ddRangeUse.GetPropertyByName(L"Font", &vFont)) return hr;
	CComDispatchDriverEx ddFont(vFont.pdispVal);	
	if (hr = ddFont.PutPropertyByName(L"Name", &vFontName)) return hr;
	if (hr = ddFont.PutPropertyByName(L"ColorIndex", &CComVariant((long)nColour, VT_I4))) return hr;
	if (hr = ddFont.PutPropertyByName(L"Bold", &CComVariant(bBold))) return hr;
	if (hr = ddFont.PutPropertyByName(L"Size", &CComVariant(nSize))) return hr;	
	return S_OK;
}
//
// set the caption according to the status of this Sirius build
//
/*static*/ HRESULT CExcelInterface::UpdateSiriusStatusCaption()
{
	HRESULT								hr;
	CComDispatchDriverEx				ddExcel;	
	CComVariant							vCaption;
	
	if (hr = _Module.GetExcel(ddExcel)) return hr;
	estring szCaption;
	if (_Module.GetDisplaySiriusStatus()){
		_Module.GetSiriusCaption(&vCaption);
	} else {		
		ddExcel.GetPropertyByName(L"Name", &vCaption);
	}
	return ddExcel.PutPropertyByName(L"Caption", &vCaption);
}
//
//	returns the key for a locked VBA project
//
/*static*/ HRESULT CExcelInterface::VbaKey(const std::string& szFileName, std::string* pszKey)
{
	char								s[128];
	char								*encrypt_key_str;

	if (!(encrypt_key_str = VbaKeyFindEncryptKey(szFileName.c_str()))){
		pszKey->assign(CStringResource(IDS_KEY_NOT_FOUND));
		return S_FALSE;
	}	
	VbaKeyDecryptKey(encrypt_key_str, s);   
	pszKey->assign(s);   
	free(encrypt_key_str);	
	return 0;
}
/*static*/ void CExcelInterface::VbaKeyDecryptKey(const char *encrypt_pass, char *s)
{
	char								str[128], ch;
	char								hs[] = { 0, 0, 0 };
	int									v1, v2, i, l;
	int									begin_found = 0;

	*s = 0;
	for (i = 0; encrypt_pass[i*2]; i++){
		hs[0] = encrypt_pass[i*2]; hs[1] = encrypt_pass[i*2+1];
		v1 = strtol(hs, 0, 16);
		hs[0] = encrypt_pass[i*2+4]; hs[1] = encrypt_pass[i*2+5];
		v2 = strtol(hs, 0, 16);
		if (!begin_found){
			if(v1 == v2) begin_found = 1;
		} else {
			if (v1 != v2){
				begin_found = 0;
			} else {
				i += 3;
				break;
			}
		}
	}
	if (!begin_found) return;

	for (ch = 0, l = 0; encrypt_pass[i*2+2]; i++, l++){
		hs[0] = encrypt_pass[(i-2)*2]; hs[1] = encrypt_pass[(i-2)*2+1];
		v1 = strtol(hs, 0, 16);
		hs[0] = encrypt_pass[i*2]; hs[1] = encrypt_pass[i*2+1];
		v2 = strtol(hs, 0, 16);
		ch = (ch + (char)v1) ^ (char)v2;
		str[l] = ch;
	}
	str[l] = 0;
	CharToOem(str, s);
}
/*static*/ LPSTR CExcelInterface::VbaKeyFindEncryptKey(const char *filename)
{
	FILE								*prot_file;
	char								*buff = 0, *s;
	size_t								rc;
	int									off;

	prot_file = fopen(filename, "rb");
	if (!prot_file) return 0;
	buff = (char*)malloc(BLOCK_SIZE);
	if(!buff) return 0;
	do {
		rc = fread(buff, 1, BLOCK_SIZE, prot_file);
		if (!rc) break;
		off = VbaKeyMemStr(buff, "DPB=\"", BLOCK_SIZE);
		if (off < 0){
			if (rc < BLOCK_SIZE) break;
			fseek(prot_file, -32, SEEK_CUR);
			continue;
		}
		fseek(prot_file, off - rc, SEEK_CUR);
		rc = fread(buff, 1, BLOCK_SIZE, prot_file);
		if (!rc) break;
		s = strchr(buff, '\"');
		*s = 0;
		return buff;
	} while (!feof(prot_file));
	free(buff);
	return 0;
}
/*static*/ int CExcelInterface::VbaKeyMemStr(const char *buf, const char *s, size_t size)
{
	int									off;
	const char*							s_temp;

	if (!buf || !s || !(*s)) return -1;
	for (off = 0; off < size; off++){
		for (s_temp = s; *s_temp != 0 && *buf == *s_temp; buf++, s_temp++){
			off++;
			if (off >= size) return -1;
		}
		if ((*(buf - 1) == *(s_temp - 1)) && *s_temp == 0) return off;
		buf++;
	}
	return -1;
}
//
//	write a value to a given Excel range constrained by a given range offset
//
/*static*/ HRESULT CExcelInterface::Write(const CComDispatchDriverEx& ddRange, const range_offset& ro, const CComVariant& vIn)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;						// range to write
	
	if (hr = ro.MapRange(ddRange, ddRangeUse)) return hr;			
	if (ddRangeUse.Invoke0(L"Clear")) ATLASSERT(false);
	return hr = ddRangeUse.PutPropertyByName(L"Value", (VARIANT*)&vIn);	
}
//
//	write a comment to the top left cell of a given Excel range constrained by a given range offset
//
/*static*/ HRESULT CExcelInterface::WriteComment(const CComDispatchDriverEx& ddRange, const range_offset& ro, const CComVariant& vIn)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse;

	if (hr = ro.TopLeft().MapRange(ddRange, ddRangeUse)) return hr;
	if (ddRangeUse.Invoke0(L"ClearComments")) ATLASSERT(false);
	return ddRangeUse.Invoke1(L"AddComment", (VARIANT*)&vIn);
}
//
//	write a value in a cell along with a combo control of possible values with weak validation
//
/*static*/ HRESULT CExcelInterface::WriteValidation(const CComDispatchDriverEx& ddRange, const range_offset& ro, const CComBSTR& sParameterName, const std::string& szList, const std::string& szInitialValue)
{
	HRESULT								hr;
	CComDispatchDriverEx				ddRangeUse, ddValidation;
	CComVariant							vParams[5] = {CComVariant(), CComBSTR(szList.c_str()), xlBetween, xlValidAlertInformation, xlValidateList};
	CComBSTR							sMessage = CStringResource(IDS_INVALID_VALUE_FOR_PARAMETER).bstr();

	if (hr = ro.TopLeft().MapRange(ddRange, ddRangeUse)) return hr;
	if (hr = ddRangeUse.GetPropertyByName(L"Validation", ddValidation)) return hr;
	if (ddValidation.Invoke0(L"Delete")) ATLASSERT(false);

	if (hr = ddValidation.InvokeN(L"Add", vParams, 5)) return hr;
	if (hr = ddValidation.PutPropertyByName(L"IgnoreBlank", &CComVariant(FALSE))) return hr;
	if (hr = ddValidation.PutPropertyByName(L"InCellDropDown", &CComVariant(TRUE))) return hr;
	if (hr = ddValidation.PutPropertyByName(L"ErrorTitle", &CComVariant(L"Sirius"))) return hr;
	sMessage += L" '";
	sMessage += sParameterName;
	sMessage += L"'.";		
	if (hr = ddValidation.PutPropertyByName(L"ErrorMessage", &CComVariant(sMessage.Copy()))) return hr;	
	if (szInitialValue.size() && (hr = ddRangeUse.PutPropertyByName(L"Value", &CComVariant(szInitialValue.c_str())))) return hr;
	return S_OK;
}

#undef BLOCK_SIZE