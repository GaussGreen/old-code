//	xmlstreamer.cpp: implementation of the CXmlStreamer class.
//
//  Author :		 David Cuin.
//
//	I would be inclined to agree with anyone who proposes the
//  hypothesis that there exists a better xml schemer. However,
//  they would have to offer proof by supplying to me a replacement
//  xmlstreamer.cpp. Otherwise, I assume that my serialisation is
//  the best!
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "xmlstreamer.h"
#include "siriusapplication.h"
#include "comobjectcollectionserialisablekey.h"

#undef min
#undef max


////////////////////////////////////////////////////////////////////////////
//	xml_attributes implementation
//
//	attempt to get a value from the attribute list
//
HRESULT CXmlStreamer::xml_attributes::GetValue(const char* szName, std::string* pszValue)
{
	for (std::vector<xml_attribute>::iterator it = begin(); it < end(); it++){
		if ((strlen(szName) == it->sz_name_end - it->sz_name_start + 1) && !strncmp(szName, it->sz_name_start, it->sz_name_end - it->sz_name_start + 1)){
			// found it
			pszValue->assign(it->sz_value_start, it->sz_value_end - it->sz_value_start + 1);
			return S_OK;
		}		
	}
	return E_FAIL;
}
HRESULT CXmlStreamer::xml_attributes::GetValue(const char* szName, long* pnValue)
{
	HRESULT								hr;
	std::string							szValue;
	
	if (hr = GetValue(szName, &szValue)) return hr;
	return estring::islong(szValue, pnValue) ? S_OK : S_FALSE;
}
HRESULT CXmlStreamer::xml_attributes::GetValue(const char* szName, unsigned short* pnValue)
{
	HRESULT								hr;
	long								nValue;
	if (hr = GetValue(szName, &nValue)) return hr;
	if (nValue < std::numeric_limits<unsigned short>::min() || nValue > std::numeric_limits<unsigned short>::max()) return E_FAIL;
	*pnValue = (unsigned short)nValue;
	return S_OK;
}


////////////////////////////////////////////////////////////////////////////
//	CXmlStreamer implementation
//
CXmlStreamer::CXmlStreamer()
{
	ATLASSERT(false);
}
CXmlStreamer::~CXmlStreamer()
{
	ATLASSERT(false);
}
//
//	Maps an input to an XML stream.
//
/*static*/ HRESULT CXmlStreamer::GetXML(const CComVariant& vIn, xmlstreamer& ssOut)
{
	return GetXML(vIn, ssOut, false, 0, NULL);
}
// All GetXML functions must call into this one. Furthermore, all GetXML calls from
// within this or other GetXML functions must be to this function. Otherwise you
// run the risk of messing up the recursion counter.
/*protected, static*/ HRESULT CXmlStreamer::GetXML(const CComVariant& vIn, xmlstreamer& ssOut, bool bNullifyDataSource, long nRecursionCount, CComPtr<IDispatch> spParent)
//	vIn - input variant to stream to XML
//	ssOut - XML stream output
//  bNullifyDataSource - true if we set all object keys to have Name.Date.NoDataSource
//	nRecursionCount - number of times this function has called itself
//	spParent - If this function was called by itself and that caller (or its caller, and so on) was processing an object, then spParent holds the pointer to that object.
{
	HRESULT								hr = S_OK;
	short int							nDims = CParameterMap::GetArrayDimensions(vIn);	
	std::string							szPrefix(nRecursionCount, '\t');

	if (!nRecursionCount){
		// Insert the XML encoding string.
		ssOut << CFastStream::s_szEncoding << "\r\n";
	}

	if (!nDims){
		// scalar
		if (vIn.vt == VT_DISPATCH){
			// Object pointer passed in. We need to decide whether we (i) decompose this object into its value
			// and stream that to XML (effectively a deep copy), or (ii) simply to record the key associated
			// with a serialisable object (effectively a shallow copy).			
			// We tend never to take a deep copy (thereby preventing data duplication) if the object is serialisable.
			// There are exceptions: (i)  We are actually trying to get the XML for a serialisable object
			//							  (indicated by the zeroth level of recursion).
			//						 (ii) If the parent object is a sirius collection.
			
			CComPtr<IDispatch>			spObject(vIn.pdispVal);
			
			if (spObject && nRecursionCount && !g_pApplication->GetObjectManager().IsSiriusCollection(spParent)){
				if (g_pApplication->GetObjectManager().IsSerialisable(spObject)){
					return GetXML(estring::GetValue(CComObjectCollectionSerialisableKey(vIn.pdispVal).GetKeyAndObjectName(bNullifyDataSource, true)), ssOut, bNullifyDataSource, nRecursionCount, spParent);
				}
			}
			
			// Write the name of the object and its index in any collection associated with spParent
			std::string					szObjectName;
			estring						szCollectionIndex;
			if (spObject){
				if (!(hr = CParameterMap::GetObjectName(spObject, &szObjectName)) && spParent){
					std::string			szCollectionName;
					g_pApplication->GetObjectManager().NameToCollectionName(szObjectName, &szCollectionName);
					// Check to see if the collection object name is a property of spParent.
					// ToDo - a more rigourous (and the only correct) way to do this is to check whether any property of spParent
					// is of type szCollectionName and use the property name instead of szCollectionName.
					g_pApplication->GetObjectManager().GetCollectionIndex(spParent, spObject, szCollectionName, &szCollectionIndex);
					CFastStream::ToXml(&szCollectionIndex, true);
				}
			} else {
				szObjectName = "Null";
			}
			if (szCollectionIndex.size()){
				ssOut << szPrefix << "<Object name=\"" << szObjectName << "\" collectionindex=\"" << szCollectionIndex << "\">\r\n";
			} else {			
				ssOut << szPrefix << "<Object name=\"" << szObjectName << "\">\r\n";
			}						
			if (hr){
				ssOut << szPrefix << "\t" << CStringResource(IDS_INVALID_OBJECT_NAME).str() << "\r\n";
			} else {
				// we use the 'Value' property of the object
				CComVariant				vValue;
				if (!spObject){
					// do nothing for the blank object
				} else if (!(hr = CComDispatchDriverEx(spObject).GetPropertyByName(L"Value", &vValue))){
					hr = GetXML(vValue, ssOut, bNullifyDataSource, nRecursionCount + 1, spObject);
				} else {
					ssOut << szPrefix << "\t" << CStringResource(IDS_COULD_NOT_GET_VALUE).str() << "\r\n";
				}
			}
			ssOut << szPrefix << "</Object>\r\n";
			return hr;
		} else {
			// scalar element
			return CFastStream::WriteScalar(szPrefix, vIn, ssOut);						
		}
	} else if (nDims == 1 || nDims == 2){
		// vector or matrix		
		CParameterMap					pm;
		long							nRow;		
		std::string						szProperty;
		std::string						szProductType;
		
		hr = pm.SetValue(vIn);		
		ssOut << szPrefix << "<Matrix rows=\"" << pm.GetRows() << "\" columns=\"" << pm.GetCols() << "\">\r\n";
		for (nRow = 0; nRow < pm.GetRows() && !hr; nRow++){
			for (long nCol = 0; nCol < pm.GetCols() && !hr; nCol++){
				CComVariant			vElement;
				pm.GetValue(nRow, nCol, &vElement);				
				// if vElement is a handle then we need to transform it to an object			
				CComPtr<IDispatch> spObject;
				if (!g_pApplication->GetObjectManager().GetObject(vElement, spObject)){
					vElement = CComVariant(spObject);
				} else {
					// Sometimes, vElement could be an array with one element in it. We decompose this into a scalar to make the XML look better.
					CParameterMap		pmElement;
					if (!pmElement.SetValue(vElement) && pmElement.IsScalar()){
						pmElement.GetValue(0, 0, &vElement);
					}
				}			
				ssOut << szPrefix << "\t" << "<Element row=\"" << nRow << "\" column=\"" << nCol << "\" type=\"" << (vElement.vt == VT_DATE ? VT_R8 : vElement.vt) << "\">\r\n";
				hr = GetXML(vElement, ssOut, bNullifyDataSource, nRecursionCount + 1, spParent);
				ssOut << szPrefix << "\t" << "</Element>\r\n";
			}
		}		
		ssOut << szPrefix << "</Matrix>\r\n";
		return hr;
	}
			
	ATLASSERT(false);	// should never reach this point
	ssOut << szPrefix << CStringResource(IDS_UNHANDLED_DATA_TYPE).str() << "\r\n";
	return E_FAIL;
}
//
//	Maps an XML stream to a variant.
//
/*static*/ HRESULT CXmlStreamer::GetVariant(const char* szIn, CComVariant& vOut)
{	
	return GetVariant(szIn, NULL, vOut);
}
/*protected, static*/ HRESULT CXmlStreamer::GetVariant(const char* szIn, DataSourceEnum* pds, CComVariant& vOut)
//	pds - If given then we modify any object keys to this data source.
{
	HRESULT								hr;
	const char*							szTagStart = NULL;				// corresponds to 'x' in "<x    >" this is NULL if we are not in a tag
	const char*							szTagEnd = NULL;				// corresponds to 'y' in "<    y>" or <"x...y attribute>"
	xml_attributes						attributes;						// set of attributes for the current tag
	xml_attribute						attribute;						// current attribute		
	std::vector<std::string>			vector_tag_stack;	
	std::vector<void*>					vector_object_stack;

	// check the encoding string
	if (!::strncmp(szIn, CFastStream::s_szEncoding.data(), CFastStream::s_szEncoding.size())){
		szIn += CFastStream::s_szEncoding.size();
	}

	vector_object_stack.push_back(reinterpret_cast<void*>(&vOut));		// this kicks the whole thing off
	while (*szIn){
		switch (*szIn){
		case '\0':
			// the end of the string has been reached
			break;
		case '<':
			// we are starting a tag
			if (szTagStart){
				ATLASSERT(false);	// ToDo - decent error
				return E_FAIL;
			}
			szTagStart = szIn + 1;
			break;
		case '>':
			// we are finishing a tag
			if (!szTagStart){
				ATLASSERT(false);	// ToDo - decent error
				return E_FAIL;
			}
			if (!szTagEnd) szTagEnd = szIn - 1;			
			if (attribute.sz_name_start){
				ATLASSERT(false);
				return E_FAIL;	// this is impossible since we process attributes in one block
			}
			if (attributes.size() && szTagStart[0] == '/'){
				// closing tags never have attributes
				ATLASSERT(false);
				return E_FAIL;
			}			
			if (hr = GetVariantProcessTag(szTagStart, szTagEnd, attributes, vector_tag_stack, vector_object_stack)){
				return hr;
			}
			attributes.clear();			
			szTagStart = szTagEnd = NULL;
			break;
		case '=':			
			if (!attribute.sz_name_start){
				ATLASSERT(false);
				return E_FAIL;
			}
			attribute.sz_name_end = szIn - 1;			
			// find the next double quote
			while (*szIn && *szIn != '\"'){
				szIn++;
			}
			if (!*szIn){
				ATLASSERT(false);
				return E_FAIL;
			}
			attribute.sz_value_start = ++szIn;
			// and the next one
			while (*szIn && *szIn != '\"'){
				szIn++;
			}
			if (!*szIn){
				ATLASSERT(false);
				return E_FAIL;
			}
			attribute.sz_value_end = szIn - 1;
			attributes.push_back(attribute);
			attribute.clear();
			break;
		case ' ': case '\t': case '\n' : case '\r' : case '\f':
			if (szTagStart){
				// We are in a tag/attribute field. We have certainly parsed the whole of the tag name now.
				// We could either be processing attributes or there could be whitespace between the tag name
				// and the final '>'.				
				if (!szTagEnd){
					szTagEnd = szIn - 1;
				}
				if (szTagStart && attribute.sz_name_start){					
					return CParameterMap::ReturnErrorR(IDS_XML_WHITESPACE);
				}
				// pass over the spaces
				while (*szIn && isspace(*szIn)){
					szIn++;
				}				
				if (!*szIn){
					ATLASSERT(false);
					return E_FAIL;
				}
				if (*szIn != '>'){
					attribute.sz_name_start = szIn;
				}				
				continue;	// this avoids a further szIn++
			}
			break;
		default:
			if (!szTagStart){
				const char* szTextStart;
				const char* szTextEnd;
				// Free text field. Pass over the initial whitespace
				while (*szIn && isspace(*szIn)){
					szIn++;
				}
				if (!*szIn){
					ATLASSERT(false);
					return E_FAIL;
				}
				szTextStart = szTextEnd = szIn;
				// Free text until the first '<'. Note that '>' is illegal.
				while (*szIn && *szIn != '<' && *szIn != '>'){
					if (!isspace(*szIn)) szTextEnd = szIn;
					szIn++;					
				}
				if (!*szIn){
					return CParameterMap::ReturnErrorR(IDS_XML_UNEXPECTED_END);
				} else if (*szIn == '>'){
					ATLASSERT(false);
					return E_FAIL;
				}
				if (hr = GetVariantProcessFreeText(szTextStart, szTextEnd, vector_tag_stack, vector_object_stack, pds)) return hr;
				continue;	// this avoids a further szIn++
			}						
		}
		szIn++;
	}
	
	if (szTagStart || szTagEnd || attribute.sz_name_start){
		ATLASSERT(false);
		return E_FAIL;
	}
	if (vector_tag_stack.size() != 0 || vector_object_stack.size() != 1){
		// this implies an unexpected end of the XML
		return CParameterMap::ReturnErrorR(IDS_XML_UNEXPECTED_END);		
	}
		
	return S_OK;
}
//
//	Called exclusively by GetVariant.
//
/*static*/ HRESULT CXmlStreamer::GetVariantProcessFreeText(const char* szTextStart, const char* szTextEnd, std::vector<std::string>& vector_tag_stack, std::vector<void*>& vector_object_stack, DataSourceEnum* pds)
//	pds - If given then we modify any object keys to this data source.
{
	estring										szText;							// text value extracted from szTextStart and szTextEnd
	std::string									szLastTag;						// this is the tag value at the top of the stack	
	CComVariant*								pvStack;						// pointer from the object stack
	unsigned short*								pnType;							// pointer from the object stack holding the variant type
	std::vector<void*>::reverse_iterator		it = vector_object_stack.rbegin();
	
	if (vector_tag_stack.size()) szLastTag.assign(*vector_tag_stack.rbegin());
	szText.assign(szTextStart, szTextEnd - szTextStart + 1);
	CFastStream::FromXml(&szText, false);	
	
	// we now process vText
	if (szLastTag == "Element"){
		pnType = reinterpret_cast<unsigned short*>(*it++);
		pvStack = reinterpret_cast<CComVariant*>(*it);
		if (pvStack->vt != VT_EMPTY){
			ATLASSERT(false);
			return E_FAIL;
		}
		// We change the type of vText to *pnType and detach it to pvStack.
		// If *pnType is VT_DISPATCH then we map to a string. This case denotes an object key.		
		if (*pnType == VT_DISPATCH){
			// The ONLY possibilities here are that we have an object key or a blank string.
			if (pds && szText.size()){								
				if (szText.right(12) == "NoDataSource"){
					szText = szText.left(szText.size() - 12) + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, *pds);
				} else if (szText.right(4) == "Last"){					
					szText = szText.left(szText.size() - 4) + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, *pds);
				} else if (szText.right(2) == "PL"){					
					szText = szText.left(szText.size() - 2) + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, *pds);
				} else {
					// Should never happen.
					throw "Impossible object key '" + szText + "' found";					
				}
			}
			return szText.GetValue(pvStack);
		} else {
			return pvStack->ChangeType(*pnType, &szText.GetValue());
		}
	} else if (szLastTag == "Object"){
		// This represents the value of an object. This is possible if the
		// value of the object is a scalar value (e.g. certian instances of the 
		// result handle).
		pvStack = reinterpret_cast<CComVariant*>(*it);
		if (pvStack->vt != VT_DISPATCH) return CParameterMap::ReturnErrorR(IDS_COULD_NOT_PARSE_XML);
		return CComDispatchDriverEx(pvStack->pdispVal).PutPropertyByName(L"Value", &szText.GetValue());
	} else {
		return CParameterMap::ReturnErrorR(IDS_COULD_NOT_PARSE_XML);	// We only handle the <Element> tag at present
	}
}
//
//	called exclusively by GetVariant
//
/*static*/ HRESULT CXmlStreamer::GetVariantProcessTag(const char* szTagStart, const char* szTagEnd, xml_attributes& attributes, std::vector<std::string>& vector_tag_stack, std::vector<void*>& vector_object_stack)
{
	std::string							szTag;
	HRESULT								hr = S_OK;
	std::string							szLastTag;						// this is the tag value at the top of the stack
	
	if (vector_tag_stack.size()) szLastTag.assign(*vector_tag_stack.rbegin());
	szTag.assign(szTagStart, szTagEnd - szTagStart + 1);			
	if (szTagStart[0] == '/'){
		// pop from stack
		if (::strcmp(szLastTag.c_str(), szTag.c_str() + 1)){
			// close tag is not identical to the corresponding opening tag
			std::stringstream ss;
			ss << "Expecting tag '/" << szLastTag << "', found '" << szTag << "'";
			return CParameterMap::ReturnErrorS(ss);
		}
		vector_tag_stack.pop_back();
		if (vector_tag_stack.size()){
			szLastTag.assign(*vector_tag_stack.rbegin());
		} else {
			szLastTag.erase();
		}
	} else {
		// push on stack
		vector_tag_stack.push_back(szTag);
	}
		
	if (szTag == "Element"){
		long nRow, nCol;
		unsigned short* pnType = new unsigned short;
		if (hr = attributes.GetValue("row", &nRow)) return hr;
		if (hr = attributes.GetValue("column", &nCol)) return hr;
		if (hr = attributes.GetValue("type", pnType)) return hr;
		// verify that the last object on the tag stack is a matrix
		if (szLastTag != "Matrix") return E_FAIL;
		// get the parameter map object from the object stack		
		CParameterMap* pm = reinterpret_cast<CParameterMap*>(*vector_object_stack.rbegin());
		// Push the parameter map element pointer on to the object stack.		
		CComVariant* pv = pm->GetElementPtr(nRow, nCol);
		vector_object_stack.push_back(pv);
		// push the variant type onto the stack
		vector_object_stack.push_back(pnType);		
	} else if (szTag == "/Element"){
		// The top element is the variant type.		
		unsigned short* pnType = reinterpret_cast<unsigned short*>(*vector_object_stack.rbegin());
		vector_object_stack.pop_back();
		delete pnType;		
		// The pointer at the top of the object stack is owned by another object lower down the stack.
		// Therefore, all we need to do is pop the value from the stack.						
		vector_object_stack.pop_back();
	} else if (szTag == "Matrix"){
		long							nRows, nCols;
		CParameterMap*					ppm;
		if (hr = attributes.GetValue("rows", &nRows)) return hr;
		if (hr = attributes.GetValue("columns", &nCols)) return hr;		
		ppm = new CParameterMap(nRows, nCols);
		vector_object_stack.push_back(ppm);			
	} else if (szTag == "/Matrix"){
		// The pointer at the top of the object stack is the parameter map object
		CParameterMap*					ppm = reinterpret_cast<CParameterMap*>(*vector_object_stack.rbegin());
		vector_object_stack.pop_back();				
		if (szLastTag == "Object" || !szLastTag.size()){		
			// The pointer at the top of the object stack is the object to which we might set the parameter map value
			CComVariant*					pv = reinterpret_cast<CComVariant*>(*vector_object_stack.rbegin());
			if (pv->vt == VT_DISPATCH){
				hr = CComDispatchDriver(pv->pdispVal).PutPropertyByName(L"Value", &ppm->GetValue());
			} else if (pv->vt == VT_EMPTY){
				hr = ppm->GetValue(pv);		
			} else {
				ATLASSERT(false);
				hr = E_FAIL;
			}
		} else if (szLastTag == "Element"){
			// The top element is the variant type
			std::vector<void*>::reverse_iterator		it = vector_object_stack.rbegin();
			unsigned short*								pnType = reinterpret_cast<unsigned short*>(*it++);
			CComVariant*								pvStack;
			if (*pnType != (VT_ARRAY | VT_VARIANT)){
				ATLASSERT(false);
				return E_FAIL;
			}
			pvStack = reinterpret_cast<CComVariant*>(*it);
			if (pvStack->vt != VT_EMPTY){
				ATLASSERT(false);
				return E_FAIL;
			}
			hr = ppm->GetValue(pvStack);
		} else {
			ATLASSERT(false);
			hr = E_FAIL;
		}
		delete ppm;
		return hr;		
	} else if (szTag == "Object"){		
		// We need to create this object, or in the null case, create an empty variant.
		estring							szObj;
		CComPtr<IDispatch>				spObject;
		CComVariant						v;
		estring							szCollectionIndex;
		if (hr = attributes.GetValue("name", &szObj)) return hr;
		attributes.GetValue("collectionindex", &szCollectionIndex);
		CFastStream::FromXml(&szCollectionIndex, true);

		if (szCollectionIndex.size()){
			std::vector<void*>::reverse_iterator		itVariant = vector_object_stack.rbegin() + 1;
			std::vector<void*>::reverse_iterator		itParameterMap = vector_object_stack.rbegin() + 2;
			CParameterMap*								ppm = reinterpret_cast<CParameterMap*>(*itParameterMap);
			CComVariant*								pv = reinterpret_cast<CComVariant*>(*itVariant);
			long										nRow, nCol;
			if (!ppm->GetElementRowCol(pv, &nRow, &nCol)){
				// we may need to augment this map to 3 columns
				if (ppm->GetCols() == 2){
					ppm->SetColumns(3);
				} else if (ppm->GetCols() != 3){
					ATLASSERT(false);
					return E_FAIL;
				}
				*itVariant = ppm->GetElementPtr(nRow, nCol);
				ppm->SetValue(nRow, 2, szCollectionIndex);
			}
		}
		if (szObj.CompareNoCaseAndSpace("Null")){
			if (hr = g_pApplication->GetObjectManager().CreateObject(szObj, spObject)) return hr;
			v = spObject;
		}
		VARIANT*						pvObj = new VARIANT();
		::VariantInit(pvObj);
		hr = v.Detach(pvObj);
		vector_object_stack.push_back(pvObj);
		return hr;
	} else if (szTag == "/Object"){		
		// The pointer at the top of the stack is now the heap-allocated variant containing an object.
		VARIANT*						pvObj = reinterpret_cast<VARIANT*>(*vector_object_stack.rbegin());
		// If pvObj is not an object not true then things could happen when we try to delete it.
		// In the VT_DISPATCH case all deletion does is decrease a reference count.
		// Note that the VT_EMPTY case is for the Null object
		ATLASSERT(pvObj->vt == VT_DISPATCH || pvObj->vt == VT_EMPTY);
		vector_object_stack.pop_back();
		// The stack now contains the variable type at the top followed by a pointer to a variant
		// to which we assign pvObj.
		std::vector<void*>::reverse_iterator	it = vector_object_stack.rbegin();
		if (vector_object_stack.size() != 1){		
			unsigned short* pnType = reinterpret_cast<unsigned short*>(*it++);
			ATLASSERT(*pnType == VT_DISPATCH);
		}
		CComVariant*							pv = reinterpret_cast<CComVariant*>(*it);
		if (pvObj->vt == VT_EMPTY){
			// Null object case
			pv->Clear();
			hr = S_OK;
		} else {
			hr = pv->Attach(pvObj);
		}
		::VariantClear(pvObj);
		delete pvObj;
		return hr;		
	} else {	
		return CParameterMap::ReturnErrorSS("Unhandled Xml tag", szTag);		
	}
	return S_OK;
}