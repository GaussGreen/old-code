#include <ARM\libarm_local\firstToBeIncluded.h>
#include "VariantTools.h"
#include <CCDate.h>
#include <ICMKernel\util\icm_macro.h>
#include <ARM\libarm_local\ARM_local_persistent.h>


//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const std::string&s, _bstr_t& ret)
{
	if (s.empty()) { ret=""; return; } 

	//	JLA: This portion of code is the A2WBSTR() function located in <atlconv.h>
	//	modified with MSDN Q241857 item.  

   USES_CONVERSION;	//JLA: declare some local variables 
   BSTR str = NULL;
   int nConvertedLen = MultiByteToWideChar(_acp, 0, s.c_str(),
     -1, NULL, NULL);
 
   // BUG FIX #1 (from Q241857): only subtract 1 from 
   // the length if the source data is nul-terminated
   nConvertedLen--;
   str = ::SysAllocStringLen(NULL, nConvertedLen);
   if(!str) 
	   ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::convert:SysAllocStringLen "<<nConvertedLen<<" bytes failed."); 
   ret = _bstr_t(str,false); 
   int output=MultiByteToWideChar(_acp, 0, s.c_str(), /*nLen**/ -1, str, nConvertedLen);
   if (output!=0  ) 
   {
		long retCode=GetLastError(); 
		std::string err ;
		switch (retCode) 
		{
		case ERROR_INSUFFICIENT_BUFFER: err ="ERROR_INSUFFICIENT_BUFFER"; break; 
		case ERROR_INVALID_FLAGS: err ="ERROR_INVALID_FLAGS"; break; 
		case ERROR_INVALID_PARAMETER: err ="ERROR_INVALID_PARAMETER"; break; 
		case ERROR_NO_UNICODE_TRANSLATION: err ="ERROR_NO_UNICODE_TRANSLATION"; break; 
		default: err ="Unknown"; break; 
		} ;
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::convert:MultiByteToWideChar failed:"<<err); 
   }
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const _variant_t &v,std::string&ret)
{
	ret=_bstr_t(v); 
}
//	------------------------------------------------------------------------------------------------
static void convert(const _variant_t & v, long &ret) {
	_variant_t V2= v;
	V2.ChangeType(VT_I4);	// might throw
	ret = V2.dblVal;
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(BSTR s,std::string&ret)
{
	if (!s) { ret=""; return; }
	ret=_bstr_t(s); 
}
//	------------------------------------------------------------------------------------------------
bool 
VariantTools::IsAnArray (const VARIANT& ref)
{
	if( ref.vt & VT_ARRAY ) return true;
	return false;
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT& vntTab,ICM_QMatrix<double>& ret)
{
	ret.Resize(0,0); 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	bool res = IsAnArray(vntTab);
	if(!res)
	{
		ICMLOG("VariantTools::convert: not an array: converting so tingle value array");
		ret.Resize(1,1);
		convert(vntTab,ret(0,0)); 
		// _variant_t cvnt (vntTab);
		// cvnt.ChangeType(VT_R8);	// might throw if not possible
		// ret(0,0)  = cvnt.dblVal;
		return;
	}
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);

	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	ret.Resize(lNumRows,lNumCols) ; 
	
	long index [2];
	index[0] = lbRows;
	VARIANT vntElt;
	for(unsigned int i =0; i< lNumRows; i++){
		index[1] = lbCols;
		for(unsigned int j =0; j< lNumCols; j++)
		{
			::SafeArrayGetElement(tmpArray,index, &vntElt) ;
			convert(vntElt,ret(i,j)); 
			// _variant_t cvnt (vntElt);
			// cvnt.ChangeType(VT_R8);	// might throw
			// ret(i,j) = cvnt.dblVal;
			index[1]++;
		}
		index[0]++;
	}
}

//	------------------------------------------------------------------------------------------------
// static
void 
VariantTools::convertXLDate(double xlJulian, ARM_Date& date)
{
	date=ARM_Date(XLDateToJulian(xlJulian)); 
}
//	------------------------------------------------------------------------------------------------
bool 
VariantTools::isXLDate(const VARIANT&ref)
{
	double ret ; 
	try { VariantTools::convert(ref,-1.,ret); }
	catch (...) { return false; }
	return isValidXLDate(ret)!=0 ;
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const std::string& item, BSTR& ret)
{
	_bstr_t tmp(item.c_str()); 
	ret=SysAllocString(tmp);
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(double item, VARIANT& ret)
{
	VariantClear(&ret); 
	ret.vt = VT_R8 ; 
	ret.dblVal = item; 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const std::string& item, VARIANT& ret)
{
	_variant_t   tmp ;
	tmp.Attach(ret) ;
	tmp.SetString(item.c_str()); 
	ret=tmp.Detach();
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const std::vector<double> & item, VARIANT& ret)
{	
	SAFEARRAYBOUND bound[1];    
	bound[0].lLbound = 0;
	bound[0].cElements = item.size();
	ret.parray = SafeArrayCreate( VT_VARIANT , 1 ,  bound ) ;
	ret.vt = VT_ARRAY | VT_VARIANT ;
 	for (long i = 0;i<item.size(); i++)
	{
		VARIANT tmp; 
		VariantInit(&tmp); 
		VariantTools::convert(item[i],tmp); 
		SafeArrayPutElement( ret.parray , &i ,&tmp) ;
	}
}
//	------------------------------------------------------------------------------------------------
void VariantTools::convert(const VARIANT& vntTab, vector<string>&  resu){
	resu.resize(0); 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	bool res = IsAnArray(vntTab);
	if(!res)
	{
		ICMLOG("VariantTools::convert: not an array: converting so tingle value array");
		_variant_t cvnt (vntTab);
		cvnt.ChangeType(VT_BSTR);	// might throw if not possible
		resu.resize(1);
		resu[0]  = cvnt.dblVal;
		return;
	}
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);

	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	long size = max(lNumRows,lNumCols);
	resu.resize(size) ; 
	
	long index [2];
	index[0] = lbRows;
	VARIANT vntElt;
	for(unsigned int i =0; i< lNumRows; i++){
		index[1] = lbCols;
		for(unsigned int j =0; j< lNumCols; j++)
		{
			::SafeArrayGetElement(tmpArray,index, &vntElt) ;
			_variant_t cvnt (vntElt);
			cvnt.ChangeType(VT_BSTR);	// might throw
			resu[i+j] = _bstr_t(cvnt.bstrVal);
			index[1]++;
		}
		index[0]++;
	}
}

void VariantTools::convert(const VARIANT& vntTab, vector<ARM_Date>&  resu){
	
	//	date=ARM_Date(XLDateToJulian(xlJulian)); 
	resu.resize(0); 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	bool res = IsAnArray(vntTab);
	if(!res)
	{
		ICMLOG("VariantTools::convert: not an array: converting so tingle value array");
		_variant_t cvnt (vntTab);
		cvnt.ChangeType(VT_BSTR);	// might throw if not possible
		resu.resize(1);
		resu[0]  = cvnt.dblVal;
		return;
	}
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);

	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	long size = max(lNumRows,lNumCols);
	resu.resize(size) ; 
	
	long index [2];
	index[0] = lbRows;
	VARIANT vntElt;
	for(unsigned int i =0; i< lNumRows; i++){
		index[1] = lbCols;
		for(unsigned int j =0; j< lNumCols; j++)
		{
			::SafeArrayGetElement(tmpArray,index, &vntElt) ;
			convert(vntElt,resu[i+j]); 
			// _variant_t cvnt (vntElt);
			// cvnt.ChangeType(VT_R8);	// long might throw
			// double d = cvnt.dblVal;
			// resu[i+j] = ARM_Date(XLDateToJulian(d));
			index[1]++;
		}
		index[0]++;
	}

}

void VariantTools::convert(const VARIANT& vntTab, vector<double>&  resu){
	resu.resize(0); 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	bool res = IsAnArray(vntTab);
	if(!res)
	{
		ICMLOG("VariantTools::convert: not an array: converting so tingle value array");
		_variant_t cvnt (vntTab);
		cvnt.ChangeType(VT_BSTR);	// might throw if not possible
		resu.resize(1);
		resu[0]  = cvnt.dblVal;
		return;
	}
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);

	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	long size = max(lNumRows,lNumCols);
	resu.resize(size) ; 
	
	long index [2];
	index[0] = lbRows;
	VARIANT vntElt;
	for(unsigned int i =0; i< lNumRows; i++){
		index[1] = lbCols;
		for(unsigned int j =0; j< lNumCols; j++)
		{
			::SafeArrayGetElement(tmpArray,index, &vntElt) ;
			convert(vntElt,resu[i+j]); 
			// _variant_t cvnt (vntElt);
			// cvnt.ChangeType(VT_R8);	// long might throw
			// resu[i+j] = cvnt.dblVal;
			index[1]++;
		}
		index[0]++;
	}

}

void VariantTools::convert(const VARIANT& vntTab, ARM_Vector&  resu){
	resu.Resize(0); 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	bool res = IsAnArray(vntTab);
	if(!res)
	{
		ICMLOG("VariantTools::convert: not an array: converting so tingle value array");
		_variant_t cvnt (vntTab);
		cvnt.ChangeType(VT_BSTR);	// might throw if not possible
		resu.Resize(1);
		resu[0]  = cvnt.dblVal;
		return;
	}
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);

	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	long size = max(lNumRows,lNumCols);
	resu.Resize(size) ; 
	
	long index [2];
	index[0] = lbRows;
	VARIANT vntElt;
	for(unsigned int i =0; i< lNumRows; i++){
		index[1] = lbCols;
		for(unsigned int j =0; j< lNumCols; j++)
		{
			::SafeArrayGetElement(tmpArray,index, &vntElt) ;
			convert(vntElt,resu[i+j]); 
			// _variant_t cvnt (vntElt);
			// cvnt.ChangeType(VT_R8);	// long might throw
			// resu[i+j] = cvnt.dblVal;
			index[1]++;
		}
		index[0]++;
	}

}

//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const std::vector<std::string> & item, VARIANT& ret)
{
	SAFEARRAYBOUND bound[1];    
	bound[0].lLbound = 0;
	bound[0].cElements = item.size();
	ret.parray = SafeArrayCreate( VT_VARIANT , 1 , bound ) ;
	ret.vt = VT_ARRAY | VT_VARIANT ;
 	for (long i = 0;i<item.size(); i++)
	{
		VARIANT tmp; 
		VariantInit(&tmp); 
		VariantTools::convert(item[i],tmp); 
		SafeArrayPutElement( ret.parray , &i ,&tmp) ;
	}
}
//	------------------------------------------------------------------------------------------------
// veteur ligne de vecteur colonne
void 
VariantTools::convert(const vector<vector<string> > & item, VARIANT& ret)
{
	SAFEARRAYBOUND bound[2];    
	// nb rows
	bound[0].lLbound = 0;
	bound[0].cElements = item.size();
	// nb cols
	bound[1].lLbound = 0;
	bound[1].cElements = item[0].size();
	ret.parray = SafeArrayCreate( VT_VARIANT , 2 , bound ) ;
	ret.vt = VT_ARRAY | VT_VARIANT ;
	for(long i=0;i<item.size();i++)
		for (long j = 0;j<item[0].size(); j++)
		{
			VARIANT tmp; 
			long l[2]; 
			VariantInit(&tmp); 
			VariantTools::convert(item[i][j],tmp); 
			l[0]=i; l[1]=j; 
			SafeArrayPutElement( ret.parray , l ,&tmp) ;
		}
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const ICM_QMatrix<double>& item, VARIANT& ret)
{
	SAFEARRAYBOUND bound[2];    
	bound[0].lLbound = 0;
	bound[0].cElements = item.Getnbrows();
	bound[1].lLbound = 0;
	bound[1].cElements = item.Getnbcols();
	ret.parray = SafeArrayCreate( VT_VARIANT , 2 ,  bound ) ;
	ret.vt = VT_ARRAY | VT_VARIANT ;
 	for (long i = 0;i<item.Getnbrows(); i++)
		for(long j=0;j<item.Getnbcols();j++)
		{
			VARIANT tmp; 
			long l[2]; 
			VariantInit(&tmp); 
			VariantTools::convert(item(i,j),tmp); 
			l[0]=i; l[1]=j; 
			SafeArrayPutElement( ret.parray , l ,&tmp) ;
		}
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const ARM_Vector& item, VARIANT& ret)
{
	SAFEARRAYBOUND bound[1];    
	bound[0].lLbound = 0;
	bound[0].cElements = item.size();
	ret.parray = SafeArrayCreate( VT_VARIANT , 1 ,  bound ) ;
	ret.vt = VT_ARRAY | VT_VARIANT ;
 	for (long i = 0;i<item.size(); i++)
	{
		VARIANT tmp; 
		long l[1]; 
		VariantInit(&tmp); 
		VariantTools::convert(item.Elt(i),tmp); 
		l[0]=i; 
		SafeArrayPutElement( ret.parray , l ,&tmp) ;
	}
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::extractColumnFromRowOfRow(const VARIANT& vntTab, unsigned int col,std::vector<string>&ret)
{
	ret.resize(0); 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	if (!IsAnArray(vntTab))
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: not an array"); 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);
	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	// 
	// if (col>=lNumCols) 
	// 	ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: out of size "<<col); 
	long size = lNumRows; 
	ret.resize(size); 

	VARIANT vntElt ;
	for(unsigned int i=0;i<lNumRows;i++)
	{
		long index[2]; 
		index[0]=lbRows+i ;
		index[1]=0; 
		::SafeArrayGetElement(tmpArray,index, &vntElt) ;
		::VariantTools::convert(vntElt,col,0,ret[i]) ; // data stored in row of row
		// _variant_t cvnt (vntElt);
		// cvnt.ChangeType(VT_BSTR);	// might throw
		// ret[i] = _bstr_t(cvnt.bstrVal);
	}
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT&vntTab,unsigned int row,unsigned int col,double&ret)
{
	// ret
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	if (!IsAnArray(vntTab))
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: not an array"); 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);
	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	// 
	if ( (col>=lNumCols) || (row>=lNumRows)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: out of size "<<col<<" " <<row); 


	long index[2]; 
	index[0] = lbRows+row;
	index[1] = lbCols+col ;
	VARIANT vntElt ;
	::SafeArrayGetElement(tmpArray,index, &vntElt) ;
	convert(vntElt,ret); 
	// _variant_t cvnt (vntElt);
	// cvnt.ChangeType(VT_R8);	// might throw
	// ret=cvnt; 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT&ref,double&ret)
{
	_variant_t cvnt (ref);
	cvnt.ChangeType(VT_R8);	// might throw
	ret=cvnt; 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT&ref,int&ret)
{
	_variant_t cvnt (ref);
	cvnt.ChangeType(VT_I4);	// might throw
	ret=cvnt.intVal; 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT&vntTab,unsigned int row,unsigned int col,std::string&ret)
{
	// ret
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	if (!IsAnArray(vntTab))
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: not an array"); 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);
	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	// 
	if ( (col>=lNumCols) || (row>=lNumRows)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: out of size "<<col<<" " <<row); 


	long index[2]; 
	index[0] = lbRows+row;
	index[1] = lbCols+col ;
	VARIANT vntElt ;
	::SafeArrayGetElement(tmpArray,index, &vntElt) ;
	VariantTools::convert(vntElt,ret); 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT&vntTab,unsigned int row,unsigned int col,int &ret)
{
	// ret
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	if (!IsAnArray(vntTab))
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: not an array"); 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);
	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	// 
	if ( (col>=lNumCols) || (row>=lNumRows)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: out of size "<<col<<" " <<row); 


	long index[2]; 
	index[0] = lbRows+row;
	index[1] = lbCols+col ;
	VARIANT vntElt ;
	::SafeArrayGetElement(tmpArray,index, &vntElt) ;
	VariantTools::convert(vntElt,ret); 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convert(const VARIANT&vntTab,unsigned int row,unsigned int col,ARM_Date &ret)
{
	// ret
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	if (!IsAnArray(vntTab))
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: not an array"); 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0, lbCols =0, ubCols =0;
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);
	long lNumRows (ubRows-lbRows+1);
	long lNumCols (ubCols-lbCols+1);
	// 
	if ( (col>=lNumCols) || (row>=lNumRows)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: out of size "<<col<<" " <<row); 


	long index[2]; 
	index[0] = lbRows+row;
	index[1] = lbCols+col ;
	VARIANT vntElt ;
	::SafeArrayGetElement(tmpArray,index, &vntElt) ;
	VariantTools::convert(vntElt,ret); 
}

//	------------------------------------------------------------------------------------------------
// static 
void VariantTools::convert(const ARM_Date&date,VARIANT&ret)
{
	SYSTEMTIME stime; 
	stime.wDay = date.GetDay() ;
	stime.wMonth = date.GetMonth()  ;	// SYSTEMTIME uses 1=January
	stime.wYear = date.GetYear() ;		// >1601 <30287  ...  
	stime.wHour= 0 ;
	stime.wMinute=0; 
	stime.wMilliseconds=0;				// unused
	stime.wSecond=0; 
	stime.wDayOfWeek= date.GetDayOfWeek() ; // SYSTEMTIME 0=Sunday 
	
	_variant_t   tmp ;
	tmp.Attach(ret) ;
	int cv1=::SystemTimeToVariantTime(&stime,&tmp.date); 
	tmp.vt=VT_DATE; 
	ret=tmp.Detach();
}
//	------------------------------------------------------------------------------------------------
// static 
void VariantTools::convert(const VARIANT&vnt,ARM_Date&ret)
{
	_variant_t   tmp(vnt) ;
	tmp.ChangeType(VT_DATE); 
	SYSTEMTIME stime; 
	int cv1=::VariantTimeToSystemTime(tmp.dblVal,&stime) ;
	ret =ARM_Date(stime.wDay,stime.wMonth,stime.wYear); 
}
//	------------------------------------------------------------------------------------------------
void 
VariantTools::convertFromRowOfRow(const VARIANT& vntTab, unsigned int row,unsigned int col,double&ret )
{
	ret=0; 
	if(vntTab.vt==VT_EMPTY) 
	{
		ICMLOG("VariantTools::convert: Variant is Empty.");
		return;
	}
	if (!IsAnArray(vntTab))
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::extractColumn: not an array"); 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0 ; 
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	long lNumRows (ubRows-lbRows+1);
	// 
	if (row>=lNumRows) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"VariantTools::convertFromRowOfRow: out of size "<<row); 
	long size = lNumRows; 


	VARIANT vntElt ;
	long index[2]; 
	index[0]=lbRows+row ;
	index[1]=0; 
	::SafeArrayGetElement(tmpArray,index, &vntElt) ;
	::VariantTools::convert(vntElt,col,0,ret) ; // data stored in row of row
	
}
//	------------------------------------------------------------------------------------------------
// static 
bool 
VariantTools::isAvailable(const VARIANT&arg)
{
	if ( arg.vt==VT_EMPTY || arg.vt==VT_NULL) return false; 
	if ( arg.vt==VT_ERROR && arg.scode==DISP_E_PARAMNOTFOUND) return false ; // MSDN Q154039
	return true; 
}
//	------------------------------------------------------------------------------------------------
// static 
void 
VariantTools::convert(const VARIANT&vAttrNames,const VARIANT&vAttrValues,const VARIANT&vAttrTypes,ICM_PropertyList&pl)
{
	if (!isAvailable(vAttrNames)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (no names)"); 
	if (!isAvailable(vAttrValues)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (no values)"); 

	std::vector<std::string> attrNames ;
	VariantTools::convert(vAttrNames,attrNames); 


	if (VariantTools::getRows(vAttrValues)!=attrNames.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (names/values size issue)"); 
	if (VariantTools::getCols(vAttrValues)!=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (values should be 1 col)"); 

	std::vector<std::string> attrTypes;
	if (isAvailable(vAttrTypes)) 
		VariantTools::convert(vAttrTypes,attrTypes); 

	for(unsigned int i=0;i<attrNames.size();i++)
	{
		std::string attrName= attrNames[i]; 
		// LPXLOPER px = xlTmp.val.array.lparray + i;
		if (i>= attrTypes.size())			{ double item; VariantTools::convert(vAttrValues,i,0,item); pl.set(attrName,item); }  // default type is double
		else if (attrTypes[i]=="")			{ double item; VariantTools::convert(vAttrValues,i,0,item) ; pl.set(attrName,item); }  // default type is double
		else if (attrTypes[i]=="string")	{ std::string item; VariantTools::convert(vAttrValues,i,0,item); pl.set(attrName,item); } 
		else if (attrTypes[i]=="double")	{ double item; VariantTools::convert(vAttrValues,i,0,item); pl.set(attrName,item); } 
		else if (attrTypes[i]=="int")		{ int item; VariantTools::convert(vAttrValues,i,0,item); pl.set(attrName,item); } 
		else if (attrTypes[i]=="object")	
		{ 
			std::string itemName; 
			VariantTools::convert(vAttrValues,i,0,itemName);
			long itemId= LocalPersistent::get().getObjectId(itemName); 
			if (itemId==-1) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't concert "<<itemName<<" to an ARM_Object"); 
			ARM_Object* item ; 
			LocalPersistent::get().convert(itemId,item); 
			pl.set(attrName,item) ;
		}  
		else if (attrTypes[i]=="date")		{ ARM_Date item; VariantTools::convert(vAttrValues,i,0,item); pl.set(attrName,item) ; }  
		else ICMTHROW(ERR_INVALID_ARGUMENT,"Can't handle argument "<<attrName); 
	}
}

//	------------------------------------------------------------------------------------------------
// static 
long 
VariantTools::getRows(const VARIANT&vntTab) 
{
	if(vntTab.vt==VT_EMPTY) return 0; 
	if (!IsAnArray(vntTab)) return 0; 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbRows = 0, ubRows= 0 ; 
	::SafeArrayGetLBound (tmpArray, 1, &lbRows);
	::SafeArrayGetUBound (tmpArray, 1, &ubRows);
	long lNumRows (ubRows-lbRows+1);
	return lNumRows ; 	
}
//	------------------------------------------------------------------------------------------------
// static 
long 
VariantTools::getCols(const VARIANT&vntTab)
{
	if(vntTab.vt==VT_EMPTY) return 0; 
	if (!IsAnArray(vntTab)) return 0; 
	SAFEARRAY* tmpArray = NULL; 
	if (vntTab.vt & VT_BYREF) tmpArray=* vntTab.pparray; 
	else tmpArray= vntTab.parray; 

	long lbCols= 0, ubCols= 0 ; 
	::SafeArrayGetLBound (tmpArray, 2, &lbCols);
	::SafeArrayGetUBound (tmpArray, 2, &ubCols);
	long lNumCols(ubCols-lbCols+1);
	return lNumCols;
}
