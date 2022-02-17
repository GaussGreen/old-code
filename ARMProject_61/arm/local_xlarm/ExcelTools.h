#ifndef _EXCEL_TOOLS_H_
#define _EXCEL_TOOLS_H_

#include <windows.h>
#include <vector>
#include <string.h>
#include <iostream>
#include <exception>
#include <string>
#include "xlcall.h"
#include "framewrk.h"

using namespace std;

class ExcelTools
{
public:
	//
	static void convert(LPXLOPER, string&);
	static void convert(LPXLOPER, double&); 
	static void convert(LPXLOPER, int&); 
	static void convert(LPXLOPER, bool&); 
	static void convert(LPXLOPER, std::vector<string>&); 
	static void convert(LPXLOPER, std::vector<double>&);

	static bool isXLDate(LPXLOPER);
	static bool isAvailable(LPXLOPER); // true if not missing, not null, not empty... 
	// 
	static void convert(LPXLOPER arg,const std::string&def,std::string&ret)
	{
		if (!isAvailable(arg)) { ret=def; return; }
		// if (!arg) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't convert: null argument"); 
		// if (arg->xltype==xltypeMissing || arg->xltype==xltypeNil) { ret=def; return; }
		ExcelTools::convert(arg,ret); 
	}

	template <class T> static void convert(LPXLOPER arg,const T& def,T&ret)
	{
		if (!isAvailable(arg)) { ret=def; return; }
		// if (!arg) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't convert: null argument"); 
		// if (arg->xltype==xltypeMissing || arg->xltype==xltypeNil) { ret=def; return; }
		ExcelTools::convert(arg,ret); 
	}
	template <class T> static void convert2(LPXLOPER arg,const T& def,T&ret)
	{
		if (!isAvailable(arg)) { ret=def; return; }
		// if (!arg) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't convert: null argument"); 
		// if (arg->xltype==xltypeMissing || arg->xltype==xltypeNil) { ret=def; return; }
		ExcelTools::convert2(arg,ret); 
	}
	//
	template <class T> static void econvert(LPXLOPER arg,const std::string&def,T&ret)
	{
		std::string s; 
		// if (arg->xltype==xltypeMissing || arg->xltype==xltypeNil) 
		if (!isAvailable(arg)) 
			s=def;  
		else 
			ExcelTools::convert(arg, s); 

		ICM_EnumsCnv::cnv(s,ret); 
	}
	static void convert(const std::string&item ,LPXLOPER ret); 
	static void convert(double item ,LPXLOPER ret); 
	static void convert(std::vector<long> item ,LPXLOPER ret); 
private:
	static void coerce(LPXLOPER from,LPXLOPER to,int xlType); 
} ; 
//	---------------------------------------------------------------------------------
//	
//	HelperClass that esaes the XL returns. 
//	
//		setError : 
//			will remove the item from corresonding object name list
//			will insert an item into error list
//
//		setObject
//			will remove the item from corresponding object name list
//			will insert an intem into object name list. 
//
class ExcelCaller
{
public:
	static ExcelCaller& get() 
	{
		return *itsInstance; 
	}

public:
	void setError(const std::string& msg) ;
	std::string setObject(long id,const std::string& className) ;
	long getObjectId() ;
	bool isCalledByWizard(); 
private:
	std::string makeObjectLabel(long id,const std::string& objectClassName); 
private:
	std::string getCaller() ; 
	void setError(const std::string & coordinates,const std::string& msg); 
	void setObject(const std::string& coordinates,const std::string& objectlabel); 
	void unsetError(const std::string& coordinates); 
	void unsetObject(const std::string& coordinates); 
private:
	static std::auto_ptr<ExcelCaller> itsInstance; 
}; 



#endif // _EXCEL_TOOLS_H_
