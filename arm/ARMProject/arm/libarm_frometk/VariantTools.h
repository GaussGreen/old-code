
#ifndef _VARIANT_TOOLS_H_
#define _VARIANT_TOOLS_H_

#include <string>
#include <atlbase.h>
#include <comdef.h>
#include <ICMKernel/util/icm_qmatrix.h>
#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel/util/icm_matrix.h"
//	
class VariantTools
{
public:
	//	----------------------------------------------------------
	//
	//		string to _bstr_t 
	//		uses this function with 'large' strings
	//
	static void convert(const std::string&s, _bstr_t& ret); 
	static void convert(const _variant_t &,std::string&ret) ; 
	static void convert(BSTR s,std::string&ret) ; 
	static void convert(const VARIANT& ret, std::vector<string>&);
	static void convert(const VARIANT& ret, std::vector<double>&);
	static void convert(const VARIANT& ret, std::vector<ARM_Date>&);
	static void convert(const VARIANT& ret, ARM_Vector&);
	static void convert(const VARIANT& var,ICM_QMatrix<double>&); 
	static void convert(const VARIANT&,ARM_Date&ret); 
	static void convert(const VARIANT&,double&ret); 
	static void convert(const VARIANT&,int&ret); 
	static void convert(const VARIANT&attrNames,const VARIANT&attrValues,const VARIANT&attrTypes, ICM_PropertyList&ref); 
	// 	
	static void convert(const std::string&, BSTR& ret);  
	static void convert(const std::string&, VARIANT& ret);  
	static void convert(double , VARIANT& ret);  
	static void convert(const ARM_Date&,VARIANT&ret); 
	static void convert(const std::vector<double> &, VARIANT& ret);  
	static void convert(const std::vector<std::string>&, VARIANT& ret); 
	static void convert(const std::vector<std::vector<std::string> > & item, VARIANT& ret);
	static void convert(const ICM_QMatrix<double>&,VARIANT&ret); 
	static void convert(const ARM_Vector&,VARIANT&ret); 
	//
	static bool IsAnArray (const VARIANT& ref) ; 
	static long getRows(const VARIANT&arg) ;
	static long getCols(const VARIANT&arg) ;
	static void extractColumnFromRowOfRow(const VARIANT& ret, unsigned int col,std::vector<string>&);
	static void convertFromRowOfRow(const VARIANT&ret,unsigned int row,unsigned int col,double&); 
	static void convert(const VARIANT&ret,unsigned int row,unsigned int col,double&); 
	static void convert(const VARIANT&ret,unsigned int row,unsigned int col,std::string&); 
	static void convert(const VARIANT&ret,unsigned int row,unsigned int col,int&); 
	static void convert(const VARIANT&ret,unsigned int row,unsigned int col,ARM_Date&); 
	// 
	static bool isXLDate(const VARIANT&); 
	static bool isAvailable(const VARIANT&) ; 
	//	----------------------------------------------------------
	//	
	//		other common conversions
	//
	static void convertXLDate(double xlJulian, ARM_Date& date); 
	//
	//
	template <class T> static void convert(BSTR arg,const std::string&def,T&ret)
	{
		std::string s; VariantTools::convert(arg,s); 
		if (s=="") s=def; 
		ICM_EnumsCnv::cnv(s,ret); 
	}
	static void convert(BSTR arg,const std::string&def,std::string&ret)
	{
		if (arg) VariantTools::convert(arg,ret); 
		else ret=def ;
	}
	template <class T> static void convert(const VARIANT& arg,const T& ref ,T&ret)
	{
		if (!isAvailable(arg)) ret=ref; 
		// if ( arg.vt==VT_EMPTY || arg.vt==VT_NULL) ret=ref;
		else convert(arg,ret); 
	}
	template <class T> static void convert(const VARIANT& arg,const std::string&def,T&ret)
	{
		std::string s; VariantTools::convert(arg,s); 
		if (s=="") s=def; 
		ICM_EnumsCnv::cnv(s,ret); 
	}
	template <class T> static void convert(const VARIANT* arg,const std::string&def,T&ret)
	{
		if (arg) { VariantTools::convert(*arg,def,ret); return ; }
		ICM_EnumsCnv::cnv(def,ret); 
	}
}; 

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
