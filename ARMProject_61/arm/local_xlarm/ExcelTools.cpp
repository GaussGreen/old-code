#include <deque>
#include "ExcelTools.h"
#include <iostream>
#include <exception>
using namespace std;

//	-------------------------------------------------------------------------
bool ExcelTools::isAvailable(LPXLOPER xloper)
{
	if (xloper==NULL) 
		return false ; 
	if (xloper->xltype==xltypeMissing || xloper->xltype==xltypeNil) 
		return false; 
	if (xloper->xltype==xltypeSRef) 
	{
		// reference to some cell. call ISBLANK ?
		XLOPER xlTmp; 
		if(xlretSuccess!=Excel(xlfIsblank ,&xlTmp,1,xloper)) return false ;
		if (xlTmp.xltype!=xltypeBool) return true;  // ISBLANK failed : #NA ?
		if (xlTmp.val.xbool) return false;  // ISBLANK=true : isAvailable=false.
		return true; 
	}
	return true; 
}

//	-------------------------------------------------------------------------
void ExcelTools::coerce(LPXLOPER from,LPXLOPER to,int xlType)
{
	/*if ( (!from) || (!to)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::xlCoerce: invalid arguments."); 
	if (xlretSuccess!=Excel (xlCoerce, to, 2, from, TempNum (xlType)))
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::xlCoerce: error while coercing to "<<xlType); 
	if (to->xltype==xltypeErr) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::xlCoerce: error while coercing to "<<xlType); */
}

//	-------------------------------------------------------------------------------------------------------
//	static 
void 
ExcelTools::convert(LPXLOPER xloper,std::string&ret)
{
	ret=""; 
	if (!isAvailable(xloper)) 
		return ;//(ERR_INVALID_ARGUMENT,"ExcelTools::convert: no string argument."); 
	
	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeStr) ; 

	char * pascalString = xlTmp.val.str ;
	int size = pascalString[0]; 
	if (size<0) size=size+256 ;
	ret.resize(size) ;
	ret.assign(pascalString+1,size); 
	ret[size]='\0'; 
	Excel (xlFree, 0, 1, &xlTmp);
	return ;
}

//	-------------------------------------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xloper,double& ret)
{
	ret=0; 
	// 
	if (!isAvailable(xloper)) 
		return;//ICMTHROW(ERR_INVALID_ARGUMENT,"Missing double argument");  
	
	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeNum) ; 

	ret= xlTmp.val.num; 
	return; 
}
//	-------------------------------------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xloper,bool& ret)
{
	ret=0; 
	// 
	if (!isAvailable(xloper)) 
		return ;//ICMTHROW(ERR_INVALID_ARGUMENT,"Missing double argument");  
	
	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeBool) ; 

	ret= xlTmp.val.xbool!=0; 
	return; 
}
//	-------------------------------------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xloper,int& ret)
{
	ret=0; 
	// 
	if (!isAvailable(xloper)) 
		return ;//(ERR_INVALID_ARGUMENT,"Missing integer argument");  

	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeNum) ; 
	ret = static_cast<int>(xlTmp.val.num); 
	
	return; 
}
//	-------------------------------------------------------------------------------------------------------
// static
void ExcelTools::convert(LPXLOPER xloper, std::vector<double>& ret ){
	ret.resize(0); 

	if (!isAvailable(xloper)) 
		return ;//(ERR_INVALID_ARGUMENT,"Missing array argument");

	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeMulti) ; 

	int rowNum = xlTmp.val.array.rows;
	int colNum = xlTmp.val.array.columns;
	ret.resize(rowNum* colNum); 
	long size(0); 
	for(long i = 0; i < xlTmp.val.array.rows  ; i++) 
	{
		for(long j=0;j<  xlTmp.val.array.columns ; j++) 
		{
			// assuming same storage scheme ?
			LPXLOPER px = xlTmp.val.array.lparray + size;
			convert(px,ret[i+j]); 
			size++;
		}
	}
	Excel (xlFree, 0, 1, &xlTmp);
}

//	-------------------------------------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xloper,vector<string>& ret)
{
	ret.resize(0); 

	if (!isAvailable(xloper)) 
		return ;//(ERR_INVALID_ARGUMENT,"Missing array argument");

	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeMulti) ; 
	
	int rowNum = xlTmp.val.array.rows;
	int colNum = xlTmp.val.array.columns;
	ret.resize(rowNum* colNum); 
	long size(0); 
	for(long i = 0; i < xlTmp.val.array.rows  ; i++) 
	{
		for(long j=0;j<  xlTmp.val.array.columns ; j++) 
		{
			LPXLOPER px = xlTmp.val.array.lparray + size;
			convert(px,ret[i+j]); 
			size++;
		}
	}
	Excel (xlFree, 0, 1, &xlTmp);
	return;
}
//	-------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(const std::string &item,LPXLOPER ret)
{
	ret->xltype = xltypeStr| xlbitDLLFree;
	int size = item.size() ; 
	char* tmp = (char*)malloc(sizeof(char)*size +2) ; 
	for(unsigned int i=0;i<size;i++) tmp[i+1]=item[i]; 
	tmp[0]=size; 
	tmp[size+1]='\0'; 
	ret->val.str = tmp; 
}
//	-------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(double item,LPXLOPER ret)
{
	ret->xltype= xltypeNum; 
	ret->val.num=item; 
}
//	--------------------------------------------------------------------------------------------------
// static
std::auto_ptr<ExcelCaller> ExcelCaller::itsInstance(new ExcelCaller()); 
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::setError(const std::string& coordinates, const std::string& msg)
{
	// 
	unsetError(coordinates); 
	unsetObject(coordinates);	
	//(*ARM_ERRORS_LIST)[strdup(coordinates.c_str())]=BuildErrorMessage(msg.c_str(),true) ;	
}
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::setObject(const std::string& coordinates, const std::string& objectlabel)
{
	// 
	unsetError(coordinates); 
	unsetObject(coordinates);	
	//(*ARM_OBJECTS_LIST)[strdup(coordinates.c_str())]=BuildErrorMessage(objectlabel.c_str(),false) ;	
}
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::unsetError(const std::string& coordinates)
{
	//ARM_ERRORS_LIST_t::iterator i = ARM_ERRORS_LIST->find(coordinates.c_str());
	//if (i==ARM_ERRORS_LIST->end()) return; 
	//delete i->second ;
	//ARM_ERRORS_LIST->erase(i); 
}
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::unsetObject(const std::string& coordinates)
{
	//ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->find(coordinates.c_str());
	//if (i==ARM_OBJECTS_LIST->end()) return; 
	//delete i->second ;
	//ARM_OBJECTS_LIST->erase(i); 
}
//	--------------------------------------------------------------------------------------------------
std::string 
ExcelCaller::getCaller()
{
	return "";// CCSTringToSTLString(XL_getCaller()) ;	
}
//	--------------------------------------------------------------------------------------------------
void 
ExcelCaller::setError(const std::string& msg) 
{
	setError(getCaller(),msg); 
}
//	--------------------------------------------------------------------------------------------------
std::string 
ExcelCaller::setObject(long id,const std::string& className) 
{
	std::string objName = makeObjectLabel(id,className);
	setObject(getCaller(),objName); 
	return objName; 
}
//	--------------------------------------------------------------------------------------------------
long 
ExcelCaller::getObjectId()
{
	/*ARM_OBJECTS_LIST_t::const_iterator i = ARM_OBJECTS_LIST->find(getCaller().c_str());
	if (i == ARM_OBJECTS_LIST->end()) return -1; 
	return LocalGetNumObjectId((CCString)(i->second->getMsg()));*/ 
	return 0;
}
//	--------------------------------------------------------------------------------------------------
std::string 
ExcelCaller::makeObjectLabel(long id,const std::string& objectClassName)
{
	return "";//CCSTringToSTLString(LocalMakeObjectId(id,objectClassName.c_str())); 
}
//	--------------------------------------------------------------------------------------------------
bool
ExcelCaller::isCalledByWizard()
{
	return true;// XL_IsCalledByFuncWiz()?true:false; 
}

