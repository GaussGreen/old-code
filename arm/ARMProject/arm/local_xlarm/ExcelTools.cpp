#include <deque>
#include <ARM\libarm_local\firstToBeIncluded.h>
#include "ExcelTools.h"
#include "ARM_local_interglob.h"
#include "ICMKernel/util/icm_matrix.h"

//	-------------------------------------------------------------------------
bool ExcelTools::isAvailable(LPXLOPER xloper)
{
	if (xloper==NULL) return false ; 
	if (xloper->xltype==xltypeMissing || xloper->xltype==xltypeNil) return false; 
	if (xloper->xltype==xltypeSRef) 
	{
		// reference to some cell. call ISBLANK ?
		XLOPER xlTmp; 
		if(xlretSuccess!=Excel(xlfIsblank ,&xlTmp,1,xloper)) return false ;
		if (xlTmp.xltype!=xltypeBool) return true;  // ISBLANK failed : #NA ?
		if (xlTmp.val.boolean) return false;  // ISBLANK=true : isAvailable=false.
		return true; 
	}
	return true; 
}

//	-------------------------------------------------------------------------
void ExcelTools::coerce(LPXLOPER from,LPXLOPER to,int xlType)
{
	if ( (!from) || (!to)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::xlCoerce: invalid arguments."); 
	if (xlretSuccess!=Excel (xlCoerce, to, 2, from, TempNum (xlType)))
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::xlCoerce: error while coercing to "<<xlType); 
	if (to->xltype==xltypeErr) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::xlCoerce: error while coercing to "<<xlType); 
}
//	-------------------------------------------------------------------------
void ExcelTools::convert(LPXLOPER xloper,ICM_QMatrix<double>&ret)
{
	ret.Resize(0,0); 
	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::convert: no array  argument") ; 

	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeMulti) ; 
	ret.Resize(xlTmp.val.array.rows ,xlTmp.val.array.columns ); 
	long size(0); 
	for(long i = 0; i < xlTmp.val.array.rows  ; i++) 
	{
		for(long j=0;j<  xlTmp.val.array.columns ; j++) 
		{
			// assuming same storage scheme ?
			LPXLOPER px = xlTmp.val.array.lparray + size;
			convert(px,ret(i,j)); 
			size++;
		}
	}
	Excel (xlFree, 0, 1, &xlTmp);
}	
//	-------------------------------------------------------------------------------------------------------
//	static 
void 
ExcelTools::convert(LPXLOPER xloper,std::string&ret)
{
	ret=""; 
	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::convert: no string argument."); 
	
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
void ExcelTools::convert(LPXLOPER xloper,ARM_Date&ret)
{
	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::convert: no date specified"); 
	double tmp; 
	convert(xloper,tmp); 
	ret=ARM_Date(XLDateToJulian(tmp)); 
}
//	-------------------------------------------------------------------------------------------------------
// static 
bool ExcelTools::isXLDate(LPXLOPER xloper)
{
	XLOPER xlTmp; 
	if (!isAvailable(xloper)) return false; 
	if (xlretSuccess !=Excel (xlCoerce, &xlTmp, 2, (LPXLOPER)xloper, TempNum(xltypeNum)))
		return false; 
	if (xlTmp.xltype==xltypeErr) return false ; 
	double tmp; 
	ExcelTools::convert(xloper,-1.,tmp); 
	return isValidXLDate(tmp)==1; 
}
//	-------------------------------------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xloper,double& ret)
{
	ret=0; 
	// 
	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing double argument");  
	
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
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing double argument");  
	
	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeBool) ; 

	ret= xlTmp.val.boolean!=0; 
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
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing integer argument");  

	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeNum) ; 
	ret= xlTmp.val.num; 
	
	return; 
}
//	-------------------------------------------------------------------------------------------------------
// static
void ExcelTools::convert(LPXLOPER xloper,vector<double>& ret ){
	ret.resize(0); 

	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing array argument");

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
void ExcelTools::convert(LPXLOPER xloper,ARM_Vector& ret )
{
	ret.Resize(0); 

	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing array argument");

	XLOPER xlTmp; 
	coerce(xloper,&xlTmp,xltypeMulti) ; 

	int rowNum = xlTmp.val.array.rows;
	int colNum = xlTmp.val.array.columns;
	ret.Resize(rowNum* colNum); 
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
/* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/
void ExcelTools::convert(LPXLOPER xloper,std::deque<bool>& ret ){
	ret.resize(0); 

	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing array argument");

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
// convert LPXLOPER in a vector of double wich represent date.
void ExcelTools::convert(LPXLOPER xloper,vector<ARM_Date>& ret)
{
	ret.resize(0); 

	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing array argument");

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

//	-------------------------------------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xloper,vector<string>& ret)
{
	ret.resize(0); 

	if (!isAvailable(xloper)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Missing array argument");

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
void ExcelTools::convert(const ICM_QMatrix<double>& mat,LPXLOPER ret)
{
	if (!ret) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ExcelTools::convert: destination null");

	ret->xltype = xltypeMulti; 
	ret->val.array.columns = mat.Getnbcols();
	ret->val.array.rows = mat.Getnbrows(); 
	ret->val.array.lparray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, ret->val.array.rows * ret->val.array.columns * sizeof (XLOPER));
	long size(0); 
	for(long i =0;i<ret->val.array.rows ;i++) 
		for (long j=0;j<ret->val.array.columns;j++) 
		{
			(ret->val.array.lparray + size)->xltype=xltypeNum; 
			(ret->val.array.lparray + size)->val.num=mat(i,j) ; 
			size++ ;
		}
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
//	-------------------------------------------------------------------------
// static 
void 
ExcelTools::convert(LPXLOPER xlAttrNames,LPXLOPER xlAttrValues,LPXLOPER xlAttrTypes,ICM_PropertyList&pl)
{
	if (!isAvailable(xlAttrNames)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (no names)"); 
	if (!isAvailable(xlAttrValues)) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (no values)"); 

	std::vector<std::string> attrNames ;
	ExcelTools::convert(xlAttrNames,attrNames); 

	XLOPER xlTmp; 
	coerce(xlAttrValues,&xlTmp,xltypeMulti) ; 

	if (xlTmp.val.array.rows!=attrNames.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (names/values size issue)"); 
	if (xlTmp.val.array.columns!=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create property list (values should be 1 col)"); 

	std::vector<std::string> attrTypes;
	if (isAvailable(xlAttrTypes)) 
		ExcelTools::convert(xlAttrTypes,attrTypes); 

	for(unsigned int i=0;i<attrNames.size();i++)
	{
		std::string attrName= attrNames[i]; 
		LPXLOPER px = xlTmp.val.array.lparray + i;
		if (i>= attrTypes.size())			{ double item; ExcelTools::convert(px,item); pl.set(attrName,item); }  // default type is double
		else if (attrTypes[i]=="")			{ double item; ExcelTools::convert(px,item); pl.set(attrName,item); }  // default type is double
		else if (attrTypes[i]=="string")	{ std::string item; ExcelTools::convert(px,item); pl.set(attrName,item); } 
		else if (attrTypes[i]=="double")	{ double item; ExcelTools::convert(px,item); pl.set(attrName,item); } 
		else if (attrTypes[i]=="int")		{ int item; ExcelTools::convert(px,item); pl.set(attrName,item); } 
		else if (attrTypes[i]=="object")	
		{ 
			std::string itemName; 
			ExcelTools::convert(px,itemName); 
			long itemId= LocalPersistent::get().getObjectId(itemName); 
			if (itemId==-1) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't concert "<<itemName<<" to an ARM_Object"); 
			ARM_Object* item ; 
			LocalPersistent::get().convert(itemId,item); 
			pl.set(attrName,item) ;
		}  
		else if (attrTypes[i]=="date")		{ ARM_Date item; ExcelTools::convert(px,item); pl.set(attrName,item) ; }  
		else ICMTHROW(ERR_INVALID_ARGUMENT,"Can't handle argument "<<attrName); 
	}
	
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
	(*ARM_ERRORS_LIST)[strdup(coordinates.c_str())]=BuildErrorMessage(msg.c_str(),true) ;	
}
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::setObject(const std::string& coordinates, const std::string& objectlabel)
{
	// 
	unsetError(coordinates); 
	unsetObject(coordinates);	
	(*ARM_OBJECTS_LIST)[strdup(coordinates.c_str())]=BuildErrorMessage(objectlabel.c_str(),false) ;	
}
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::unsetError(const std::string& coordinates)
{
	ARM_ERRORS_LIST_t::iterator i = ARM_ERRORS_LIST->find(coordinates.c_str());
	if (i==ARM_ERRORS_LIST->end()) return; 
	delete i->second ;
	ARM_ERRORS_LIST->erase(i); 
}
//	--------------------------------------------------------------------------------------------------
void
ExcelCaller::unsetObject(const std::string& coordinates)
{
	ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->find(coordinates.c_str());
	if (i==ARM_OBJECTS_LIST->end()) return; 
	delete i->second ;
	ARM_OBJECTS_LIST->erase(i); 
}
//	--------------------------------------------------------------------------------------------------
std::string 
ExcelCaller::getCaller()
{
	return CCSTringToSTLString(XL_getCaller()) ;	
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
	ARM_OBJECTS_LIST_t::const_iterator i = ARM_OBJECTS_LIST->find(getCaller().c_str());
	if (i == ARM_OBJECTS_LIST->end()) return -1; 
	return LocalGetNumObjectId((CCString)(i->second->getMsg())); 
}
//	--------------------------------------------------------------------------------------------------
std::string 
ExcelCaller::makeObjectLabel(long id,const std::string& objectClassName)
{
	return CCSTringToSTLString(LocalMakeObjectId(id,objectClassName.c_str())); 
}
//	--------------------------------------------------------------------------------------------------
bool
ExcelCaller::isCalledByWizard()
{
	return XL_IsCalledByFuncWiz()?true:false; 
}

