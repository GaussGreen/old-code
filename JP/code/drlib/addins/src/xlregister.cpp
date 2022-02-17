/**********************************************************************
*    Group       : EDG DR
*
*    Filename    : xlregister.c
*
*    Description : handle registration with Excel
*
*    Author      : Mark A Robson
*
*    Date        : 19 Feb 2001
*
*
**********************************************************************/



#include "edginc/config.hpp"
#include "edginc/xlapi.hpp"
#include "edginc/xlapi.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/XLConvert.hpp"


USING_DRLIB_NAMESPACE
extern "C" 
{
    EDG_XLFUNC(void) xlAutoFree (XL_OPER *x){
        if (x) {
            XLConvert::operFreeSimple(*x);
            FREE(x);
        }
    }
}

// don't bother on unix
#ifndef UNIX

#include "edginc/XLAddin.hpp"
#include "edginc/Library.hpp"
#include "edginc/Format.hpp"
#include "edginc/XLTest.hpp"
#include "edginc/AddinLib.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include "windows.h"

DRLIB_BEGIN_NAMESPACE

#define EDG_MAX_NUM_ADDIN_FNCS 800
#define EDG_MAX_NUM_ADDIN_ZP_FNCS 50

#define EXCEL_ALERT(buf) Excel(xlcAlert,0,2,TempStr(buf),TempInt(2))

/* Function wizard dll loading */
extern "C"
{
	typedef HINSTANCE TLibHandle;
}
static const char fwLibName[] = "funcwiz"; 
static const char registerFuncName[] = "QlibRegisterFunctionInfo";
typedef int
(*RegisterFunc)(const char*, const char*);
static TLibHandle funcWizHandle = LoadLibrary(fwLibName);

/* structure just used in registration */
static XLOPER   xDLL;
static bool haveValidxDLL = false;
static XL_OPER staticOperIntxltypeMulti; // for performance in coerce

static void FunctionRegister(
    const XLAddin*     addin,     /* (I) identifies function */
    const string&      addinName,  /* (I) function name */
    bool               returnHandle, /* (I) TRUE: addin func returns 
                                        a handle */ 
    const string&      category,   /* (I) category */
    const CFieldArray& paramDesc); /* (I) text description of each parameter */

static XLAddin::TCallerIsVB              callerIsVB;
static Handle::SpreadSheetCellNameMethod ExcelLoc;
static XLConvert::XLCoerceToMulti        coerceToMulti;
static XLConvert::XLFreeCoerce           xlFreeCoerce;
static SpreadSheetMode::SpreadSheetAbort xlCancelCalc;


static void GetExcel4vFuncPtr();
/* methods from excel developer's kit (at the end of the file). Note not
   exported as these are not multithread safe */
static void InitFramework(void);
static int Excel(int xlfn, LPXLOPER pxResult, int count, ...);
LPXLOPER TempNum(double d);
LPXLOPER TempStr(const char* lpstr);
LPXLOPER TempInt(short int i);
LPXLOPER TempBool(int b);
#if 0
/* these aren't actually needed - so removed to stop compiler warnings */
LPXLOPER TempErr(WORD i);
LPXLOPER TempMissing(void);
#endif


/* join two strings together */
static char * stringJoin(const char *s1, const char *s2){
    static const string routine = "stringJoin";
    if (!s1 || !s2)
    {
        throw ModelException("stringJoin", "Either/both strings are NULL");
    }
    int s1Len = strlen(s1);
    char        *joined  = NEW_ARRAY(char, s1Len + strlen(s2) + 1);
    strcpy(joined, s1);
    strcpy(joined+s1Len, s2);
    return joined;
}



static void ExcelDateOriginInitialise()
{
    static const string routine = "ExcelDateOriginInitialise";
    XLOPER      x;
    int         xlStatus;
    long        xlDateOffset = 109205L; /* Corresponds to 1/1/1901 */

    /* With thanks to ALIB.....
    ** We get the date for 1/1/1904.
    ** If date format is 1904, then the result will be zero.
    ** Otherwise the result will be 4 years worth. */    
    if ((xlStatus = Excel(xlfDate, (LPXLOPER)&x, 3, TempNum(1904),
                          TempNum(1), TempNum(1))) != xlretSuccess) {
        throw ModelException("ExcelDateOriginInitialise",
                             "Failed to get Excel date origin\n"
                             "Excel returned error code"+ 
                             Format::toString(xlStatus));
    }

    if (x.type & xltypeNum) {
        if (x.val.num <= 1) 
        { /* 1904 date system in use */
            xlDateOffset += 1462;
        }
    }
    XLConvert::setExcelDateOffsetSet(xlDateOffset);
}


/* implements kit register method in class Addin */
static void xlRegisterWrapper(
    const string&      addinName,
    const string&      category,
    CClassConstSP      dataClass,   // identifies class representing addin
    const CFieldArray& params,      // description for each parameter
    bool               handleNameParam, /* true: extra input for handle
                                           name */
    Addin::ReturnStyle returnStyle,  // how to return the output
    bool              returnIsNative, /* true - parameter is returned
                                         as an atomic type */
    Addin::MethodType  methodType,   /* constructor/instanceMethod/
                                        classMethod */
    Addin::Method      method,       // holds function pointer
    CClassConstSP      returnType)  // the type of the return parameter
{
    /* turn  addin into our XL equivalent */
    const XLAddin* addin = XLAddin::create(addinName,
                                           dataClass,
                                           params,
                                           handleNameParam,
                                           returnStyle,
                                           returnIsNative,
                                           methodType,
                                           method,
                                           returnType);

    FunctionRegister(addin, addinName, handleNameParam, category, params);
}
DRLIB_END_NAMESPACE

extern "C" 
{	
    USING_DRLIB_NAMESPACE
	/*
	*   Function wizard suppport
	*/
    EDG_XLFUNC(short) xlAutoOpen(void)
    {
        /* see comment above */
        LPXLOPER         pxDLL = &xDLL;
        int              xlStatus;
        static const string  routine = "xlAutoOpen";
        static bool      initialised = false;

        try{
            /* need to copy these strings as TempStr modifies given
               string and this causes problems on optimised build where the
               compiler pools duplicate strings */
            char     *tempStr1 = stringJoin(" 1", "");
            char     *tempStr0 = stringJoin(" ", "");
            if (!initialised)
            {
                // write errors to output file etc - but very careful 
                // that we don't fail if we can't open the file etc
                try{
                    SpreadSheetMode::on();
                    Library::configureErrorLog();
                } catch (exception){
                    // do nothing
                }

                GetExcel4vFuncPtr();
                InitFramework();
                if ((xlStatus = Excel(xlGetName, pxDLL, 0)) != xlretSuccess){
                    throw ModelException(routine, 
                                         "Excel returned error code: "+
                                         Format::toString(xlStatus)+"\n"
                                         "Failed to get xll's name\n");
                }
                haveValidxDLL = true;
                // initialise staticOperIntxltypeMulti
                staticOperIntxltypeMulti.type = xltypeInt;
                staticOperIntxltypeMulti.val.w = xltypeMulti;
                // set up date offset
                ExcelDateOriginInitialise();
                // ensure all symbols are linked in
                AddinLib::linkInAllLibs();
                /* set up other call backs to allow function to identify
                   the current location in the sheet */
                Handle::setSpreadSheetCellNameMethod(ExcelLoc);
                XLAddin::setSpreadSheetCellNameMethod(ExcelLoc);
                SpreadSheetMode::setAbortMethod(xlCancelCalc);
                /* invoke type registration */
                Library::startup();
                /* register all addin functions with excel */
                Addin::registerKit(xlRegisterWrapper);
                /* set up other call backs to allow function to identify
                   the current location in the sheet */
                XLAddin::setCallerIsVBMethod(callerIsVB);
                XLConvert::setXLCoerceMethods(coerceToMulti, xlFreeCoerce);
                // add our extra variable argument length addin functions
                /* register generic addin function */
                string xlPrefix1(Library::XL_PREFIX);
                xlPrefix1.erase(xlPrefix1.end()-1); // remove "_" - hacky
                string xlPrefix2(XLAddin::ADDIN_PREFIX_ALT);
                xlPrefix2.erase(xlPrefix2.end()-1); // remove "_" - hacky
                if (((xlStatus = 
                      Excel(xlfRegister,0,8,&xDLL,
                            TempStr(string(" EdrExcelGeneric").c_str()),
                            // note addin name is of type P
                            TempStr(string(" RPRRRRRRRRRRRRRRRRRRRRRRRRRRRR").
                                    c_str()),
                            /*123456789012345678901234567890*/
                            TempStr(string(" "+xlPrefix1).c_str()),
                            TempStr(string(" Addin function name,"
                                           "param1, "
                                           "param2,...").c_str()),
                            TempStr(tempStr1),
                            TempStr(string(" "+string(Addin::XL_TESTS)).c_str()),
                            TempStr(tempStr0))) != xlretSuccess) ||
                    /* register special addin regression tester */
                    ((xlStatus = 
                      Excel(xlfRegister,0,8,&xDLL,
                            TempStr(string(" EdrExcelCreateRegressionTest").c_str()),
                            // note addin name and file name are of type P
                            TempStr(string(" RPPRRRRRRRRRRRRRRRRRRRRRRRRRRR").c_str()),
                            /*123456789012345678901234567890*/
                            TempStr(string(" "+string(Library::XL_PREFIX)+
                                           "TEST").c_str()),
                            TempStr(string(" Addin function name,"
                                           "Filename for test, param1, "
                                           "param2,...").c_str()),
                            TempStr(tempStr1),
                            TempStr(string(" "+string(Addin::XL_TESTS)).c_str()),
                            TempStr(tempStr0))) != xlretSuccess) ||
                    /* register special low level excel kit tester */
                    ((xlStatus = 
                      Excel(xlfRegister,0,8,&xDLL,
                            TempStr(string(" EdrExcelCreateTest").c_str()),
                            // note addin name and file name are of type P
                            TempStr(string(" RPPRRRRRRRRRRRRRRRRRRRRRRRRRRR").c_str()),
                            /*123456789012345678901234567890*/
                            /* MAX_PARAMS + 1 */
                            TempStr(string(" "+string(Library::XL_PREFIX)+
                                           "XL_TEST").c_str()),
                            TempStr(string(" Addin function name,"
                                           "Filename for test, param1, "
                                           "param2,...").c_str()),
                            TempStr(tempStr1),
                            TempStr(string(" "+string(Addin::XL_TESTS)).c_str()),
                            TempStr(tempStr0))) != xlretSuccess) ||
((xlStatus = 
                      Excel(xlfRegister,0,8,&xDLL,
                            TempStr(string(" EdrExcelGeneric").c_str()),
                            // note addin name is of type P
                            TempStr(string(" RPRRRRRRRRRRRRRRRRRRRRRRRRRRRR").
                                    c_str()),
                            /*123456789012345678901234567890*/
                            TempStr(string(" "+ xlPrefix2).c_str()),
                            TempStr(string(" Addin function name,"
                                           "param1, "
                                           "param2,...").c_str()),
                            TempStr(tempStr1),
                            TempStr(string(" "+string(Addin::XL_TESTS)).c_str()),
                            TempStr(tempStr0))) != xlretSuccess) ||
                    /* register special addin regression tester */
                    ((xlStatus = 
                      Excel(xlfRegister,0,8,&xDLL,
                            TempStr(string(" EdrExcelCreateRegressionTest").c_str()),
                            // note addin name and file name are of type P
                            TempStr(string(" RPPRRRRRRRRRRRRRRRRRRRRRRRRRRR").c_str()),
                            /*123456789012345678901234567890*/
                            TempStr(string(" "+
                                           string(XLAddin::ADDIN_PREFIX_ALT)+
                                           "TEST").c_str()),
                            TempStr(string(" Addin function name,"
                                           "Filename for test, param1, "
                                           "param2,...").c_str()),
                            TempStr(tempStr1),
                            TempStr(string(" "+string(Addin::XL_TESTS)).c_str()),
                            TempStr(tempStr0))) != xlretSuccess) ||
                    /* register special low level excel kit tester */
                    ((xlStatus = 
                      Excel(xlfRegister,0,8,&xDLL,
                            TempStr(string(" EdrExcelCreateTest").c_str()),
                            // note addin name and file name are of type P
                            TempStr(string(" RPPRRRRRRRRRRRRRRRRRRRRRRRRRRR").c_str()),
                            /*123456789012345678901234567890*/
                            /* MAX_PARAMS + 1 */
                            TempStr(string(" "+string(XLAddin::ADDIN_PREFIX_ALT)+
                                           "XL_TEST").c_str()),
                            TempStr(string(" Addin function name,"
                                           "Filename for test, param1, "
                                           "param2,...").c_str()),
                            TempStr(tempStr1),
                            TempStr(string(" "+string(Addin::XL_TESTS)).c_str()),
                            TempStr(tempStr0))) != xlretSuccess))
                {
                    throw ModelException(routine, 
                                         "Excel returned error code "+
                                         Format::toString(xlStatus));
                }
                /* and then return */
                Excel(xlFree,0,1,pxDLL);
                haveValidxDLL = false;
            }
        } catch (exception& e){
            try{
                ModelException e2(e, routine);
                e2.errorLog();
            } catch (exception){
                // again ignore
            }
            EXCEL_ALERT(string(" Excel registration failed for "+
                               string(Library::XL_PREFIX)+
                               " functions").c_str());
        }
        initialised = true;
        return 1;
    }


    EDG_XLFUNC(short) xlAutoClose(void){
        //Free the function wizard dll handle
		// (decrease the count of number of times the library is loaded
		//  when count reaches 0 the library is unloaded )
		if(funcWizHandle != NULL)
			FreeLibrary(funcWizHandle);

        return 1;
    }

}
DRLIB_BEGIN_NAMESPACE

static const XLAddin* addinTypes[EDG_MAX_NUM_ADDIN_FNCS] = {0};
static int    numAddinTypes = 0;
/* for functions that take no input parameters */
static const XLAddin* addinZPTypes[EDG_MAX_NUM_ADDIN_ZP_FNCS] = {0};
static int    numAddinZPTypes = 0;
/* how many characters to turn numAddinTypes into a string */
#define EXCEL_FUNC_PREFIX "EdrAddinFunc_"
#define EXCEL_FUNC_ZP_PREFIX "EdrAddinFuncZP_"
#define MAX_LEN_ADDIN_NAME (sizeof(EXCEL_FUNC_PREFIX) -1 + 5)
#define EXCEL_MAX_STR_LEN  255

/* register a single function with Excel */
static void FunctionRegister(
    const XLAddin*     addin,     /* (I) identifies function */
    const string&      addinName,  /* (I) function name */
    bool               handleNameParam, /* (I) TRUE: addin func returns 
                                           a handle */ 
    const string&      category,   /* (I) category */
    const CFieldArray& params)     /* (I) description of each parameter */
{
    static char       routine[] = "FunctionRegister";
    int               numParams = params.size();
    char             *name      = stringJoin(" ", (XLAddin::ADDIN_PREFIX +
                                                   addinName).c_str());
    char             *nameAlt   = stringJoin(" ", (XLAddin::ADDIN_PREFIX_ALT +
                                                   addinName).c_str());
    char             *categoryL = stringJoin(" ", category.c_str());
    char              funcNum[MAX_LEN_ADDIN_NAME +1] = "";
    char              descL[EXCEL_MAX_STR_LEN +1]   = "";
    /* need one for null termination, one for XL count, one for return type
       and optionally one for handle name */
    char             *paramsL   = NEW_ARRAY(char, numParams + 4);

    /* need to copy these strings as TempStr modifies given string and this
       causes problems on optimised build where the compiler pools
       duplicate strings */
    char             *tempStr1 = stringJoin(" 1", "");
    char             *tempStr0 = stringJoin(" ", "");

	/* Funtion wizard library */
	RegisterFunc funcWizRegister=NULL;
	if(funcWizHandle != NULL)
		funcWizRegister = (RegisterFunc) GetProcAddress(funcWizHandle,registerFuncName);
		

    int xlNumParams = numParams + (handleNameParam? 1: 0);
    if ((xlNumParams > 0 && numAddinTypes == EDG_MAX_NUM_ADDIN_FNCS) ||
        (xlNumParams == 0 && numAddinZPTypes == EDG_MAX_NUM_ADDIN_ZP_FNCS)){
        throw ModelException(routine, string("max number of ")+ 
                             (xlNumParams == 0? "zero parameter ": "")+
                             "addins exceeded - this file "
                             "needs to be updated\n");
    }
    if (handleNameParam){
        strcpy(paramsL, " RR");
        strcpy(descL, " Handle Name");
    } else {
        strcpy(paramsL, " R");
    }
    if (xlNumParams == 0){
        strcpy(descL, " ");  // TempStr is a bit weak
    }
    int   idx, jdx, descLen;
    for (idx = 0, jdx = strlen(paramsL+1), descLen = strlen(descL);
         idx < numParams; idx++)
    {
        const string myParamDesc = params[idx]->getName();
        int len = myParamDesc.size() +1; /* +1 for comma */
        jdx++;
        paramsL[jdx] = 'R';
        if (len + descLen > EXCEL_MAX_STR_LEN +1) {
            len = EXCEL_MAX_STR_LEN - descLen +1;
        }
        if (len > 0) {
            descL[descLen] = ','; /* desc[0] is overwritten by TempStr */
            strncpy(descL +1 + descLen, myParamDesc.c_str(), len-1);
            descLen += len;
        }
    }
        
    if (xlNumParams > EDG_EXCEL_MAX_PARAMS) {
        throw ModelException(routine, "Excel cannot cope with more than "+
                             Format::toString(EDG_EXCEL_MAX_PARAMS)+
                             " parameters\n"+addinName+" has "+
                             Format::toString(xlNumParams));
    }
    if (xlNumParams == 0){
        sprintf(funcNum, " " EXCEL_FUNC_ZP_PREFIX "%d", numAddinZPTypes);
        addinZPTypes[numAddinZPTypes] = addin;
        numAddinZPTypes++;
    } else {
        sprintf(funcNum, " " EXCEL_FUNC_PREFIX "%d", numAddinTypes);
        addinTypes[numAddinTypes] = addin;
        numAddinTypes++;
    }
    
    int    xlStatus;
    if (haveValidxDLL && 
        (((xlStatus = Excel(xlfRegister,0,8,&xDLL,
                            TempStr(funcNum),
                            TempStr(paramsL),
                            TempStr(name),
                            TempStr(descL),
                            TempStr(tempStr1),
                            TempStr(categoryL),
                            TempStr(tempStr0))) != xlretSuccess) ||
         (xlStatus = Excel(xlfRegister,0,8,&xDLL,
                           TempStr(funcNum),
                           TempStr(paramsL),
                           TempStr(nameAlt),
                           TempStr(descL),
                           TempStr(tempStr1),
                           TempStr(categoryL),
                           TempStr(tempStr0))) != xlretSuccess)) {
        // so it might leak - but it really doesn't matter here
        throw ModelException(routine, "Excel returned error code "+
                             Format::toString(xlStatus));
    }
	
    /* Function wizard registration */
    const char* cName;
	if(name[0]==' ')
		cName = name+1;
	else
		cName = name;
	const char* cCat = category.c_str();
	//Function wizard is an option so if library or function are not
	// found just do nothing
	if(funcWizHandle != NULL && funcWizRegister != NULL)
		funcWizRegister(cName,cCat);


    FREE(paramsL);
    FREE(name);
    FREE(nameAlt);
    FREE(categoryL);
    FREE(tempStr0);
    FREE(tempStr1);
}

/** Requests conversion of XLOPER to type xlmulti. If successful, output
    must be passed to XLFree */
static int coerceToMulti(const XLOPER& input, XLOPER& output){
    return Excel(xlCoerce, &output, 2, &input, &staticOperIntxltypeMulti);
}

/** Free the contents of the OPER created in XLCoerceToMulti method */
static void xlFreeCoerce (XLOPER& toFreeContentsOf){
    Excel(xlFree, 0, 1, &toFreeContentsOf);
}

/* are we being called from VB */
static bool callerIsVB(void)
{
    bool inVB = false; /* default */
    XLOPER   ref;
    if (Excel(xlfCaller, &ref, 0) == xlretSuccess &&
        !(ref.type & xltypeSRef))
    {
        /* called from vb */
        inVB = true;
    }
    Excel(xlFree,0,1,&ref);         
    return inVB;
}

/* returns string identifying current position */
static string ExcelLoc()
{
    static const string routine = "ExcelLoc";
    string cellName;
    bool   freeRef     = false;
    bool   freeSheet   = false;   
    bool   freeAddress = false;
    XLOPER ref, sheet, address;
    try{

        int  sheetNameLen, addressLen;
    
        if (Excel(xlfCaller, &ref, 0) != xlretSuccess){              
            throw ModelException(routine, "Excel(xlfCaller) failed");
        }
        freeRef = true;

        if (!(ref.type & xltypeSRef)) {
            /* called from vb */
            cellName = "vb";
        } else {
            if (Excel(xlSheetNm, &sheet, 1, &ref) != xlretSuccess ||
                (freeSheet = true, !(sheet.type & xltypeStr))) {         
                throw ModelException(routine, "Excel(xlSheetNm) failed");
            }
            
            if (Excel(xlfAddress,
                      &address,
                      3,
                      TempNum(ref.val.sref.ref.rwFirst + 1),
                      TempNum(ref.val.sref.ref.colFirst + 1),
                      TempNum(4)) != xlretSuccess ||
                (freeAddress = true, !(address.type & xltypeStr))){
                throw ModelException(routine, "Excel(xlfAddress) failed");
            }
            
            sheetNameLen = sheet.val.str[0];
            addressLen   = address.val.str[0];
            char *name = NEW_ARRAY(char, sheetNameLen + addressLen + 2);
            strncpy (name, sheet.val.str + 1, sheetNameLen);
            name[sheetNameLen] = '!';
            strncpy (name + sheetNameLen + 1, address.val.str + 1, addressLen);
            name[sheetNameLen + addressLen + 1] = '\0';
            cellName = string(name);
            FREE(name);
        }
    } catch (ModelException& e){
        e.errorLog();
    }
    catch (exception& e){
        ModelException e2(e, routine);
        e2.errorLog();
    }

    if (freeRef)
    {
        Excel(xlFree,0,1,&ref);         
    }
    if (freeSheet)
    {
        Excel(xlFree,0,1,&sheet);
    }
    if (freeAddress)
    {
        Excel(xlFree,0,1,&address);
    }                       
    return cellName;
}

// CALLBACK fn, used to test if user want to cancel
static bool xlCancelCalc()
{
    XLOPER xAbortOper;
    if (Excel(xlAbort, &xAbortOper, 0) == xlretSuccess) {
        if (xAbortOper.val.boolVal) {
            /* this dialog box does not display if function called
               from spreadsheet but wroks from VBA ??? */
            XLOPER xConfirmOper;
            Excel(xlcAlert, &xConfirmOper, 2, 
                  TempStr(" Stop calculation?"), TempNum(1));
            if (xConfirmOper.val.boolVal) {
                Excel(xlcMessage, 0, 1, TempBool(0));	// clean sthing ?
                return true;
            } else {
                Excel(xlAbort, 0, 1, TempBool(0)); /* clean abort if we do not
                                                      want to stop */
                return false;
            }
        }
    }
    return false;
}

extern "C" 
{
    /* Generic way of calling any addin function via EDR(...) or QLIB(...) */
    EDG_XLFUNC(XL_OPER *) EdrExcelGeneric(
        XL_OPER *addinNameOper, /* (I) */
        XL_OPER *a0,            /* (I) first parameter */
        ...)
    {
        va_list  args;
        XL_OPER    *xlOutput;
        XL_OPER     fileNameOper;
        fileNameOper.type = xltypeNil; // as long it's not a string
        /* get hold of all the parameters */
        va_start(args, a0);
        /* and call generic function */
        xlOutput = XLTest::generic(addinNameOper, &fileNameOper, a0, args);
        va_end(args);
        return xlOutput;
    }
    
    /* Creates regression test for all standard addin functions. The
       input and the output 'interface structure' are serialised out */
    EDG_XLFUNC(XL_OPER *) EdrExcelCreateRegressionTest(
        XL_OPER *addinNameOper, /* (I) */
        XL_OPER *fileNameOper,  /* (I) */
        XL_OPER *a0,            /* (I) first parameter */
        ...)
    {
        va_list  args;
        XL_OPER    *xlOutput;
        /* get hold of all the parameters */
        va_start(args, a0);
        /* and call generic function */
        xlOutput = XLTest::generic(addinNameOper, fileNameOper, a0, args);
        va_end(args);
        return xlOutput;
    }
    
    /* Create platform interface regression test file - captures explicit 
       parameters passed from and to excel (ie serialises opers) */
    EDG_XLFUNC(XL_OPER *) EdrExcelCreateTest(
        XL_OPER *addinNameOper, /* (I) */
        XL_OPER *fileNameOper,  /* (I) */
        ...)
    {
        va_list  args;
        XL_OPER    *xlOutput;
        /* get hold of all the parameters */
        va_start(args, fileNameOper);
        /* and call generic low level XL function */
        xlOutput = XLTest::genericXL(callerIsVB(), addinNameOper, 
                                     fileNameOper, args);
        va_end(args);
        return xlOutput;
    }
    
    
    /* generic function addin function 0 parameters */
#define EDG_ADDIN_FUNC_ZP_MACRO(addinNum)  \
EDG_XLFUNC(XL_OPER *) EdrAddinFuncZP_##addinNum() \
{ \
    return addinZPTypes[addinNum]->execute(NULL, NULL); \
}
    
    /* generic function addin function (for > 0 parameters) */
#define EDG_ADDIN_FUNC_MACRO(addinNum)  \
EDG_XLFUNC(XL_OPER *) EdrAddinFunc_##addinNum(XL_OPER *a0, ...) \
{ \
    XL_OPER       *output; \
    va_list     args; \
    va_start(args, a0); \
    output = addinTypes[addinNum]->execute(a0, args); \
    va_end(args); \
    return output; \
}

/* there must be at least EDG_MAX_NUM_ADDIN_ZP_FNCS number of these */
EDG_ADDIN_FUNC_ZP_MACRO(0);
EDG_ADDIN_FUNC_ZP_MACRO(1);
EDG_ADDIN_FUNC_ZP_MACRO(2);
EDG_ADDIN_FUNC_ZP_MACRO(3);
EDG_ADDIN_FUNC_ZP_MACRO(4);
EDG_ADDIN_FUNC_ZP_MACRO(5);
EDG_ADDIN_FUNC_ZP_MACRO(6);
EDG_ADDIN_FUNC_ZP_MACRO(7);
EDG_ADDIN_FUNC_ZP_MACRO(8);
EDG_ADDIN_FUNC_ZP_MACRO(9);
EDG_ADDIN_FUNC_ZP_MACRO(10);
EDG_ADDIN_FUNC_ZP_MACRO(11);
EDG_ADDIN_FUNC_ZP_MACRO(12);
EDG_ADDIN_FUNC_ZP_MACRO(13);
EDG_ADDIN_FUNC_ZP_MACRO(14);
EDG_ADDIN_FUNC_ZP_MACRO(15);
EDG_ADDIN_FUNC_ZP_MACRO(16);
EDG_ADDIN_FUNC_ZP_MACRO(17);
EDG_ADDIN_FUNC_ZP_MACRO(18);
EDG_ADDIN_FUNC_ZP_MACRO(19);
EDG_ADDIN_FUNC_ZP_MACRO(20);
EDG_ADDIN_FUNC_ZP_MACRO(21);
EDG_ADDIN_FUNC_ZP_MACRO(22);
EDG_ADDIN_FUNC_ZP_MACRO(23);
EDG_ADDIN_FUNC_ZP_MACRO(24);
EDG_ADDIN_FUNC_ZP_MACRO(25);
EDG_ADDIN_FUNC_ZP_MACRO(26);
EDG_ADDIN_FUNC_ZP_MACRO(27);
EDG_ADDIN_FUNC_ZP_MACRO(28);
EDG_ADDIN_FUNC_ZP_MACRO(29);
EDG_ADDIN_FUNC_ZP_MACRO(30);
EDG_ADDIN_FUNC_ZP_MACRO(31);
EDG_ADDIN_FUNC_ZP_MACRO(32);
EDG_ADDIN_FUNC_ZP_MACRO(33);
EDG_ADDIN_FUNC_ZP_MACRO(34);
EDG_ADDIN_FUNC_ZP_MACRO(35);
EDG_ADDIN_FUNC_ZP_MACRO(36);
EDG_ADDIN_FUNC_ZP_MACRO(37);
EDG_ADDIN_FUNC_ZP_MACRO(38);
EDG_ADDIN_FUNC_ZP_MACRO(39);
EDG_ADDIN_FUNC_ZP_MACRO(40);
EDG_ADDIN_FUNC_ZP_MACRO(41);
EDG_ADDIN_FUNC_ZP_MACRO(42);
EDG_ADDIN_FUNC_ZP_MACRO(43);
EDG_ADDIN_FUNC_ZP_MACRO(44);
EDG_ADDIN_FUNC_ZP_MACRO(45);
EDG_ADDIN_FUNC_ZP_MACRO(46);
EDG_ADDIN_FUNC_ZP_MACRO(47);
EDG_ADDIN_FUNC_ZP_MACRO(48);
EDG_ADDIN_FUNC_ZP_MACRO(49);

/* there must be at least EDG_MAX_NUM_ADDIN_FNCS number of these */
EDG_ADDIN_FUNC_MACRO(0);
EDG_ADDIN_FUNC_MACRO(1);
EDG_ADDIN_FUNC_MACRO(2);
EDG_ADDIN_FUNC_MACRO(3);
EDG_ADDIN_FUNC_MACRO(4);
EDG_ADDIN_FUNC_MACRO(5);
EDG_ADDIN_FUNC_MACRO(6);
EDG_ADDIN_FUNC_MACRO(7);
EDG_ADDIN_FUNC_MACRO(8);
EDG_ADDIN_FUNC_MACRO(9);
EDG_ADDIN_FUNC_MACRO(10);
EDG_ADDIN_FUNC_MACRO(11);
EDG_ADDIN_FUNC_MACRO(12);
EDG_ADDIN_FUNC_MACRO(13);
EDG_ADDIN_FUNC_MACRO(14);
EDG_ADDIN_FUNC_MACRO(15);
EDG_ADDIN_FUNC_MACRO(16);
EDG_ADDIN_FUNC_MACRO(17);
EDG_ADDIN_FUNC_MACRO(18);
EDG_ADDIN_FUNC_MACRO(19);
EDG_ADDIN_FUNC_MACRO(20);
EDG_ADDIN_FUNC_MACRO(21);
EDG_ADDIN_FUNC_MACRO(22);
EDG_ADDIN_FUNC_MACRO(23);
EDG_ADDIN_FUNC_MACRO(24);
EDG_ADDIN_FUNC_MACRO(25);
EDG_ADDIN_FUNC_MACRO(26);
EDG_ADDIN_FUNC_MACRO(27);
EDG_ADDIN_FUNC_MACRO(28);
EDG_ADDIN_FUNC_MACRO(29);
EDG_ADDIN_FUNC_MACRO(30);
EDG_ADDIN_FUNC_MACRO(31);
EDG_ADDIN_FUNC_MACRO(32);
EDG_ADDIN_FUNC_MACRO(33);
EDG_ADDIN_FUNC_MACRO(34);
EDG_ADDIN_FUNC_MACRO(35);
EDG_ADDIN_FUNC_MACRO(36);
EDG_ADDIN_FUNC_MACRO(37);
EDG_ADDIN_FUNC_MACRO(38);
EDG_ADDIN_FUNC_MACRO(39);
EDG_ADDIN_FUNC_MACRO(40);
EDG_ADDIN_FUNC_MACRO(41);
EDG_ADDIN_FUNC_MACRO(42);
EDG_ADDIN_FUNC_MACRO(43);
EDG_ADDIN_FUNC_MACRO(44);
EDG_ADDIN_FUNC_MACRO(45);
EDG_ADDIN_FUNC_MACRO(46);
EDG_ADDIN_FUNC_MACRO(47);
EDG_ADDIN_FUNC_MACRO(48);
EDG_ADDIN_FUNC_MACRO(49);
EDG_ADDIN_FUNC_MACRO(50);
EDG_ADDIN_FUNC_MACRO(51);
EDG_ADDIN_FUNC_MACRO(52);
EDG_ADDIN_FUNC_MACRO(53);
EDG_ADDIN_FUNC_MACRO(54);
EDG_ADDIN_FUNC_MACRO(55);
EDG_ADDIN_FUNC_MACRO(56);
EDG_ADDIN_FUNC_MACRO(57);
EDG_ADDIN_FUNC_MACRO(58);
EDG_ADDIN_FUNC_MACRO(59);
EDG_ADDIN_FUNC_MACRO(60);
EDG_ADDIN_FUNC_MACRO(61);
EDG_ADDIN_FUNC_MACRO(62);
EDG_ADDIN_FUNC_MACRO(63);
EDG_ADDIN_FUNC_MACRO(64);
EDG_ADDIN_FUNC_MACRO(65);
EDG_ADDIN_FUNC_MACRO(66);
EDG_ADDIN_FUNC_MACRO(67);
EDG_ADDIN_FUNC_MACRO(68);
EDG_ADDIN_FUNC_MACRO(69);
EDG_ADDIN_FUNC_MACRO(70);
EDG_ADDIN_FUNC_MACRO(71);
EDG_ADDIN_FUNC_MACRO(72);
EDG_ADDIN_FUNC_MACRO(73);
EDG_ADDIN_FUNC_MACRO(74);
EDG_ADDIN_FUNC_MACRO(75);
EDG_ADDIN_FUNC_MACRO(76);
EDG_ADDIN_FUNC_MACRO(77);
EDG_ADDIN_FUNC_MACRO(78);
EDG_ADDIN_FUNC_MACRO(79);
EDG_ADDIN_FUNC_MACRO(80);
EDG_ADDIN_FUNC_MACRO(81);
EDG_ADDIN_FUNC_MACRO(82);
EDG_ADDIN_FUNC_MACRO(83);
EDG_ADDIN_FUNC_MACRO(84);
EDG_ADDIN_FUNC_MACRO(85);
EDG_ADDIN_FUNC_MACRO(86);
EDG_ADDIN_FUNC_MACRO(87);
EDG_ADDIN_FUNC_MACRO(88);
EDG_ADDIN_FUNC_MACRO(89);
EDG_ADDIN_FUNC_MACRO(90);
EDG_ADDIN_FUNC_MACRO(91);
EDG_ADDIN_FUNC_MACRO(92);
EDG_ADDIN_FUNC_MACRO(93);
EDG_ADDIN_FUNC_MACRO(94);
EDG_ADDIN_FUNC_MACRO(95);
EDG_ADDIN_FUNC_MACRO(96);
EDG_ADDIN_FUNC_MACRO(97);
EDG_ADDIN_FUNC_MACRO(98);
EDG_ADDIN_FUNC_MACRO(99);
EDG_ADDIN_FUNC_MACRO(100);
EDG_ADDIN_FUNC_MACRO(101);
EDG_ADDIN_FUNC_MACRO(102);
EDG_ADDIN_FUNC_MACRO(103);
EDG_ADDIN_FUNC_MACRO(104);
EDG_ADDIN_FUNC_MACRO(105);
EDG_ADDIN_FUNC_MACRO(106);
EDG_ADDIN_FUNC_MACRO(107);
EDG_ADDIN_FUNC_MACRO(108);
EDG_ADDIN_FUNC_MACRO(109);
EDG_ADDIN_FUNC_MACRO(110);
EDG_ADDIN_FUNC_MACRO(111);
EDG_ADDIN_FUNC_MACRO(112);
EDG_ADDIN_FUNC_MACRO(113);
EDG_ADDIN_FUNC_MACRO(114);
EDG_ADDIN_FUNC_MACRO(115);
EDG_ADDIN_FUNC_MACRO(116);
EDG_ADDIN_FUNC_MACRO(117);
EDG_ADDIN_FUNC_MACRO(118);
EDG_ADDIN_FUNC_MACRO(119);
EDG_ADDIN_FUNC_MACRO(120);
EDG_ADDIN_FUNC_MACRO(121);
EDG_ADDIN_FUNC_MACRO(122);
EDG_ADDIN_FUNC_MACRO(123);
EDG_ADDIN_FUNC_MACRO(124);
EDG_ADDIN_FUNC_MACRO(125);
EDG_ADDIN_FUNC_MACRO(126);
EDG_ADDIN_FUNC_MACRO(127);
EDG_ADDIN_FUNC_MACRO(128);
EDG_ADDIN_FUNC_MACRO(129);
EDG_ADDIN_FUNC_MACRO(130);
EDG_ADDIN_FUNC_MACRO(131);
EDG_ADDIN_FUNC_MACRO(132);
EDG_ADDIN_FUNC_MACRO(133);
EDG_ADDIN_FUNC_MACRO(134);
EDG_ADDIN_FUNC_MACRO(135);
EDG_ADDIN_FUNC_MACRO(136);
EDG_ADDIN_FUNC_MACRO(137);
EDG_ADDIN_FUNC_MACRO(138);
EDG_ADDIN_FUNC_MACRO(139);
EDG_ADDIN_FUNC_MACRO(140);
EDG_ADDIN_FUNC_MACRO(141);
EDG_ADDIN_FUNC_MACRO(142);
EDG_ADDIN_FUNC_MACRO(143);
EDG_ADDIN_FUNC_MACRO(144);
EDG_ADDIN_FUNC_MACRO(145);
EDG_ADDIN_FUNC_MACRO(146);
EDG_ADDIN_FUNC_MACRO(147);
EDG_ADDIN_FUNC_MACRO(148);
EDG_ADDIN_FUNC_MACRO(149);
EDG_ADDIN_FUNC_MACRO(150);
EDG_ADDIN_FUNC_MACRO(151);
EDG_ADDIN_FUNC_MACRO(152);
EDG_ADDIN_FUNC_MACRO(153);
EDG_ADDIN_FUNC_MACRO(154);
EDG_ADDIN_FUNC_MACRO(155);
EDG_ADDIN_FUNC_MACRO(156);
EDG_ADDIN_FUNC_MACRO(157);
EDG_ADDIN_FUNC_MACRO(158);
EDG_ADDIN_FUNC_MACRO(159);
EDG_ADDIN_FUNC_MACRO(160);
EDG_ADDIN_FUNC_MACRO(161);
EDG_ADDIN_FUNC_MACRO(162);
EDG_ADDIN_FUNC_MACRO(163);
EDG_ADDIN_FUNC_MACRO(164);
EDG_ADDIN_FUNC_MACRO(165);
EDG_ADDIN_FUNC_MACRO(166);
EDG_ADDIN_FUNC_MACRO(167);
EDG_ADDIN_FUNC_MACRO(168);
EDG_ADDIN_FUNC_MACRO(169);
EDG_ADDIN_FUNC_MACRO(170);
EDG_ADDIN_FUNC_MACRO(171);
EDG_ADDIN_FUNC_MACRO(172);
EDG_ADDIN_FUNC_MACRO(173);
EDG_ADDIN_FUNC_MACRO(174);
EDG_ADDIN_FUNC_MACRO(175);
EDG_ADDIN_FUNC_MACRO(176);
EDG_ADDIN_FUNC_MACRO(177);
EDG_ADDIN_FUNC_MACRO(178);
EDG_ADDIN_FUNC_MACRO(179);
EDG_ADDIN_FUNC_MACRO(180);
EDG_ADDIN_FUNC_MACRO(181);
EDG_ADDIN_FUNC_MACRO(182);
EDG_ADDIN_FUNC_MACRO(183);
EDG_ADDIN_FUNC_MACRO(184);
EDG_ADDIN_FUNC_MACRO(185);
EDG_ADDIN_FUNC_MACRO(186);
EDG_ADDIN_FUNC_MACRO(187);
EDG_ADDIN_FUNC_MACRO(188);
EDG_ADDIN_FUNC_MACRO(189);
EDG_ADDIN_FUNC_MACRO(190);
EDG_ADDIN_FUNC_MACRO(191);
EDG_ADDIN_FUNC_MACRO(192);
EDG_ADDIN_FUNC_MACRO(193);
EDG_ADDIN_FUNC_MACRO(194);
EDG_ADDIN_FUNC_MACRO(195);
EDG_ADDIN_FUNC_MACRO(196);
EDG_ADDIN_FUNC_MACRO(197);
EDG_ADDIN_FUNC_MACRO(198);
EDG_ADDIN_FUNC_MACRO(199);
EDG_ADDIN_FUNC_MACRO(200);
EDG_ADDIN_FUNC_MACRO(201);
EDG_ADDIN_FUNC_MACRO(202);
EDG_ADDIN_FUNC_MACRO(203);
EDG_ADDIN_FUNC_MACRO(204);
EDG_ADDIN_FUNC_MACRO(205);
EDG_ADDIN_FUNC_MACRO(206);
EDG_ADDIN_FUNC_MACRO(207);
EDG_ADDIN_FUNC_MACRO(208);
EDG_ADDIN_FUNC_MACRO(209);
EDG_ADDIN_FUNC_MACRO(210);
EDG_ADDIN_FUNC_MACRO(211);
EDG_ADDIN_FUNC_MACRO(212);
EDG_ADDIN_FUNC_MACRO(213);
EDG_ADDIN_FUNC_MACRO(214);
EDG_ADDIN_FUNC_MACRO(215);
EDG_ADDIN_FUNC_MACRO(216);
EDG_ADDIN_FUNC_MACRO(217);
EDG_ADDIN_FUNC_MACRO(218);
EDG_ADDIN_FUNC_MACRO(219);
EDG_ADDIN_FUNC_MACRO(220);
EDG_ADDIN_FUNC_MACRO(221);
EDG_ADDIN_FUNC_MACRO(222);
EDG_ADDIN_FUNC_MACRO(223);
EDG_ADDIN_FUNC_MACRO(224);
EDG_ADDIN_FUNC_MACRO(225);
EDG_ADDIN_FUNC_MACRO(226);
EDG_ADDIN_FUNC_MACRO(227);
EDG_ADDIN_FUNC_MACRO(228);
EDG_ADDIN_FUNC_MACRO(229);
EDG_ADDIN_FUNC_MACRO(230);
EDG_ADDIN_FUNC_MACRO(231);
EDG_ADDIN_FUNC_MACRO(232);
EDG_ADDIN_FUNC_MACRO(233);
EDG_ADDIN_FUNC_MACRO(234);
EDG_ADDIN_FUNC_MACRO(235);
EDG_ADDIN_FUNC_MACRO(236);
EDG_ADDIN_FUNC_MACRO(237);
EDG_ADDIN_FUNC_MACRO(238);
EDG_ADDIN_FUNC_MACRO(239);
EDG_ADDIN_FUNC_MACRO(240);
EDG_ADDIN_FUNC_MACRO(241);
EDG_ADDIN_FUNC_MACRO(242);
EDG_ADDIN_FUNC_MACRO(243);
EDG_ADDIN_FUNC_MACRO(244);
EDG_ADDIN_FUNC_MACRO(245);
EDG_ADDIN_FUNC_MACRO(246);
EDG_ADDIN_FUNC_MACRO(247);
EDG_ADDIN_FUNC_MACRO(248);
EDG_ADDIN_FUNC_MACRO(249);
EDG_ADDIN_FUNC_MACRO(250);
EDG_ADDIN_FUNC_MACRO(251);
EDG_ADDIN_FUNC_MACRO(252);
EDG_ADDIN_FUNC_MACRO(253);
EDG_ADDIN_FUNC_MACRO(254);
EDG_ADDIN_FUNC_MACRO(255);
EDG_ADDIN_FUNC_MACRO(256);
EDG_ADDIN_FUNC_MACRO(257);
EDG_ADDIN_FUNC_MACRO(258);
EDG_ADDIN_FUNC_MACRO(259);
EDG_ADDIN_FUNC_MACRO(260);
EDG_ADDIN_FUNC_MACRO(261);
EDG_ADDIN_FUNC_MACRO(262);
EDG_ADDIN_FUNC_MACRO(263);
EDG_ADDIN_FUNC_MACRO(264);
EDG_ADDIN_FUNC_MACRO(265);
EDG_ADDIN_FUNC_MACRO(266);
EDG_ADDIN_FUNC_MACRO(267);
EDG_ADDIN_FUNC_MACRO(268);
EDG_ADDIN_FUNC_MACRO(269);
EDG_ADDIN_FUNC_MACRO(270);
EDG_ADDIN_FUNC_MACRO(271);
EDG_ADDIN_FUNC_MACRO(272);
EDG_ADDIN_FUNC_MACRO(273);
EDG_ADDIN_FUNC_MACRO(274);
EDG_ADDIN_FUNC_MACRO(275);
EDG_ADDIN_FUNC_MACRO(276);
EDG_ADDIN_FUNC_MACRO(277);
EDG_ADDIN_FUNC_MACRO(278);
EDG_ADDIN_FUNC_MACRO(279);
EDG_ADDIN_FUNC_MACRO(280);
EDG_ADDIN_FUNC_MACRO(281);
EDG_ADDIN_FUNC_MACRO(282);
EDG_ADDIN_FUNC_MACRO(283);
EDG_ADDIN_FUNC_MACRO(284);
EDG_ADDIN_FUNC_MACRO(285);
EDG_ADDIN_FUNC_MACRO(286);
EDG_ADDIN_FUNC_MACRO(287);
EDG_ADDIN_FUNC_MACRO(288);
EDG_ADDIN_FUNC_MACRO(289);
EDG_ADDIN_FUNC_MACRO(290);
EDG_ADDIN_FUNC_MACRO(291);
EDG_ADDIN_FUNC_MACRO(292);
EDG_ADDIN_FUNC_MACRO(293);
EDG_ADDIN_FUNC_MACRO(294);
EDG_ADDIN_FUNC_MACRO(295);
EDG_ADDIN_FUNC_MACRO(296);
EDG_ADDIN_FUNC_MACRO(297);
EDG_ADDIN_FUNC_MACRO(298);
EDG_ADDIN_FUNC_MACRO(299);
EDG_ADDIN_FUNC_MACRO(300);
EDG_ADDIN_FUNC_MACRO(301);
EDG_ADDIN_FUNC_MACRO(302);
EDG_ADDIN_FUNC_MACRO(303);
EDG_ADDIN_FUNC_MACRO(304);
EDG_ADDIN_FUNC_MACRO(305);
EDG_ADDIN_FUNC_MACRO(306);
EDG_ADDIN_FUNC_MACRO(307);
EDG_ADDIN_FUNC_MACRO(308);
EDG_ADDIN_FUNC_MACRO(309);
EDG_ADDIN_FUNC_MACRO(310);
EDG_ADDIN_FUNC_MACRO(311);
EDG_ADDIN_FUNC_MACRO(312);
EDG_ADDIN_FUNC_MACRO(313);
EDG_ADDIN_FUNC_MACRO(314);
EDG_ADDIN_FUNC_MACRO(315);
EDG_ADDIN_FUNC_MACRO(316);
EDG_ADDIN_FUNC_MACRO(317);
EDG_ADDIN_FUNC_MACRO(318);
EDG_ADDIN_FUNC_MACRO(319);
EDG_ADDIN_FUNC_MACRO(320);
EDG_ADDIN_FUNC_MACRO(321);
EDG_ADDIN_FUNC_MACRO(322);
EDG_ADDIN_FUNC_MACRO(323);
EDG_ADDIN_FUNC_MACRO(324);
EDG_ADDIN_FUNC_MACRO(325);
EDG_ADDIN_FUNC_MACRO(326);
EDG_ADDIN_FUNC_MACRO(327);
EDG_ADDIN_FUNC_MACRO(328);
EDG_ADDIN_FUNC_MACRO(329);
EDG_ADDIN_FUNC_MACRO(330);
EDG_ADDIN_FUNC_MACRO(331);
EDG_ADDIN_FUNC_MACRO(332);
EDG_ADDIN_FUNC_MACRO(333);
EDG_ADDIN_FUNC_MACRO(334);
EDG_ADDIN_FUNC_MACRO(335);
EDG_ADDIN_FUNC_MACRO(336);
EDG_ADDIN_FUNC_MACRO(337);
EDG_ADDIN_FUNC_MACRO(338);
EDG_ADDIN_FUNC_MACRO(339);
EDG_ADDIN_FUNC_MACRO(340);
EDG_ADDIN_FUNC_MACRO(341);
EDG_ADDIN_FUNC_MACRO(342);
EDG_ADDIN_FUNC_MACRO(343);
EDG_ADDIN_FUNC_MACRO(344);
EDG_ADDIN_FUNC_MACRO(345);
EDG_ADDIN_FUNC_MACRO(346);
EDG_ADDIN_FUNC_MACRO(347);
EDG_ADDIN_FUNC_MACRO(348);
EDG_ADDIN_FUNC_MACRO(349);
EDG_ADDIN_FUNC_MACRO(350);
EDG_ADDIN_FUNC_MACRO(351);
EDG_ADDIN_FUNC_MACRO(352);
EDG_ADDIN_FUNC_MACRO(353);
EDG_ADDIN_FUNC_MACRO(354);
EDG_ADDIN_FUNC_MACRO(355);
EDG_ADDIN_FUNC_MACRO(356);
EDG_ADDIN_FUNC_MACRO(357);
EDG_ADDIN_FUNC_MACRO(358);
EDG_ADDIN_FUNC_MACRO(359);
EDG_ADDIN_FUNC_MACRO(360);
EDG_ADDIN_FUNC_MACRO(361);
EDG_ADDIN_FUNC_MACRO(362);
EDG_ADDIN_FUNC_MACRO(363);
EDG_ADDIN_FUNC_MACRO(364);
EDG_ADDIN_FUNC_MACRO(365);
EDG_ADDIN_FUNC_MACRO(366);
EDG_ADDIN_FUNC_MACRO(367);
EDG_ADDIN_FUNC_MACRO(368);
EDG_ADDIN_FUNC_MACRO(369);
EDG_ADDIN_FUNC_MACRO(370);
EDG_ADDIN_FUNC_MACRO(371);
EDG_ADDIN_FUNC_MACRO(372);
EDG_ADDIN_FUNC_MACRO(373);
EDG_ADDIN_FUNC_MACRO(374);
EDG_ADDIN_FUNC_MACRO(375);
EDG_ADDIN_FUNC_MACRO(376);
EDG_ADDIN_FUNC_MACRO(377);
EDG_ADDIN_FUNC_MACRO(378);
EDG_ADDIN_FUNC_MACRO(379);
EDG_ADDIN_FUNC_MACRO(380);
EDG_ADDIN_FUNC_MACRO(381);
EDG_ADDIN_FUNC_MACRO(382);
EDG_ADDIN_FUNC_MACRO(383);
EDG_ADDIN_FUNC_MACRO(384);
EDG_ADDIN_FUNC_MACRO(385);
EDG_ADDIN_FUNC_MACRO(386);
EDG_ADDIN_FUNC_MACRO(387);
EDG_ADDIN_FUNC_MACRO(388);
EDG_ADDIN_FUNC_MACRO(389);
EDG_ADDIN_FUNC_MACRO(390);
EDG_ADDIN_FUNC_MACRO(391);
EDG_ADDIN_FUNC_MACRO(392);
EDG_ADDIN_FUNC_MACRO(393);
EDG_ADDIN_FUNC_MACRO(394);
EDG_ADDIN_FUNC_MACRO(395);
EDG_ADDIN_FUNC_MACRO(396);
EDG_ADDIN_FUNC_MACRO(397);
EDG_ADDIN_FUNC_MACRO(398);
EDG_ADDIN_FUNC_MACRO(399);
EDG_ADDIN_FUNC_MACRO(400);
EDG_ADDIN_FUNC_MACRO(401);
EDG_ADDIN_FUNC_MACRO(402);
EDG_ADDIN_FUNC_MACRO(403);
EDG_ADDIN_FUNC_MACRO(404);
EDG_ADDIN_FUNC_MACRO(405);
EDG_ADDIN_FUNC_MACRO(406);
EDG_ADDIN_FUNC_MACRO(407);
EDG_ADDIN_FUNC_MACRO(408);
EDG_ADDIN_FUNC_MACRO(409);
EDG_ADDIN_FUNC_MACRO(410);
EDG_ADDIN_FUNC_MACRO(411);
EDG_ADDIN_FUNC_MACRO(412);
EDG_ADDIN_FUNC_MACRO(413);
EDG_ADDIN_FUNC_MACRO(414);
EDG_ADDIN_FUNC_MACRO(415);
EDG_ADDIN_FUNC_MACRO(416);
EDG_ADDIN_FUNC_MACRO(417);
EDG_ADDIN_FUNC_MACRO(418);
EDG_ADDIN_FUNC_MACRO(419);
EDG_ADDIN_FUNC_MACRO(420);
EDG_ADDIN_FUNC_MACRO(421);
EDG_ADDIN_FUNC_MACRO(422);
EDG_ADDIN_FUNC_MACRO(423);
EDG_ADDIN_FUNC_MACRO(424);
EDG_ADDIN_FUNC_MACRO(425);
EDG_ADDIN_FUNC_MACRO(426);
EDG_ADDIN_FUNC_MACRO(427);
EDG_ADDIN_FUNC_MACRO(428);
EDG_ADDIN_FUNC_MACRO(429);
EDG_ADDIN_FUNC_MACRO(430);
EDG_ADDIN_FUNC_MACRO(431);
EDG_ADDIN_FUNC_MACRO(432);
EDG_ADDIN_FUNC_MACRO(433);
EDG_ADDIN_FUNC_MACRO(434);
EDG_ADDIN_FUNC_MACRO(435);
EDG_ADDIN_FUNC_MACRO(436);
EDG_ADDIN_FUNC_MACRO(437);
EDG_ADDIN_FUNC_MACRO(438);
EDG_ADDIN_FUNC_MACRO(439);
EDG_ADDIN_FUNC_MACRO(440);
EDG_ADDIN_FUNC_MACRO(441);
EDG_ADDIN_FUNC_MACRO(442);
EDG_ADDIN_FUNC_MACRO(443);
EDG_ADDIN_FUNC_MACRO(444);
EDG_ADDIN_FUNC_MACRO(445);
EDG_ADDIN_FUNC_MACRO(446);
EDG_ADDIN_FUNC_MACRO(447);
EDG_ADDIN_FUNC_MACRO(448);
EDG_ADDIN_FUNC_MACRO(449);
EDG_ADDIN_FUNC_MACRO(450);
EDG_ADDIN_FUNC_MACRO(451);
EDG_ADDIN_FUNC_MACRO(452);
EDG_ADDIN_FUNC_MACRO(453);
EDG_ADDIN_FUNC_MACRO(454);
EDG_ADDIN_FUNC_MACRO(455);
EDG_ADDIN_FUNC_MACRO(456);
EDG_ADDIN_FUNC_MACRO(457);
EDG_ADDIN_FUNC_MACRO(458);
EDG_ADDIN_FUNC_MACRO(459);
EDG_ADDIN_FUNC_MACRO(460);
EDG_ADDIN_FUNC_MACRO(461);
EDG_ADDIN_FUNC_MACRO(462);
EDG_ADDIN_FUNC_MACRO(463);
EDG_ADDIN_FUNC_MACRO(464);
EDG_ADDIN_FUNC_MACRO(465);
EDG_ADDIN_FUNC_MACRO(466);
EDG_ADDIN_FUNC_MACRO(467);
EDG_ADDIN_FUNC_MACRO(468);
EDG_ADDIN_FUNC_MACRO(469);
EDG_ADDIN_FUNC_MACRO(470);
EDG_ADDIN_FUNC_MACRO(471);
EDG_ADDIN_FUNC_MACRO(472);
EDG_ADDIN_FUNC_MACRO(473);
EDG_ADDIN_FUNC_MACRO(474);
EDG_ADDIN_FUNC_MACRO(475);
EDG_ADDIN_FUNC_MACRO(476);
EDG_ADDIN_FUNC_MACRO(477);
EDG_ADDIN_FUNC_MACRO(478);
EDG_ADDIN_FUNC_MACRO(479);
EDG_ADDIN_FUNC_MACRO(480);
EDG_ADDIN_FUNC_MACRO(481);
EDG_ADDIN_FUNC_MACRO(482);
EDG_ADDIN_FUNC_MACRO(483);
EDG_ADDIN_FUNC_MACRO(484);
EDG_ADDIN_FUNC_MACRO(485);
EDG_ADDIN_FUNC_MACRO(486);
EDG_ADDIN_FUNC_MACRO(487);
EDG_ADDIN_FUNC_MACRO(488);
EDG_ADDIN_FUNC_MACRO(489);
EDG_ADDIN_FUNC_MACRO(490);
EDG_ADDIN_FUNC_MACRO(491);
EDG_ADDIN_FUNC_MACRO(492);
EDG_ADDIN_FUNC_MACRO(493);
EDG_ADDIN_FUNC_MACRO(494);
EDG_ADDIN_FUNC_MACRO(495);
EDG_ADDIN_FUNC_MACRO(496);
EDG_ADDIN_FUNC_MACRO(497);
EDG_ADDIN_FUNC_MACRO(498);
EDG_ADDIN_FUNC_MACRO(499);
EDG_ADDIN_FUNC_MACRO(500);
EDG_ADDIN_FUNC_MACRO(501);
EDG_ADDIN_FUNC_MACRO(502);
EDG_ADDIN_FUNC_MACRO(503);
EDG_ADDIN_FUNC_MACRO(504);
EDG_ADDIN_FUNC_MACRO(505);
EDG_ADDIN_FUNC_MACRO(506);
EDG_ADDIN_FUNC_MACRO(507);
EDG_ADDIN_FUNC_MACRO(508);
EDG_ADDIN_FUNC_MACRO(509);
EDG_ADDIN_FUNC_MACRO(510);
EDG_ADDIN_FUNC_MACRO(511);
EDG_ADDIN_FUNC_MACRO(512);
EDG_ADDIN_FUNC_MACRO(513);
EDG_ADDIN_FUNC_MACRO(514);
EDG_ADDIN_FUNC_MACRO(515);
EDG_ADDIN_FUNC_MACRO(516);
EDG_ADDIN_FUNC_MACRO(517);
EDG_ADDIN_FUNC_MACRO(518);
EDG_ADDIN_FUNC_MACRO(519);
EDG_ADDIN_FUNC_MACRO(520);
EDG_ADDIN_FUNC_MACRO(521);
EDG_ADDIN_FUNC_MACRO(522);
EDG_ADDIN_FUNC_MACRO(523);
EDG_ADDIN_FUNC_MACRO(524);
EDG_ADDIN_FUNC_MACRO(525);
EDG_ADDIN_FUNC_MACRO(526);
EDG_ADDIN_FUNC_MACRO(527);
EDG_ADDIN_FUNC_MACRO(528);
EDG_ADDIN_FUNC_MACRO(529);
EDG_ADDIN_FUNC_MACRO(530);
EDG_ADDIN_FUNC_MACRO(531);
EDG_ADDIN_FUNC_MACRO(532);
EDG_ADDIN_FUNC_MACRO(533);
EDG_ADDIN_FUNC_MACRO(534);
EDG_ADDIN_FUNC_MACRO(535);
EDG_ADDIN_FUNC_MACRO(536);
EDG_ADDIN_FUNC_MACRO(537);
EDG_ADDIN_FUNC_MACRO(538);
EDG_ADDIN_FUNC_MACRO(539);
EDG_ADDIN_FUNC_MACRO(540);
EDG_ADDIN_FUNC_MACRO(541);
EDG_ADDIN_FUNC_MACRO(542);
EDG_ADDIN_FUNC_MACRO(543);
EDG_ADDIN_FUNC_MACRO(544);
EDG_ADDIN_FUNC_MACRO(545);
EDG_ADDIN_FUNC_MACRO(546);
EDG_ADDIN_FUNC_MACRO(547);
EDG_ADDIN_FUNC_MACRO(548);
EDG_ADDIN_FUNC_MACRO(549);
EDG_ADDIN_FUNC_MACRO(550);
EDG_ADDIN_FUNC_MACRO(551);
EDG_ADDIN_FUNC_MACRO(552);
EDG_ADDIN_FUNC_MACRO(553);
EDG_ADDIN_FUNC_MACRO(554);
EDG_ADDIN_FUNC_MACRO(555);
EDG_ADDIN_FUNC_MACRO(556);
EDG_ADDIN_FUNC_MACRO(557);
EDG_ADDIN_FUNC_MACRO(558);
EDG_ADDIN_FUNC_MACRO(559);
EDG_ADDIN_FUNC_MACRO(560);
EDG_ADDIN_FUNC_MACRO(561);
EDG_ADDIN_FUNC_MACRO(562);
EDG_ADDIN_FUNC_MACRO(563);
EDG_ADDIN_FUNC_MACRO(564);
EDG_ADDIN_FUNC_MACRO(565);
EDG_ADDIN_FUNC_MACRO(566);
EDG_ADDIN_FUNC_MACRO(567);
EDG_ADDIN_FUNC_MACRO(568);
EDG_ADDIN_FUNC_MACRO(569);
EDG_ADDIN_FUNC_MACRO(570);
EDG_ADDIN_FUNC_MACRO(571);
EDG_ADDIN_FUNC_MACRO(572);
EDG_ADDIN_FUNC_MACRO(573);
EDG_ADDIN_FUNC_MACRO(574);
EDG_ADDIN_FUNC_MACRO(575);
EDG_ADDIN_FUNC_MACRO(576);
EDG_ADDIN_FUNC_MACRO(577);
EDG_ADDIN_FUNC_MACRO(578);
EDG_ADDIN_FUNC_MACRO(579);
EDG_ADDIN_FUNC_MACRO(580);
EDG_ADDIN_FUNC_MACRO(581);
EDG_ADDIN_FUNC_MACRO(582);
EDG_ADDIN_FUNC_MACRO(583);
EDG_ADDIN_FUNC_MACRO(584);
EDG_ADDIN_FUNC_MACRO(585);
EDG_ADDIN_FUNC_MACRO(586);
EDG_ADDIN_FUNC_MACRO(587);
EDG_ADDIN_FUNC_MACRO(588);
EDG_ADDIN_FUNC_MACRO(589);
EDG_ADDIN_FUNC_MACRO(590);
EDG_ADDIN_FUNC_MACRO(591);
EDG_ADDIN_FUNC_MACRO(592);
EDG_ADDIN_FUNC_MACRO(593);
EDG_ADDIN_FUNC_MACRO(594);
EDG_ADDIN_FUNC_MACRO(595);
EDG_ADDIN_FUNC_MACRO(596);
EDG_ADDIN_FUNC_MACRO(597);
EDG_ADDIN_FUNC_MACRO(598);
EDG_ADDIN_FUNC_MACRO(599);
EDG_ADDIN_FUNC_MACRO(600);
EDG_ADDIN_FUNC_MACRO(601);
EDG_ADDIN_FUNC_MACRO(602);
EDG_ADDIN_FUNC_MACRO(603);
EDG_ADDIN_FUNC_MACRO(604);
EDG_ADDIN_FUNC_MACRO(605);
EDG_ADDIN_FUNC_MACRO(606);
EDG_ADDIN_FUNC_MACRO(607);
EDG_ADDIN_FUNC_MACRO(608);
EDG_ADDIN_FUNC_MACRO(609);
EDG_ADDIN_FUNC_MACRO(610);
EDG_ADDIN_FUNC_MACRO(611);
EDG_ADDIN_FUNC_MACRO(612);
EDG_ADDIN_FUNC_MACRO(613);
EDG_ADDIN_FUNC_MACRO(614);
EDG_ADDIN_FUNC_MACRO(615);
EDG_ADDIN_FUNC_MACRO(616);
EDG_ADDIN_FUNC_MACRO(617);
EDG_ADDIN_FUNC_MACRO(618);
EDG_ADDIN_FUNC_MACRO(619);
EDG_ADDIN_FUNC_MACRO(620);
EDG_ADDIN_FUNC_MACRO(621);
EDG_ADDIN_FUNC_MACRO(622);
EDG_ADDIN_FUNC_MACRO(623);
EDG_ADDIN_FUNC_MACRO(624);
EDG_ADDIN_FUNC_MACRO(625);
EDG_ADDIN_FUNC_MACRO(626);
EDG_ADDIN_FUNC_MACRO(627);
EDG_ADDIN_FUNC_MACRO(628);
EDG_ADDIN_FUNC_MACRO(629);
EDG_ADDIN_FUNC_MACRO(630);
EDG_ADDIN_FUNC_MACRO(631);
EDG_ADDIN_FUNC_MACRO(632);
EDG_ADDIN_FUNC_MACRO(633);
EDG_ADDIN_FUNC_MACRO(634);
EDG_ADDIN_FUNC_MACRO(635);
EDG_ADDIN_FUNC_MACRO(636);
EDG_ADDIN_FUNC_MACRO(637);
EDG_ADDIN_FUNC_MACRO(638);
EDG_ADDIN_FUNC_MACRO(639);
EDG_ADDIN_FUNC_MACRO(640);
EDG_ADDIN_FUNC_MACRO(641);
EDG_ADDIN_FUNC_MACRO(642);
EDG_ADDIN_FUNC_MACRO(643);
EDG_ADDIN_FUNC_MACRO(644);
EDG_ADDIN_FUNC_MACRO(645);
EDG_ADDIN_FUNC_MACRO(646);
EDG_ADDIN_FUNC_MACRO(647);
EDG_ADDIN_FUNC_MACRO(648);
EDG_ADDIN_FUNC_MACRO(649);
EDG_ADDIN_FUNC_MACRO(650);
EDG_ADDIN_FUNC_MACRO(651);
EDG_ADDIN_FUNC_MACRO(652);
EDG_ADDIN_FUNC_MACRO(653);
EDG_ADDIN_FUNC_MACRO(654);
EDG_ADDIN_FUNC_MACRO(655);
EDG_ADDIN_FUNC_MACRO(656);
EDG_ADDIN_FUNC_MACRO(657);
EDG_ADDIN_FUNC_MACRO(658);
EDG_ADDIN_FUNC_MACRO(659);
EDG_ADDIN_FUNC_MACRO(660);
EDG_ADDIN_FUNC_MACRO(661);
EDG_ADDIN_FUNC_MACRO(662);
EDG_ADDIN_FUNC_MACRO(663);
EDG_ADDIN_FUNC_MACRO(664);
EDG_ADDIN_FUNC_MACRO(665);
EDG_ADDIN_FUNC_MACRO(666);
EDG_ADDIN_FUNC_MACRO(667);
EDG_ADDIN_FUNC_MACRO(668);
EDG_ADDIN_FUNC_MACRO(669);
EDG_ADDIN_FUNC_MACRO(670);
EDG_ADDIN_FUNC_MACRO(671);
EDG_ADDIN_FUNC_MACRO(672);
EDG_ADDIN_FUNC_MACRO(673);
EDG_ADDIN_FUNC_MACRO(674);
EDG_ADDIN_FUNC_MACRO(675);
EDG_ADDIN_FUNC_MACRO(676);
EDG_ADDIN_FUNC_MACRO(677);
EDG_ADDIN_FUNC_MACRO(678);
EDG_ADDIN_FUNC_MACRO(679);
EDG_ADDIN_FUNC_MACRO(680);
EDG_ADDIN_FUNC_MACRO(681);
EDG_ADDIN_FUNC_MACRO(682);
EDG_ADDIN_FUNC_MACRO(683);
EDG_ADDIN_FUNC_MACRO(684);
EDG_ADDIN_FUNC_MACRO(685);
EDG_ADDIN_FUNC_MACRO(686);
EDG_ADDIN_FUNC_MACRO(687);
EDG_ADDIN_FUNC_MACRO(688);
EDG_ADDIN_FUNC_MACRO(689);
EDG_ADDIN_FUNC_MACRO(690);
EDG_ADDIN_FUNC_MACRO(691);
EDG_ADDIN_FUNC_MACRO(692);
EDG_ADDIN_FUNC_MACRO(693);
EDG_ADDIN_FUNC_MACRO(694);
EDG_ADDIN_FUNC_MACRO(695);
EDG_ADDIN_FUNC_MACRO(696);
EDG_ADDIN_FUNC_MACRO(697);
EDG_ADDIN_FUNC_MACRO(698);
EDG_ADDIN_FUNC_MACRO(699);
EDG_ADDIN_FUNC_MACRO(700);
EDG_ADDIN_FUNC_MACRO(701);
EDG_ADDIN_FUNC_MACRO(702);
EDG_ADDIN_FUNC_MACRO(703);
EDG_ADDIN_FUNC_MACRO(704);
EDG_ADDIN_FUNC_MACRO(705);
EDG_ADDIN_FUNC_MACRO(706);
EDG_ADDIN_FUNC_MACRO(707);
EDG_ADDIN_FUNC_MACRO(708);
EDG_ADDIN_FUNC_MACRO(709);
EDG_ADDIN_FUNC_MACRO(710);
EDG_ADDIN_FUNC_MACRO(711);
EDG_ADDIN_FUNC_MACRO(712);
EDG_ADDIN_FUNC_MACRO(713);
EDG_ADDIN_FUNC_MACRO(714);
EDG_ADDIN_FUNC_MACRO(715);
EDG_ADDIN_FUNC_MACRO(716);
EDG_ADDIN_FUNC_MACRO(717);
EDG_ADDIN_FUNC_MACRO(718);
EDG_ADDIN_FUNC_MACRO(719);
EDG_ADDIN_FUNC_MACRO(720);
EDG_ADDIN_FUNC_MACRO(721);
EDG_ADDIN_FUNC_MACRO(722);
EDG_ADDIN_FUNC_MACRO(723);
EDG_ADDIN_FUNC_MACRO(724);
EDG_ADDIN_FUNC_MACRO(725);
EDG_ADDIN_FUNC_MACRO(726);
EDG_ADDIN_FUNC_MACRO(727);
EDG_ADDIN_FUNC_MACRO(728);
EDG_ADDIN_FUNC_MACRO(729);
EDG_ADDIN_FUNC_MACRO(730);
EDG_ADDIN_FUNC_MACRO(731);
EDG_ADDIN_FUNC_MACRO(732);
EDG_ADDIN_FUNC_MACRO(733);
EDG_ADDIN_FUNC_MACRO(734);
EDG_ADDIN_FUNC_MACRO(735);
EDG_ADDIN_FUNC_MACRO(736);
EDG_ADDIN_FUNC_MACRO(737);
EDG_ADDIN_FUNC_MACRO(738);
EDG_ADDIN_FUNC_MACRO(739);
EDG_ADDIN_FUNC_MACRO(740);
EDG_ADDIN_FUNC_MACRO(741);
EDG_ADDIN_FUNC_MACRO(742);
EDG_ADDIN_FUNC_MACRO(743);
EDG_ADDIN_FUNC_MACRO(744);
EDG_ADDIN_FUNC_MACRO(745);
EDG_ADDIN_FUNC_MACRO(746);
EDG_ADDIN_FUNC_MACRO(747);
EDG_ADDIN_FUNC_MACRO(748);
EDG_ADDIN_FUNC_MACRO(749);
EDG_ADDIN_FUNC_MACRO(750);
EDG_ADDIN_FUNC_MACRO(751);
EDG_ADDIN_FUNC_MACRO(752);
EDG_ADDIN_FUNC_MACRO(753);
EDG_ADDIN_FUNC_MACRO(754);
EDG_ADDIN_FUNC_MACRO(755);
EDG_ADDIN_FUNC_MACRO(756);
EDG_ADDIN_FUNC_MACRO(757);
EDG_ADDIN_FUNC_MACRO(758);
EDG_ADDIN_FUNC_MACRO(759);
EDG_ADDIN_FUNC_MACRO(760);
EDG_ADDIN_FUNC_MACRO(761);
EDG_ADDIN_FUNC_MACRO(762);
EDG_ADDIN_FUNC_MACRO(763);
EDG_ADDIN_FUNC_MACRO(764);
EDG_ADDIN_FUNC_MACRO(765);
EDG_ADDIN_FUNC_MACRO(766);
EDG_ADDIN_FUNC_MACRO(767);
EDG_ADDIN_FUNC_MACRO(768);
EDG_ADDIN_FUNC_MACRO(769);
EDG_ADDIN_FUNC_MACRO(770);
EDG_ADDIN_FUNC_MACRO(771);
EDG_ADDIN_FUNC_MACRO(772);
EDG_ADDIN_FUNC_MACRO(773);
EDG_ADDIN_FUNC_MACRO(774);
EDG_ADDIN_FUNC_MACRO(775);
EDG_ADDIN_FUNC_MACRO(776);
EDG_ADDIN_FUNC_MACRO(777);
EDG_ADDIN_FUNC_MACRO(778);
EDG_ADDIN_FUNC_MACRO(779);
EDG_ADDIN_FUNC_MACRO(780);
EDG_ADDIN_FUNC_MACRO(781);
EDG_ADDIN_FUNC_MACRO(782);
EDG_ADDIN_FUNC_MACRO(783);
EDG_ADDIN_FUNC_MACRO(784);
EDG_ADDIN_FUNC_MACRO(785);
EDG_ADDIN_FUNC_MACRO(786);
EDG_ADDIN_FUNC_MACRO(787);
EDG_ADDIN_FUNC_MACRO(788);
EDG_ADDIN_FUNC_MACRO(789);
EDG_ADDIN_FUNC_MACRO(790);
EDG_ADDIN_FUNC_MACRO(791);
EDG_ADDIN_FUNC_MACRO(792);
EDG_ADDIN_FUNC_MACRO(793);
EDG_ADDIN_FUNC_MACRO(794);
EDG_ADDIN_FUNC_MACRO(795);
EDG_ADDIN_FUNC_MACRO(796);
EDG_ADDIN_FUNC_MACRO(797);
EDG_ADDIN_FUNC_MACRO(798);
EDG_ADDIN_FUNC_MACRO(799);
}


/**********************************************************************
 * from Microsoft Excel Developer's Kit - put into this file to minimise
 to just one file for building the xll from the library 
**********************************************************************/

/*
**  Microsoft Excel Developer's Kit
**
**  File:           SAMPLE\FRAMEWRK\FRAMEWRK.C
**  Description:    Framework library for Microsoft Excel
**  Platform:       Microsoft Windows
**
**  This library provides some basic functions
**  that help you write Excel DLLs. It includes
**  simple functions for managing memory with XLOPERs,
**  creating temporary XLOPERs, robustly calling
**  Excel4(), and printing debugging strings on
**  a terminal attached to COM1.
**
**  The main purpose of this library is to help
**  you to write cleaner C code for calling Excel.
**  For example, using the framework library you
**  can write
**
**      Excel(xlcDisplay, 0, 2, TempMissing(), TempBool(0));
**
**  instead of the more verbose
**
**      XLOPER xMissing, xBool;
**      xMissing.xltype = xltypeMissing;
**      xBool.xltype = xltypeBool;
**      xBool.val.bool = 0;
**      Excel4(xlcDisplay, 0, 2, (LPXLOPER) &xMissing, (LPXLOPER) &xBool);
**
**
**  The library is non-reentrant: it assumes that
**  if there are multiple copies of Excel using the
**  DLL, they will not be preempted. This is
**  acceptable under Windows/DOS, but may present
**  some problems under Windows/NT. In particular,
**  the function Excel() frees all the temporary
**  memory, which means that this DLL must not be
**  reentered while temporary memory is in use.
**
**  Define DEBUG to use the debugging functions.
**
**  Source code is provided so that you may
**  enhance this library or optimize it for your
**  own application.
**
*/

/* avoid explicit linking to excel dll by using LoadLibrary to open dll and
   from there get the function address of Excel4v. The reason for this is
   to be able to use the xll as a regular dll for other apps (eg java). */
typedef int (pascal MS_Excel4v)(
    int xlfn, LPXLOPER operRes, int count, LPXLOPER far opers[]);
static MS_Excel4v *excel4vFuncPtr = NULL;
#define EXCEL_DLL "xlcall32.dll"
#define EXCEL_DLL_FUNC "Excel4v"

static void GetExcel4vFuncPtr()
{
    static const string routine = "GetExcel4vFuncPtr";
    HINSTANCE dllHandle = LoadLibrary(EXCEL_DLL);
    if (!dllHandle) {
        throw ModelException(routine, "Failed to load " EXCEL_DLL +
                             Format::toString((int)(GetLastError())));
    }
    if (!(excel4vFuncPtr = (MS_Excel4v *)
          GetProcAddress(dllHandle, EXCEL_DLL_FUNC))) {
        throw ModelException(routine, "Failed to lookup " EXCEL_DLL_FUNC+
                             Format::toString((int)(GetLastError())));
    }
}

/*
** Total amount of memory to allocate for all temporary XLOPERs
*/

#define MEMORYSIZE 1024
#define PRIVATE static

PRIVATE LPSTR GetTempMemory(int cBytes);
PRIVATE void FreeAllTempMemory(void);

/*
** Globals (see the comment about reentrancy 
** earlier in this file)
*/

PRIVATE char vMemBlock[MEMORYSIZE]; /* Memory for temporary XLOPERs */
PRIVATE int vOffsetMemBlock=0;    /* Offset of next memory block to allocate */


/*
** GetTempMemory
**
** Allocates temporary memory. Temporary memory
** can only be freed in one chunk, by calling
** FreeAllTempMemory(). This is done by Excel().
**
** Arguments:
**
**      int cBytes      How many bytes to allocate
**
** Returns:
**
**      LPSTR           A pointer to the allocated memory,
**                      or 0 if more memory cannot be
**                      allocated. If this fails,
**                      check that you are initializing
**                      vOffsetMemBlock to 0, and check that
**                      MEMORYSIZE is big enough.
**
** Algorithm:
**
**      The memory allocation algorithm is extremely
**      simple: on each call, allocate the next cBytes
**      bytes of a static memory buffer. If the buffer
**      becomes too full, simply fail. To free memory,
**      simply reset the pointer (vOffsetMemBlock)
**      back to zero. This memory scheme is very fast
**      and is optimized for the assumption that the
**      only thing you are using temporary memory
**      for is to hold arguments while you call Excel().
**      We rely on the fact that you will free all the
**      temporary memory at the same time. We also
**      assume you will not need more memory than
**      the amount required to hold a few arguments
**      to Excel().
*/

PRIVATE LPSTR GetTempMemory(int cBytes)
{
    LPSTR lpMemory;

    if (vOffsetMemBlock + cBytes > MEMORYSIZE)
    {
        return 0;
    }
    else
    {
        lpMemory = (LPSTR) &vMemBlock + vOffsetMemBlock;
        vOffsetMemBlock += cBytes;

        /* Prevent odd pointers */
        if (vOffsetMemBlock & 1) vOffsetMemBlock++;
        return lpMemory;
    }
}


/*
** FreeAllTempMemory
**
** Frees all temporary memory that has been allocated.
**
** Arguments:
**
**      None.
**
** Return value:
**
**      None.
*/

PRIVATE void FreeAllTempMemory(void)
{
    vOffsetMemBlock = 0;
}


/*
** Excel
**
** A fancy wrapper for the Excel4() function. It also
** does the following:
**
**  (1) Checks that none of the LPXLOPER arguments are 0,
**      which would indicate that creating a temporary XLOPER
**      has failed. In this case, it doesn't call Excel
**      but it does print a debug message.
**  (2) If an error occurs while calling Excel,
**      print a useful debug message.
**  (3) When done, free all temporary memory.
**
**  #1 and #2 require DEBUG to be defined.
**
** Arguments (same as Excel4()):
**
**      int xlfn            Function number (xl...) to call
**      LPXLOPER pxResult   Pointer to a place to stuff the result,
**                          or 0 if you don't care about the result.
**      int count           Number of arguments
**      ...                 (all LPXLOPERs) - the arguments.
**
** Return value:
**
**      A return code (Some of the xlret... values, as defined
**      in XLCALL.H, OR'ed together).
**
** Note:
**
**      Be sure to cast all the arguments after the third
**      to LPXLOPERs. If you accidentally pass a near pointer
**      instead of a far pointer, you will probably crash Excel.
*/
PRIVATE int Excel(int xlfn, LPXLOPER pxResult, int count, ...)
{
    int xlret;
    LPXLOPER xArray[30]; /* Cannot have more than 30 arguments */
    int maxCount = sizeof(xArray) / sizeof(LPXLOPER);
    int pos;
    va_list ap;

/*
** In the published version of this function, this walking through the variable
** argument list was not done "properly". This is an attempt to do so...
*/
    va_start (ap, count);
    if (count > maxCount) count = maxCount;
    for (pos = 0; pos < count; ++pos)
    {
        xArray[pos] = va_arg(ap, LPXLOPER);
    }
    va_end (ap);

    xlret = excel4vFuncPtr?
        excel4vFuncPtr(xlfn,pxResult,count,xArray): xlretFailed;

    FreeAllTempMemory();

    return xlret;
}



/*
** TempNum
**
** Creates a temporary numeric (IEEE floating point) XLOPER.
**
** Arguments:
**
**      double d        The value
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/

LPXLOPER TempNum(double d)
{
    LPXLOPER lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));
    if (lpx)
    {
        lpx->type = xltypeNum;
        lpx->val.num = d;
    }
    return lpx;
}


/*
** TempStr
**
** Creates a temporary string XLOPER.
**
** Arguments:
**
**      LPSTR lpstr     The string, as a null-terminated
**                      C string, with the first byte
**                      undefined. This function will
**                      count the bytes of the string
**                      and insert that count in the
**                      first byte of lpstr. Excel cannot
**                      handle strings longer than 255
**                      characters.
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
** Notes:
**
**      (1) This function has the side effect of inserting
**          the byte count as the first character of
**          the created string.
**
**      (2) For highest speed, with constant strings,
**          you may want to manually count the length of
**          the string before compiling, and then avoid
**          using this function.
**
**      (3) Behavior is undefined for non-null terminated
**          input or strings longer than 255 characters.
**
*/

LPXLOPER TempStr(const char* constLpstr)
{
    LPXLOPER lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (lpx)
    {
        int size = lstrlen (constLpstr);
        LPSTR lpstr = (LPSTR) GetTempMemory(size+1);
        strcpy(lpstr, constLpstr);
        lpstr[0] = (BYTE) (size-1);
        lpx->type = xltypeStr;
        lpx->val.str = lpstr;
    }
    return lpx;
}

/*
** TempBool
**
** Creates a temporary logical (true/false) XLOPER.
**
** Arguments:
**
**      int b           0 - for a FALSE XLOPER
**                      Anything else - for a TRUE XLOPER
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/

LPXLOPER TempBool(int b)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->type = xltypeBool;
    lpx->val.boolVal = b?1:0;

    return lpx;
}

/*
** TempInt
**
** Creates a temporary integer XLOPER.
**
** Arguments:
**
**      short int i          The integer
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/

LPXLOPER TempInt(short int i)
{
    LPXLOPER lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (lpx)
    {
        lpx->type = xltypeInt;
        lpx->val.w = i;
    }
    return lpx;
}

#if 0
/*
** TempErr
**
** Creates a temporary error XLOPER.
**
** Arguments:
**
**      WORD err        The error code. One of the xlerr...
**                      constants, as defined in XLCALL.H.
**                      See the Excel user manual for
**                      descriptions about the interpretation
**                      of various error codes.
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/

LPXLOPER TempErr(WORD err)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeErr;
    lpx->val.err = err;

    return lpx;
}



/*
** TempMissing
**
** This is used to simulate a missing argument when
** calling Excel(). It creates a temporary
** "missing" XLOPER.
**
** Arguments:
**
**      none.
**
** Returns:
**
**      LPXLOPER        The temporary XLOPER, or 0
**                      if GetTempMemory() failed.
**
*/

LPXLOPER TempMissing(void)
{
    LPXLOPER lpx;

    lpx = (LPXLOPER) GetTempMemory(sizeof(XLOPER));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeMissing;

    return lpx;
}
#endif


/*
** InitFramework
**
** Initializes all the framework functions.
**
** Arguments:
**
**      None.
**
** Return value:
**
**      None.
*/

void InitFramework(void)
{
    vOffsetMemBlock = 0;
}


DRLIB_END_NAMESPACE
#endif
