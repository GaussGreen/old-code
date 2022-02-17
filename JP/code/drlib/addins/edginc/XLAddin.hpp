//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLAddin.hpp
//
//   Description : Class for representing addins
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_XLADDIN_HPP
#define EDG_XLADDIN_HPP

#include "edginc/XLConvert.hpp"
#include "edginc/Handle.hpp"
#include <stdarg.h>

DRLIB_BEGIN_NAMESPACE

/** Class for representing addins. Each instance
    of this class corresponds to one excel addin function */
class ADDINS_DLL XLAddin{
public:
    /** global prefix for Excel */
    static const string ADDIN_PREFIX;
    static const string ADDIN_PREFIX_ALT; /* alternative prefix for backwards
                                             compatibility */

    /** Builds an addin - passes all the data needed to handle an addin
        function */
    static const XLAddin* create(
        const string&       addinName,   // excel name
        CClassConstSP       addinClass,  // identifies class representing addin
        const CFieldArray&  params,      // describes each parameter
        bool                handleNameParam, /* true: extra input for handle
                                                name */
        Addin::ReturnStyle  returnStyle, // how to return the output
        bool                returnIsNative, /* true - parameter is returned
                                               as an atomic type */
        Addin::MethodType   methodType,  /* constructor/instanceMethod/
                                            classMethod */
        Addin::Method       method,      // holds function pointer
        CClassConstSP       returnType); // the type of the return parameter
    
    /* implements default register method needed for regression tester */
    static void defaultRegister(
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
        CClassConstSP      returnType);  // the type of the return parameter

    /** execute the addin function using the given parameters as arguments */
    XL_OPER* execute(const XL_OPER*   firstArg,
                     va_list       args) const;

    /** execute the addin function using the given parameters as arguments.
        Additionally, if inpFileName is non null, creates regression input file
        together with output file if outFileName is non null */
    XL_OPER* execute(const string*   inpFileName, // to create .inp file
                     const string*   outFileName, // to create .out file
                     const XL_OPER*     firstArg,
                     va_list         args) const;

    /** Sets the method used for finding the current cell in the spreadsheet.
        Method can be null (default), in which case "NULL" is used */
    static void setSpreadSheetCellNameMethod(
        Handle::SpreadSheetCellNameMethod* sSheetMethod);

    /** Returns the XLAddin* instance for the addin with the current name.
        Will try given name, and then again having removed any ADDIN_PREFIX
        from addinName */
    static const XLAddin* lookUp(const string& addinName);

    /** Strips the ADDIN_PREFIX, if any, from a string */
    static string stripAddinPrefix(const string& addinName);

    /** Returns the number of XL_OPERs this addin expects */
    int numParams() const;

    /** Return the number of optional parameters this addin has */
    int numOptionalParams()const;

    /* method which determines whether we are being called from VB or not */
    typedef bool (TCallerIsVB)(void);
 
    /* register function which determines whether we are being called from VB
       or not */
    static void setCallerIsVBMethod(TCallerIsVB *method);

    /* uses registered method to determine whether we are being called
       from VB */
    static bool callerIsVB();

private:
    /** Builds an addin - passes all the data needed to handle an addin
        function */
    XLAddin(
        const string&       addinName,   // excel name
        CClassConstSP       addinClass,  // identifies class representing addin
        const CFieldArray&  params,      // describes each parameter
        bool                handleNameParam, /* true: extra input for handle
                                                name */
        Addin::ReturnStyle  returnStyle, // how to return the output
        bool                returnIsNative, /* true - parameter is returned
                                               as an atomic type */
        Addin::MethodType   methodType,  /* constructor/instanceMethod/
                                            classMethod */
        Addin::Method       method,      // holds function pointer
        CClassConstSP       returnType); // the type of the return parameter

    XLAddin(const XLAddin &rhs);
    XLAddin& operator=(const XLAddin& rhs);

    /** 'Expand' the given object in a single column */
    void expand(jobject&         jOutput,
                const string&    handleName,
                XL_OPER&            oper) const;

    void expandArrayAcrossColumns(IArray&             arr,
                                  bool                canSimplify,
                                  const CStringArray& names,
                                  const string&       handleName,
                                  XL_OPER&               output) const;

    void expandArrayDownRows(IArray&             arr,
                             const CStringArray& names,
                             const string&       handleName,
                             XL_OPER&               output) const;

    static TCallerIsVB*   callerIsVBMethod;
    static Handle::SpreadSheetCellNameMethod* sheetNameMethod;
    static const XLConvert* handleConvert;   // how to convert handle
    static const XLConvert* objectConvert;   // how to convert generic objects
    static const XLConvert* objArrayConvert; // how to convert generic arrays

    /** Write useful message to error log */
    void reportError(
        ModelException& e, 
        const string&   handleName, // if returning a handle
        XL_OPER*           output) const;

    string             addinName;    // excel name
    CClassConstSP      addinClass;   // identifies class representing addin
    bool               handleNameParam; // true: extra input for handle name
    Addin::ReturnStyle returnStyle;  // how to return the output
    Addin::MethodType  methodType;   // constructor/instanceMethod/classMethod
    Addin::Method      method;       // holds function pointer
    CFieldArray        params;       // the type of each parameter
    XLConvertArray     paramConvert; // how to convert each parameter
    CClassConstSP      returnType;   // the type of the return parameter
    bool               returnIsNative; /* true - method returns a 
                                          native type */
    const XLConvert*   returnConvert;// how to convert return parameter
};


DRLIB_END_NAMESPACE
#endif
