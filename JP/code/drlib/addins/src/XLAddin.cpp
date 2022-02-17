//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLAddin.cpp
//
//   Description : Class for representing addins for XL
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLAddin.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/AddinTest.hpp"
#include "edginc/Format.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/Library.hpp"
#include ext_hash_map
#include "edginc/SystemException.hpp"       // this includes <windows.h>

DRLIB_BEGIN_NAMESPACE


/** Hash function needed for strings */
struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};


typedef hash_map <string, XLAddin*, MyStringHash> XLAddinHash;

/** Holder for hash table - avoid slowing header file */
class XLAddinHelper{
public:
    static XLAddinHash addinHash;
};

XLAddinHash XLAddinHelper::addinHash;


const XLConvert* XLAddin::handleConvert = 0; // how to convert handle
const XLConvert* XLAddin::objectConvert = 0; // how to convert generic objects
const XLConvert* XLAddin::objArrayConvert= 0; // how to convert generic arrays
Handle::SpreadSheetCellNameMethod* XLAddin::sheetNameMethod = 0;
XLAddin::TCallerIsVB*   XLAddin::callerIsVBMethod = 0;

const string XLAddin::ADDIN_PREFIX(Library::XL_PREFIX);
const string XLAddin::ADDIN_PREFIX_ALT("EDR_");

/** Sets the method used for finding the current cell in the spreadsheet.
    Method can be null (default), in which case "???" is used for error
    reporting purposes */
void XLAddin::setSpreadSheetCellNameMethod(
    Handle::SpreadSheetCellNameMethod* sSheetMethod){
    sheetNameMethod = sSheetMethod;
}

/** Builds an XLAddin - passes all the data needed to handle an XLAddin
    function */
const XLAddin* XLAddin::create(
    const string&       addinName,   // excel name
    CClassConstSP       addinClass,  // identifies class representing addin
    const CFieldArray&  params,      // describes each parameter
    bool                handleNameParam, // true: extra input for handle name
    Addin::ReturnStyle  returnStyle, // how to return the output
    bool                returnIsNative, /* true - parameter is returned
                                           as an atomic type */
    Addin::MethodType   methodType,  /* constructor/instanceMethod/
                                      classMethod */
    Addin::Method       method,      // holds function pointer
    CClassConstSP       returnType) // the type of the return parameter
{
    XLAddin* addin = new XLAddin(addinName, addinClass, params,
                                 handleNameParam,
                                 returnStyle,returnIsNative, methodType,
                                 method, returnType);
    // add to hash table
    XLAddinHelper::addinHash[addinName] = addin;
    return addin;
}

/** Returns the XLAddin* instance for the addin with the current name.
    Will try given name, and then again having removed any ADDIN_PREFIX
    from addinName */
const XLAddin* XLAddin::lookUp(const string& addinName){
    XLAddinHash::const_iterator iter =
        XLAddinHelper::addinHash.find(addinName);
    if (iter == XLAddinHelper::addinHash.end()){
        // no match, see if string starts with ADDIN_PREFIX
        // note: string::compare not standard on gcc
        iter = XLAddinHelper::addinHash.find(stripAddinPrefix(addinName));
        if (iter == XLAddinHelper::addinHash.end()){
            throw ModelException("XLAddin::lookUp", "Unknown addin function "+
                                 addinName);
        }
    }
    return iter->second;
}

/** Builds an XLAddin - passes all the data needed to handle an XLAddin
    function */
XLAddin::XLAddin(
    const string&       addinName,   // excel name
    CClassConstSP       addinClass,  // identifies class representing addin
    const CFieldArray&  params,      // describes each parameter
    bool                handleNameParam, // true: extra input for handle name
    Addin::ReturnStyle  returnStyle, // how to return the output
    bool                returnIsNative, /* true - parameter is returned
                                           as an atomic type */
    Addin::MethodType   methodType,  /* constructor/instanceMethod/
                                      classMethod */
    Addin::Method       method,      // holds function pointer
    CClassConstSP       returnType): // the type of the return parameter
    addinName(addinName), addinClass(addinClass),
    handleNameParam(handleNameParam),
    returnStyle(returnStyle), methodType(methodType),
    method(method),
    params(params),
    returnType(returnType), returnIsNative(returnIsNative)
{
    if (!handleConvert){
        // initialise static member;
        handleConvert = XLConvertFactory::create(CString::TYPE);
        objectConvert = XLConvertFactory::create(IObject::TYPE);
        objArrayConvert = XLConvertFactory::create(ObjectArray::TYPE);
    }

    // populate paramConvert, returnConvert
    paramConvert = XLConvertArray(params.size());
    for (unsigned int idx = 0; idx < params.size(); idx++){
        paramConvert[idx] = XLConvertFactory::create(params[idx]->getType());
    }
    returnConvert = XLConvertFactory::create(returnType);
}

/** Write useful message to error log */
void XLAddin::reportError(
    ModelException& e,
    const string&   handleName, // if returning a handle
    XL_OPER*        output) const{
    if (returnStyle == Addin::returnHandle){
        e.addMsg("Failed whilst trying to create handle "+handleName);
    }
    string cellName = sheetNameMethod? sheetNameMethod(): "???";
    e.addMsg(ADDIN_PREFIX+addinName+ " failed at "+cellName+"\n");
    e.errorLog();
    if (output){
        XLConvert::setToError(*output);
    }
}

/** execute the addin function using the given parameters as arguments */
XL_OPER* XLAddin::execute(const XL_OPER  *firstArg,
                          va_list         args) const
{
    return execute(0, 0, firstArg, args);
}

/** generate handle name to use */
static void generateHandleName(string&        handleName,
                               bool           handleNameParam,
                               const jobject& paramObj){
    if (!handleNameParam || handleName.empty()){
        Handle::IName* hName = dynamic_cast<Handle::IName*>(paramObj.get());
        if (hName){
            handleName = hName->defaultHandleName();
        }
    }
    if (handleName.empty()){
        handleName = "NULL";
    }
}

/** execute the addin function using the given parameters as arguments.
    Additionally, if inpFileName is non null, creates regression input file
    together with output file if outFileName is non null */
XL_OPER* XLAddin::execute(const string*    inpFileName, // to create .inp file
                          const string*    outFileName, // to create .out file
                          const XL_OPER*   firstArg,
                          va_list          args) const{
    static const string  routine  = "XLAddin::execute";
    XL_OPER             *output  = 0;
    string               handleName; // may be used in error message
    XL_OPER              coerceOutput;
    bool                 useCoerceOutput = false;
    static int           debugCount = 0;

#if defined(WIN32)
    Win32StructuredException::HookGuard hookGuard;  // win32 structured exception handler
#endif

    try{
        bool           skip   = false; /* default */
        jvalue         jOutput;

        // check type of first oper
        if (handleNameParam && (firstArg->type == xltypeRef ||
                                firstArg->type == xltypeSRef)){
            if (XLConvert::coerceToMulti(*firstArg,
                                         coerceOutput)==xlretUncalced){
                return 0; // what MS tell us to do
            }
            useCoerceOutput = true;
        }
        /* create oper to hold output - also used for error codes on failure */
        output  = NEW(XL_OPER);

        if (handleNameParam){
            bool useDefaultHandleName = true;
            try{
                const XL_OPER& operToUse = XLConvert::coerceOperFromMulti(
                    useCoerceOutput?
                    (const XL_OPER&)coerceOutput: *firstArg);
                if (handleConvert->
                    validateInput(operToUse,
                                  CString::TYPE,
                                  1,
                                  *output,
                                  useDefaultHandleName)){
                    if (!useDefaultHandleName){
                        handleName = XLConvert::operToString(operToUse);
                    }
                } else {
                    skip = true; // terminate calculation
                }
            } catch (exception& e){
                throw ModelException(e, routine," Failed whilst trying to read"
                                     " handle name from first parameter");
            }
            if (useCoerceOutput){
                XLConvert::freeCoerce(coerceOutput);
                useCoerceOutput = false;
            }
        }

        /* create empty object to fill with with addin parameters ie we convert
           each of the XL_OPERs into our equivalent types and populate this
           structure with them */
        jobject paramObj(addinClass->newInstance());

        for (unsigned idx = 0; idx < params.size() && !skip; idx++) {
            CFieldConstSP field = params[idx];
            bool useDefault = field->isOptional();
            /* get oper */
            const XL_OPER *oper = (idx == 0 && !handleNameParam)?
                firstArg: va_arg(args, const XL_OPER *);
            if ((oper->type & xltypeRef) || (oper->type & xltypeSRef)){
                if (XLConvert::coerceToMulti(*oper,
                                             coerceOutput) == xlretUncalced){
                    debugCount++;
                    XLConvert::operFreeSimple(*output);
                    FREE(output);
                    // avoid calling xlAutoFree(output); - EAS link issues
                    return 0; // what MS tell us to do
                }
                useCoerceOutput = true;
            }
            const XL_OPER& operToUse = XLConvert::coerceOperFromMulti(
                useCoerceOutput? coerceOutput: *oper);
            // get hold of the required type
            CClassConstSP requiredType = field->getType();
            // validate the oper
            try{
                skip = !(paramConvert[idx]->
                         validateInput(operToUse,
                                       requiredType,
                                       idx + 1 + (handleNameParam? 1: 0),
                                       *output,
                                       useDefault));


                /* and then convert it (if we're not to use the default) */
                if (!skip && !useDefault){
                    paramConvert[idx]->convertGenericInput(requiredType,
                                                           operToUse,
                                                           params[idx],
                                                           paramObj);
                }
            } catch (exception& e){
                string m = "Failed to read '"+field->getDescription()+"'\n"
                    "(parameter name '"+field->getName()+"' "+
                    "parameter number "+
                    Format::toString((int)(idx+1+(handleNameParam?1:0)))+
                    ", type "+ field->getType()->getName()+")";
                throw ModelException(e, routine, m);
            }
            if (useCoerceOutput){
                XLConvert::freeCoerce(coerceOutput);
                useCoerceOutput = false;
            }
        }

        if (!skip){
            std::vector<std::string>  exceptionStack;
            ErrorHandler* origErrorHandler = 0;
            // write out regression file - note do this before any validation
            // - this allows us to create cases to test for failure
            if (inpFileName){
                // if running as an addin, prevent overwriting existing files
                if (SpreadSheetMode::isOn()) {
                    FILE *fp = fopen(inpFileName->c_str(), "r");
                    if (fp) {
                        fclose(fp);
                        throw ModelException(routine,
                                             "file (" + *inpFileName +
                                             ") already exists");
                    }
                }

                XMLWriter    xml(*inpFileName);
                AddinTest    addinTest(addinName, paramObj);

                addinTest.write("ADDIN-TEST", &xml);
                origErrorHandler = new ErrorHandler();
                *origErrorHandler =
                    ErrorHandler::set(ErrorHandler(exceptionStack), false);
            }
            try{
                // invoke any validation needed for this object
                paramObj->validatePop2Object();
                // then clone in order to protect the inputs from being
                // changed as this can really screw up a spreadsheet.
                // DO NOT REMOVE THIS
                paramObj = IObjectSP(paramObj->clone());                
                /* now execute method */
                returnConvert->invoke(method, paramObj, jOutput);
            } catch (exception& e){
                if (outFileName){
                    ModelException modelEx(e);
                    modelEx.errorLog();
                    // restore error handler
                    ErrorHandler::set(*origErrorHandler, true);
                    delete origErrorHandler;
                    OutputFile fileOut(*outFileName);
                    fileOut.write(0, CStringArray(exceptionStack.begin(), exceptionStack.end()));
                }
                throw;
            }
            if (outFileName){
                // restore error handler
                ErrorHandler::set(*origErrorHandler, true);
                delete origErrorHandler;
                // record output to file
                if (returnIsNative){
                    // need to wrap any native types
                    returnConvert->nativeToObject(jOutput);
                }
                OutputFile fileOut(*outFileName);
                fileOut.write(jOutput.l.get());
            }
            // sort out handle name to pass down
            generateHandleName(handleName, handleNameParam, paramObj);
            // and then switch depending on how the data is to be returned
            switch (returnStyle){
            case Addin::returnHandle:
                if (returnIsNative){
                    // need to wrap any native types
                    returnConvert->nativeToObject(jOutput);
                }
                XLConvert::createHandle(handleName, "", // additional name
                                        jOutput.l,
                                        0,  // one and only parameter
                                        -1, // don't append to userName
                                        *output);
                break;
            case Addin::expandSimple:
            case Addin::expandMulti:
                if (returnIsNative){
                    int width = returnConvert->width();
                    if (width > 1){ // set num of cols if need xltypeMulti
                        XLConvert::setNumColsInOper(width, *output);
                    }
                    returnConvert->convertNativeOutput(
                        handleNameParam?
                        handleName: addinClass->getName(),
                        jOutput,
                        0, // column index
                        *output);
                } else {
                    if (!jOutput.l){
                        throw ModelException(routine, "Addin returned null!");
                    }
                    expand(jOutput.l, handleName, *output);
                    /* if called from the spreadsheet (cf from VB)
                       then ensure that the number of rows returned is
                       at least 2 (otherwise excel repeats the same
                       number in all cells) */
                    if ((output->type & xltypeMulti) &&
                        (output->val.xlarray.columns < 1 ||
                         output->val.xlarray.rows < 2) &&
                        !callerIsVB()){
                        XLConvert::ensureOperIsMinimumSize(*output);
                    } else if ((output->type & xltypeMulti) &&
                               output->val.xlarray.columns *
                               output->val.xlarray.rows == 0){
                        if ( !callerIsVB()) {
                            // XL doesn't like an XL_OPER like this (it cores)
                            output->val.xlarray.lparray = NEW_ARRAY(XL_OPER,1);
                            output->val.xlarray.lparray[0].type = xltypeErr;
                            output->val.xlarray.lparray[0].val.err = xlerrNA;
                            output->val.xlarray.columns = output->val.xlarray.rows = 1;
                        } else {
                            output->type = xltypeErr;
                            output->val.err = xlerrNA;
                        }
                    }
                }
                break;
        default:
            throw ModelException(routine, "Unknown return style");
            }
        }
    }
    catch (exception& e)
    {
        if (useCoerceOutput){
            XLConvert::freeCoerce(coerceOutput);
        }
        ModelException e2(e, routine);
        reportError(e2, handleName, output);
    }
    catch (...)     // Don't let Excel crash out completely!
    {
        if (callerIsVB())
            throw;      // Let VB handle it 
        if (useCoerceOutput)
            XLConvert::freeCoerce(coerceOutput);
        ModelException e2(routine, "FATAL ERROR - please contact QLib Support Team!");
        reportError(e2, handleName, output);
    }

#if 0
    /* sort out calling from vb bit */
#endif

    if (output)
    {
        output->type |= xlbitDLLFree;
    }
    return output;
}


/** 'Expand' the given object in a single column */
void XLAddin::expand(jobject&         obj,
                     const string&    handleName,
                     XL_OPER&         output) const {
    try{
        // ensure we show an object's public representation
        obj = CObject::convertToPublicRep(obj);
        // look up how to convert this type of object
        const XLConvert* convertMethod =
            XLConvertFactory::create(obj->getClass());

        // use provided its not the default object/default array method
        if (convertMethod == objArrayConvert){
            // this block is for handling generic arrays
            IArray& array = *DYNAMIC_CAST(IArray, obj.get());
            CStringArray names(array.getLength());
            // need array of names for each component
            for (int i = 0; i < array.getLength(); i++){
                names[i] = "List"+Format::toString(i);
            }
            if (array.getLength() == 0){
                // set to multitype with numRows = 0
                XLConvert::setNumColsInOper(0, output);
                XLConvert::setNumberOfRowsInOper(0, 0, output);
            } else if (returnStyle == Addin::expandMulti){
                // here we expand the array horizontally across columns
                expandArrayAcrossColumns(array, true, names,
                                         handleName, output);
            } else {
                expandArrayDownRows(array, names, handleName, output);
           }
        } else if (convertMethod != objectConvert){
            /* this block is for handling objects (or arrays of
               objects) with specialised methods */
            int width = convertMethod->width();
            if (width > 1){
                // set cols if we need a multi type oper
                XLConvert::setNumColsInOper(width, output);
            }
            convertMethod->convertObjectOutput(handleName,
                                               "",
                                               obj,
                                               0, // column index
                                               output);
        } else {
            // this block is for handling generic objects which aren't arrays
            // get the list of fields for this object (includes parents)
            CFieldArray fields = Addin::getDataClassFields(obj->getClass());
            // turn into ObjectArray and use expandArrayAcrossColumns
            // or expandArrayDownRows as appropriate
            ObjectArray  objArray(fields.size());

            CStringArray names(fields.size());
            for (unsigned int i = 0; i < fields.size(); i++){
                objArray[i] = fields[i]->get(obj);
                names[i] = fields[i]->getName();
            }
            if (returnStyle == Addin::expandMulti){
                expandArrayAcrossColumns(objArray, false, names,
                                         handleName, output);
            } else {
                expandArrayDownRows(objArray, names, handleName, output);
            }
        }
    } catch (exception& e){
        throw ModelException(e, "XLAddin::expand");
    }
}

/** Specialised version of expand. If canSimplify is true, then if it
    makes more sense to display this array in a single column then that is
    what is done. */
void XLAddin::expandArrayAcrossColumns(IArray&             arr,
                                       bool                canSimplify,
                                       const CStringArray& names,
                                       const string&       handleName,
                                       XL_OPER&               output) const {
    try{
        int numElts = arr.getLength();
        // find out total width required
        int numCols = 0;
        for (int i = 0; i < arr.getLength(); i++){
            IObjectSP elt(arr.get(i));
            if (!elt){
                numCols++;
            } else {
                const XLConvert* convertMethod
                    = XLConvertFactory::create(elt->getClass());
                int width = convertMethod->width();
                numCols += width > 0? width: 1;
                if (canSimplify && convertMethod != objectConvert){
                    // if we are just going to return an array of handles
                    // it is nicer to return them in a column
                    canSimplify = false;
                }
            }
        }
        if (canSimplify){
            expandArrayDownRows(arr, names, handleName, output);
        } else {
            // ensure our output oper is big enough
            XLConvert::setNumColsInOper(numCols, output);
            int colIdx = 0;
            for (int jdx = 0; jdx < numElts; jdx++){
                // get object out of array
                IObjectSP elt(arr.get(jdx));
                if (!elt){
                    XLConvert::blankColumn(0, colIdx, output);
                    colIdx++;
                } else {
                    if (elt->getRefCount() == 0){
                        // need to be careful if we get reference
                        elt = IObjectSP(elt->clone());
                    }
                    // look up method to deal with this type of object
                    const XLConvert* convertMethod =
                        XLConvertFactory::create(elt->getClass());
                    string name1 = handleNameParam? handleName: names[jdx];
                    string name2 = handleNameParam? names[jdx]:
                        elt->getClass()->getName();
                    if (convertMethod->width() < 1){
                        //  return a handle for variable width
                        XLConvert::makeRoomForSingleOper(colIdx, output);
                        XLConvert::createHandle(name1, name2, elt,
                                                colIdx, // column index
                                                -1,     // ignore
                                                output);
                        XLConvert::blankColumn(1, colIdx, output);
                        colIdx++;
                    } else {
                        // populate oper
                        convertMethod->convertObjectOutput(name1, name2, elt,
                                                           colIdx,
                                                           output);
                        colIdx += convertMethod->width();
                    }
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, "XLAddin::expandArrayAcrossColumns");
    }
}

/** Specialised version of expand */
void XLAddin::expandArrayDownRows(IArray&             arr,
                                  const CStringArray& names,
                                  const string&       handleName,
                                  XL_OPER&            output) const {
    try{
        int numElts = arr.getLength();
        // create an array the right length
        XLConvert::setNumColsInOper(1, output);
        XLConvert::setNumberOfRowsInOper(numElts, 0, output);
        for (int jdx = 0; jdx < numElts; jdx++){
            // get object out of array
            IObjectSP elt(arr.get(jdx));
            if (!elt){
                // could create handle to CNull object instead?
                output.val.xlarray.lparray[jdx].val.str = NEW(char);
                output.val.xlarray.lparray[jdx].type = xltypeStr;
            } else {
                // if the element is a atomic then it's nice just to return
                // it directly
                if (CDouble::TYPE->isInstance(elt) ||
                    CString::TYPE->isInstance(elt) ||
                    CInt::TYPE->isInstance(elt)    ||
                    CBool::TYPE->isInstance(elt)){
                    // get hold of XLConvert object
                    const XLConvert* convertMethod =
                        XLConvertFactory::create(elt->getClass());
                    XL_OPER& operToPop = XLConvert::offset(0, jdx, output);
                    convertMethod->convertObjectOutput(handleName, names[jdx],
                                                       elt, 0,
                                                       operToPop);
                } else {
                    if (elt->getRefCount() == 0){
                        // need to be careful if we get reference
                        elt = IObjectSP(elt->clone());
                    }
                    // and create a handle to it
                    XLConvert::createHandle(
                        handleNameParam? handleName:names[jdx],
                        handleNameParam? names[jdx]:elt->getClass()->getName(),
                        elt,
                        0, // all in 0'th column
                        jdx, // row index
                        output);
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, "XLAddin::expandArrayDownRows");
    }
}


/** Strips the ADDIN_PREFIX, if any, from a string */
string XLAddin::stripAddinPrefix(const string& addinName){
    if (strncmp(addinName.c_str(), ADDIN_PREFIX.c_str(),
                ADDIN_PREFIX.size())==0){
        return addinName.substr(ADDIN_PREFIX.size());
    }
    if (strncmp(addinName.c_str(), ADDIN_PREFIX_ALT.c_str(),
                ADDIN_PREFIX_ALT.size())==0){
        return addinName.substr(ADDIN_PREFIX_ALT.size());
    }
    return addinName;
}

int XLAddin::numParams() const{
    return (params.size() + (handleNameParam? 1: 0));
}

/** Return the number of optional parameters this addin has */
int XLAddin::numOptionalParams()const
{
    int numOptionalParams = 0;

    for (int i = 0; i < (int)params.size(); i++)
    {
        if (params[i]->isOptional())
        {
            numOptionalParams++;
        }
    }
    return numOptionalParams;
}

/* implements default register method needed for regression tester */
void XLAddin::defaultRegister(
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
    /* just need to create our instance */
    XLAddin::create(addinName, dataClass, params,
                    handleNameParam, returnStyle, returnIsNative,
                    methodType, method, returnType);
}

/* register function which determines whether we are being called from VB
   or not */
void XLAddin::setCallerIsVBMethod(TCallerIsVB *method){
    callerIsVBMethod = method;
}

/* uses registered method to determine whether we are being called
   from VB */
bool XLAddin::callerIsVB(){
    return callerIsVBMethod();
}


DRLIB_END_NAMESPACE
