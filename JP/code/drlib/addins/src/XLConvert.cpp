//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvert.cpp
//
//   Description : Base Class for converting to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/Handle.hpp"
#include "edginc/XLConvertFactory.hpp"
DRLIB_BEGIN_NAMESPACE

const string XLConvert::EXCEL_TRUE = "TRUE";
const string XLConvert::EXCEL_FALSE = "FALSE";
int XLConvert::xlAddinDateOffset = 109205; /* default value */;

const string VB_PREFIX = "vb_";

XLConvert::~XLConvert(){}

void XLConvert::setExcelDateOffsetSet(int xlDateOffset){
    xlAddinDateOffset = xlDateOffset;
}

/** Provides methods for converting data to/from excel implementation.
    Specialisations of this class exist for dealing with special types of
    objects (eg doubles, dates etc) */

/** invoke addin function which has a signature specific to the instance
    of this class. Result is passed via the jvalue structure */
void XLConvert::invoke(
    const Addin::Method& method,        // the method to invoke
    const jobject&       addin,         // contains addin arguments
    jvalue&              output) const // the return value
{
    try{
        IObject*   addinObj = addin.get();
        if (method.objMethod){
            output.l = method.objMethod->execute(addinObj);
        } else {
            /* no method - this could be a constructor or something like 
               object viewer */
            // no need to take copy here since XLAddin copies the input object
            // which here is "addin"
            output.l = CObject::convertToPrivateRep(addin);
        }
    } catch (exception& e){
        throw ModelException(e, "XLConvert::invoke");
    }
}

/** convert XLstring in XL_OPER to C++ style string  - no
    checking of type of XL_OPER */
string XLConvert::operStringToCPPString(const XL_OPER& oper){
    char buffer[EXCEL_MAX_STR_LEN];
    XLConvert::operStringToMemBuff(oper, buffer);
    return string(buffer);
}

/** Returns handle with given name. objectClass, if non null, is used
    to validate the type of the object */
IObjectSP XLConvert::getObjectFromHandle(
    const XL_OPER&    oper,
    CClassConstSP     objectClass) // can be null 
{
    string handleName(operStringToCPPString(oper));
    return Handle::fetch(handleName, objectClass);
}


/** Parse OPER and place data in given field in output object. This is
    the same as convertInput except that if the OPER is a single cell
    containing a handle which matches the desired type then that
    handle is used otherwise convertInput is invoked */
void XLConvert::convertGenericInput(
    CClassConstSP     desiredType,   // what type we want
    const XL_OPER&    oper,          // oper to convert
    const CFieldConstSP& field,   // where to store the converted data
    jobject&          output) const  // store value in here using field
{
    bool converted = false;
    if (operContainsHandle(oper)){
        // get object from handle
        IObjectSP obj = getObjectFromHandle(oper, 0);
        if (desiredType->isInstance(obj)){
             field->set(output, obj);
             converted = true;
        } 
    }
    if (!converted){
        convertInput(desiredType, oper, field, output);
    }
}

IObjectSP XLConvert::operToObject(
    CClassConstSP     desiredType,   // what type we want
    const XL_OPER&    oper)          // oper to convert
{
    IObjectSP obj;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit
    {
    case xltypeNum:
        obj = IObjectSP(CDouble::create(oper.val.num));
        break;
    case xltypeBool:
        obj = IObjectSP(CBool::create(oper.val.boolVal? true: false));
        break;
    case xltypeInt:
        obj = IObjectSP(CInt::create(oper.val.w));
        break;
    case xltypeStr:
    {
        char  buffer[EXCEL_MAX_STR_LEN];
        XLConvert::operStringToMemBuff(oper, buffer);
        obj = IObjectSP(CString::create(buffer));
        break;
    }
    default:
        throw ModelException("XLConvert::operToObject",
                             "Cell contains \""+operToString(oper)+"\"");
    }
    CObject::checkType(obj, desiredType);
    return obj;
}

/** Parse OPER and place data in given field in output object */
void XLConvert::convertInput(
    CClassConstSP     desiredType,   // what type we want
    const XL_OPER&    oper,          // oper to convert
    const CFieldConstSP& field,   // where to store the converted data
    jobject&          output) const  // store value in here using field
{
    static const string routine = "XLConvert::convertInput";
    IObjectSP obj;
    // check we've got a handle
    if (!operContainsHandle(oper)){
        // try converting cell to an object
        try{
            if (oper.type & xltypeMulti){
                throw ModelException(routine, "(A range of more"
                                     " than one cell was supplied)");
            }
            obj = operToObject(desiredType, oper);
        } catch (exception& e){
            throw ModelException(e, routine, "For parameters of type "+
                                 desiredType->getName()+
                                 " a handle must be supplied. ");
        }
    } else {
        // get object from handle
        obj = getObjectFromHandle(oper, desiredType);
    }
    // set in output
    field->set(output, obj);
}

/** Take relevant value from data and fill column idx in output 
    accordingly. Default is to create a handle */
void XLConvert::convertObjectOutput(
    const string&  handleName,    // for handles
    const string&  name,          /* name for this bit of data
                                     (for handles/error messages) */
    const jobject& data,          // populate OPER with the object
    int            idx,           // column index number */
    XL_OPER&       output) const  // converted output
{
    // ensure our output oper is big enough
    makeRoomForSingleOper(idx, output);
    createHandle(handleName, name, data, idx, -1, output); 
    blankColumn(1, idx, output);
}


/** Take relevant value from data and fill column idx in output 
    accordingly. Default is to create a handle */
void XLConvert::convertNativeOutput(
    const string&  name,          /* name for this bit of data
                                     (for handles/error messages) */
    const jvalue&  data,          // value to be put into OPER 
    int            idx,           // column index number */
    XL_OPER&       output) const  // converted output
{
    convertObjectOutput(name, "", data.l, idx, output); 
}

/** convert value held in data to an object. Default implementation
    does nothing since it is already an object */
void XLConvert::nativeToObject(jvalue&  data) const{}

/** Returns the width needed to dsiplay this type of object in
    XL. A return value of -1 indicates a variable width. This is
    used for addins that use Addin:expandMulti */
int XLConvert::width() const{
    // default for most types (including all complex objects/lists thereof)
    return 1;
}


/** Returns true if the oper contains a handle */
bool XLConvert::operContainsHandle(const XL_OPER& oper){
    return (oper.type & xltypeStr && Handle::valid(oper.val.str));
}

/* throw exception for unrecognised oper type */
void XLConvert::unknownOper(const string& routine,  /* (I) */
                            const XL_OPER&   oper){    /* (I) */
    throw ModelException(routine.empty()? "XLConvert::unknownOper": routine,
                         "Unknown excel oper type ("+
                         Format::toString(oper.type)+")!!");
}

/* blank out the remaining cells in a column */
void XLConvert::blankColumn(
    int      nextCell,  /* (I) cell to start from  */
    int      idx,       /* (I) column index number */
    XL_OPER& oper)      /* (M) array gets altered */
{
    if (oper.type & xltypeMulti) {
        int jdx;
        for (jdx = nextCell; jdx < oper.val.xlarray.rows; jdx++) {
            XL_OPER *cell = oper.val.xlarray.lparray +
                offset(idx, jdx, oper.val.xlarray.columns, oper.val.xlarray.rows);
            cell->type = xltypeErr;
            cell->val.err = xlerrNA;
        }
    }
}

/* convert XLstring in OPER to C style string by filling up a
   charbuffer - no checking of type of OPER */
void XLConvert::operStringToMemBuff(const XL_OPER& oper, 
                                    char buffer[EXCEL_MAX_STR_LEN])
{
    char *xlString = oper.val.str;
    memcpy(buffer, xlString + 1, (size_t) ((unsigned char) xlString[0]));
    buffer[(unsigned char)xlString[0]] = '\0';
}

/** Populate oper with supplied string */
void XLConvert::populateOperFromString(const string& value,
                                       XL_OPER&      oper){
    static const string routine = "XLConvert::populateOperFromString";
    int   length = value.size();
    if (length > EXCEL_MAX_STR_LEN){
        throw ModelException(routine, "The following string was too long "
                             "for excel\n<------ start of string ------>\n"
                             +value+"<------ end of string ------>");
    }
    char *xlString = NEW_ARRAY(char, length+1);
    memcpy(xlString+1, value.c_str(), length); /* CANNOT use strcpy */
    xlString[0] = (char) length;
    oper.val.str = xlString;
    oper.type = xltypeStr;
}
  

/** Populate oper given by supplied offsets with supplied string */
void XLConvert::populateOperFromString(const string& value,
                                       int           xPos,
                                       int           yPos,
                                       XL_OPER&      oper){
    static const string routine = "XLConvert::populateOperFromString";
    XL_OPER *cell;

    if (!(oper.type & xltypeMulti)) {
        if (xPos != 0 || yPos != 0)  {
            throw ModelException(routine, "OPER must be of type xltypeMulti"
                                 " for non-origin coords");
        }
        cell = &oper;
    } else {
        cell = oper.val.xlarray.lparray + 
            offset(xPos, yPos, oper.val.xlarray.columns, oper.val.xlarray.rows);
    }
    populateOperFromString(value, *cell);
}
    

/** build handle to object and populate OPER with string representing
    handle name. The paramIdx is used as a x coord for writing to the
    oper array */
void XLConvert::createHandle(
    const string&   name1,    /* (I) ignored if "" */
    const string&   name2,    /* (I) ignored if "" */
    const jobject   object,   /* (I) */
    int             paramIdx, /* (I) x coord for writing to oper */
    int             listIdx,  /* (I) >=0:append to handle name and use 
                                 as y coord*/
    XL_OPER&        oper)      /* (M) */
{
    string userName;
    if (!name1.empty() && !name2.empty()) {
        userName = name1+"_"+name2;
    } else if (name1.empty()) {
        userName = !name2.empty()? name2: "NULL";
    } else {
        userName = name1;
    }

    /* return handle - first need handle name */
    if ( name1 == VB_PREFIX) {
        static int handleCount = 0;
        userName += Format::toString(handleCount);
        handleCount = (handleCount+1) % Handle::MAX_ID;
        /* return handle - first need handle name */
        if (listIdx >= 0) {
            userName += Format::toString((listIdx + 1) % Handle::MAX_ID);
        }
    } else if (listIdx >= 0) {
        userName += "_"+Format::toString((listIdx + 1) % Handle::MAX_ID);
    }
    string handleName(Handle::create(userName, object));
    populateOperFromString(handleName, 
                           paramIdx,
                           listIdx >= 0? listIdx: 0,
                           oper);
}

/** If the oper is not already signifying an error, alter so that it
    is - do not throw an exception */
void XLConvert::setToError(XL_OPER&  oper) throw(){
    if (!(oper.type & xltypeErr)){
        if (oper.type & xltypeMulti){
            size_t pos;    /* Iterator for array */     
            size_t size = oper.val.xlarray.rows * oper.val.xlarray.columns;
            for (pos = 0; pos < size; ++pos)
            {
                XL_OPER *xpos = oper.val.xlarray.lparray + pos;
                if (xpos->type & xltypeStr)
                {
                    FREE(xpos->val.str);
                }
            }
            FREE(oper.val.xlarray.lparray);
        } else if (oper.type & xltypeStr){
            FREE(oper.val.str);
        }
        oper.type = xltypeErr;
        oper.val.err = xlerrValue;
    }
}

/** Check the validity of the supplied OPER. The return value
    indicates whether processing should continue for this object ie
    true everything is okay, false: silently return N/A to the
    spreadsheet */
bool XLConvert::validateOPER(
    const XL_OPER&   oper,         /* (I) */
    CClassConstSP    desiredType,  // (I) required type
    int              paramNum,     /* (I) which parameter we're on */
    XL_OPER&         output)       /* (M) error code written here on failure */
{
    static const string routine = "XLConvert::validateOPER";
    bool         processOper    = true; // default
    // refuse to use opers which are in error, missing or blank (ie type nil)
    // unless they are for arrays
    if (oper.type & xltypeErr || oper.type & xltypeMissing
        || (oper.type & xltypeNil && !desiredType->isArray()))
    {
        setToError(output);
        if (oper.type & xltypeMissing || oper.type & xltypeNil) {
            output.val.err = xlerrValue;
            string m = (oper.type & xltypeMissing)? "missing": 
                "blank (ie supplied from an empty cell)";
            throw ModelException(routine, "Non optional Parameter "+
                                 Format::toString(paramNum)+" is "+m);
        }

        /* suppress error messages - dependent cell in error */
        processOper = false;
        /* report all other errors as NA - stop cascade of error messages*/
        output.val.err = xlerrNA;
    }
    else if (oper.type & xltypeMulti)
    {
        int numElts = oper.val.xlarray.rows * oper.val.xlarray.columns;
        XL_OPER *subOper = oper.val.xlarray.lparray;
        for (int idx = 0; idx < numElts && processOper; idx++, subOper++)
        {
            if (subOper->type & xltypeErr || subOper->type & xltypeMissing){
                // generate exception or fail silently
                validateOPER(*subOper, desiredType, paramNum, output);
                processOper = false;
            }
        }
    }
    return processOper;
}

   

/* Preliminary examination of input data. Determine whether we should
   immediately fail or not (and if so whether to do so
   quietly). Additionally 'useDefault' flag indicates, on input,
   whether a default is available. On output, its value reflects
   whether the default value for the parameter should be used.
   The return value indicates whether processing should continue for
   this object ie true everything is okay, false: silently return N/A
   to the spreadsheet */
bool XLConvert::validateInput(
    const XL_OPER&   oper,           /* (I) oper to validate */
    CClassConstSP    desiredType,    // (I) required type
    int              paramNum,       // (I) which parameter (for error message)
    XL_OPER&         output,         /* (M) error code written 
                                        here on failure */
    bool&            useDefault) const  /* (M) see above */
{
    bool processOper = true; // the return value
    // if multiple cells blank and not an array then use default
    if (useDefault && (oper.type & xltypeMulti) && !desiredType->isArray()){
        int numCells = oper.val.xlarray.rows * oper.val.xlarray.columns;
        for (int i = 0; i < numCells && useDefault; i++){
            if (!(oper.val.xlarray.lparray[i].type & xltypeNil)){
                useDefault = false;
            }
        }
    }
    /* else if parameter is optional then use default value if 
       type missing/nil or (xltype string and string = "" and required type not
       string) */
    else if (useDefault){
        if ((oper.type & xltypeMissing) || 
            (oper.type & xltypeNil && !desiredType->isArray())) {
            // use default value
        } else {
            if (!(oper.type & xltypeStr && *oper.val.str == 0 &&
                  desiredType != CString::TYPE)) {
                // in general don't use default value
                useDefault = false;
            }
        }
    }
    if (!useDefault) {
        /* we have some sort of data so use it */
        /* check that the parameter is valid */
        processOper = validateOPER(oper, desiredType, paramNum, output);
    }
    return processOper;
}

/** Sets the number of cols in an oper (and manages memory accordingly) and
    the number of rows to 1. Oper must not be of type xltypeMulti already
    New opers have all bits set to 0 */
void XLConvert::setNumColsInOper(int         numCols,
                                 XL_OPER&    oper){
    if (oper.type & xltypeMulti){
        throw ModelException("XLConvert::setNumColsInOper", 
                             "oper is already of type xltypeMulti");
    }
    oper.val.xlarray.lparray = 0;
    oper.type = xltypeMulti;
    oper.val.xlarray.rows = 0;
    oper.val.xlarray.columns = (unsigned short) numCols;
}

/** Increase the number of rows in an oper by the number given. New opers
    have all bits set to 0 */
void XLConvert::increaseRowsInOper(int       numRows,
                                   XL_OPER&  oper){
    if (oper.type & xltypeMulti){
        if (numRows > 0){
            oper.val.xlarray.lparray = REALLOC(
                oper.val.xlarray.lparray, XL_OPER, 
                (oper.val.xlarray.rows + numRows)* oper.val.xlarray.columns);
            // clear new memory allocated
            memset(oper.val.xlarray.lparray + 
                   oper.val.xlarray.rows * oper.val.xlarray.columns, 
                   0, sizeof(XL_OPER) * numRows * oper.val.xlarray.columns);
        }
        oper.val.xlarray.rows += numRows;
    } else {
        XL_OPER* block = 0;
        if (numRows > 0){
            block = NEW_ARRAY(XL_OPER, numRows);
            block[0] = oper; // structure copy
        } 
        oper.type = xltypeMulti;
        oper.val.xlarray.rows = (unsigned short) numRows;
        oper.val.xlarray.columns = (unsigned short) 1;
        oper.val.xlarray.lparray = block;
    }
}

/** Increases the number of rows in an oper to the number
    given. Any opers created for columns to the left (ie <) of
    colIdx are set to #NA otherwise new opers have all bits set to
    0. Does nothing if there are already given number of rows */
void XLConvert::setNumberOfRowsInOper(int       numRows,
                                      int       colIdx,
                                      XL_OPER&  oper){
    if (!(oper.type & xltypeMulti)){
        increaseRowsInOper(numRows, oper);
    } else if (oper.val.xlarray.rows < numRows){
        int origNumRows = oper.val.xlarray.rows;
        increaseRowsInOper(numRows - oper.val.xlarray.rows, oper);
        for (int i = 0; i < colIdx; i++){
            for (int j = origNumRows; j < numRows; j++){
                XL_OPER *cell = oper.val.xlarray.lparray +
                    offset(i, j, oper.val.xlarray.columns, oper.val.xlarray.rows);
                cell->type = xltypeErr;
                cell->val.err = xlerrNA;
            }
        }
    }
}

/** Ensure oper is big enough to contain a single oper output. If the oper is
    xltypeMulti then setNumberOfRowsInOper() is used to ensure the oper 
    contains at least one row otherwise no action is taked */
void XLConvert::makeRoomForSingleOper(int        colIdx,
                                      XL_OPER&   oper){
    if (oper.type & xltypeMulti){
        setNumberOfRowsInOper(1, colIdx, oper);
    }
}

/** Increase size of multi oper so at least 2 rows tall and 1 col wide */
void XLConvert::ensureOperIsMinimumSize(XL_OPER&   oper){
    if (oper.type & xltypeMulti){
        int origNumRows = oper.val.xlarray.rows;
        if (oper.val.xlarray.columns < 1){
            if (origNumRows > 0){
                oper.val.xlarray.lparray = NEW_ARRAY(XL_OPER, origNumRows);
            }
            oper.val.xlarray.columns = 1;
        }
        if (origNumRows < 2) {
            XLConvert::increaseRowsInOper(2 - origNumRows, oper);
            for (int i = 0; i < oper.val.xlarray.columns; i++){
                XLConvert::blankColumn(origNumRows, i, oper);
            }
        }
    }
}

/** ensure xlType is multi - fail otherwise */
void XLConvert::ensureOperIsMultiType(const XL_OPER& oper){
    if (!(oper.type & xltypeMulti)){
        throw ModelException("XLConvert::ensureOperIsMultiType",
                             "Internal error - expected oper of type multi");
    }
}


/** Converts an oper to a double */
double XLConvert::operToDouble(const XL_OPER& oper){
    static const string routine = "XLConvert::operToDouble";
    double val;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit
    {
    case xltypeNum:
        val = oper.val.num;
        break;
    case xltypeInt:
        val = oper.val.w;
        break;
    case xltypeNil:
    case xltypeMissing:
        val = 0.0;
        break;
    case xltypeBool:
        val = (double) oper.val.boolVal;
        break;
    case xltypeStr:
    {
        if (operContainsHandle(oper)){
            IObjectSP object = getObjectFromHandle(oper, CDouble::TYPE);
            CDoubleSP doubleObj = CDoubleSP::dynamicCast(object);
            val = doubleObj->doubleValue();
        } else {
            char *endPtr;
            char  buffer[EXCEL_MAX_STR_LEN];
            operStringToMemBuff(oper, buffer);
            val = strtod(buffer, &endPtr);
            if (*endPtr != '\0') {
                throw ModelException(routine, "Failed to convert string "+
                                     string(buffer)+ " to a double");
            }
        }
        break;
    }
    case xltypeMulti:
        throw ModelException(routine, "A single double was expected but a "
                             "range was suppplied");
    default:
        unknownOper(routine, oper);
        val = 0.0; // shut compiler up
    }
    return val;
}

/* Returns the length of the array of opers. The count finishes on
   the cell before the first nil cell. Counting is performed by
   going down vertically */
int XLConvert::arrayLength(
    bool           emptyStringAsNil, /* (I) treat "" as xltypeNil */
    const XL_OPER& oper)            /* (I) of type xltypeMulti */
{
    int   actualLen = oper.val.xlarray.rows * oper.val.xlarray.columns;
    int   idx;
    for (idx = 0; idx < actualLen &&
             !(oper.val.xlarray.lparray[idx].type & xltypeNil) &&
             !(emptyStringAsNil &&
               (oper.val.xlarray.lparray[idx].type & xltypeStr) &&
               *oper.val.xlarray.lparray[idx].val.str == 0);
         idx++); /* empty */
    return idx;
}    

/* Returns the number of rows in the grid of opers. The count
   finishes on the cell before the first nil cell. Counting is
   performed by going down vertically on the LHS column. */
int XLConvert::gridLength(
    bool            emptyStringAsNil,  /* (I) treat "" as xltypeNil */
    const XL_OPER&  oper)              /* (I) of type xltypeMulti */
{
    int actualLen = oper.val.xlarray.rows;
    int idx;
    int cols;
    for (idx = 0, cols = 0; idx < actualLen &&
             !(oper.val.xlarray.lparray[cols].type & xltypeNil) &&
             !(emptyStringAsNil && 
               (oper.val.xlarray.lparray[cols].type & xltypeStr) &&
               *oper.val.xlarray.lparray[cols].val.str == 0);
         idx++, cols += oper.val.xlarray.columns); /* empty */
    return idx;
}    

/* fill in tmpOper to be of type xltypeMulti and point to origOper */
void XLConvert::coerceOperToMulti(const XL_OPER& origOper, /* (I) */
                                  XL_OPER&       tmpOper){  /* (M) */
    /* coerce to xltypeMulti */
    tmpOper.type = xltypeMulti;
    tmpOper.val.xlarray.rows = 1;
    tmpOper.val.xlarray.columns = 1;
    // const_cast seems the only way out here
    tmpOper.val.xlarray.lparray = &const_cast<XL_OPER&>(origOper);
}

/* if xltypeMulti and num cols * num rows == 1, returns oper.val.xlarray
   else returns origOper */
const XL_OPER& XLConvert::coerceOperFromMulti(const XL_OPER& origOper){
    if ((origOper.type & xltypeMulti) && 
        origOper.val.xlarray.rows * origOper.val.xlarray.columns == 1){
        return *origOper.val.xlarray.lparray;
    }
    return origOper;
}

/** Converts an oper to an int */
int XLConvert::operToInt(const XL_OPER& oper){
    static const string routine = "XLConvert::operToInt";
    int  val;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit
    {
    case xltypeNum:
        val = (int)oper.val.num;
        if (!Maths::equals((double)val, oper.val.num)) {
            throw ModelException(routine, "Cannot convert "+
                                 Format::toString(oper.val.num)+ " to an int");
        }
        break;
    case xltypeInt:
        val = oper.val.w;
        break;
    case xltypeNil:
    case xltypeMissing:
        val = 0;
        break;
    case xltypeBool:
        val = (int) oper.val.boolVal;
        break;
    case xltypeStr:
    {
        if (operContainsHandle(oper)){
            IObjectSP object = getObjectFromHandle(oper, CInt::TYPE);
            CIntSP intObj = CIntSP::dynamicCast(object);
            val = intObj->intValue();
        } else {
            char *endPtr;
            char  buffer[EXCEL_MAX_STR_LEN];
            operStringToMemBuff(oper, buffer);
            val = strtol(buffer, &endPtr, 10 /* base 10! */);
            if (*endPtr != '\0') {
                throw ModelException(routine, "Failed to convert string "+
                                     string(buffer)+ " to a int");
            }
        }
        break;
    }
    case xltypeMulti:
        throw ModelException(routine, "A single int was expected but a "
                             "range was suppplied");
    default:
        unknownOper(routine, oper);
        val = 0; // shut compiler up
    }
    return val;
}

/** Converts an oper to a string */
string XLConvert::operToString(const XL_OPER& oper){
    return operToString(oper, true);
}

/** Converts an oper to a string */
string XLConvert::operToString(const XL_OPER& oper, bool checkForHandles){
    static const string routine = "XLConvert::operToString";
    string  val;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit
    {
    case xltypeNum:
        val = Format::toString(oper.val.num);
        break;
    case xltypeInt:
        val = Format::toString(oper.val.num);
        break;
    case xltypeNil:
    case xltypeMissing:
        val = "";
        break;
    case xltypeBool:
        val = oper.val.boolVal? XLConvert::EXCEL_TRUE: 
            XLConvert::EXCEL_FALSE;
        break;
    case xltypeStr:
    {
        if (checkForHandles && operContainsHandle(oper)){
            IObjectSP object = getObjectFromHandle(oper, CString::TYPE);
            CStringSP stringObj = CStringSP::dynamicCast(object);
            val = stringObj->stringValue();
        } else {
            char  buffer[EXCEL_MAX_STR_LEN];
            operStringToMemBuff(oper, buffer);
            val = string(buffer);
        }
        break;
    }
    case xltypeMulti:
        throw ModelException(routine, "A single string was expected but a "
                             "range was suppplied");
    default:
        unknownOper(routine, oper);
    }
    return val;
}

/** Converts an oper to a bool */
bool XLConvert::operToBool(const XL_OPER& oper){
    static const string routine = "XLConvert::operToBool";
    bool  val;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit
    {
    case xltypeNum:
        if (Maths::isZero(oper.val.num)){
            val = false;
        } else if (Maths::equals(oper.val.num, 1.0)){
            val = true;
        } else {
            throw ModelException(routine, "Cannot convert "+
                                 Format::toString(oper.val.num)+" to a"
                                 " boolean");
        }
        break;
    case xltypeInt:
        if (oper.val.w == 0){
            val = false;
        } else if (oper.val.w == 1){
            val = true;
        } else {
            throw ModelException(routine, "Cannot convert "+
                                 Format::toString(oper.val.w)+" to a"
                                 " boolean");
        }
        break;
    case xltypeNil:
    case xltypeMissing:
        val = false;
        break;
    case xltypeBool:
        val = oper.val.boolVal? true: false;
        break;
    case xltypeStr:
    {
        if (operContainsHandle(oper)){
            IObjectSP object = getObjectFromHandle(oper, CBool::TYPE);
            CBoolSP boolObj = CBoolSP::dynamicCast(object);
            val = boolObj->boolValue();
        } else {
            char *pos;
            char  buffer[EXCEL_MAX_STR_LEN];
            operStringToMemBuff(oper, buffer);
            /* convert to upper case */
            for (pos = buffer; (*pos = toupper(*pos));
                 pos++); /* empty loop */
            
            /* supporting TRUE, T, YES, Y, 1 and FALSE, F, NO, F, 0 */
            if (!strcmp(buffer, EXCEL_TRUE.c_str()) ||
                !strcmp(buffer, "T") ||
                !strcmp(buffer, "YES") || !strcmp(buffer, "Y") ||
                !strcmp(buffer, "1")) {
                val = true;
            } 
            else if (!strcmp(buffer, EXCEL_FALSE.c_str()) ||
                     !strcmp(buffer, "F")  || 
                     !strcmp(buffer, "NO") || !strcmp(buffer, "N") ||
                     !strcmp(buffer, "0")) {
                val = false;
            } else {
                throw ModelException(routine, "Failed to convert "+
                                     string(buffer)+" to a boolean");
            }
        }
        break;
    case xltypeMulti:
        throw ModelException(routine, "A single boolean was expected but a "
                             "range was suppplied");
    }
    default:
        unknownOper(routine, oper);
        val = false; // shut compiler up
    }
    return val;
}

/** same as xlAutoFree but does not free container memory */
void XLConvert::operFreeSimple(XL_OPER& x){
    if (x.type & xltypeMulti) {   
        size_t pos;    /* Iterator for array */     
        size_t size = x.val.xlarray.rows * x.val.xlarray.columns;
        for (pos = 0; pos < size; ++pos) {
            XL_OPER *xpos = x.val.xlarray.lparray + pos;
            if (xpos->type & xltypeStr) {
                FREE(xpos->val.str);
            }
        }
        FREE(x.val.xlarray.lparray);
    } else if (x.type & xltypeStr) {
        FREE(x.val.str);
    } else if (x.type & xltypeRef){
        FREE(x.val.mref.lpmref);
    }
}

bool  XLConvert::isSimpleType() const {
    return false;
}
    
/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDouble::TYPE is initialised. Outside of class to avoid necessity
    of header file - in this case consistency with other XLConvert classes */
void XLConvertRegister(){
    XLConvertFactory::registerXLObjectConvert(IObject::TYPE, new XLConvert());
}

XLConvert::XLCoerceToMulti*   XLConvert::xlCoerceToMulti = 0;
XLConvert::XLFreeCoerce*      XLConvert::xlFreeCoerce = 0;

/** register a method which can be used for XLCoerceToMulti */
void XLConvert::setXLCoerceMethods(XLCoerceToMulti* coerceMethod,
                                   XLFreeCoerce*    freeMethod){
    xlCoerceToMulti = coerceMethod;
    xlFreeCoerce = freeMethod;
}

/** Requests conversion of XLOPER to type xlmulti. If successful, output
    must be passed to Free */
int XLConvert::coerceToMulti(const XLOPER& input, XLOPER& output){
    /* in theory can just do 
       "return Excel(xlCoerce, &output, 2, &input, TempInt(xltypeMulti))"
       However, with massive arrays Excel can't cope. The solution is to
       manually request each individual cell. This works but is painfully slow
       when doing a global recalc. The solution to this is to try and start
       calculating cells from where you were last time. This is what the
       previousRowOffset is all about */
    static int        previousRowOffset = -1;
    static const int  CALC_THRESHOLD = 10; // arbitrary really
    bool   isSRef = input.type & xltypeSRef? true: false;
    if (isSRef || (input.type & xltypeRef)){
        const XLREF& xlRef = isSRef? 
            input.val.sref.ref: *input.val.mref.lpmref->reftbl;
        if ((xlRef.colFirst+1-xlRef.colLast)* 
            (xlRef.rwLast+1-xlRef.rwFirst) > CALC_THRESHOLD){
            XLOPER theCell;
            XLMREF xlmRef;
            theCell.type = input.type;
            if (isSRef){
                theCell.val.sref.count = 1;
            } else {
                theCell.val.mref.idSheet = input.val.mref.idSheet;
                theCell.val.mref.lpmref = &xlmRef;
                xlmRef.count = 1;
            }
            XLREF& newRef = isSRef? theCell.val.sref.ref: *xlmRef.reftbl;
            for (int col = xlRef.colFirst; col <= xlRef.colLast; col++){
                newRef.colFirst = col;
                newRef.colLast = col;
                // reset previousRow if out of bounds
                int previousRow = previousRowOffset + xlRef.rwFirst;
                if (previousRow < xlRef.rwFirst){
                    previousRow = xlRef.rwFirst-1;
                } else if (previousRow > xlRef.rwLast){
                    previousRow = xlRef.rwLast;
                }
                int numRows = xlRef.rwLast - xlRef.rwFirst + 1;
                for (int offset = 0; offset < numRows; offset++){
                    // split into two parts: from previousRow to rwLast and
                    // rwFirst to previousRow
                    int row = previousRow+1+offset;
                    if (row > xlRef.rwLast){
                        row -= numRows;
                    }
                    newRef.rwFirst = row;
                    newRef.rwLast = row;
                    if ((*xlCoerceToMulti)(theCell, output) == xlretUncalced){
                        previousRowOffset = row - xlRef.rwFirst;
                        return xlretUncalced;
                    }
                    xlFreeCoerce(output);
                }
            }
        }
    }       
    return ((*xlCoerceToMulti)(input, output));
}

/** Free the contents of the OPER created in CoerceToMulti method */
void XLConvert::freeCoerce(XLOPER& toFreeContentsOf){
    (*xlFreeCoerce)(toFreeContentsOf);
}


DRLIB_END_NAMESPACE

