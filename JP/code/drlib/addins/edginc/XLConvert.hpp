//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvert.hpp
//
//   Description : Base Class for converting to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_XLCONVERT_HPP
#define EDG_XLCONVERT_HPP

#include "edginc/Object.hpp"
#include "edginc/Field.hpp"
#include "edginc/xlapi.hpp"
#include "edginc/Addin.hpp"
DRLIB_BEGIN_NAMESPACE


typedef bool	  jboolean;
typedef int       jint;
typedef double    jdouble;
typedef IObjectSP jobject;

/* nice way of dealing with different signatures - except that can't have
   jobject inside union as its got constructors/destructors*/
typedef union jAtomicValue{
    jboolean    z;
    // jbyte    b;
    // jchar    c;
    // jshort   s;
    jint        i;
    // jlong    j;
    // jfloat   f;
    jdouble     d;
} jAtomicValue;

typedef struct jvalue {
    jAtomicValue a;
    jobject      l;
} jvalue;


// MSVC can't handle this being a constant static field
#define EXCEL_MAX_STR_LEN 255
/* max number of input parameters which can go through excel */
#define EDG_EXCEL_MAX_PARAMS 29

/** Provides methods for converting data to/from excel implementation.
    Specialisations of this class exist for dealing with special types of
    objects (eg doubles, dates etc). The virtual methods here deal with
    data types of XL_OPER. If XL_OPER is actually XLOPER then it can be
    taken that any reference types have already been converted to 
    xltype multi. */
class ADDINS_DLL XLConvert{
public:
    static const string EXCEL_TRUE;
    static const string EXCEL_FALSE;

    virtual ~XLConvert();

    /** Requests conversion of XLOPER to type xlmulti. If successful, output
        must be passed to Free */
    static int coerceToMulti(const XLOPER& input, XLOPER& output);

    /** Free the contents of the OPER created in CoerceToMulti method */
    static void freeCoerce(XLOPER& toFreeContentsOf);

    /* how to access (x, y) position in XL_OPER * for a 
       rectangle width * height.
       height is used artificially to stop compiler warnings */
    inline static int offset(int x, int y, int width, int height){
        return ((y * width) + x);
    }

    /** Returns XL_OPER offset by given coordinates in oper (which must be of
        type xltypeMulti */
    inline static XL_OPER& offset(int x, int y, XL_OPER& oper){
        return *(oper.val.xlarray.lparray + offset(x, y, 
                                                 oper.val.xlarray.columns,
                                                 oper.val.xlarray.rows));
    }

    /** Returns XL_OPER offset by given coordinates in oper (which must be of
        type xltypeMulti */
    inline static const XL_OPER& offset(int x, int y, const XL_OPER& oper){
        return *(oper.val.xlarray.lparray + offset(x, y, 
                                                 oper.val.xlarray.columns,
                                                 oper.val.xlarray.rows));
    }

    /** invoke addin function which has a signature specific to the instance
        of this class. Result is passed via the jvalue structure */
    virtual void invoke(
        const Addin::Method& method,        // the method to invoke
        const jobject&       addin,         // contains addin arguments
        jvalue&              output) const; // the return value

    /** Parse XL_OPER and place data in given field in output object. This is
        the same as convertInput except that if the XL_OPER is a single cell
        containing a handle which matches the desired type then that
        handle is used otherwise convertInput is invoked */
    virtual void convertGenericInput(
        CClassConstSP  desiredType,   // what type we want
        const XL_OPER& oper,          // oper to convert
        const CFieldConstSP& field,   // where to store the converted data
        jobject&       output) const; // store value in here using field

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP        desiredType,   // what type we want
        const XL_OPER&       oper,          // oper to convert
        const CFieldConstSP& field,     // where to store the converted data
        jobject&             output) const; // store value in here using field

    /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertNativeOutput(
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jvalue&  data,          // value to be put into XL_OPER 
        int            idx,           // column index number */
        XL_OPER&       output) const; // converted output

    /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& data,          // populate XL_OPER with the object
        int            idx,           // column index number */
        XL_OPER&       output) const; // converted output

    /* Preliminary examination of input data. Determine whether we should
       immediately fail or not (and if so whether to do so
       quietly). Additionally 'useDefault' flag indicates, on input,
       whether a default is available. On output, its value reflects
       whether the default value for the parameter should be used.
       The return value indicates whether processing should continue for
       this object ie true everything is okay, false: silently return N/A
       to the spreadsheet */
    virtual bool validateInput(
        const XL_OPER& oper,          /* (I) oper to validate */
        CClassConstSP  desiredType,   // (I) required type
        int            paramNum,     // (I) which parameter (for error message)
        XL_OPER&       output,        /* (M) error code written 
                                           here on failure */
        bool&          useDefault) const; /* (M) see above */

    /** convert value held in data to an object. Default implementation
        does nothing since it is already an object */
    virtual void nativeToObject(jvalue&  data) const;

    /** Returns the width needed to dsiplay this type of object in
        XL. A return value of -1 indicates a variable width. This is
        used for addins that use Addin:expandMulti */
    virtual int width() const;

    /** returns true if object is of simple type */
    virtual bool isSimpleType() const; 

    /** build handle to object and populate XL_OPER with string representing
        handle name. The paramIdx is used as a x coord for writing to the
        oper array */
    static void createHandle(
        const string&   name1,    /* (I) ignored if "" */
        const string&   name2,    /* (I) ignored if "" */
        const jobject   object,   /* (I) */
        int             paramIdx, /* (I) x coord for writing to oper */
        int             listIdx,  /* (I) >=0:append to handle name and use 
                                     as y coord*/
        XL_OPER&        oper);     /* (M) */

    /** Return handle name from oper */
    static string convertHandleName(XL_OPER&  oper);

    /** If the oper is not already signifying an error, alter so that it
        is - do not throw an exception */
    static void setToError(XL_OPER&  oper) throw();

    /** convert XLstring in OPER to C++ style string  - no
        checking of type of OPER */
    static string operStringToCPPString(const XL_OPER& oper);

    /** Increases the number of rows in an oper to the number
        given. Any opers created for columns to the left (ie <) of
        colIdx are set to #NA otherwise new opers have all bits set to
        0. Does nothing if there are already given number of rows */
    static void setNumberOfRowsInOper(int    numRows,
                                      int    colIdx,
                                      XL_OPER&  oper);
                                 
    /** Sets the number of cols in an oper (and manages memory
        accordingly) and the number of rows to 1. Oper must not be of type
        xltypeMulti already. New opers have all bits set to 0 */
    static void setNumColsInOper(int      numCols,
                                 XL_OPER& oper);
    
    /** Ensure oper is big enough to contain a single oper output. If
        the oper is xltypeMulti then setNumberOfRowsInOper() is used
        to ensure the oper contains at least one row otherwise no
        action is taked */
    static void makeRoomForSingleOper(int        colIdx,
                                      XL_OPER&   oper);

    /** Increase size of multi oper so at least 2 rows tall and 1 col wide */
    static void ensureOperIsMinimumSize(XL_OPER&   oper);

    /** ensure xlType is multi - fail otherwise */
    static void ensureOperIsMultiType(const XL_OPER& oper);

    /** Converts an oper to a double */
    static double operToDouble(const XL_OPER& oper);

    /** Converts an oper to an int */
    static int operToInt(const XL_OPER& oper);

    /** Converts an oper to a string. Any strings which are handles are
        converted to the object in the handle */
    static string operToString(const XL_OPER& oper);

    /** Converts an oper to a string. checkForHandles drives whether strings
     which are handles are converted to the wrapped object or not */
    static string operToString(const XL_OPER& oper, bool checkForHandles);

    /** Converts an oper to a bool */
    static bool operToBool(const XL_OPER& oper);

    /** same as xlAutoFree but does not free container memory */
    static void operFreeSimple(XL_OPER& oper);

    /** Populate oper with supplied string */
    static void populateOperFromString(const string& value,
                                       XL_OPER&      oper);

    /** Returns true if the oper contains a handle */
    static bool operContainsHandle(const XL_OPER& oper);

    /** Returns handle with given name. objectClass, if non null, is used
        to validate the type of the object */
    static IObjectSP getObjectFromHandle(
        const XL_OPER& oper,
        CClassConstSP  objectClass); // can be null 

    /** blank out the remaining cells in a column */
    static void blankColumn(
        int       nextCell,  /* (I) cell to start from  */
        int       idx,       /* (I) column index number */
        XL_OPER&  oper);     /* (M) array gets altered */

    /** Requests conversion of XLOPER to type xlmulti. If successful, output
        must be passed to XLFree. Return code as defined by XL */
    typedef int (XLCoerceToMulti)(const XLOPER& input, XLOPER& output);

    /** Free the contents of the OPER created in XLCoerceToMulti method */
    typedef void (XLFreeCoerce)(XLOPER& toFreeContentsOf);

    /** register a method which can be used for XLCoerceToMulti */
    static void setXLCoerceMethods(XLCoerceToMulti* coerceMethod,
                                   XLFreeCoerce*    freeMethod);

    /* if xltypeMulti and num cols * num rows == 1, returns oper.val.xlarray
       else returns origOper */
    static const XL_OPER& coerceOperFromMulti(const XL_OPER& origOper);

    static void setExcelDateOffsetSet(int xlDateOffset);

protected:

    /** convert XLstring in OPER to C style string by filling up a
        charbuffer - no checking of type of OPER */
    static void operStringToMemBuff(const XL_OPER&  oper,
                                    char  buffer[EXCEL_MAX_STR_LEN]);

    /* throw exception for unrecognised oper type */
    static void unknownOper(const string&   routine,  /* (I) */
                            const XL_OPER&  oper);    /* (I) */

    /** Populate oper given by supplied offsets with supplied string */
    static void populateOperFromString(const string& value,
                                       int           xPos,
                                       int           yPos,
                                       XL_OPER&      oper);

    /** Check the validity of the supplied OPER. The return value
        indicates whether processing should continue for this object ie
        true everything is okay, false: silently return N/A to the
        spreadsheet */
    static bool validateOPER(
        const XL_OPER&   oper,        /* (I) */
        CClassConstSP    desiredType, // (I) required type
        int              paramNum,    /* (I) which parameter we're on */
        XL_OPER&         output);     /* (M) error code written here on failure */

    /* Returns the length of the array of opers. The count finishes on
       the cell before the first nil cell. Counting is performed by
       going down vertically */
    static int arrayLength(
        bool           emptyStringAsNil, /* (I) treat "" as xltypeNil */
        const XL_OPER& oper);            /* (I) of type xltypeMulti */
    
    /* Returns the number of rows in the grid of opers. The count
       finishes on the cell before the first nil cell. Counting is
       performed by going down vertically on the LHS column. */
    static int gridLength(
        bool            emptyStringAsNil,  /* (I) treat "" as xltypeNil */
        const XL_OPER&  oper);             /* (I) of type xltypeMulti */

    /* fill in tmpOper to be of type xltypeMulti and point to origOper */
    static void coerceOperToMulti(const XL_OPER& origOper, /* (I) */
                                  XL_OPER&       tmpOper);  /* (M) */

    static int xlAddinDateOffset;
private:
    /** Increase the number of rows in an oper by the number given. New
        opers have all bits set to 0 */
    static void increaseRowsInOper(int       numRows,
                                   XL_OPER&  oper);
                                 
    static IObjectSP operToObject(
        CClassConstSP  desiredType,   // what type we want
        const XL_OPER& oper);         // oper to convert

    static XLCoerceToMulti*   xlCoerceToMulti;
    static XLFreeCoerce*      xlFreeCoerce;

};

typedef vector<const XLConvert*> XLConvertArray;

DRLIB_END_NAMESPACE

/** Free function for XL Oper - must be a global static method. Note:
    defined in xlregister.cpp to allow EAS to link old and new libraries */
extern "C" 
{
    EDG_XLFUNC(void) xlAutoFree (XL_OPER *x);
}


#endif
