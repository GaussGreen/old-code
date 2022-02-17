//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertStringArray.cpp
//
//   Description : Class for converting string arrays to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
// keep these three includes above in this order otherwise the compiler gets
// confused by the array template and the array struct in xlapi.hpp
#include "edginc/Handle.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling string arrays */
class XLConvertStringArray: public XLConvert{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{// store value in here using field
        
        XL_OPER           tmpOper;
        const XL_OPER*    operToUse;
        if (!(oper.type & xltypeMulti)){
            coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
            operToUse = &tmpOper;
        } else {
            operToUse = &oper;
        }

        int numElts = arrayLength(true, /* stop at first blank cell */
                                  *operToUse);
        CStringArraySP  array(new StringArray(numElts));
        StringArray&   stringArray = *array;
        for (int idx = 0; idx < numElts; idx++){
            stringArray[idx] = operToString(operToUse->val.xlarray.lparray[idx],
                                            true);
        }
        // set in output
        field->set(output, array);
    }

    /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       oper) const{   // converted output

        CStringArraySP array = CStringArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        StringArray&   stringArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER* currentOper = arrayStart + 
                offset(idx, jdx, cols, rows);
            populateOperFromString(stringArray[jdx], *currentOper);
        }
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};

/** Specialisation of XLConvert for handling 'raw' string arrays ie
    strings that will not be processed as handles even if they are */
class XLConvertRawStringArray: public XLConvert{
public:
    // need to override this to avoid handle look up
    virtual void convertGenericInput(
        CClassConstSP  desiredType,   // what type we want
        const XL_OPER& oper,          // oper to convert
        const CFieldConstSP& field,   // where to store the converted data
        jobject&       output) const{ // store value in here using field
        
        convertInput(desiredType, oper, field, output);
    }

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        XL_OPER           tmpOper;
        const XL_OPER*    operToUse;
        if (!(oper.type & xltypeMulti)){
            coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
            operToUse = &tmpOper;
        } else {
            operToUse = &oper;
        }

        int numElts = arrayLength(true, /* stop at first blank cell */
                                  *operToUse);
        Handle::RawStringArraySP  arr(new Handle::RawStringArray(numElts));
        Handle::RawStringArray&   stringArray = *arr;
        for (int idx = 0; idx < numElts; idx++){
            string val = operToString(operToUse->val.xlarray.lparray[idx],
                                      false);
            stringArray[idx] = Handle::RawStringSP(new Handle::RawString(val));
        }
        // set in output
        field->set(output, arr);
    }

   /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       oper) const{   // converted output

        Handle::RawStringArraySP arr = 
            Handle::RawStringArraySP::dynamicCast(object);
        int numElts = arr.get()? arr->size(): 0;
        Handle::RawStringArray&   stringArray = *arr;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER* currentOper = arrayStart + 
                offset(idx, jdx, cols, rows);
            populateOperFromString(stringArray[jdx]->getString(), 
                                   *currentOper);
        }
        blankColumn(numElts, idx, oper);
     }

};


/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CString::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertStringArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        CString::TYPE, new XLConvertStringArray());
    XLConvertFactory::registerXLArrayConvert(
        Handle::RawString::TYPE, new XLConvertRawStringArray());
}

DRLIB_END_NAMESPACE
