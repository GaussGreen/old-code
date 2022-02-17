//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertBoolArray.cpp
//
//   Description : Class for converting bool arrays to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling bool arrays */
class XLConvertBoolArray: public XLConvert{
public:

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
        CBoolArraySP  array(new CBoolArray(numElts));
        CBoolArray&   boolArray = *array;
        for (int idx = 0; idx < numElts; idx++){
            boolArray[idx] = operToBool(operToUse->val.xlarray.lparray[idx]);
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

        CBoolArraySP array = CBoolArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        CBoolArray&   boolArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER  *currentOper = arrayStart + offset(idx, jdx, cols, rows);
            currentOper->val.boolVal = boolArray[jdx];
            currentOper->type = xltypeBool;
        }
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CBool::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertBoolArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        CBool::TYPE, new XLConvertBoolArray());
}

DRLIB_END_NAMESPACE
