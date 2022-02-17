//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertDoubleArray.cpp
//
//   Description : Class for converting double arrays to/from Excel
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


/** Specialisation of XLConvert for handling double arrays */
class XLConvertDoubleArray: public XLConvert{
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
        DoubleArraySP  array(new DoubleArray(numElts));
        DoubleArray&   doubleArray = *array;
        for (int idx = 0; idx < numElts; idx++){
            doubleArray[idx] = operToDouble(operToUse->val.xlarray.lparray[idx]);
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

        DoubleArraySP array = DoubleArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        DoubleArray&   dbArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER  *currentOper = arrayStart + offset(idx, jdx, cols, rows);
            currentOper->val.num = dbArray[jdx];
            currentOper->type = xltypeNum;
        }
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }

};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDouble::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertDoubleArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        CDouble::TYPE, new XLConvertDoubleArray());
}

DRLIB_END_NAMESPACE
