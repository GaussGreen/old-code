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
#include "edginc/XLConvertUtils.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling expiry arrays */
class XLConvertExpiryArray: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertExpiryArray::convertInput";
        // count the number of opers in the array - stop at first nil cell
        int numElts = oper.type & xltypeMulti? gridLength(true, oper): 1;
        ExpiryArraySP  array(new ExpiryArray(numElts));
        ExpiryArray&   expiryArray = *array;
        if (!(oper.type & xltypeMulti)){
            XL_OPER oper2;
            oper2.type = xltypeNil;
            expiryArray[0] = opersToExpiry(oper, oper2);
        } else {
            // for ease, cache some variables
            unsigned short cols = oper.val.xlarray.columns;
            unsigned short rows = oper.val.xlarray.rows;
            XL_OPER *arrayStart = oper.val.xlarray.lparray;
            if (cols == 1){
                // if a single column treat as blank 2nd column
                XL_OPER oper2;
                oper2.type = xltypeNil;
                for (int idx = 0; idx < numElts; idx++){
                    XL_OPER  *oper1 = arrayStart + idx;
                    expiryArray[idx] = opersToExpiry(*oper1, oper2);
                }
            } else if (cols != 2){
                static const string m("Expiries must be supplied as "
                                      "in a column of width 1 or 2");
                throw ModelException(routine, m);
            } else {
                for (int idx = 0; idx < numElts; idx++){
                    XL_OPER  *oper1 = arrayStart + offset(0, idx, cols, rows);
                    XL_OPER  *oper2 = arrayStart + offset(1, idx, cols, rows);
                    expiryArray[idx] = opersToExpiry(*oper1, *oper2);
                }
            }
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

        ExpiryArraySP array = ExpiryArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        ExpiryArray&   expiryArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER  *oper1 = arrayStart + offset(idx, jdx, cols, rows);
            XL_OPER  *oper2 = arrayStart + offset(idx+1, jdx, cols, rows);
            opersFromExpiry(expiryArray[jdx], *oper1, *oper2);
        }
        blankColumn(numElts, idx, oper);
        blankColumn(numElts, idx+1, oper);
    }

    /** Returns the width needed to dsiplay this type of object in
        XL. A return value of -1 indicates a variable width. This is
        used for addins that use Addin:expandMulti */
    virtual int width() const{
        // display expiries using 2 columns
        return 2;
    }

    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CExpiry::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertExpiryArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        Expiry::TYPE, new XLConvertExpiryArray());
}

DRLIB_END_NAMESPACE
