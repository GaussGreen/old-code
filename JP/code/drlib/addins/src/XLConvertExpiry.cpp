//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertExpiry.cpp
//
//   Description : Class for converting expiries to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 1 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/XLConvertUtils.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling expiries */
class XLConvertExpiry: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertExpiry::convertInput";
        IObjectSP  expiry;
        if (oper.type != xltypeMulti){
            XL_OPER oper2;
            oper2.type = xltypeNil;
            expiry = opersToExpiry(oper, oper2);
        } else {
            if (oper.val.xlarray.rows * oper.val.xlarray.columns != 2){
                throw ModelException(routine, "Expirys must be supplied in a "
                                     "single cell or 2 x 1 range");
            }
            expiry = opersToExpiry(oper.val.xlarray.lparray[0],
                                   oper.val.xlarray.lparray[1]);
        }
        field->set(output, expiry);
    }

    /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       output) const{ // converted output

        // set number of rows in output to at least one
        setNumberOfRowsInOper(1, idx, output);
        opersFromExpiry(object, offset(idx, 0, output), 
                        offset(idx+1, 0, output));
        blankColumn(1, idx, output);   // blank out any cells below us
        blankColumn(1, idx+1, output); // blank out any cells below us
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
void XLConvertExpiryRegister(){
    XLConvertFactory::registerXLObjectConvert(Expiry::TYPE, 
                                              new XLConvertExpiry());
}

DRLIB_END_NAMESPACE
