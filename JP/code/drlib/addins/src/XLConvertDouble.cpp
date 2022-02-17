//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertDouble.cpp
//
//   Description : Class for converting doubles to/from Excel
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


/** Specialisation of XLConvert for handling doubles */
class XLConvertDouble: public XLConvert{
public:

    /** invoke addin function which has a signature specific to the instance
        of this class. Result is passed via the jvalue structure */
    virtual void invoke(
        const Addin::Method& method,        // the method to invoke
        const jobject&       addin,         // contains addin arguments
        jvalue&              output) const{ // the return value

        try{
            output.a.d = method.doubleMethod->execute(addin.get());
        } catch (exception& e){
            throw ModelException(e, "XLConvertDouble::invoke");
        }
    }
    
    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        field->setDouble(output, operToDouble(oper));
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

        CDoubleSP db = CDoubleSP::dynamicCast(object);
        jvalue  data;
        data.a.d = db->doubleValue();
        convertNativeOutput(handleName, data, idx, output);
    }


    /** Take relevant value from data and fill column idx in output 
        accordingly. Default is to create a handle */
    virtual void convertNativeOutput(
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jvalue&  data,          // value to be put into OPER 
        int            idx,           // column index number */
        XL_OPER&       oper) const {  // converted output

        makeRoomForSingleOper(idx, oper);
        XL_OPER *output = oper.type == xltypeMulti? oper.val.xlarray.lparray + 
            offset(idx, 0,  oper.val.xlarray.columns, oper.val.xlarray.rows):
            &oper;
        output->val.num = data.a.d;
        output->type = xltypeNum;
        blankColumn(1, idx, oper);
    }

    /** convert value held in data to an object. */
    virtual void nativeToObject(jvalue&  data) const{
        data.l = IObjectSP(CDouble::create(data.a.d));
    }
   
    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDouble::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertDoubleRegister(){
    XLConvertFactory::registerXLObjectConvert(CDouble::TYPE, 
                                              new XLConvertDouble());
}

DRLIB_END_NAMESPACE
