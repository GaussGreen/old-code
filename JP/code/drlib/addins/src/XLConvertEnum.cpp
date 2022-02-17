//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Class for converting enums to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 14 Sep 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/Handle.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling strings */
class XLConvertEnum: public XLConvert{
public:

    /** invoke addin function which has a signature specific to the instance
        of this class. Result is passed via the jvalue structure */
    virtual void invoke(
        const Addin::Method& method,        // the method to invoke
        const jobject&       addin,         // contains addin arguments
        jvalue&              output) const{ // the return value

        try{
#if 0
            string str = method.stringMethod->execute(addin.get());
            output.l = IObjectSP(CString::create(str));
#endif
            // work needed here in Addin.hpp to flag the return as an Enum
            // but also in FunctorMember to convert the enum into a string
            // (or into an Enum - and then don't bother having this method?)
            throw ModelException("XLConvertEnum::invoke", 
                                 "Not yet implemented");
        } catch (exception& e){
            throw ModelException(e, "XLConvertEnum::invoke");
        }
    }
    
    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        string val = operToString(oper);
        // then build Enum
        IObjectSP enumAsObj(Enum::create(desiredType, val));
        // and just use generic field set method, could optimise this if we
        // had a Field::setEnum(int) method
        field->set(output, enumAsObj);
    }

    /** Take relevant value from data and fill column idx in output 
        accordingly. */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       oper) const{  // converted output

        makeRoomForSingleOper(idx, oper);
        EnumSP e = EnumSP::dynamicCast(object);
        const string& enumAsString = e->enumValueAsString();
        XL_OPER *output = oper.type == xltypeMulti? oper.val.xlarray.lparray + 
            offset(idx, 0,  oper.val.xlarray.columns, oper.val.xlarray.rows):
            &oper;
        XLConvert::populateOperFromString(enumAsString, *output);
        blankColumn(1, idx, oper);
    }


    /** Take relevant value from data and fill column idx in output 
        accordingly. */
    virtual void convertNativeOutput(
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jvalue&  data,          // value to be put stringo OPER 
        int            idx,           // column index number */
        XL_OPER&       oper) const {  // converted output
        
        convertObjectOutput(name, "", data.l, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    Enum::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertEnumRegister(){
    XLConvertFactory::registerXLObjectConvert(Enum::TYPE, new XLConvertEnum());
}


DRLIB_END_NAMESPACE
