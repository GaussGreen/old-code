//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertString.cpp
//
//   Description : Class for converting strings to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 24 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/Handle.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling strings */
class XLConvertString: public XLConvert{
public:

    /** invoke addin function which has a signature specific to the instance
        of this class. Result is passed via the jvalue structure */
    virtual void invoke(
        const Addin::Method& method,        // the method to invoke
        const jobject&       addin,         // contains addin arguments
        jvalue&              output) const{ // the return value

        try{
            string str = method.stringMethod->execute(addin.get());
            output.l = IObjectSP(CString::create(str));
        } catch (exception& e){
            throw ModelException(e, "XLConvertString::invoke");
        }
    }
    
    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        string val = operToString(oper);
        field->setString(output, val);
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
        CStringSP s = CStringSP::dynamicCast(object);
        const string& str = s->stringValue();
        XL_OPER *output = oper.type == xltypeMulti? oper.val.xlarray.lparray + 
            offset(idx, 0,  oper.val.xlarray.columns, oper.val.xlarray.rows):
            &oper;
        XLConvert::populateOperFromString(str, *output);
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

/** Specialisation of XLConvert for handling 'raw' strings ie strings that
    will not be processed as handles even if they are */
class XLConvertRawString: public XLConvert{
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

        string val = operToString(oper, false);
        IObjectSP rawString(new Handle::RawString(val));
        field->set(output, rawString);
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
        Handle::RawString& s = *DYNAMIC_CAST(Handle::RawString, object.get());
        const string& str = s.getString();
        XL_OPER *output = oper.type == xltypeMulti? oper.val.xlarray.lparray + 
            offset(idx, 0,  oper.val.xlarray.columns, oper.val.xlarray.rows):
            &oper;
        XLConvert::populateOperFromString(str, *output);
        blankColumn(1, idx, oper);
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDouble::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertStringRegister(){
    XLConvertFactory::registerXLObjectConvert(
        CString::TYPE, new XLConvertString());
    XLConvertFactory::registerXLObjectConvert(Handle::RawString::TYPE, 
                                              new XLConvertRawString());
}


DRLIB_END_NAMESPACE
