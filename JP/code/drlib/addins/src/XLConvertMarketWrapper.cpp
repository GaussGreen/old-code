//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertMarketWrapper.cpp
//
//   Description : Class for converting to/from MarketObjectWrappers
//
//   Author      : Mark A Robson
//
//   Date        : 29 Mar 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/MarketObject.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling market object wrappers */
class XLConvertMarketWrapper: public XLConvert{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertMarketWrapper::convertInput";
        try{
            IObjectSP    object;
            if (oper.type & xltypeMulti){
                throw ModelException(routine,
                                     "Market Object must be supplied in "
                                     "a single cell");
            }
            if (operContainsHandle(oper)){
                // get the market object from the handle - or rather rely
                // on CObject::checkType to get from what is actually supplied
                // to what we want
                object = getObjectFromHandle(oper, desiredType);
            } else {
                // get string in oper
                string name(operToString(oper));
                // need to convert to appropriate market object wrapper type
                object = IObjectSP(MarketObjectWrapper::
                                   createFromString(desiredType, name));
            }
            field->set(output, object);
        } catch (exception& e){
            throw ModelException(e, routine);
        }
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

        makeRoomForSingleOper(idx, output);
        MarketObjectWrapper& wrapper =
            *DYNAMIC_CAST(MarketObjectWrapper, object.get());
        if (wrapper.usingCache()){
            // just return name of market object to sheet
            string name = wrapper.getName();
            populateOperFromString(name, idx, 0, output);
        } else {
            // return a handle to the pointer
            IObjectSP marketObject(wrapper.getMO());
            createHandle(handleName, name, marketObject, idx, -1, output); 
        }
        blankColumn(1, idx, output);
    }

    bool isSimpleType() const {
        return false;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    MarketObjectWrapper::TYPE is initialised. Outside of class to
    avoid necessity of header file */
void XLConvertMarketWrapperRegister(){
    XLConvertFactory::registerXLObjectConvert(MarketObjectWrapper::TYPE, 
                                              new XLConvertMarketWrapper());
}

DRLIB_END_NAMESPACE
