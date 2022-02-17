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
class XLConvertMarketWrapperArray: public XLConvert{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = 
            "XLConvertMarketWrapperArray::convertInput";
        try{
            XL_OPER         tmpOper;
            const XL_OPER*  operToUse;
            if (!(oper.type & xltypeMulti)){
                coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
                operToUse = &tmpOper;
            } else {
                operToUse = &oper;
            }
            CClassConstSP eltType = desiredType->getComponentType();
            int numElts = arrayLength(true, /* stop at first blank cell */
                                      *operToUse);
            if (numElts == 1){
                // support people passing a StringArray via a handle. You can
                // of course just pass in the strings directly. Note the case
                // when a single handle of the right type is supplied is
                // already supported by the calling function
                XL_OPER& cell = operToUse->val.xlarray.lparray[0];
                if (operContainsHandle(cell)){
                    IObjectSP object = getObjectFromHandle(cell,
                                                           0 /* any type */);
                    if (object->getClass()->isArray()){
                        CObject::checkType(object, desiredType);
                        field->set(output, object);
                        return; // all done
                    }
                }
            }

            IArraySP arr(desiredType->newArrayInstance(numElts));
            for (int idx = 0; idx < numElts; idx++){
                IObjectSP    object;
                XL_OPER& cell = operToUse->val.xlarray.lparray[idx];
                if (operContainsHandle(cell)){
                    // get the market object from the handle (rely on
                    // CObject::checkType to do any required conversions)
                    object =
                        getObjectFromHandle(cell,
                                            desiredType->getComponentType());
                } else {
                    // get string in oper
                    string name(operToString(cell));
                    // need to convert to appropriate marketObject wrapper type
                    object = IObjectSP(MarketObjectWrapper::
                                       createFromString(eltType, name));
                }
                arr->set(idx, object); // place in array
            }
            // the put array into output
            field->set(output, arr);
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
        XL_OPER&       oper) const{ // converted output

        static const string routine = 
            "XLConvertMarketWrapperArray::convertObjectOutput";
        try{
            IArray& arr = *DYNAMIC_CAST(IArray, object.get());
            int numElts = object.get()? arr.getLength(): 0;
            // ensure our output oper is big enough
            setNumberOfRowsInOper(numElts, idx, oper);
            for (int jdx = 0; jdx < numElts; jdx++){
                // pull out the current element out of the array
                IObjectSP elt =  arr.get(jdx);
                MarketObjectWrapper& wrapper =
                    *DYNAMIC_CAST(MarketObjectWrapper, elt.get());
                if (wrapper.usingCache()){
                    // just return name of market object to sheet
                    string name = wrapper.getName();
                    populateOperFromString(name, idx, jdx, oper);
                } else {
                    // return a handle to the pointer
                    IObjectSP marketObject(wrapper.getMO());
                    createHandle(handleName, name, marketObject, idx, 
                                 jdx, oper); 
                }
            }
            blankColumn(numElts, idx, oper);
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    bool isSimpleType() const {
        return false;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    MarketObjectWrapper::TYPE is initialised. Outside of class to
    avoid necessity of header file */
void XLConvertMarketWrapperArrayRegister(){
    /* need one of the function calls for each wrapper array (above code is
       generic) */
    XLConvertFactory::registerXLArrayConvert(MarketObjectWrapper::TYPE, 
                                             new XLConvertMarketWrapperArray());
}

DRLIB_END_NAMESPACE
