//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertObjArray.cpp
//
//   Description : Class for converting [object] arrays to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include "edginc/Format.hpp"
#include "edginc/Malloc.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling arrays of objects */
class XLConvertObjArray: public XLConvert{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertObjArray::convertInput";
        try{
            XL_OPER           tmpOper;
            const XL_OPER*    operToUse;
            if (!(oper.type & xltypeMulti)){
                coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
                operToUse = &tmpOper;
            } else {
                operToUse = &oper;
            }
            
            int numElts = arrayLength(false, /* don't stop at "" */
                                      *operToUse);
            CClassConstSP componentType = desiredType->getComponentType();
            // create empty instance of array
            IArraySP  objArray(desiredType->newArrayInstance(numElts));
            IArray&   array = *objArray;
            // then fill in each value
            for (int idx = 0; idx < numElts; idx++){
                // check we've got a handle
                XL_OPER& currentOper = operToUse->val.xlarray.lparray[idx];
                // skip if a string and it's ""
                if (!(currentOper.type & xltypeStr) || 
                    currentOper.val.str[0] != 0)
                {
                    // skip over cells containing "" (will become null)
                    if (!operContainsHandle(currentOper)){
                        // could consider converting cell to an object - to do
                        throw ModelException(routine, "Couldn't read handle "
                                             "from cell "+
                                             Format::toString(idx)+ " in range"
                                             ". Cell contains "+
                                             operToString(currentOper)+
                                             ". For parameters of type "+
                                             componentType->getName()+
                                             " a handle must be supplied");
                    }
                    // get object from handle
                    IObjectSP obj(getObjectFromHandle(currentOper, 
                                                      componentType));
                    // set in array
                    array.set(idx, obj);
                }
            }
            // set in output
            field->set(output, objArray);
        } catch (exception& e){
            throw ModelException(e, routine, "Failed whilst trying to read in"
                                 " array of type "+desiredType->getName());
        }
    }

    /** Take relevant value from data and fill column idx in output 
        accordingly. */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       oper) const{   // converted output

        IArraySP array = IArraySP::dynamicCast(object);
        int numElts = array.get()? array->getLength(): 0;
        if (numElts > 0){
            // ensure our output oper is big enough
            setNumberOfRowsInOper(numElts, idx, oper);
            for (int jdx = 0; jdx < numElts; jdx++){
                // get object out of array
                IObjectSP obj(array->get(jdx));
                if (!obj){
                    // could create handle to CNull object instead?
                    XL_OPER& cell = offset(idx, jdx, oper);
                    cell.val.str = NEW(char);
                    cell.type = xltypeStr;
                } else {
                    if (obj->getRefCount() == 0){
                        // need to be careful if we get reference
                        obj = IObjectSP(obj->clone());
                    }
                    // then create handle to it(and return handle name in oper)
                    createHandle(handleName, name, obj, idx, jdx, oper);
                }
            }
        }
        // then blank out any remaining cells in this column
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return false;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    C::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertObjArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        IObject::TYPE, new XLConvertObjArray());
}

DRLIB_END_NAMESPACE
