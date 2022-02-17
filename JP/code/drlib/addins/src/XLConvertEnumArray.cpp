//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Class for converting enum arrays to/from Excel
//
//   Author      : Mark A Robson
//
//   Date        : 14 Sep 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
// keep these three includes above in this order otherwise the compiler gets
// confused by the array template and the array struct in xlapi.hpp
#include "edginc/Handle.hpp"
DRLIB_BEGIN_NAMESPACE


/** Specialisation of XLConvert for handling enum arrays */
class XLConvertEnumArray: public XLConvert{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{// store value in here using field
        
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
        // type varies depending on particular flavour of enum, so need to
        // use generic methods
        IArraySP  enumArray(desiredType->newArrayInstance(numElts));
        CClassConstSP enumType = desiredType->getComponentType();
        for (int idx = 0; idx < numElts; idx++){
            const string& enumAsString = 
                operToString(operToUse->val.xlarray.lparray[idx], true);
            IObjectSP enumAsObject(Enum::create(enumType, enumAsString));
            enumArray->set(idx, enumAsObject);
        }
        // set in output
        field->set(output, enumArray);
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

        IArraySP array = IArraySP::dynamicCast(object);
        int numElts = array.get()? array->getLength(): 0;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER* currentOper = arrayStart + 
                offset(idx, jdx, cols, rows);
            EnumSP elt(EnumSP::dynamicCast(array->get(jdx)));
            const string& enumAsString = elt->enumValueAsString();
            populateOperFromString(enumAsString, *currentOper);
        }
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};


/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CString::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertEnumArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        Enum::TYPE, new XLConvertEnumArray());
}

DRLIB_END_NAMESPACE
