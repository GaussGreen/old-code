//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertDateTimeArray.cpp
//
//   Description : Class for converting DateTime arrays to/from Excel
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


/** Specialisation of XLConvert for handling dateTime arrays */
class XLConvertDateTimeArray: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine("XLConvertDateTimeArray::convertInput");
        static const string m("Date-times must be supplied as "
                              "handles or in a column of width 2");
        XL_OPER             tmpOper;
        const XL_OPER*      operToUse = &oper; // default
        DateTimeArraySP  array;

        if (!(oper.type & xltypeMulti)){
            coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
            operToUse = &tmpOper;
        }

        // count the number of opers in the array - stop at first nil cell
        int numElts = gridLength(true, *operToUse);
        array = DateTimeArraySP(new DateTimeArray(numElts));
        DateTimeArray&   dateTimeArray = *array;  // for ease

        if (operToUse->val.xlarray.columns == 1){
            /* support for array of handles to date-times */
            for (int i = 0; i < numElts; i++){
                XL_OPER& cell = operToUse->val.xlarray.lparray[i];
                if (!(cell.type & xltypeStr)){
                    throw ModelException(routine, m);
                }
                // get object from handle
                IObjectSP obj = getObjectFromHandle(cell, DateTime::TYPE);
                // cast to DateTime
                DateTime& theDate = *DYNAMIC_CAST(DateTime, obj.get());
                // put into array
                dateTimeArray[i] = theDate;
            }
        } else if (operToUse->val.xlarray.columns != 2){
            throw ModelException(routine, m);
        } else {
            // for ease, cache some variables
            unsigned short cols = operToUse->val.xlarray.columns;
            unsigned short rows = operToUse->val.xlarray.rows;
            XL_OPER *arrayStart = operToUse->val.xlarray.lparray;
            for (int idx = 0; idx < numElts; idx++){
                XL_OPER  *oper1 = arrayStart + offset(0, idx, cols, rows);
                XL_OPER  *oper2 = arrayStart + offset(1, idx, cols, rows);
                dateTimeArray[idx] = opersToDateTime(*oper1, *oper2);
            }
        }
        // set in output
        field->set(output, array);
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

        DateTimeArraySP array = DateTimeArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        DateTimeArray&   dateTimeArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER  *oper1 = arrayStart + offset(idx, jdx, cols, rows);
            XL_OPER  *oper2 = arrayStart + offset(idx+1, jdx, cols, rows);
                xlDateTimeFromDateTime(dateTimeArray[jdx], *oper1, *oper2);
        }
        blankColumn(numElts, idx, oper);
        blankColumn(numElts, idx+1, oper);
    }

    /** Returns the width needed to dsiplay this type of object in
        XL. A return value of -1 indicates a variable width. This is
        used for addins that use Addin:expandMulti */
    virtual int width() const{
        // display dates using 2 columns
        return 2;
    }

    bool isSimpleType() const {
        return true;
    }

};

/** Specialisation of XLConvert for handling dateTime arrays */
class XLConvertDateArray: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        XL_OPER                tmpOper;
        const XL_OPER*         operToUse = &oper; // default
        DateTime::DateArraySP  array;

        if (!(oper.type & xltypeMulti)){
            coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
            operToUse = &tmpOper;
        }

        // count the number of opers in the array - stop at first nil cell
        int numElts = arrayLength(true, *operToUse);
        array = DateTime::DateArraySP(new DateTime::DateArray(numElts));
        DateTime::DateArray&   dateArray = *array;  // for ease

        for (int i = 0; i < numElts; i++){
            XL_OPER& cell = operToUse->val.xlarray.lparray[i];
            int dateAsInt = operToDate(cell);
            DateTime::Date* date = new DateTime::Date(dateAsInt);
            dateArray[i] = DateTime::DateSP(date);
        }
        // set in output
        field->set(output, array);
    }

    /** Take relevant value from data and fill column idx in output 
        accordingly.  */
    virtual void convertObjectOutput(
        const string&  handleName,    // for handles
        const string&  name,          /* name for this bit of data
                                         (for handles/error messages) */
        const jobject& object,        // populate OPER with the object
        int            idx,           // column index number */
        XL_OPER&       oper) const{   // converted output

        DateTime::DateArraySP array = 
            DateTime::DateArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        DateTime::DateArray&   dateArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER  *currentOper = arrayStart + 
                offset(idx, jdx, cols, rows);
            xlDateFromDate(dateArray[jdx]->getDate(), *currentOper);
        }
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};

/** Specialisation of XLConvert for handling dateTime arrays */
class XLConvertTimeArray: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        XL_OPER                tmpOper;
        const XL_OPER*         operToUse = &oper; // default
        DateTime::TimeArraySP  array;

        if (!(oper.type & xltypeMulti)){
            coerceOperToMulti(oper, tmpOper); /* coerce to xltypeMulti */
            operToUse = &tmpOper;
        }

        // count the number of opers in the array - stop at first nil cell
        int numElts = arrayLength(true, *operToUse);
        array = DateTime::TimeArraySP(new DateTime::TimeArray(numElts));
        DateTime::TimeArray&   timeArray = *array;  // for ease

        for (int i = 0; i < numElts; i++){
            XL_OPER& cell = operToUse->val.xlarray.lparray[i];
            int timeAsInt = operToTime(cell);
            DateTime::Time* time = new DateTime::Time(timeAsInt);
            timeArray[i] = DateTime::TimeSP(time);
        }
        // set in output
        field->set(output, array);
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

        DateTime::TimeArraySP array = 
            DateTime::TimeArraySP::dynamicCast(object);
        int numElts = array.get()? array->size(): 0;
        DateTime::TimeArray&   timeArray = *array;
        // ensure our output oper is big enough
        setNumberOfRowsInOper(numElts, idx, oper);
        // for ease, cache some variables
        unsigned short cols = oper.val.xlarray.columns;
        unsigned short rows = oper.val.xlarray.rows;
        XL_OPER *arrayStart = oper.val.xlarray.lparray;
        for (int jdx = 0; jdx < numElts; jdx++){
            XL_OPER  *currentOper = arrayStart + 
                offset(idx, jdx, cols, rows);
            xlTimeFromTime(timeArray[jdx]->getTime(), *currentOper);
        }
        blankColumn(numElts, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDateTime::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertDateTimeArrayRegister(){
    XLConvertFactory::registerXLArrayConvert(
        DateTime::TYPE, new XLConvertDateTimeArray());
    XLConvertFactory::registerXLArrayConvert(
        DateTime::Date::TYPE, new XLConvertDateArray());
    XLConvertFactory::registerXLArrayConvert(
        DateTime::Time::TYPE, new XLConvertTimeArray());
}

DRLIB_END_NAMESPACE
