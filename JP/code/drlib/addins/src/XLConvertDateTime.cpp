//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertDateTime.cpp
//
//   Description : Class for converting dateTimes to/from Excel
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


/** Specialisation of XLConvert for handling dateTimes */
class XLConvertDateTime: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertDateTime::convertInput";
        bool         invalidSize = false;
        DateTimeSP   dateTime;
        if (!(oper.type & xltypeMulti)){
            if (operContainsHandle(oper)){
                // get object from handle
                IObjectSP obj = getObjectFromHandle(oper, DateTime::TYPE);
                dateTime = DateTimeSP::dynamicCast(obj);
            } else {
                invalidSize = true;
            }
        } else if (oper.val.xlarray.rows * oper.val.xlarray.columns != 2){
            invalidSize = true;
        } else {
            dateTime = DateTimeSP(
                new DateTime(opersToDateTime(oper.val.xlarray.lparray[0], 
                                             oper.val.xlarray.lparray[1])));
        }
        if (invalidSize){
            throw ModelException(routine, "DateTimes must be supplied in a "
                                 "2 x 1 range");
        }
        field->set(output, dateTime);
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
        DateTimeSP dateTime = DateTimeSP::dynamicCast(object);
        xlDateTimeFromDateTime(*dateTime, offset(idx, 0, output), 
                               offset(idx+1, 0, output));
        blankColumn(1, idx, output);   // blank out any cells below us
        blankColumn(1, idx+1, output); // blank out any cells below us
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

/** Specialisation of XLConvert for handling dates */
class XLConvertDate: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertDate::convertInput";
        int dateAsInt = operToDate(oper);
        DateTime::Date* date = new DateTime::Date(dateAsInt);
        IObjectSP dateObj(date);
        field->set(output, dateObj);
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

        DateTime::Date& dateObj = *DYNAMIC_CAST(DateTime::Date, object.get());
        makeRoomForSingleOper(idx, oper);
        XL_OPER *output = oper.type == xltypeMulti? oper.val.xlarray.lparray + 
            offset(idx, 0,  oper.val.xlarray.columns, oper.val.xlarray.rows):
            &oper;
        xlDateFromDate(dateObj.getDate(), *output);
        blankColumn(1, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }    
    
};

/** Specialisation of XLConvert for handling times */
class XLConvertTime: public XLConvertUtils{
public:

    /** Parse OPER and place data in given field in output object */
    virtual void convertInput(
        CClassConstSP  desiredType,        // what type we want
        const XL_OPER&       oper,         // oper to convert
        const CFieldConstSP& field,        // where to store the converted data
        jobject&             output) const{ // store value in here using field

        static const string routine = "XLConvertTime::convertInput";
        int timeAsInt = operToTime(oper);
        DateTime::Time* time = new DateTime::Time(timeAsInt);
        IObjectSP timeObj(time);
        field->set(output, timeObj);
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

        makeRoomForSingleOper(idx, oper);
        DateTime::Time& timeObj = *DYNAMIC_CAST(DateTime::Time, object.get());
        XL_OPER *output = oper.type == xltypeMulti? oper.val.xlarray.lparray + 
            offset(idx, 0,  oper.val.xlarray.columns, oper.val.xlarray.rows):
            &oper;
        xlTimeFromTime(timeObj.getTime(), *output);
        blankColumn(1, idx, oper);
    }

    bool isSimpleType() const {
        return true;
    }
};

/** non static method which is used by XLConvertFactory to ensure that
    this class gets linked in and the method is called after
    CDateTime::TYPE is initialised. Outside of class to avoid necessity
    of header file */
void XLConvertDateTimeRegister(){
    XLConvertFactory::registerXLObjectConvert(DateTime::TYPE, 
                                              new XLConvertDateTime());
    XLConvertFactory::registerXLObjectConvert(DateTime::Date::TYPE, 
                                              new XLConvertDate());
    XLConvertFactory::registerXLObjectConvert(DateTime::Time::TYPE, 
                                              new XLConvertTime());
}

DRLIB_END_NAMESPACE
