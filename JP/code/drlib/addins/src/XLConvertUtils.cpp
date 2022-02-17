//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertUtils.cpp
//
//   Description : Provides some utility methods for dealing with 
//                 EDR specific classes (formerly part of XLConvert)
//
//   Author      : Mark A Robson
//
//   Date        : 25 Nov 2003
//
//
//   $Log: XLConvertUtils.cpp,v $
//   Revision 1.1  2003/11/26 17:26:59  mrobson
//   EDR specific version of XLConvert
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XLConvertUtils.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/EnergyContractLabel.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
DRLIB_BEGIN_NAMESPACE
// for derived classes
XLConvertUtils::XLConvertUtils(){}


int XLConvertUtils::fracOfDayToTime(double  timeFrac){ /* (I) */
    if (timeFrac < 0.0 || timeFrac > 1.0) {
        throw ModelException("XLConvertUtils::FracOfDayToTime",
                             "Fraction of day "+Format::toString(timeFrac)+
                             " must be in [0,1]");
    }
    /* round to nearest time unit */
    return (DateTime::START_OF_DAY_TIME +
            (int)((double)(DateTime::END_OF_DAY_TIME - 
                           DateTime::START_OF_DAY_TIME) * timeFrac + 0.5));
}

/** Converts an oper to a 'time of day' */
int XLConvertUtils::operToTime(const XL_OPER& oper){
    static const string routine = "XLConvertUtils::operToTime";
    int val;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit)
    {
    case xltypeNum:
    case xltypeInt:
    {
        double number = (oper.type & (~xlbitDLLFree)) == xltypeInt?
            oper.val.w: oper.val.num;
        val = fracOfDayToTime(number);
        break;
    }
    case xltypeNil:
    case xltypeMissing:
        val = DateTime::START_OF_DAY_TIME; /* default to start of day */
        break;
    case xltypeBool:
        val = oper.val.boolVal? 
            DateTime::END_OF_DAY_TIME: DateTime::START_OF_DAY_TIME;
        break;
    case xltypeStr:
    {
        double    timeFrac; 
        char     *pos;
        char      buffer[EXCEL_MAX_STR_LEN];
        /* read string into buffer */
        operStringToMemBuff(oper, buffer);
        /* turn to upper case */
        for (pos = buffer; (*pos = toupper(*pos)); pos++); /* empty loop*/
        /* need to replace these with hash defines */
        if (!strcmp(buffer, "AM") || !strcmp(buffer, "S") ||
            !strcmp(buffer, DateTime::START_OF_DAY.c_str())){
            val = DateTime::START_OF_DAY_TIME;
        }
        else if (!strcmp(buffer, "PM") || !strcmp(buffer, "E") ||
                 !strcmp(buffer, DateTime::END_OF_DAY.c_str())) {
            val = DateTime::END_OF_DAY_TIME;
        }
        else if (!strcmp(buffer, DateTime::BEFORE_EX_DIV.c_str())) {
            val = DateTime::BEFORE_EX_DIV_TIME;
        } 
        else if (!strcmp(buffer, DateTime::EX_DIV_BEFORE_START.c_str())) {
            val = DateTime::EX_DIV_BEFORE_START_TIME;
        } else {
            char *endPtr;
            timeFrac = strtod(buffer, &endPtr);
            if (*endPtr != '\0') {
                throw ModelException(routine, "Unknown time of day ("+
                                     string(buffer)+")");
            }
            val = fracOfDayToTime(timeFrac);
        }
        break;
    }
    case xltypeMulti:
        throw ModelException(routine, "A single cell indicating time of day"
                             " was expected but a range was suppplied");
    default:
        unknownOper(routine, oper);
        val = 0; // shut compiler up
    }
    return val;
}

/** Converts an oper to an integer date */
int XLConvertUtils::operToDate(const XL_OPER& oper){
    static const string routine = "XLConvertUtils::operToDate";
    int          val = 0;;
    switch (oper.type & (~xlbitDLLFree)) // remove xlbitDLLFree bit
    {
    case xltypeNum:
    case xltypeInt:
        val = (oper.type & (~xlbitDLLFree)) == xltypeInt?
            oper.val.w: (int) oper.val.num; 
        if (val < 0){
            throw ModelException(routine, "Date "+
                                 Format::toString(val)+" is < 0");
        }
        if (!Maths::equals((double)val, oper.val.num)) {
            throw ModelException(routine, "Cannot convert "+
                                 Format::toString(oper.val.num)+ " to a date");
        }
        val += xlAddinDateOffset;
        break;
    case xltypeNil:
    case xltypeMissing:
        throw ModelException(routine, "Date is missing");
    case xltypeBool:
        throw ModelException(routine, "Incorrect type (Boolean) "
                             "supplied for a date");
    case xltypeStr:
    {
        char      buffer[EXCEL_MAX_STR_LEN];
        /* read string into buffer */
        operStringToMemBuff(oper, buffer);
        DateTime dateTime = DateTime(buffer, DateTime::START_OF_DAY);
        val = dateTime.getDate();
        break;
    }
    case xltypeMulti:
        throw ModelException(routine, "A single date was expected but a "
                             "range was suppplied");
    default:
        unknownOper(routine, oper);
    }

    return val;
}

/** Converts an oper to a DateTime */
DateTime XLConvertUtils::opersToDateTime(const XL_OPER& dateOper,
                                         const XL_OPER& timeOper){
    static const string routine = "XLConvertUtils::opersToDateTime";
    int date = operToDate(dateOper);
    int time = operToTime(timeOper);
    return DateTime(date, time);
}

/* construct the xl representation of the date part of a DateTime */
void XLConvertUtils::xlDateFromDate(int    date,
                                    XL_OPER&  oper){
    double xlDate = (double)(date - xlAddinDateOffset);
    /* don't return dates before excel date origin - mucks up vb etc */
    oper.val.num = xlDate < 0.0? 0.0: xlDate;
    oper.type = xltypeNum;
}

/* construct the xl representation of the time part of a DateTime */
void XLConvertUtils::xlTimeFromTime(int       time, /* (I) */
                                    XL_OPER&     oper) /* (O) xlType is set too */
{
    if (time == 0){ // needs to be hash define
        populateOperFromString(DateTime::START_OF_DAY, oper);
    } else if (time == DateTime::START_OF_DAY_TIME){
        populateOperFromString(DateTime::START_OF_DAY, oper);
    } else if (time == DateTime::END_OF_DAY_TIME){
        populateOperFromString(DateTime::END_OF_DAY, oper);
    } else if (time == DateTime::BEFORE_EX_DIV_TIME){
        populateOperFromString(DateTime::BEFORE_EX_DIV, oper);
    } else if (time == DateTime::EX_DIV_BEFORE_START_TIME){
        populateOperFromString(DateTime::EX_DIV_BEFORE_START, oper);
    } else {
        oper.type = xltypeNum;
        oper.val.num = (double)(time - DateTime::START_OF_DAY_TIME)/
            (double)(DateTime::END_OF_DAY_TIME - DateTime::START_OF_DAY_TIME);
    }
}

/* construct the xl representation of a DateTime */
void XLConvertUtils::xlDateTimeFromDateTime(const DateTime& dateTime,
                                            XL_OPER&        dateOper,
                                            XL_OPER&        timeOper){
    xlDateFromDate(dateTime.getDate(), dateOper);
    xlTimeFromTime(dateTime.getTime(), timeOper);
}

/** indicates if the oper contains a maturity rather than a date say. Here
    a matruity is defined as <numbers><letter> eg 3M, 10Y etc. A return
    value of true guarantees that the oper type is xltypeStr */
bool XLConvertUtils::operIsMaturity(const XL_OPER& oper){
    if (oper.type & xltypeStr){
        // may want to make this more sophisticated - need to add
        // relevant function to MaturityPeriod.cpp (which should check
        // upon creation)
        return true;
    }
    return false;
}

/** Constructs expiry object from supplied opers */
ExpirySP XLConvertUtils::opersToExpiry(const XL_OPER& oper1,
                                       const XL_OPER& oper2){
    static const string routine = "XLConvertUtils::opersToExpiry";
    ExpirySP expiry;
    try{
        // if the first cell is something like 3M then we've got some kind of
        // maturity
        if (operIsMaturity(oper1)){
            char  buffer[EXCEL_MAX_STR_LEN];
            operStringToMemBuff(oper1, buffer);
            string matPeriod = buffer;
            if (oper2.type & xltypeNil){
                // MaturityPeriod type
                expiry = ExpirySP(new MaturityPeriod(matPeriod));
            } else {
                int time = operToTime(oper2);
                expiry = ExpirySP(new MaturityTimePeriod(matPeriod, time));
            }
        } else {
            // must be a date
            DateTime date = opersToDateTime(oper1, oper2);
            expiry = ExpirySP(new BenchmarkDate(date));
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
    return expiry;
}

/** populate the given opers using the supplied expiryObj. Will fail
    if expiryObj is not a recognised expiry */
void XLConvertUtils::opersFromExpiry(const IObjectConstSP& expiryObj,
                                     XL_OPER&              oper1,
                                     XL_OPER&              oper2){
    CClassConstSP expiryType = expiryObj->getClass();
    // switch on type of supplied expiry - not very oo but would require
    // an extra class for each one
    if (expiryType == MaturityPeriod::TYPE){
        MaturityPeriodConstSP expiry = 
            MaturityPeriodConstSP::dynamicCast(expiryObj);
        populateOperFromString(expiry->toString(), oper1);
        populateOperFromString("", oper2);
    } else if (expiryType == MaturityTimePeriod::TYPE){
        MaturityTimePeriodConstSP expiry = 
            MaturityTimePeriodConstSP::dynamicCast(expiryObj);
        populateOperFromString(expiry->getMaturity(), oper1);
        xlTimeFromTime(expiry->getTime(), oper2);
    } else if (expiryType == BenchmarkDate::TYPE){
        BenchmarkDateConstSP expiry = 
            BenchmarkDateConstSP::dynamicCast(expiryObj);
        DateTime date = expiry->toDate();
        xlDateTimeFromDateTime(date, oper1, oper2);
    } else if (expiryType == EnergyContractLabel::TYPE){
        EnergyContractLabelConstSP expiry = 
            EnergyContractLabelConstSP::dynamicCast(expiryObj);
        populateOperFromString(expiry->toString(), oper1);
        populateOperFromString("", oper2);
	} else {
        throw ModelException("XLConvertUtils::opersFromExpiry",
                             "Unkwown expiry type: " + expiryType->getName());
    }
}
    


DRLIB_END_NAMESPACE

