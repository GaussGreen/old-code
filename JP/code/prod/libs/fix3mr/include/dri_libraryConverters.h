// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
//
// $Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/include/dri_libraryConverters.h,v 1.1 2004/04/19 15:50:39 markss Exp $
//

#ifndef _IR_LIBRARY_CONVERTERS_
#define _IR_LIBRARY_CONVERTERS_


namespace IR {

//////////////////////////////
//
//  E N U M E R A T I O N S 
//
//////////////////////////////

//
//  frequency('M', 'Q', 'S', 'A') and converter
//

LG_ENUM4(Frequency,
    IR_ANNUAL,   "A",   "Annual",
    IR_SEMI,     "S",   "Semi",
    IR_QUARTER,  "Q",   "Quarterly",
    IR_MONTHLY,  "M",   "Monthly")

LG_ENUM_BEGIN(Frequency)
LG_DESCRIPTION("Specify the curve frequency")
LG_END()

struct FrequencyEnumAsChar 
{
    typedef Frequency ExportAsType;
    static void toPublic(char src, ExportAsType &value ) 
    {
        switch (src) 
        {
            case 'A': value = IR_ANNUAL; break;
            case 'S': value = IR_SEMI; break;
            case 'Q': value = IR_QUARTER; break;
            case 'M': value = IR_MONTHLY; break;
            default: 
                throw Error("Unhandled enum character in FrequencyEnumAsChar converter");
        }
    }
    static void fromPublic(char &dst, ExportAsType value ) 
    {
        switch (value)
        {
            case IR_ANNUAL:   dst = 'A'; break;
            case IR_SEMI:     dst = 'S'; break;
            case IR_QUARTER:  dst = 'Q'; break;
            case IR_MONTHLY:  dst = 'M'; break;
            default: 
                throw Error("Unhandled enum value in FrequencyEnumAsChar converter");
        }
    }
};


//
//  True/False and converter
//

LG_ENUM2(TrueFalse,
    IR_TRUE,  "T",   "True",
    IR_FALSE, "F",   "False")

LG_ENUM_BEGIN(TrueFalse)
LG_DESCRIPTION("Specify if true or false condition")
LG_END()

struct TrueFalseEnumAsInt
{
    typedef TrueFalse ExportAsType;
    static void toPublic(int src, ExportAsType &value)
    {
        // true may be any non zero value in library
        if (src)  
            value = IR_TRUE;
        else
            value = IR_FALSE;
    }
    static void fromPublic(int &src, ExportAsType value)
    {
        if (value == IR_TRUE)
            src = 1;
        else if (value == IR_FALSE)
            src = 0;
        else
            throw Error("Unhandled enum value in TrueFalseAsInt converter");
    }
};


//
//  Yes/No and converter
//

LG_ENUM2(YesNo,
    IR_YES,  "Y",   "Yes",
    IR_NO,   "N",   "No")

LG_ENUM_BEGIN(YesNo)
LG_DESCRIPTION("Specify if Yes or No condition")
LG_END()

struct YesNoEnumAsChar
{
    typedef YesNo ExportAsType;
    static void toPublic(char src, ExportAsType &value)
    {
        if (src == 'Y')
            value = IR_YES;
        else if (src == 'N')
            value = IR_NO;
        else
            throw Error("Unhandled enum character in YesNoEnumAsChar converter");
    }
    static void fromPublic(char &src, ExportAsType value)
    {
        if (value == IR_YES)
            src = 'Y';
        else if (value == IR_NO)
            src = 'N';
        else
            throw Error("Unhandled enum value in YesNoEnumAsChar converter");
    }
};

//
//  Long/Short and converter
//

LG_ENUM2(LongShort,
    IR_LONG,  "L",   "Long",
    IR_SHORT, "S",   "Short")

LG_ENUM_BEGIN(LongShort)
LG_DESCRIPTION("Specify if Morgan is Long or Short the Trade")
LG_END()

struct LongShortEnumAsChar
{
    typedef LongShort ExportAsType;
    static void toPublic(char src, ExportAsType &value)
    {
        if (src == 'L')
            value = IR_LONG;
        else if (src == 'S')
            value = IR_SHORT;
        else
            throw Error("Unhandled enum character in LongShortEnumAsChar converter");
    }
    static void fromPublic(char &src, ExportAsType value)
    {
        if (value == IR_LONG)
            src = 'L';
        else if (value == IR_SHORT)
            src = 'S';
        else
            throw Error("Unhandled enum value in LongShortEnumAsChar converter");
    }
};

/*
//
//  Index label converter - most likely a temporary solution until char[N]
//  can be exported as string rather than as an array of 1 character strings
//

LG_ENUM19(IndexLabel,
    IR_1M,  "1m",   "One Month",
    IR_3M,  "3m",   "Three Months",
    IR_6M,  "6m",   "Six Months",
    IR_12M, "12m",  "Twelve Months",
    IR_1Y,  "1Y",   "One Year",
    IR_2Y,  "2Y",   "Two Years",
    IR_3Y,  "3Y",   "Three Years",
    IR_4Y,  "4Y",   "Four Years",
    IR_5Y,  "5Y",   "Five Years",
    IR_6Y,  "6Y",   "Six Years",
    IR_7Y,  "7Y",   "Seven Years",
    IR_8Y,  "8Y",   "Eight Years",
    IR_9Y,  "9Y",   "Nine Years",
    IR_10Y, "10Y",  "Ten Years",
    IR_12Y, "12Y",  "Twelve Years",
    IR_15Y, "15Y",  "Fifteen Years",
    IR_20Y, "20Y",  "Twenty Years",
    IR_25Y, "25Y",  "Twenty Five Years",
    IR_30Y, "30Y",  "Thirty Years")

LG_ENUM_BEGIN(IndexLabel)
LG_DESCRIPTION("String label representing an index date")
LG_END()

struct IndexLabelEnumAsString
{
    typedef IndexLabel ExportAsType;
    static void toPublic(char *src, ExportAsType &value)
    {
        if (src[0] == '1' && src[1] == 'M')
            value = IR_1M;
        else if (src[0] == '3' && src[1] == 'M')
            value = IR_3M;
        else if (src[0] == '6' && src[1] == 'M')
            value = IR_6M;
        else if (src[0] == '1' && src[1] == '2' && src[1] == 'M')
            value = IR_12M;
        else if (src[0] == '1' && src[1] == 'Y')
            value = IR_1Y;
        else if (src[0] == '2' && src[1] == 'Y')
            value = IR_2Y;
        else if (src[0] == '3' && src[1] == 'Y')
            value = IR_3Y;
        else if (src[0] == '4' && src[1] == 'Y')
            value = IR_4Y;
        else if (src[0] == '5' && src[1] == 'Y')
            value = IR_5Y;
        else if (src[0] == '6' && src[1] == 'Y')
            value = IR_6Y;
        else if (src[0] == '7' && src[1] == 'Y')
            value = IR_7Y;
        else if (src[0] == '8' && src[1] == 'Y')
            value = IR_8Y;
        else if (src[0] == '9' && src[1] == 'Y')
            value = IR_9Y;
        else if (src[0] == '1' && src[1] == '0' && src[2] == 'Y')
            value = IR_10Y;
        else if (src[0] == '1' && src[1] == '2' && src[2] == 'Y')
            value = IR_12Y;
        else if (src[0] == '1' && src[1] == '5' && src[2] == 'Y')
            value = IR_15Y;
        else if (src[0] == '2' && src[1] == '0' && src[2] == 'Y')
            value = IR_20Y;
        else if (src[0] == '2' && src[1] == '5' && src[2] == 'Y')
            value = IR_25Y;
        else if (src[0] == '3' && src[1] == '0' && src[2] == 'Y')
            value = IR_30Y;
        else
            throw Error("Invalid value for Index string");
    }
    static void fromPublic(char *src, ExportAsType value)
    {
        if (sizeof(src) < 4)
            throw Error("C index string must be at least 4 characters long");

        switch (value)
        {
            case IR_1M:
                strcpy(src, "1m"); break;
            case IR_3M:
                strcpy(src, "3m"); break;
            case IR_6M:
                strcpy(src, "6m"); break;
            case IR_12M:
                strcpy(src, "12m"); break;
            case IR_1Y:
                strcpy(src, "1y"); break;
            case IR_2Y:
                strcpy(src, "2y"); break;
            case IR_3Y:
                strcpy(src, "3y"); break;
            case IR_4Y:
                strcpy(src, "4y"); break;
            case IR_5Y:
                strcpy(src, "5y"); break;
            case IR_6Y:
                strcpy(src, "6y"); break;
            case IR_7Y:
                strcpy(src, "7y"); break;
            case IR_8Y:
                strcpy(src, "8y"); break;
            case IR_9Y:
                strcpy(src, "9y"); break;
            case IR_10Y:
                strcpy(src, "10y"); break;
            case IR_12Y:
                strcpy(src, "12y"); break;
            case IR_15Y:
                strcpy(src, "15y"); break;
            case IR_20Y:
                strcpy(src, "20y"); break;
            case IR_25Y:
                strcpy(src, "25y"); break;
            case IR_30Y:
                strcpy(src, "30y"); break;
            default:
                throw Error ("Unhandled enumeration in index string");
        }
    }
};
*/
//
//  Year Basis and converter
//

LG_ENUM3(YearBasis,
    IR_ACT, "ACT",   "Actual number of days in year",
    IR_360, "360",   "360 days in year",
    IR_365, "365",   "365 days in year")
    

LG_ENUM_BEGIN(YearBasis)
LG_DESCRIPTION("Specify year basis used in calculations")
LG_END()

// converter not to be used to export type, but may be used within the DRI layer
// to manually convert between the enums and the string
struct YearBasisEnumAsString
{
    
    static void toPublic(const char *_src, YearBasis &value)
    {
        // reference comes in as reference to array - convert to pointer in order
        // to use normal array operators etc. on value
        //char* _src = &src;  
        if (_src[0] == 'A' && _src[1] == 'C' && _src[2] == 'T')
            value = IR_ACT;
        else if (_src[0] == '3' && _src[1] == '6' && _src[2] == '0')
            value = IR_360;
        else if (_src[0] == '3' && _src[1] == '6' && _src[2] == '5')
            value = IR_365;
        else
            throw Error(std::string("Invalid value for Year Basis") <<
                _src);
    }
    static void fromPublic(char* src, const YearBasis value)
    {

       /*if (sizeof(*src) < 4)
            throw Error("C index string must be at least 4 characters long"); */

        char* _src = src;
        switch (value)
        {
            case IR_ACT:
                strcpy(_src, "ACT"); break;
            case IR_360:
                strcpy(_src, "360"); break;
            case IR_365:
                strcpy(_src, "365"); break;
            default:
                throw Error("Unhandled enum value for Year Basis");
        }
    }
};




//////////////////////////////
//
//  C O N V E R T E R S
//
//////////////////////////////

// convert to and from yyyymmdd
struct DateAsLong 
{
    typedef lego::Date ExportAsType;

    static void toPublic( long src, ExportAsType &value ) 
    {
        value.year  =  src / 10000L;
        value.month = (src % 10000L) / 100L;
        value.day   = (src % 100L  );
    }
    static void fromPublic( long &dst, ExportAsType value ) 
    {
        dst = value.year * 10000L + value.month * 100L + value.day;
    }
};


}
#endif
