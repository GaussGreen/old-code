// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
//
// $Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/include/fix123components.h,v 1.1 2004/04/19 15:50:39 markss Exp $
//


// collection of components that may be used by the fix3 products


#ifndef _fix123components_
#define _fix123components_

#include <vector>
#include <string>
#include "IRTraits.h"


using namespace std;

namespace IR {


//
//   O p t i o n D a t a
//

LG_CLASS_BEGIN(OptionData)
LG_DESCRIPTION("This structure stores the option information")
LG_END()

class OptionData {
public:

    LG_CLASS_MEMBERS8(
        "Call or Put on Underlying",
        char, CoP,
        "Exercise Frequency for Option (American or European)",
        char, ExerFreq,

        "Number of Exercise Dates",
        int, NbExer,
        "Notification dates",
        vector<Date>, NotifDate,         
        "Effective Notification Dates",
        vector<Date>, NotifEffDate,      
        "Exercise dates",
        vector<Date>, ExerDate,
        "Strike levels (as decimal)",
        vector<double>, Strikes,

        "Linear or staircase strike interpolation",
        char, ProfileInterp)
        
};

//
//   F i x e d L e g D a t a
//

LG_CLASS_BEGIN(FixedLegData)
LG_DESCRIPTION("This structure stores the fixed leg data for a swaption")
LG_END()

class FixedLegData {

public:
    LG_CLASS_MEMBERS9(
        // fixed dates
        "Nunber of fixed pay and accrue dates",
        int, NbFixDates,
        "Fix accrual start dates",
        vector<Date>, FixAccStart,
        "Fix accrual end dates",
        vector<Date>, FixAccEnd,
        "Pay date for fix",
        vector<Date>, FixPayDate,
        "Cpn day count for fix",
        vector<double>, FixDCF,
        "Cpn notional for fix",
        vector<double>, FixNotional,
        "Cpn rate for fix (as decimal)",
        vector<double>, FixRate,

        "fix DayCount for exer midway",
        char, DayCountFixed,
        "Fix Stub conv for exer midway",
        char, FixOptionStubConv)

};

//
//   F l o a t L e g D a t a
//

LG_CLASS_BEGIN(FloatLegData)
LG_DESCRIPTION("This structure stores the floating leg data for a swaption")
LG_END()

class FloatLegData {

public:
    LG_CLASS_MEMBERS17(        
        // floating dates
        "Number of floating pay and accrue dates",
        int, NbFltDates,
        "Reset notif date for flt",
        vector<Date>, FltResetDate,      
        "Reset effective date for flt",
        vector<Date>, FltResetEffDate,   
        "Flt acc start dates",
        vector<Date>, FltAccStart,
        "Flt acc end dates",
        vector<Date>, FltAccEnd,
        "Pay date for flt",
        vector<Date>,  FltPayDate,
        "Cpn day count for flt",
        vector<double>, FltDCF,
        "Cpn notional for flt",
        vector<double>, FltNotional,
        "Spread to float index (as decimal)",
        vector<double>, FltSpread,

        "Floating leg DayCount convention (for part-way exercise)",
        char, DayCountFloat,
        "Floating leg Stub convention (for part-way exercise)",
        char, FltOptionStubConv,
        "Simple or compound floating cpn",
        char, CompFloat,

        // floating leg index
        "Maturity of floating index (in months)",
        int, PayIndexMat,
        "Frequency of floating index",
        char, PayIndexFreq,
        "Day count convention of index",
        char, PayIndexDayCount,
        "Curve for floating index",
        int, PayIndexIoD,
        
        "Estimate previous floating fixing from curve, or set rate (as decimal)",
        std::string, LastRefixSetting)

};





}  // end namespace fix3










#endif

