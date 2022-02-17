//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : AsianASWStrike.hpp
//
//   Description : Asset swap strike object for Asian asset swaps
//
//   Author      : André Segger
//
//   Date        : 03 June 2003
//
//
//----------------------------------------------------------------------------

#ifndef ASIAN_ASW_STRIKE_HPP
#define ASIAN_ASW_STRIKE_HPP

#include "edginc/AssetSwapStrike.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL AsianASWStrike: public AssetSwapStrike {
public:
    static CClassConstSP const TYPE;
    friend class AsianASWStrikeHelper;

private:

    DateTime             valueDate;
    // YieldCurveWrapper    discount;

    DateTimeArraySP      refixDates;
    DateTimeArraySP      accrualDates;
    DateTimeArraySP      payDates;
    DoubleArraySP        fixings;

    double               swapNotional;
    double               spread;
    string               swapInterval; // $unregistered
    string               swapDCC;
    string               rateBadDayConv; // $unregistered
    HolidayWrapper       hol; // $unregistered


    ScheduleSP           exerSched;
    double               swapBackEndFee;

    /** calculate the strike of the asset swap */
    virtual double getStrike(const DateTime&    baseDate,
                             ConvBondConstSP    cvb,
                             YieldCurveConstSP  discountCurve) const;

    /** returns true if the asset swap is exercisable on the baseDate */
    virtual bool   isExercisable(const DateTime& baseDate);

    virtual void validatePop2Object(); // initialize a few params
};

typedef smartPtr<AsianASWStrike>        AsianASWStrikeSP;
typedef smartConstPtr<AsianASWStrike>   AsianASWStrikeConstSP;

DRLIB_END_NAMESPACE
#endif
