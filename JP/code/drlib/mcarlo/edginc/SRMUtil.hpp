

//----------------------------------------------------------------------------
//
//   Group       : QR cross asset
//
//   Filename    : SRMUtil.hpp
//
//   Description : 
//
//   Author      : Henrik Rasmussen
//
//   Date        : 16 Dec 2005
//
//
//----------------------------------------------------------------------------

#ifndef SRM_UTIL_HPP
#define SRM_UTIL_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

// round to nearest int - slow
int SRMRound(double x);

// from Mark's SRMRatesUtil::yearFrac - rewrite/fixme?
double SRMYearFrac(const DateTime &dateFrom,
                   const DateTime &dateTo);
    
double SRMExp(double exponent);

//////////////////////////////////////////////////////////////////////////////
//  From eqvol.c:SimpleFwdCurve
//
//   Given NbTP points, it outputs (NbTP - 1) "Forward rates * TimePeriod",
//     where FwdRate[i] spans (TPDate[i],TPDate[i+1]).
//
//  Note : this is a general utility, so we should remove versions from other
//         fiels as soon as convenient
//////////////////////////////////////////////////////////////////////////////

MCARLO_DLL vector<double> simpleFwdCurve(
    IYieldCurveConstSP   yc,       // (I)
    const DateTimeArray& TPDate);  // (I) date of time point 


struct MCARLO_DLL SRMUtil // namespace would be better, but beware of VC6
{

/*****  volExtend **************************************************/
/*
*      Extends a curve using the flat interp method and modifies
*      the date indexing of the extended curve. It copes with cases where
*      a period in the extended curve covers multiple periods in the source
*      curve.
*
*      On input:
*      ---------
*      SrcCurve[t] is applicable between InputDates[t-1] to
*      InputDates[t]. (when t=0, it will be between Today and InputDates[0])
*
*      Note: InputDates must be in strictly ascending order
*
*      On output:
*      ----------
*      ExtCurve[t] is the spot vol applicable between ExtendDates[t] and 
*      ExtendDates[t+1].
*
*      Note: 1) We assume TimePt[0] to be Today
*            2) Size of ExtSpotVol must be the same as that of TimePts
*            3) We assume ExtSpotVol[NbTimePts-1] is never used because
*               it is meaningless.
*            4) TimePts must be in strictly ascending order
*
*  From util_s::ExtendSpotVol 
*
*/


    static vector<double> extendVol(
                const DateTimeArray&  InputDates, // vol expiry dates
                const vector<double>& SrcCurve,
                const DateTimeArray&  ExtendDates);

    static bool SRMExplosion(double logIr);

    //// the integral of exp(-b*t) between t and t+dt
    static double GFAC(double t, double dt, double b);
	static double GFACSimple(double t, double b);

    // Return today and all future dates inside allDates;
    static DateTimeArray calcDiffusionDates(const DateTime& today, const DateTimeArray& allDates);
    // Return the size of DateTimeArray returned by calcDiffusionDates()
    static int           getNumSimDates(const DateTime& today, const DateTimeArray& allDates);

    static vector<double> computeSqrtYearFrac(const DateTimeArray& dates);
    
	// year fractions, offset from today
	static vector<double> computeYearFrac(const DateTime& today, const DateTimeArray& dates);

    /** From the list of dates requested add additional points as specified
        by timePointsPerYear */
    static DateTimeArraySP fillInTimeLine(
        const DateTime&            start,      // must be <= inputDates[0]
        const DateTimeArray&       inputDates, // dates in order
        double                     timePointsPerYear);


};

//Deprecated: QLIB_VERIFY now automatically deduces the method, using __FUNCTION__
#define SRM_VERIFY(A,B,C) QLIB_VERIFY(A,B)

DRLIB_END_NAMESPACE

#endif

