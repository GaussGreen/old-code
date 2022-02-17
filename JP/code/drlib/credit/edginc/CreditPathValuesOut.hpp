//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditValueOut.cpp
//
//   Description : Credit engine
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_PATH_VALUES_OUT_HPP
#define CREDIT_PATH_VALUES_OUT_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/MarginAcct.hpp"
#include "edginc/CreditDebug.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/CreditSpreadCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/** simulated path values */
class CREDIT_DLL CreditPathValuesOut: public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
	
    static IObject* defaultCreditPathValuesOut(){
        return new CreditPathValuesOut();
    }

    /** add another CreditPathValuesOut object's values to this one*/
    void add(const CreditPathValuesOut* rhs);

    /** calculate peak and average exposures from path values */
    void calcPeakAverageStdAndDRE(const CMarketDataSP&        market,
	 						 const MarginAcctSP&        marginAcct,
	 						 double		                percentile,     
							 const string&              creditCcyName,
							 const CreditSpreadCurveSP&	riskyCurve,
                             const IModelSP              model = IModelSP(   ));

	/** calculate profile of DRE */
	void computeDreProfile(
		int							numExposureDates,
		const YieldCurveSP			&discount,
		const CreditSpreadCurveSP	&creditSpread);

    /** access methods */
    DateTimeArraySP&  getPathDates();
    CIntArraySP&  getPathDateTypes();
    DoubleArrayArraySP&  getPathValues();
	CreditPathValuesOut();

	/** remove pathValues from output */
	void flushGrids();

    static inline int numHistoricalDates(const DateTimeArraySP& dates, 
                                     const DateTime& theDate, 
                                     bool todayIsHistory = false)
    {
        DateTimeArray& datesArr = *dates;
        int count = 0;
        while (count < datesArr.size() && 
                ( todayIsHistory ? 
                    datesArr[count].getDate() <= theDate.getDate() : 
                    datesArr[count].getDate() < theDate.getDate() )
               )
        {
            count++;
        }
        return count;
    }
                  
private:
	void arraySum(DoubleArrayArraySP target,
				  const DoubleArrayArraySP source);

    friend class CreditEngine;

    DateTimeArraySP     pathDates;
    CIntArraySP         pathDateTypes; // 0=exposure date, 1=mtm date, 2=both
    DoubleArrayArraySP  pathValues; // values on each path date
	
    DateTimeArraySP     pathMTMDates; // not registered.  contains all mtm dates. $unregistered
	DateTimeArraySP     pathExposureDates;
	DoubleArraySP		peakProfile;
	DoubleArraySP		aveProfile;
	DoubleArraySP		dreProfile;
	DoubleArraySP		stdProfile;

    DateTimeArray       gemBandDates;  /** dates for GEMS21 banding (input) */
    DoubleArraySP       gemPeak;       /** GEMS21 peak values */
    DoubleArraySP       gemAve;        /** GEMS21 avg values */
    DoubleArraySP       gemDRE;        /** GEMS21 dre values */
    DoubleArraySP       gemStd;        /** GEMS21 dre values */

    double              peakExposure;
    double              DRE;
    double              aveExposure;
	double              cvr;

    /** only dates on or before this date are included in the average computation. */
    DateTime            lastExposureDate; 

	// contains pathSpots, pathDates, pathValues and calcTime
	CreditDebugSP		creditDebug;

    // transient
    
	/**  exposure ccy, copied from input*/
    string              exposureCcy;   

    /**  Buckets an exposure profile into a new set of dates. */
	void calcExposureProfileBandForGEMS(const YieldCurveSP& discount);

    /** Computes peak and ave Exposure members, given an exposureWithMargin profile for each exposure date and path */
    void computeExposureProfiles(
        int numExposureDates,
        double percentileCount,
        const DoubleArrayArraySP& exposureWithMargin);

    void resetLastExposureDate(const CreditPathValuesOut* rhs); 

	void EdgCreditBandedStatistics(
   YieldCurveSP        discount,     /* (I) Currency to use for discounting  */
   CreditSpreadCurveSP creditSpread); /* (I) Credit spreads  (NULL if no CVR) */
//   double            *expected,     /* (O) Expected exposure                */
//   double            *peak,         /* (O) Peak loan equivalent             */
//   double            *cvr);          /* (O) Capital Valuation Reserve        */

   int numNonZeroExposureDates() const;

};

typedef smartConstPtr<CreditPathValuesOut> CreditPathValuesOutConstSP;
typedef smartPtr<CreditPathValuesOut> CreditPathValuesOutSP;

typedef array<CreditPathValuesOutSP, CreditPathValuesOut> CreditPathValuesOutArray;

DRLIB_END_NAMESPACE

#endif
