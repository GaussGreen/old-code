/*
 * $Log: infcurv.h,v $
 * Revision 1.23  2004/05/11 10:45:40  jpriaudel
 * function added by mingzhi
 *
 * Revision 1.22  2003/11/25 10:13:16  rguillemot
 * Inflation Season Manager
 *
 * Revision 1.21  2003/09/23 17:32:38  ebenhamou
 * accessor to MonthyInterpType
 *
 * Revision 1.20  2003/09/22 13:50:30  ebenhamou
 * more strict using
 *
 * Revision 1.19  2003/09/11 10:42:34  ebenhamou
 * GetCPIIndexValue
 *
 * Revision 1.18  2003/09/05 07:20:15  ebenhamou
 * added more constant
 *
 * Revision 1.17  2003/08/27 07:51:06  ebenhamou
 * using namespace
 *
 * Revision 1.16  2003/08/21 08:43:50  ebenhamou
 * inherit from ARM_ZERO_CURVE
 *
 * Revision 1.15  2003/08/20 08:43:12  ebenhamou
 * added virtual for easier reading
 *
 * Revision 1.14  2003/08/18 11:07:17  ebenhamou
 * remove namespace as not supported by unix yet
 *
 * Revision 1.13  2003/08/15 07:07:55  ebenhamou
 * move constant to a namespace
 *
 * Revision 1.12  2003/08/12 08:11:51  ebenhamou
 * added cutoffdate
 *
 * Revision 1.11  2003/08/06 11:31:21  ebenhamou
 * method is public now
 *
 * Revision 1.10  2003/08/05 08:32:20  ebenhamou
 * keep inf curve derived from arm_object ... use ARM_Currency
 *
 * Revision 1.7  2003/07/18 08:59:50  ebenhamou
 * more explicit terminology
 *
 * Revision 1.5  2003/07/16 07:01:17  ebenhamou
 * version with monthly and daily interp
 *
 * Revision 1.4  2003/06/30 17:14:03  ebenhamou
 * dos2unix
 *
 * Revision 1.3  2003/06/30 16:11:13  ebenhamou
 * discounting version
 *
 */

/*!--------------------------------------------------------------------------*
    Declaration of the ARM_InfCurv class
	an inflation curve is defined by historical data and a set 
		- of dates and HCPI
		- a reference date
		- an index name
	Ability to interpolate according to various type the curve
*----------------------------------------------------------------------------*/ 

#ifndef _INGPINFLATION_INFCRV_H
#define _INGPINFLATION_INFCRV_H

#include <glob/firsttoinc.h>

#include <string>
#include <map>
#include <vector>

#include <glob/dates.h>		// because of the use of date
#include <crv/zerocurv.h>
#include "gpbase/port.h"

CC_USING_NS(std,less)
CC_USING_NS(std,map)
CC_USING_NS(std,pair)

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ResetManager;
class ARM_InfIdx;
class ARM_SeasonalityManager;



// The argument less is not necessary, and VC is OK without it.
// Unfortunately not Sun 4.2 ...
typedef map< double, double, less< double > > doubleDoubleMap;


/*
 * class for an inflation curve
 * this is really a forward inflation curve
 * that can computes either Zero Coupon Inflation Rate
 * or CPI forward rate
 * 
 * It is a curve and therefore inherit from a ZEROCURVE
 */
class ARM_InfCurv : public ARM_Object
{
private:
	/// Starting date of the inflation crv
	ARM_Date itsCPIIndexDate;
	double itsCPIIndexValue;

	/// name of the index
	string itsIndexName;
	ARM_Date itsAsOfDate;
	ARM_Date itsCutOffDate;

	/// to allow to generate shift curve, we store 
	/// the mkt data
	vector<string> itsMktTerms;			/// mkt terms for the curve construction
	vector<double> itsMktValues;		/// corresponding mkt data

	/// the curve is stored as a vector of expiry terms
	/// and their corresponding CPI and ZC values
	vector<double> itsExpiryTermsVec;	/// vector of double of expiries
	vector<double> itsCPIValuesVec;		/// vector of CPI values
	vector<double> itsZCValuesVec;		/// vector of ZCValues

	/// these are stored in table
	long itsMonthlyInterpType;	/// type of interpolation for monthly data
	long itsDailyInterpType;	/// type of interpolation for daily data
	long itsDCFMonthly;			/// Monthly daycount
	long itsDCFDaily;			/// Daily daycount
	long itsExtrapolType;		/// weight for the extrapolation

	/// this is stored to avoid many recomputations
	string itsDCFLag;
	string itsResetLag;
	string itsCalendar;
	ARM_Currency* itsCurrency;

	/// gestion des resets
	ARM_ResetManager* itsResetManager; // pointor to the reset manager
	/// seasonality managment
	ARM_SeasonalityManager* itsSeasonalityManager; // pointor to the seasonality manager

	/// construction of the monthly terms
	/// keep growing as more and more is computed
	doubleDoubleMap itsMonthlyCPIValue;
	
	// various cleanup and copy functions
	void CopyNoCleanUp(const ARM_InfCurv& srcobj);
	void CleanUp();

	/// dates function
	double ConvertString2Double( const string& s, const char* calendar ) const;
	void BuildCurve(const vector<string>& MktTerms, const vector<double>& MktValues );

	/// function for the CPI Interpolation
	/// compute CPI Monthly
	double CPIOneDateInterpolateAndStore( const ARM_Date& Date,
		long MonthlyInterpType = -1 );

	/// function to get bounding month dates
	pair<ARM_Date,ARM_Date> GetMonthBoundingDates( const ARM_Date& CPIDate ) const;

	/// to make the difference between a CPI value
	/// we assume that any CPI Value is always above MINCPIVal
	static const double minCPIVal;

	/// default constructor to be compliant with ARM
	ARM_InfCurv();
	void Init();

	double itsBucketStartPeriod, itsBucketEndPeriod;

public:
	/// when we use -1... it will use the default
	/// data in infdata
	/// constructor based on expiry as string
	ARM_InfCurv(
		const ARM_Date& asOfDate, 
		const string&			indexName,
		double					CPIIndexValue,
		const ARM_Date&			CPIIndexDate,
		const vector<string>&	ExpiryTerms,
		const vector<double>&	CPIValues,
		long					MonthlyInterpType	= -1,
		long					DailyInterpType		= -1,
		long					DCFMonthly			= -1,
		long					DCFDaily			= -1,
		long					ExtrapolType		= -1,
		ARM_ResetManager*		ResetManager		= NULL,
		ARM_SeasonalityManager*		SeasonalityManager	= NULL);
	
	/// copy constructor
	ARM_InfCurv( const ARM_InfCurv& srcInfCurv  );

	/// note that the ARM_Object destructor is virtual
	/// hence virtual
	virtual ~ARM_InfCurv(); 

	ARM_InfCurv& operator = ( const ARM_InfCurv& srcInfCurv);
	
	/// function for ARM_Object compatibility
	virtual ARM_Object* Clone();

	/// CPIInterpolation
	/// API function : version with date
	double CPIInterpolate( const ARM_Date& date, 
		const string& DCFLag	= "-1",
		long DailyInterpType	= -1,
		const string& ResetLag	= "-1", 
		double weight			= 0 );		

	/// CPIInterpolation
	/// API function : version with time from AsOf
	inline double CPIInterpolate( double CPItime, 
		const string& DCFLag	= "-1",
		long DailyInterpType	= -1,
		const string& ResetLag	= "-1", 
		double weight			= 0 ) { return CPIInterpolate( ARM_Date( CPItime + itsAsOfDate.GetJulian() ), DCFLag, DailyInterpType, ResetLag, weight ); }

	/// CPIInterpolation
	/// API function : version with time from AsOf
	inline double CPIInterpolate( double CPItime, 
		double DCFtime,
		long DailyInterpType= -1,
		double weight		= 0.0) { double AsOfDate = itsAsOfDate.GetJulian(); return CPIInterpolate( ARM_Date(CPItime+AsOfDate), ARM_Date( DCFtime+AsOfDate), DailyInterpType, weight ); }

	/// version for ZCRate
	double ZCRateInterpolate( const ARM_Date& date, 
		const string& DCFLag	= "-1",
		long DailyInterpType	= -1,
		const string& ZCRateLag	= "-1", 
		double weight			= 0 );		

	/// fast function with no DCF ...
	/// general function
	double CPIInterpolate( 	const ARM_Date& CPIDate, 
		const ARM_Date& DCFDate,
		long DailyInterpType= -1,
		double weight		= 0.0);

	// function to convert CPI to ZC and vice versa
	double ConvertZCRate2CPI( double ZCRate, double expiry ) const;
	double ConvertCPI2ZCRate( double CPIVal, double expiry ) const;

    /// view function in the ARM Menu
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	/// overwritte the function of ARM_ZEROCURVE to allow pricing
	// these two methods are declared virtual in the base class
	virtual double DiscountFunction(double yearTerm ); /// returns discount price
	virtual double D1DiscountFunction(double yearTerm);/// returns derivatives of discount price

	/// we use rather vector of string here
	/// than the infamous char**

	ARM_InfCurv* GenerateShiftCurve(	
		const vector<string>& mktTerms, 
		const vector<double>& epsilon );

	double GetCPIIndexValue() const;
	long GetMonthlyInterpType() const;
	long GetDailyInterpType() const;
	long GetDCFMonthly() const;
	long GetDCFDaily() const;
	long GetExtrapolType() const;
	
	const vector<string>& GetMktTerms()		const	{	return itsMktTerms;		}
	const vector<double>& GetMktValues()	const	{	return itsMktValues;	}

	void SetCurve(const vector<string>& MktTerms, const vector<double>& MktValues ) { 
		itsCPIIndexValue		= MktValues[0];
		itsMktValues			= MktValues;
		itsMonthlyCPIValue.clear();
		BuildCurve( itsMktTerms, itsMktValues );
	}


	// accessor to the reset manager
	ARM_ResetManager* GetResetManager() const {return itsResetManager;};
	ARM_SeasonalityManager* GetSeasonalityManager() const {return itsSeasonalityManager;};

	void SetResetManager(ARM_ResetManager* ResetManager);
	void SetSeasonalityManager(ARM_SeasonalityManager* SeasonalityManager);

	inline ARM_Date GetAsOf() const{ return itsAsOfDate;}; 
	inline ARM_Date GetCPIIndexDate() const{ return itsCPIIndexDate;}; 

	ARM_InfIdx* GetInfIdx() const;
	inline string GetInfIdxName() const { return itsIndexName;}
	inline string GetCalendar() const	{ return itsCalendar;  }


    inline double GetBucketStartPeriod(void) const { return(itsBucketStartPeriod); }
    void SetBucketEndPeriod(double endPeriod) { itsBucketEndPeriod=endPeriod; }
    void SetBucketStartPeriod(double startPeriod) { itsBucketStartPeriod=startPeriod; }
    inline double GetBucketEndPeriod(void) const { return(itsBucketEndPeriod); }

	inline  string			GetIndexName()		const{ return itsIndexName;		}
	inline 	vector<double>	GetExpiryTermsVec()	const{ return itsExpiryTermsVec;}
	inline 	vector<double>	GetCPIValuesVec()	const{ return itsCPIValuesVec;	}
	inline 	vector<double>	GetZCValuesVec()	const{ return itsZCValuesVec;	}

};

CC_END_NAMESPACE()

#endif
