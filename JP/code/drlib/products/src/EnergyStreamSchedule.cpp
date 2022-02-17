//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyStreamSchedule.cpp
//
//   Description : Energy Swap stream schedule
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
//#include "edginc/ErrorHandler.hpp"

#include "edginc/EnergyStreamSchedule.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE

EnergyStreamSchedule::EnergyStreamSchedule() : CObject(TYPE)
{    
}

EnergyStreamSchedule::EnergyStreamSchedule(CClassConstSP clazz) : CObject(clazz)
{
}

EnergyStreamSchedule::~EnergyStreamSchedule()
{    
}

EnergyUnderlyerSP EnergyStreamSchedule::getEnergyUnderlyer() const
{    
    return energyUnderlyer;
}
 

// Coupon start/end dates are derived from the good day schedule (ensuring coupons end on good days).
// Notional start/end dates are derived from a combination of good day
// and bad day scheduling as deliveries can occur on bad days.

EnergyStreamSchedule::EnergyStreamSchedule (
    const EnergyUnderlyerSP&   energyUnderlyer,
	const DateTime&            today,
    const DateTime&            startDate,
    const DateTime&            endDate,
	int                        settleDays,
	const EnergyDWMYPeriod&    paymentFrequency,  
	int                        avgDaysB4End, // 0, 1, ...days from END
	const EnergyDWMYPeriod&    averagingPeriod,
	const EnergyDWMYPeriod&    averagingFrequency,
    int                        rollDay,
	const string&              notionalSetting)
    : CObject(TYPE),
    energyUnderlyer(energyUnderlyer),today(today), startDate(startDate),
	endDate(endDate), settleDays(settleDays),paymentFrequency(paymentFrequency),
    avgDaysB4End(avgDaysB4End),averagingPeriod(averagingPeriod),
    averagingFrequency(averagingFrequency),rollDay(rollDay),
	notionalSetting(notionalSetting),
	rollDate(startDate), alignDate(startDate),hasStub(0)
{
    //buildStreamSchedule();
}

// a version for children
EnergyStreamSchedule::EnergyStreamSchedule (
	CClassConstSP clazz,
    const EnergyUnderlyerSP&   energyUnderlyer,
	const DateTime&            today,
    const DateTime&            startDate,
    const DateTime&            endDate,
	int                        settleDays,
	const EnergyDWMYPeriod&    paymentFrequency,  
	int                        avgDaysB4End, // 0, 1, ...days from END
	const EnergyDWMYPeriod&    averagingPeriod,
	const EnergyDWMYPeriod&    averagingFrequency,
    int                        rollDay,
	const string&              notionalSetting)
    : CObject(clazz),
    energyUnderlyer(energyUnderlyer),today(today), startDate(startDate),
	endDate(endDate), settleDays(settleDays),paymentFrequency(paymentFrequency),
    avgDaysB4End(avgDaysB4End),averagingPeriod(averagingPeriod),
    averagingFrequency(averagingFrequency),rollDay(rollDay),
	notionalSetting(notionalSetting),
	rollDate(startDate), alignDate(startDate),hasStub(0)
{
}

// Would not add this if vector<A> thing were working in FIELD for child class, since we always prefer
// clone() from default object implimentation
IObject* EnergyStreamSchedule::clone() const
{
    static const string routine = "EnergyStreamSchedule::clone";
    try {
        IObject* Copy = CObject::clone(); // call parent
        EnergyStreamSchedule* copy = dynamic_cast<EnergyStreamSchedule*>(Copy);
        if(!copy) {
                throw ModelException("Clone method failed");
        }
       
        copy->energyUnderlyer = energyUnderlyer;
        copy->today = today;
        copy->startDate = startDate;
	    copy->endDate = endDate;
        copy->settleDays = settleDays;
        copy->paymentFrequency = paymentFrequency;
        copy->avgDaysB4End = avgDaysB4End;
  	    copy->averagingPeriod = averagingPeriod;
        copy->averagingFrequency = averagingFrequency;
        copy->rollDay = rollDay; 
        copy->notionalSetting = notionalSetting;
	    copy->rollDate = rollDate;
	    copy->alignDate = alignDate;
	    copy->hasStub = hasStub;
	    copy->currentPeriodStart = currentPeriodStart;

	    copy->coupons = coupons;
        return copy;
    } catch(exception& e) {
            throw ModelException(e, routine);
    }
}


void EnergyStreamSchedule::buildStreamSchedule()
{
    // this is void DRIRStream::BuildCommoditySchedule(  in drirstreamcpp...

    if (paymentFrequency.getInterval()=="M")
    {
        buildMonthlySchedule();
		// roll and align dates have been determined. 
    }
    else if (paymentFrequency.getInterval()=="I" || paymentFrequency.getInterval()=="O")
	{
		// if coming to here, FutureSwapSchedule has been used and hence startDate 
		// is the fixing date of previous contranct that is used for alignment, 
		// but not the start date of the stream.
		alignDate = startDate;
		startDate = energyUnderlyer->addDaysToDateAndAdjust(startDate,1);
	}
	else
    {
		// note, interval has been converted to either "M" or "D". So, "D" hear
        alignDate = findRollDate();  // rollDate has been figured out inside the function. 
    }
    
    // ensure last date of swap is a good day
    DateTime finalBoundary = energyUnderlyer->busDayAdjust(endDate);

    // Make 2 lists of coupon boundary dates - a "good" boundary dates list and a "bad" boundary list
    // Bad boundary dates are coupon end dates that may fall on bad days.  These are required
    // for working out the notional start/end dates which may fall on bad days
    DateTimeArray couponDates;
	DateTimeArray notionalDates; // "bad" list

    // coupon starts on startDate if it's a good day or next good day
    couponDates.push_back(energyUnderlyer->busDayAdjustFollowing(startDate));

    // For notional dates
    notionalDates.push_back(startDate);

    int i = 1;
	int j;
	DateTime currentDate;
	DateTime currentBadDate;
	int payPeriod = paymentFrequency.getCount();

	if (paymentFrequency.getInterval()=="M")
    {
		// interval is a month
        DateTime::MonthDayYear alignMDY = alignDate.toMDY();
        DateTime::MonthDayYear couponMDY = alignMDY;
		
		do
		{
            // adding a coupon period, first by adding coupon/notional end date
			
            currentBadDate = energyUnderlyer->addIBORMonthsToDateNoAdj(alignDate,(i++)*payPeriod);

			currentDate = energyUnderlyer->busDayAdjust(currentBadDate);
		
            couponDates.push_back(currentDate);
            notionalDates.push_back(currentBadDate);
		}
        while (finalBoundary.daysDiff(currentDate) > 0);
    }
	else if (paymentFrequency.getInterval()=="I" || paymentFrequency.getInterval()=="O")
	{
		int type = (paymentFrequency.getInterval()=="I") ? 0 : 1;
		do
		{
            // should always get a good day compared to above
            currentDate = energyUnderlyer->addContractPeriod(alignDate,i++, type);

			// should never get a bad day, so commented out
			//currentDate = energyUnderlyer->busDayAdjust(currentBadDate);
		
            couponDates.push_back(currentDate);
            notionalDates.push_back(currentDate);
		}
        while (finalBoundary.daysDiff(currentDate) > 0);
	}
    else
    {	
        do
		{
            // adding a coupon period (e.g. 5D)
   
			currentDate = energyUnderlyer->addDaysToDateAndAdjust(alignDate, i*payPeriod);
            currentBadDate = alignDate.rollDate((i++)*payPeriod);
  
            couponDates.push_back(currentDate);
            notionalDates.push_back(currentBadDate);
		}
        while (finalBoundary.daysDiff(currentDate) > 0);
	}
 
    // if overshoot ...
    if (couponDates.back().daysDiff(finalBoundary) >0 )
        couponDates.back() = finalBoundary;
    if (notionalDates.back().daysDiff(finalBoundary) > 0 )
        notionalDates.back() = endDate;
        
    
    // Calculate corresponding payment dates
    DateTimeArray payDates(couponDates.size());

    for (i=0; i<couponDates.size(); ++i)
    {
		// note, payDates[0], ie startDate is also adjusted.
        payDates[i] = energyUnderlyer->busDayAdjustFollowing(couponDates[i]);
    }

    if (settleDays!=0)
    {
        // adding good days to coupon
        //DRPeriodGoodDay paymentLagAdjuster(fPaymentOffset, fIndex->GetBadDays());
        for (i=0; i<couponDates.size(); ++i)
            payDates[i] = energyUnderlyer->addBusDaysToDate(payDates[i],settleDays);
    }

    // nextCoupon is the next coupon to pay
    // Stub may have already paid, so we keep track of whether we're actually creating a stub payment
    
    j = 0;
    DateTime notionalStart;
	DateTime notionalEnd;
	DateTime couponStart;
    DateTime couponBadStart;
    DateTime couponEnd;
    DateTime couponBadEnd;
	DateTime paymentDate;

    for (i=1; i<payDates.size(); ++i, ++j)
    {
        // Build schedule, ignoring "obviously" realised payments
        // note payDates[0] is start date, not relavant. 
        if (payDates[i].daysDiff(today)>=0)
        {
            if (j > 0)
            {
                // coupon start = first good/bad day after end of previous coupon
				// Make sure its end of previous coupon??????????????????????????
                couponStart = energyUnderlyer->addBusDaysToDate(couponDates[i-1],1);
                couponBadStart = notionalDates[i-1].rollDate(1);
            }
            else 
            {
                // if first coupon 
                couponStart = couponDates[i-1];
                couponBadStart = notionalDates[i-1];
            }
            couponEnd = couponDates[i];
            couponBadEnd = notionalDates[i];
            paymentDate = payDates[i];

            // if notional dates are for delivery month (eg. gas swaps) find the start and end dates
            // of the delivery month.  This is only valid for monthly swaps
            if (notionalSetting == "MATCHDELIVERY")
            {
                energyUnderlyer->calculateContractMonthBoundaryDates(
                    couponEnd, EnergyUnderlyer::INDEX, &notionalStart, &notionalEnd);
            }
            else
            {
                // notinonal start is the first good or bad coupon start date that falls after
                // the notional end of the previous coupon
                notionalStart = (couponStart < couponBadStart)? couponStart: couponBadStart;
                notionalStart = (notionalStart > notionalEnd.rollDate(1))? notionalStart: notionalEnd.rollDate(1);
                notionalEnd = (couponEnd > couponBadEnd)? couponEnd: couponBadEnd;
            }


            EnergySwapCouponSP newCoupon(
				createCouponPayment(
                    today,
                    energyUnderlyer,
                    couponStart, 
                    couponEnd, 
                    notionalStart,
                    notionalEnd,
                    paymentDate,
					averagingPeriod,
	                averagingFrequency,
	                avgDaysB4End )   );
            
            newCoupon->setNumTotalCoupons(couponDates.size() - 1);
            coupons.push_back(newCoupon);
        }
    }

}

void EnergyStreamSchedule::buildMonthlySchedule()
{
    static const string method = "EnergyStreamSchedule::buildMonthlySchedule";
    
    int payMonths = paymentFrequency.getCount(); // num of Months in payment period
    
    if (payMonths<=0)
        throw ModelException(method,"Bad payment interval months");


	DateTime::MonthDayYear MDY1 = startDate.toMDY();
    DateTime::MonthDayYear MDY2 = endDate.toMDY();
	int d1 = MDY1.day;
	int d2 = MDY2.day;

	// num of months between start and end ofthe deal
    int numMonths = (MDY2.month-MDY1.month) + 12*(MDY2.year-MDY1.year);
	
	DateTime streamEnd = energyUnderlyer->addIBORMonthsToDate(startDate, numMonths);

	 
    int fwdAdjMonths = numMonths%payMonths;

	// We always schedule from end
    bool couponsMatch = fwdAdjMonths==0 && (d1==d2 || streamEnd==energyUnderlyer->busDayAdjust(endDate));
    
	// if coupons match, return, otherwise, continue...
    
	if (!couponsMatch)
    {
		// either total num of months is not a multiple of pay frequency or more/less days than the multiple.
        hasStub = 1;
       
        // Front stub is forced, need to find roll date( from which coupons will match rest of deal)
        // A roll period will take us forward to the right day of a month forward.
        
        if (fwdAdjMonths==0 && d2<d1)
		{
			// less days. Then, roll this num of months forward
            fwdAdjMonths = payMonths;
		}

		//DRPeriodRollDay::Apply() logic
		
		MDY1.month = MDY1.month + fwdAdjMonths;
        rollDate = MDY1.toDateTime();
		MDY1 = rollDate.toMDY();
		MDY2 = rollDate.returnEndOfMonth(false).toMDY();

		// last day for the month
		////////////////////////////////////////////////////////////////////////////////////
		// this is the logic inside DRPeriodRollDay::Apply(), which is wrong to use.
		// We should make sure that roll date is end of month, if endDate is end of month.
		//MDY1.day = (d2<MDY2.day)? d2: MDY2.day; -- commented out.
		DateTime endDay = endDate.returnEndOfMonth(false); // new code starts
	    DateTime endGoodDay = energyUnderlyer->busDayAdjust(endDate.rollDate(endDay.daysDiff(endDate)) );
	    bool isEndOfMonth = (endDate.daysDiff(endGoodDay)>=0) ? true : false;

		if (isEndOfMonth)
			rollDate = rollDate.returnEndOfMonth(false);
		else
		{
			MDY1.day = (d2<MDY2.day)? d2: MDY2.day;
            rollDate = MDY1.toDateTime();
		}
        ////////////////////////////////////////////////////////////////////////////////////

        // implement override roll date
		
        if (!rollDay)
            //MDY2.day = d2; this is wrong also. out.
            alignDate = rollDate;
        else
		{
            MDY2.day = rollDay;
            alignDate = MDY2.toDateTime();
		}             
    }
}

DateTime EnergyStreamSchedule::findRollDate()
{
	static const string method = "EnergyStreamSchedule::findRollDate";

	// Calculate coupons from the start date until we equal or exceed the end date
	
	int payDays = paymentFrequency.getCount();
	int numRegularPeriods = 0;
	DateTime finalCoupon = startDate;
    DateTime oldFinalCoupon;

	// determining how many days to endDate
	while (endDate.daysDiff(finalCoupon) > 0)
	{
		oldFinalCoupon = finalCoupon;
		// counting how many coupons
		finalCoupon = energyUnderlyer->addDaysToDateAndAdjust(rollDate,++numRegularPeriods*payDays);		
	}
	
	// If finalCoupon!=finalAccrueDate, do we truncate the final coupon, or does the swap
	// have a stub?  To find out, we repeatedly decrease rollDate and recalculate the corresponding
	// finalCoupon until finalCoupon decreases.  If finalCoupon<finalAccrueDate, we truncate the
	// final period.  Otherwise, the swap has a stub.
	
	if ( finalCoupon.daysDiff(endDate) != 0 )
	{
		// finalCoupon > endDate.
		DateTime rollDate1 = rollDate;
		DateTime finalCoupon1 = finalCoupon;
		while (finalCoupon1.daysDiff(finalCoupon) == 0 )
		{
			rollDate1 = rollDate1.rollDate(-1);
			finalCoupon1 = energyUnderlyer->addDaysToDateAndAdjust(rollDate1, numRegularPeriods*payDays);	
		}

		hasStub = (finalCoupon1.daysDiff(endDate) >= 0 ); // endDate<finalCoupon1<finalCoupon, has a front stub.
	}
	
	// 	If the swap has a front stub, find the right rollDate
	if (hasStub)
	{
		// We're replacing a regular period with a special stub period	
		--numRegularPeriods;
		
		if (numRegularPeriods==0)
		{
			rollDate = endDate;
		}
		else
		{
			// for binary search, invariant is:
			//   rollDate==left    =>   finalCoupon<finalAccrueDate
			//   rollDate==right   =>   finalCoupon>=finalAccrueDate
			// so right==left+1    =>   right is the first sync date for which final coupon date equals or exceeds swap end date
			DateTime left = startDate;
			DateTime right = energyUnderlyer->addDaysToDateAndAdjust(startDate,payDays);
			DateTime mid, endDateForMid;
			while (right.rollDate(-1).daysDiff(left) > 0)
			{
				int middle = right.daysDiff(left) >> 1;
				mid = left.rollDate(middle);
				endDateForMid = energyUnderlyer->addDaysToDateAndAdjust(mid, numRegularPeriods*payDays);

				if (endDateForMid.daysDiff(endDate)<0)
					left = mid;
				else
					right = mid;
			}
			rollDate = right;
		}
	}
	
	return rollDate;
}

void EnergyStreamSchedule::getScheduleDetails(

		StringArray& thePaymentDates,
		StringArray& theCouponStartDates,
		StringArray& theCouponEndDates,
		StringArray& theNotionalStartDates,
		StringArray& theNotionalEndDates,
		IntArray& theNumFixings,
		StringArray& theFixingDates,
		StringArray& theFixingLabels,
		DoubleArray& theFixingRates ) const
{
    static const string method = "EnergyStreamSchedule::getScheduleDetails";
	throw ModelException(method, "Should not come to here");
}

vector<EnergySwapCouponSP> EnergyStreamSchedule::getCoupons() const
{
	return coupons;
}

class EnergyStreamScheduleHelper
{
  
public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergyStreamSchedule, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEnergyStreamSchedule);

		FIELD(today, "");
		FIELD_MAKE_TRANSIENT(today);
        FIELD(energyUnderlyer, "");
		FIELD_MAKE_TRANSIENT(energyUnderlyer);
        FIELD(paymentFrequency,  ""); 
        FIELD_MAKE_TRANSIENT(paymentFrequency);
		FIELD(averagingPeriod,       "");
        FIELD_MAKE_TRANSIENT(averagingPeriod); 
		FIELD(averagingFrequency, ""); 
        FIELD_MAKE_TRANSIENT(averagingFrequency); 
        FIELD(avgDaysB4End, "");
        FIELD_MAKE_TRANSIENT(avgDaysB4End); 
        FIELD(settleDays,"");
        FIELD_MAKE_TRANSIENT(settleDays);
        FIELD(rollDay,  "");
        FIELD_MAKE_TRANSIENT(rollDay);
        FIELD(hasStub, ""); 
        FIELD_MAKE_TRANSIENT(hasStub); 

		FIELD(startDate,           "");
		FIELD_MAKE_TRANSIENT(startDate);
        FIELD(endDate,             "");
		FIELD_MAKE_TRANSIENT(endDate);
	    FIELD(rollDate,            "");
		FIELD_MAKE_TRANSIENT(rollDate);
		FIELD(alignDate,            "");
		FIELD_MAKE_TRANSIENT(alignDate);
		FIELD(notionalSetting,     "");
        FIELD_MAKE_TRANSIENT(notionalSetting);  
		FIELD(currentPeriodStart,"");
		FIELD_MAKE_TRANSIENT(currentPeriodStart);
		/**** below not working, so the clone overwrite in this class
		FIELD(coupons,"");
		FIELD_MAKE_TRANSIENT(coupons);
		***/

    }

    static IObject* defaultEnergyStreamSchedule()
    {
        return new EnergyStreamSchedule();
    }
};

CClassConstSP const EnergyStreamSchedule::TYPE = CClass::registerClassLoadMethod(
	 "EnergyStreamSchedule", typeid(EnergyStreamSchedule), EnergyStreamScheduleHelper::load);


// ********************************************************************************
// FIXED
//*********************************************************************************

EnergyStreamScheduleFixed::EnergyStreamScheduleFixed(
    const EnergyUnderlyerSP&   energyUnderlyer,
	const DateTime&            today,
    const DateTime&            startDate,
    const DateTime&            endDate,
	int                        settleDays,
	const EnergyDWMYPeriod&    paymentFrequency,  // swaplet period
	int                        avgDaysB4End, // 0, 1, ...days from END
	const EnergyDWMYPeriod&    averagingPeriod,
	const EnergyDWMYPeriod&    averagingFrequency,
    int                        rollDay,
	const string&              notionalSetting)
	: 
    EnergyStreamSchedule(TYPE,energyUnderlyer, today, startDate, endDate,settleDays,
                         paymentFrequency,avgDaysB4End,averagingPeriod,
						 averagingFrequency, rollDay, notionalSetting), rate(0.)
{
}


EnergySwapCoupon* EnergyStreamScheduleFixed::createCouponPayment(
	const DateTime& /*today*/,
	const EnergyUnderlyerSP& energyUnderlyer,
	const DateTime& couponStart, 
	const DateTime& couponEnd, 
	const DateTime& notionalStart,
	const DateTime& notionalEnd,
	const DateTime& cashflowDate,
	const EnergyDWMYPeriod& averagingPeriod,
	const EnergyDWMYPeriod& averagingFrequency,
	int   avgDaysB4End )
{
	const static string method = "EnergyStreamScheduleFixed::createCouponPayment";

	EnergySwapCoupon* out = new EnergySwapCouponFixed(
	                             energyUnderlyer, couponStart, couponEnd, notionalStart, 
		                         notionalEnd, cashflowDate, averagingPeriod,
								 averagingFrequency, avgDaysB4End );

	if (!out)
		throw ModelException(method,"Error building energy swap coupon");

	return out;
}

double EnergyStreamScheduleFixed::getPV(const EnergyFuturesCurveSP& energyFuturesCurve,
		                      const YieldCurveSP& yieldCurve, 
			                  double fixedRate,
		                      double notionalAmount, 
							  const string& notionalType,
							  const string& dealDirection)
{
	static const string method = "EnergyStreamScheduleFixed::getPV";

	string underlyerCurrency = energyUnderlyer->getPricingCurrency();
    if (underlyerCurrency.size() == 0 )
        underlyerCurrency = yieldCurve->getCcy();
	double localNotional = notionalAmount * (dealDirection == "PAYFIXED"?-1.0:1.0); 
	double result = 0;
	double couponAmount;
	
	DateTime valueDate = today;
	DateTime cashflowDate;

    rate = fixedRate;


	if (Maths::isZero(fixedRate))
		throw ModelException(method,"Energy swap must have fixed rate in order to value instrument");
		
	double df;
	for (unsigned i=0; i<coupons.size(); ++i)
	{
		cashflowDate = coupons[i]->getCashflowDate();

		if ( cashflowDate.daysDiff(valueDate)>=0)
		{
			couponAmount = coupons[i]->getCashflowAmount(fixedRate,
				                                         localNotional, notionalType);
			
			// cashflows will always be in pricing currencys.
			// DF = yeildCurve->pv(valueDate, cashflowDate)
			df = yieldCurve->pv(valueDate, cashflowDate);
		
			couponAmount = couponAmount* df;
			
			result += couponAmount;
		}
	}

	return result;
}


EnergyStreamScheduleFixed::EnergyStreamScheduleFixed() : EnergyStreamSchedule(TYPE), rate(0.0)
{
}


EnergyStreamScheduleFixed::~EnergyStreamScheduleFixed()
{
}

void EnergyStreamScheduleFixed::getScheduleDetails(
		StringArray& thePaymentDates,
		StringArray& theCouponStartDates,
		StringArray& theCouponEndDates,
		StringArray& theNotionalStartDates,
		StringArray& theNotionalEndDates,
		IntArray& theNumFixings,
		StringArray& theFixingDates,
		StringArray& theFixingLabels,
		DoubleArray& theFixingRates ) const
{
    static const string method = "EnergyStreamScheduleFixed::getScheduleDetails";

	DateTime d1, d2, d3, d4, d5;

	for (unsigned i=0; i<coupons.size(); ++i)
	{
		d1 = coupons[i]->getCashflowDate();
		d2 = coupons[i]->getCouponStartDate();
		d3 = coupons[i]->getCouponEndDate();
		d4 = coupons[i]->getNotionalStartDate();
		d5 = coupons[i]->getNotionalEndDate();
	
		thePaymentDates.push_back(d1.dateFormat(d1.getDate()));
		theCouponStartDates.push_back(d2.dateFormat(d2.getDate()));
		theCouponEndDates.push_back(d3.dateFormat(d3.getDate()));
		theNotionalStartDates.push_back(d4.dateFormat(d4.getDate()));
		theNotionalEndDates.push_back(d5.dateFormat(d5.getDate()));
		theFixingRates.push_back(rate);
	}

}

class EnergyStreamScheduleFixedHelper
{
  
public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergyStreamScheduleFixed, clazz);
        SUPERCLASS(EnergyStreamSchedule);
        EMPTY_SHELL_METHOD(defaultEnergyStreamScheduleFixed);

    }

    static IObject* defaultEnergyStreamScheduleFixed()
    {
        return new EnergyStreamScheduleFixed();
    }
};

CClassConstSP const EnergyStreamScheduleFixed::TYPE = CClass::registerClassLoadMethod(
	 "EnergyStreamScheduleFixed", typeid(EnergyStreamScheduleFixed), EnergyStreamScheduleFixedHelper::load);

//****************************************************************************
// FLOATING STREAM .......
//***************************************************************************

EnergyStreamScheduleFloating::EnergyStreamScheduleFloating(
    const EnergyUnderlyerSP&   energyUnderlyer,
    const DateTime&            today,
    const DateTime&            startDate,
    const DateTime&            endDate,
    int                        settleDays,
	const EnergyDWMYPeriod&    paymentFrequency,  // swaplet period
	int                        avgDaysB4End, // 0,1,2... days before END
	const EnergyDWMYPeriod&    averagingPeriod,
	const EnergyDWMYPeriod&    averagingFrequency,
    int                        rollDay,
	int                        nearbyRel,
    const EnergyContractLabel& nearbyAbsLabel,
	const string&              notionalSetting)
	: 
    EnergyStreamSchedule(TYPE, energyUnderlyer, today, startDate, endDate,settleDays,
                         paymentFrequency,avgDaysB4End,averagingPeriod,
						 averagingFrequency, rollDay,notionalSetting),
		   nearbyRel(nearbyRel), nearbyAbsLabel(nearbyAbsLabel)
{    
}

IObject* EnergyStreamScheduleFloating::clone() const
{
    static const string routine = "EnergyStreamScheduleFloating::clone";
    try {
        IObject* Copy = EnergyStreamSchedule::clone(); // call parent
        EnergyStreamScheduleFloating* copy = dynamic_cast<EnergyStreamScheduleFloating*>(Copy);
        if(!copy) {
                throw ModelException("Clone method failed");
        }

        copy->nearbyRel = nearbyRel;
        copy->nearbyAbsLabel = nearbyAbsLabel;
 
        return copy;
    } catch(exception& e) {
            throw ModelException(e, routine);
    }
}

EnergySwapCoupon* EnergyStreamScheduleFloating::createCouponPayment(
	const DateTime& today,  // this is the today date of when the swap was created
	const EnergyUnderlyerSP& energyUnderlyer,
	const DateTime& couponStart, 
	const DateTime& couponEnd, 
	const DateTime& notionalStart,
	const DateTime& notionalEnd,
	const DateTime& cashflowDate,
	const EnergyDWMYPeriod& averagingPeriod,
	const EnergyDWMYPeriod& averagingFrequency,
	int   avgDaysB4End )
{
	const static string method = "EnergyStreamScheduleFloating::createCouponPayment";

	EnergySwapCoupon* out;

	EnergyFixing energyFixing;
	energyFixing.setNearbyRel(nearbyRel);
	energyFixing.setEnergyUnderlyer(energyUnderlyer);

	if (nearbyAbsLabel.isValid())
	{
		// Absolute contract label
		energyFixing.setContractLabelForFixing(nearbyAbsLabel);
	}
	
	// 0 means no nearby
	out = new EnergySwapCouponFloating(energyUnderlyer, couponStart, couponEnd, notionalStart, notionalEnd,
		cashflowDate, averagingPeriod, averagingFrequency, avgDaysB4End, energyFixing );

	if (!out)
		throw ModelException(method,"Error building Energy swap coupon");

	return out;
}


double EnergyStreamScheduleFloating::getPV(const EnergyFuturesCurveSP& energyFuturesCurve,
		                                           const YieldCurveSP& yieldCurve, 
								                   double notionalAmount, 
												   const string& notionalType,
								                   const string& dealDirection,
								                   double pastAverage, 
								                   bool pastAverageInclToday)
{
	const static string method = "EnergyStreamScheduleFloating::GetCommoidtyPV";

	string underlyerCurrency = energyUnderlyer->getPricingCurrency();
	double localNotional = notionalAmount * (dealDirection == "PAYFIXED"?1.0:-1.0); 
	double result = 0;
	double couponAmount;
	double couponRate;

    DateTime valueDate = today;
	DateTime cashflowDate;

	double df;
	
	for (unsigned i=0; i<coupons.size(); ++i)
	{
		cashflowDate = coupons[i]->getCashflowDate();
		if ( cashflowDate.daysDiff(valueDate)>=0)
		{
		    couponRate = coupons[i]->getRate(valueDate, energyFuturesCurve, underlyerCurrency, 
			                                       pastAverage, pastAverageInclToday);

			couponAmount = coupons[i]->getCashflowAmount(couponRate, localNotional,
				                                         notionalType);

			// cashflows will always be in pricing currencys.
			// DF = yeildCurve->pv(valueDate, cashflowDate)
			df = yieldCurve->pv(valueDate, cashflowDate);
			
			couponAmount = couponAmount* df;
			

			result += couponAmount;
		}
	}

    return result;
}


EnergyStreamScheduleFloating::EnergyStreamScheduleFloating() : EnergyStreamSchedule(TYPE)
{
}

EnergyStreamScheduleFloating::~EnergyStreamScheduleFloating()
{
}

void EnergyStreamScheduleFloating::getScheduleDetails(
		StringArray& thePaymentDates,
		StringArray& theCouponStartDates,
		StringArray& theCouponEndDates,
		StringArray& theNotionalStartDates,
		StringArray& theNotionalEndDates,
		IntArray& theNumFixings,
		StringArray& theFixingDates,
		StringArray& theFixingLabels,
		DoubleArray& theFixingRates) const
{
    static const string method = "EnergyStreamScheduleFloating::getScheduleDetails";

	unsigned numFix;
	DateTime d1, d2, d3, d4, d5, d6;
	for (unsigned i=0; i<coupons.size(); ++i)
	{
		d1 = coupons[i]->getCashflowDate();
		d2 = coupons[i]->getCouponStartDate();
		d3 = coupons[i]->getCouponEndDate();
		d4 = coupons[i]->getNotionalStartDate();
		d5 = coupons[i]->getNotionalEndDate();

		numFix = coupons[i]->getNumFixings();
		for (unsigned j=0; j<numFix; j++)
		{
			d6 = coupons[i]->getFixingDate(j);

		    thePaymentDates.push_back(d1.dateFormat(d1.getDate()));
		    theCouponStartDates.push_back(d2.dateFormat(d2.getDate()));
		    theCouponEndDates.push_back(d3.dateFormat(d3.getDate()));
		    theNotionalStartDates.push_back(d4.dateFormat(d4.getDate()));
		    theNotionalEndDates.push_back(d5.dateFormat(d5.getDate()));
		    theNumFixings.push_back(numFix);
		
		    theFixingDates.push_back(d6.dateFormat(d6.getDate()));
			theFixingLabels.push_back(coupons[i]->getFixingLabel(j));
			theFixingRates.push_back(coupons[i]->getFixingRate(j));
		}

	}	

}

class EnergyStreamScheduleFloatingHelper
{
  
public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergyStreamScheduleFloating, clazz);
        SUPERCLASS(EnergyStreamSchedule);
        EMPTY_SHELL_METHOD(defaultEnergyStreamScheduleFloating);

		FIELD(nearbyRel, "");
		FIELD_MAKE_TRANSIENT(nearbyRel);
        FIELD(nearbyAbsLabel, "");
		FIELD_MAKE_TRANSIENT(nearbyAbsLabel);

    }

    static IObject* defaultEnergyStreamScheduleFloating()
    {
        return new EnergyStreamScheduleFloating();
    }
};

CClassConstSP const EnergyStreamScheduleFloating::TYPE = CClass::registerClassLoadMethod(
 "EnergyStreamScheduleFloating", typeid(EnergyStreamScheduleFloating), EnergyStreamScheduleFloatingHelper::load);

DRLIB_END_NAMESPACE



