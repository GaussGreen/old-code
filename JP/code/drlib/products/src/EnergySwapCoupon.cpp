//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwapCoupon.cpp
//
//   Description : Energy Swap Coupon
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//---------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EnergySwapCoupon.hpp"

DRLIB_BEGIN_NAMESPACE

EnergySwapCoupon::EnergySwapCoupon(
    CClassConstSP clazz,
    const EnergyUnderlyerSP& energyUnderlyer,
    const DateTime& couponStartDate,
    const DateTime& couponEndDate,
    const DateTime& notionalStartDate,
    const DateTime& notionalEndDate,
    const DateTime& cashflowDate,
	const EnergyDWMYPeriod& averagingPeriod, 
    const EnergyDWMYPeriod& averagingFrequency,
    int avgDaysB4End)
    : CObject(clazz),
    energyUnderlyer(energyUnderlyer), cashflowDate(cashflowDate),
    notionalStartDate(notionalStartDate), notionalEndDate(notionalEndDate),
    couponStartDate(couponStartDate), couponEndDate(couponEndDate),
    averagingPeriod(averagingPeriod),averagingFrequency(averagingFrequency),
    avgDaysB4End(avgDaysB4End)
{
}

EnergySwapCoupon::EnergySwapCoupon(CClassConstSP clazz) : CObject(clazz)
{
}

EnergySwapCoupon::EnergySwapCoupon() : CObject(TYPE)
{
}



EnergySwapCoupon::~EnergySwapCoupon()
{
}


// Would not add this if vector<A> thing working in FIELD for child class, since we always prefer
// using clone() with default implementation.
IObject* EnergySwapCoupon::clone() const {
    static const string routine = "EnergySwapCoupon::clone";
    try {
        IObject* Copy = CObject::clone(); // call parent
        EnergySwapCoupon* copy = dynamic_cast<EnergySwapCoupon*>(Copy);
        if(!copy) {
                throw ModelException("Clone method failed");
        }
		copy->energyUnderlyer = energyUnderlyer;
        copy->cashflowDate = cashflowDate;
        copy->couponStartDate = couponStartDate;
        copy->couponEndDate = couponEndDate;
        copy->notionalStartDate = notionalStartDate;
        copy->notionalEndDate = notionalEndDate;
        copy->numFixings = numFixings;

        return copy;
    } catch(exception& e) {
            throw ModelException(e, routine);
    }
}

const EnergyUnderlyerSP& EnergySwapCoupon::getEnergyUnderlyer() const
{
    return energyUnderlyer;
}

DateTime EnergySwapCoupon::getCashflowDate() const
{
    return cashflowDate;
}

DateTime EnergySwapCoupon::getNotionalStartDate() const
{
    return notionalStartDate;
}

DateTime EnergySwapCoupon::getNotionalEndDate() const
{
    return notionalEndDate;
}

DateTime EnergySwapCoupon::getCouponStartDate() const
{
    return couponStartDate;
}

DateTime EnergySwapCoupon::getCouponEndDate() const
{
    return couponEndDate;
}

int EnergySwapCoupon::getNumFixings() const
{
    return numFixings;
}

int EnergySwapCoupon::getNumTotalCoupons() const
{
    return numTotalCoupons;
}


// Fixed coupon will come to here, while floating one will go to the child class

double EnergySwapCoupon::getCashflowAmount(
    double couponRate,
    double notional,
    const string& notionalType) const
{
    
	static const string method = "EnergySwapCoupon::GetCashflowAmount";

    if ( notionalType =="COUPON")
        couponRate *= notional;
    else if (notionalType =="DAILY")
	{
        // +1 because notional start is inclusive
        DateTime startDate = getNotionalStartDate();
        DateTime endDate = getNotionalEndDate();
            
        int numDays = endDate.daysDiff(startDate);
        numDays += 1;  // because 
        couponRate *= numDays * notional;
    }
	else if (notionalType =="FIXING")
	{
        couponRate *= getNumFixings() * notional;
    }
    else if (notionalType =="TOTAL")
	{
        if ( numTotalCoupons < 1)
            throw ModelException(method,"Total number of coupons is incorrect");
        couponRate *= (double) notional/numTotalCoupons; 
    }
    
    return couponRate;
}


void EnergySwapCoupon::setNumFixings(int theNumFixings)
{
    numFixings = theNumFixings;
}

void EnergySwapCoupon::setNumTotalCoupons(int theTotalCoupons)
{
    numTotalCoupons = theTotalCoupons;
}

class EnergySwapCouponHelper
{
  
public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergySwapCoupon, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEnergySwapCoupon);

        FIELD(energyUnderlyer, "");
		FIELD_MAKE_TRANSIENT(energyUnderlyer);

        FIELD(couponStartDate, "");
        FIELD_MAKE_TRANSIENT(couponStartDate); 
        FIELD(couponEndDate,"");
        FIELD_MAKE_TRANSIENT(couponEndDate);
        FIELD(notionalStartDate,  "");
        FIELD_MAKE_TRANSIENT(notionalStartDate);
        FIELD(notionalEndDate, ""); 
        FIELD_MAKE_TRANSIENT(notionalEndDate); 
		FIELD(cashflowDate,           "");
		FIELD_MAKE_TRANSIENT(cashflowDate);
        FIELD(numFixings,             "");
		FIELD_MAKE_TRANSIENT(numFixings);
	    FIELD(numTotalCoupons,            "");
		FIELD_MAKE_TRANSIENT(numTotalCoupons);
		FIELD(averagingPeriod,            "");
		FIELD_MAKE_TRANSIENT(averagingPeriod);
		FIELD(averagingFrequency,     "");
        FIELD_MAKE_TRANSIENT(averagingFrequency);  
		FIELD(avgDaysB4End,"");
		FIELD_MAKE_TRANSIENT(avgDaysB4End);
	
    }

    static IObject* defaultEnergySwapCoupon()
    {
        return new EnergySwapCoupon();
    }
};

CClassConstSP const EnergySwapCoupon::TYPE = CClass::registerClassLoadMethod(
			"EnergySwapCoupon", typeid(EnergySwapCoupon), EnergySwapCouponHelper::load);

//**************************************************************************************
// FIXED LEG
//**************************************************************************************

EnergySwapCouponFixed::EnergySwapCouponFixed(
    const EnergyUnderlyerSP& energyUnderlyer,
    const DateTime& couponStart,
    const DateTime& couponEnd,
    const DateTime& notionalStart, 
    const DateTime& notionalEnd, 
    const DateTime& cashflowDate,
    const EnergyDWMYPeriod& averagingPeriod,
    const EnergyDWMYPeriod& averagingFrequency,
    int avgDaysB4End)
    :
    EnergySwapCoupon(TYPE, energyUnderlyer, couponStart, couponEnd, notionalStart, 
                     notionalEnd, cashflowDate,averagingPeriod,
					 averagingFrequency,avgDaysB4End )
{
   // calculate numfixings up front as notional calculations may require notional/Fixing
    calculateNumFixings();
}

EnergySwapCouponFixed::EnergySwapCouponFixed() : EnergySwapCoupon(TYPE)
{
}

EnergySwapCouponFixed::~EnergySwapCouponFixed()
{
}


void EnergySwapCouponFixed::calculateNumFixings()
{
    DateTime startDate = getCouponStartDate();
    DateTime endDate = getCouponEndDate();
    
    int theNumFixings = 1;

	// Note, averaging window ends 'avgDaysB4End' before end of the coupon.
    
    if (averagingFrequency.getCount() == 1)
	{
		// if count=1, convention is averaging only on last day
        setNumFixings(theNumFixings);
        return;
    }
            
    DateTime toDate;
    DateTime fromDate = 
		    energyUnderlyer->busDayAdjustPrior(endDate.rollDate(-1*avgDaysB4End));

    // if averaging Period not set (i.e. count=0), assume averaging over full coupon
	int count = averagingPeriod.getCount();
    if (count)
	{
        if(averagingPeriod.getInterval() == "M")
            toDate = energyUnderlyer->addIBORMonthsToDate(fromDate,-1*count);
		else
            toDate = energyUnderlyer->addDaysToDateAndAdjust(fromDate,-1*count);
	}
    else
    {
        toDate = startDate;  
		toDate = toDate.rollDate(-1);
    }
 
	

    DateTime fixingDate = fromDate;
	DateTime nextFixingDate;
	int i = 1;
    // work backwards from the end of the coupon until the start of the averaging window is reached
    for (;;)
    {
		count = averagingFrequency.getCount()*i++;
		if(averagingFrequency.getInterval() == "M")
            nextFixingDate = energyUnderlyer->addIBORMonthsToDate(fromDate,-1*count);
		else
            nextFixingDate = energyUnderlyer->addDaysToDateAndAdjust(fromDate,-1*count);
           
        if (nextFixingDate.daysDiff(toDate)>0 && nextFixingDate.daysDiff(startDate)>=0)
            theNumFixings++;
        else
            break;  // exceeded averaging window, so leave loop
    }
    
    setNumFixings(theNumFixings);
}

class EnergySwapCouponFixedHelper
{
  
public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergySwapCouponFixed, clazz);
        SUPERCLASS(EnergySwapCoupon);
        EMPTY_SHELL_METHOD(defaultEnergySwapCouponFixed);
    }

    static IObject* defaultEnergySwapCouponFixed()
    {
        return new EnergySwapCouponFixed();
    }
};

CClassConstSP const EnergySwapCouponFixed::TYPE = CClass::registerClassLoadMethod(
			"EnergySwapCouponFixed", typeid(EnergySwapCouponFixed), EnergySwapCouponFixedHelper::load);


//***********************************************************
// FLOATING
//***********************************************************

EnergySwapCouponFloating::EnergySwapCouponFloating(
    const EnergyUnderlyerSP& energyUnderlyer,
    const DateTime& couponStart, 
    const DateTime& couponEnd,
    const DateTime& notionalStart,
    const DateTime& notionalEnd,
    const DateTime& cashflowDate,
    const EnergyDWMYPeriod& averagingPeriod, 
    const EnergyDWMYPeriod& averagingFrequency,
    int avgDaysB4End,
	const EnergyFixing& fixing)
    : 
    EnergySwapCoupon(TYPE, energyUnderlyer, couponStart, couponEnd, notionalStart, 
                     notionalEnd, cashflowDate,averagingPeriod,
					 averagingFrequency,avgDaysB4End ), numPastFixings(0), actualPastAverage(0.0)
{
    
    calculateFixingDates(fixing);
    calculateNumFixings();
}

EnergySwapCouponFloating::EnergySwapCouponFloating() : EnergySwapCoupon(TYPE), 
                                                       numPastFixings(0), actualPastAverage(0.0)
{
}

EnergySwapCouponFloating::~EnergySwapCouponFloating()
{
}

int EnergySwapCouponFloating::getNumPastFixings() const
{
    return numPastFixings;
}

double EnergySwapCouponFloating::getActualPastAverage() const
{
    return actualPastAverage;
}

IObject* EnergySwapCouponFloating::clone() const
{
    static const string routine = "EnergySwapCouponFloating::clone";
    try {
        IObject* Copy = EnergySwapCoupon::clone(); // call parent
        EnergySwapCouponFloating* copy = dynamic_cast<EnergySwapCouponFloating*>(Copy);
        if(!copy) {
                throw ModelException("Clone method failed");
        }
		copy->fixings = fixings;
        
        return copy;
    } catch(exception& e) {
            throw ModelException(e, routine);
    }
}

void EnergySwapCouponFloating::calculateFixingDates(const EnergyFixing& fixing)
{
  
	DateTime startDate = getCouponStartDate();
    DateTime endDate = getCouponEndDate();
    DateTimeArray fixingDates;  // local scope

    int theNumFixings = 1;


	// Note, averaging window ends 'avgDaysB4End' before end of the coupon.
	DateTime fromDate = 
		    energyUnderlyer->busDayAdjustPrior(endDate.rollDate(-1*avgDaysB4End));

	// traversing the averaging period from the END
    fixingDates.push_back(fromDate);

    if (averagingFrequency.getCount() == 0)
	{
		// if count=0, convention is averaging only on last day
		fixings.push_back(createFixing(fixing, fixingDates[0]));
        return;
    }
            
    DateTime toDate;
    
    // if averaging Period not set (i.e. count=0), assume averaging over full coupon
	int count = averagingPeriod.getCount();

    if (count)
	{
		
        if(averagingPeriod.getInterval() == "M")
            toDate = energyUnderlyer->addIBORMonthsToDate(fromDate,-1*count);
		else
            toDate = energyUnderlyer->addDaysToDateAndAdjust(fromDate,-1*count);
	}
    else
    {
        toDate = startDate;
        toDate = toDate.rollDate(-1);
    }
 

    DateTime fixingDate = fromDate;
	DateTime nextFixingDate = fromDate;
	count = averagingFrequency.getCount();
	int mCount;
	int i = 1;

    // work backwards from the end of the coupon until the start of the averaging window is reached
    for (;;)
    {
		mCount = count*i++;

		if(averagingFrequency.getInterval() == "M")
            nextFixingDate = energyUnderlyer->addIBORMonthsToDate(fromDate,-1*mCount);
		else
            nextFixingDate = energyUnderlyer->addDaysToDateAndAdjust(nextFixingDate,-1*count);
           
        if (nextFixingDate.daysDiff(toDate)>0 && nextFixingDate.daysDiff(startDate)>=0)
            fixingDates.push_back(nextFixingDate);
        else
            break;  // exceeded averaging window, so leave loop
    }
   
	// construct energy fixing, using input as a model and adding fixing dates
    // dates in fixingDates array are in reverse order, so put them into chronological order as well
    for (int j = fixingDates.size()-1; j >= 0; j--)
    {
        fixings.push_back(createFixing(fixing, fixingDates[j]));
    }
}


EnergyFixing EnergySwapCouponFloating::createFixing(
    const EnergyFixing& fixing, 
    const DateTime& fixingDate) const
{
    EnergyFixing localFixing = fixing;
    localFixing.setFixingDate(fixingDate);

	if (localFixing.getContractLabelForFixing().isValid())
	{
		// each fixing in this coupon observes rate for this future contract
        return localFixing;
    }


    if (localFixing.getNearbyRel())
	{
		// 'nearby'th Label relative to current contract represented by fixingDate
		localFixing.setContractLabelForFixing(
			           energyUnderlyer->calculateContractLabel(
					          fixingDate,
							  localFixing.getNearbyRel(),
							  EnergyUnderlyer::INDEX));  
	}
	else
	{
		throw ModelException("EnergySwapCouponFloating::createFixing", "Nearby must be >= 1");
	}
    

    return localFixing;
}

void EnergySwapCouponFloating::calculateNumFixings()
{
    setNumFixings(fixings.size());
}

double EnergySwapCouponFloating::getRate(const DateTime today,
										 const EnergyFuturesCurveSP& energyFuturesCurve,
		                                 const string& currency, 
								         double pastAverage, 
								         bool pastAverageInclToday )
{
	static const string method = "EnergySwapCouponFloating::getRate";


	// convert rates to pricing currency and average them
	DateTime fixingDate;
	double rate = 0, tmpRate;
	bool isPastAverage = false;

	for (unsigned i = 0; i < fixings.size(); i++)
	{
		fixingDate = fixings[i].getFixingDate();
		tmpRate = determineFixingRate(today, fixings[i], energyFuturesCurve,
			                         pastAverage, pastAverageInclToday, isPastAverage);
		if(isPastAverage)
        {
			numPastFixings += 1;
	        actualPastAverage += tmpRate;
		}
		
		rate += tmpRate;
        ((EnergyFixing&)fixings[i]).setFixingRate(tmpRate);
	}
	
	rate /= fixings.size();  // calc average in the pricing currency

	return rate;  // rate
}

double EnergySwapCouponFloating::determineFixingRate(
	const DateTime& today,
    const EnergyFixing& fixing,
    const EnergyFuturesCurveSP& futuresCurve,
    const double pastAverage,
	bool pastAverageInclToday,
	bool& isPastAverage) const
{
    const string& energyCurrency = getEnergyUnderlyer()->getPricingCurrency();
    
    double fixingRate;

    if (Maths::isZero(pastAverage))
    {
		// no past average to use
		if ( !today.isGreater(fixing.getFixingDate()))
        {
	        if (fixing.getContractLabelForFixing().isValid())
                fixingRate = futuresCurve->fixing(fixing.getContractLabelForFixing().asString());
	        else
                fixingRate = futuresCurve->fixing(fixing.getFixingDate());
		}
        else
		{
            isPastAverage = true;
			// get fixed price from historical curves
			// to be implemented
        }
    }
    else
    {            
        // only use past average if coupon has not ended and payment has not been made yet
		if ( !today.isGreater(getCouponEndDate()))
        {
	        if ( today.isGreater(fixing.getFixingDate()))
            {
	            fixingRate = pastAverage;
				isPastAverage = true;
			}
            else if ( today.daysDiff(fixing.getFixingDate()) == 0 &&  pastAverageInclToday )
			{
				fixingRate = pastAverage;
				isPastAverage = true;
			}
		    else
            {
				// fixing is after today ...
                if (fixing.getContractLabelForFixing().isValid())
                    fixingRate = futuresCurve->fixing(fixing.getContractLabelForFixing().asString());
		        else
                    fixingRate = futuresCurve->fixing(fixing.getFixingDate());
            }
        }
        else
		{
            isPastAverage = true;
            // get fixed price from historical curves
			// to be implemented
		}
    }

    return fixingRate;
}

DateTime EnergySwapCouponFloating::getFixingDate(int index) const
{
	return fixings[index].getFixingDate();
}

double EnergySwapCouponFloating::getFixingRate(int index) const
{
	return fixings[index].getFixingRate();
}

string EnergySwapCouponFloating::getFixingLabel(int index) const
{
	if (fixings[index].getContractLabelForFixing().isValid())
	{
		return fixings[index].getContractLabelForFixing().asString();
	}

	return "";
}

class EnergySwapCouponFloatingHelper
{
  
public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergySwapCouponFloating, clazz);
        SUPERCLASS(EnergySwapCoupon);
        EMPTY_SHELL_METHOD(defaultEnergySwapCouponFloating);
// not working, that's why I added clone() in this class
//		FIELD(fixings, "");
//        FIELD_MAKE_TRANSIENT(fixings); 
    }

    static IObject* defaultEnergySwapCouponFloating()
    {
        return new EnergySwapCouponFloating();
    }
};

CClassConstSP const EnergySwapCouponFloating::TYPE = CClass::registerClassLoadMethod(
			"EnergySwapCouponFloating", typeid(EnergySwapCouponFloating), EnergySwapCouponFloatingHelper::load);

bool  EnergySwapCouponLoad() { return (EnergySwapCoupon::TYPE != 0);   }


DRLIB_END_NAMESPACE
