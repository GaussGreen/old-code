//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwap.cpp
//
//   Description : Energy Swap instrument
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/EnergySwap.hpp"
#include "edginc/Results.hpp"
#include "edginc/Control.hpp"

#include <cctype>
#include <algorithm>
using namespace std;


DRLIB_BEGIN_NAMESPACE


void EnergySwap::validatePop2Object()
{
    static const string method("EnergySwap::validatePop2Object");

    // notional is 1 is not supplied
    // now the number of underlying can be supplied or the number of contracts

    transform(notionalType.begin(),notionalType.end(),notionalType.begin(),(int(*)(int))toupper);
    if (!notionalType.size())
       notionalType = "TOTAL";
    else if (notionalType != "COUPON" && notionalType != "TOTAL" &&
             notionalType != "FIXING" && notionalType != "DAILY")
        throw ModelException(method, "Notional type must be COUPON, TOTAL, FIXING or DAILY");

	transform(dealDirection.begin(),dealDirection.end(),dealDirection.begin(),(int(*)(int))toupper);
    if (!dealDirection.size())
       dealDirection = "PAYFIXED";
    else if (dealDirection !="PAYFIXED" && dealDirection != "PAYFLOATING")
        throw ModelException(method, "Deal direction must be PATFIXED or PAYFLOATING");

	if ( !Maths::isZero(pastAverage))
	{
		if(pastAverageInclToday != true && pastAverageInclToday !=false)
			throw ModelException(method, "must be true or false for including today's fixing field");
	}
}


void EnergySwap::GetMarket(const IModel*   model, 
                           const CMarketDataSP  market)
{
    static const string method("EnergySwap::GetMarket");

  
    //swapScheduleSP->getMarket(model, market);

    futuresCurve.getData(model, market);

    yieldCurve.getData(model, market);

}

DateTime EnergySwap::getValueDate() const
{
    return swapScheduleSP->getValueDate();
}

/** when to stop tweaking */
DateTime EnergySwap::endDate(const Sensitivity* sensControl) const 
{ 
    return swapScheduleSP->getEndDate();
}

/** when to begin tweaking */
DateTime EnergySwap::beginDate(const SensControl* sensControl) const 
{
    return swapScheduleSP->getStartDate();
}

EnergySwapScheduleBaseSP EnergySwap::getSwapSchedule() const 
{
    return swapScheduleSP;
}

string EnergySwap::getNotionalType() const 
{
    return notionalType;
}

EnergyFuturesCurveConstSP EnergySwap::getFuturesCurve() const 
{
    return futuresCurve.getSP();
}

YieldCurveConstSP EnergySwap::getYieldCurve() const 
{
    return yieldCurve.getSP();
}

EnergySwap::EnergySwap() : CInstrument(TYPE), 
              rate(0.0), notionalAmount(1.0), pastAverage(0.0),pastAverageInclToday(true)
{    
}

EnergySwap::~EnergySwap()
{    
}


double EnergySwap::price()
{
	double value = swapScheduleSP->calculatePV(futuresCurve.getSP(), yieldCurve.getSP(),
		                                       rate, notionalAmount, notionalType, 
											   dealDirection, pastAverage,
											   pastAverageInclToday );
	return value;
}

string EnergySwap::getCcy() const
{
	return yieldCurve.getSP()->getCcy();
}

double EnergySwap::getStrike() const
{
	return rate;
}

double EnergySwap::getNotional() const 
{
	return notionalAmount;
}

void EnergySwap::Validate()
{
}

string EnergySwap::discountYieldCurveName() const
{
	return yieldCurve.getSP().get()?yieldCurve.getSP()->getName():"";
}

void EnergySwap::addMoreToResults(Control* control, Results* results)
{
    // take care of additional outputs. 
	// ObjectArraySP is the type of outputFixed and outputFloating

	outputFixed =  swapScheduleSP->getFixedLegDetails();
	outputFloating =  swapScheduleSP->getFloatingLegDetails();

    if ( results && outputFixed.get() && outputFloating.get() )
    {
        // FIXED_LEG_DETAILS
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::FIXED_LEG_DETAILS, request)) 
        {
            results->storeRequestResult(request, outputFixed); 
		}

		// FLOATING_LEG_DETAILS
		if (control->requestsOutput(OutputRequest::FLOATING_LEG_DETAILS, request)) 
        {
            results->storeRequestResult(request, outputFloating); 
		}
	}
}



/** private class */
class EnergySwapClosedForm : public ClosedFormEnergy::IProduct
{
public:

    EnergySwapClosedForm(const EnergySwap* swap): swap(swap){}

    void price(ClosedFormEnergy*   model,
               Control*        control, 
               CResults*       results) const;
    
private:

    double priceFixedStream() const;
    double priceFloatStream() const;

    const EnergySwap*  swap;

};

/** Implementation of ClosedFormEnergy::IntoProduct interface */
ClosedFormEnergy::IProduct* EnergySwap::createProduct(
                                     ClosedFormEnergy* model) const
{
    return new EnergySwapClosedForm(this);
}


void EnergySwapClosedForm::price(ClosedFormEnergy* model,
                                 Control*       control, 
                                 CResults*      results) const
{
    static const string method = "EnergySwapClosedForm::price";
    
    double value;

    try
    {
        value = ((EnergySwap*)swap)->price();
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }

    results->storePrice(value, swap->getCcy());

	if ( control->isPricing() )  
	    ((EnergySwap*)swap)->addMoreToResults(control, results);

}
    
class EnergySwapHelper
{

public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergySwap, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedFormEnergy::IIntoProduct);
        IMPLEMENTS(LastSensDate);
		IMPLEMENTS(FirstSensDate);
        EMPTY_SHELL_METHOD(defaultSwap);

        FIELD(swapScheduleSP,       "Schedule Handle");
        FIELD(futuresCurve,       "Futures Curve Wrapper");
        FIELD(yieldCurve,         "Yield Curve Wrapper");
        
        FIELD(rate,               "Rate");
        //FIELD_MAKE_OPTIONAL(rate);
        FIELD(notionalAmount,     "Notional Amount");
        //FIELD_MAKE_OPTIONAL(notionalAmount);
        // COUPON TOTAL FIXING DAILY
        FIELD(notionalType,       "Notional Type");
        FIELD_MAKE_OPTIONAL(notionalType);
        // PAYFIXED PAYFLOATING
        FIELD(dealDirection,      "Deal Direction");
        FIELD_MAKE_OPTIONAL(dealDirection);
		FIELD(pastAverage,         "Past Average Rate"); // forced on coupons > today
        FIELD_MAKE_OPTIONAL(pastAverage);
		FIELD(pastAverageInclToday,"Includes Today Fixing"); // above includes today's fixing?
        FIELD_MAKE_OPTIONAL(pastAverageInclToday); // Default is true. Mandatory if pastAverage above is given


    }

    static IObject* defaultSwap()
    {
        return new EnergySwap();
    }
};

CClassConstSP const EnergySwap::TYPE = CClass::registerClassLoadMethod(
    "EnergySwap", typeid(EnergySwap), EnergySwapHelper::load);

bool  EnergySwapLoad() { return (EnergySwap::TYPE != 0);}


DRLIB_END_NAMESPACE

