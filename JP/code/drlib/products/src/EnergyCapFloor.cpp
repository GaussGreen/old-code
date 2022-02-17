//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyCapFloor.cpp
//
//   Description : Energy Cap/Floor instrument
//
//   Author      : Sean Chen
//
//   Date        : Dec. 23, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/EnergyCapFloor.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Black.hpp"


DRLIB_BEGIN_NAMESPACE


void EnergyCapFloor::validatePop2Object()
{
    static const string method("EnergyCapFloor::validatePop2Object");

	if (Maths::isZero(strike))
		strike = energySwap->getStrike();
	if (Maths::isZero(notional))
		notional = energySwap->getNotional();
}

void EnergyCapFloor::GetMarket(const IModel*   model, 
                              const CMarketDataSP  market)
{
    static const string method("EnergyCapFloor::GetMarket");

    market->GetReferenceDate(valueDate);
	energySwap->GetMarket(model, market);

    volData.getData(model, market);

}

DateTime EnergyCapFloor::getValueDate() const
{
    return valueDate;
}

/** when to stop tweaking */
DateTime EnergyCapFloor::endDate(const Sensitivity* sensControl) const 
{ 
    return energySwap->endDate(sensControl);
}

/** when to begin tweaking */
DateTime EnergyCapFloor::beginDate(const SensControl* sensControl) const 
{
    return valueDate;
}

string EnergyCapFloor::discountYieldCurveName() const
{
	return "USD";
	//return energySwap->getYieldCurve()->getName();
}

EnergyCapFloor::EnergyCapFloor(): CInstrument(TYPE), strike(0.0), notional(0.0)
{    
}

EnergyCapFloor::~EnergyCapFloor()
{    
}

void EnergyCapFloor::Validate()
{    
}

void EnergyCapFloor::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        // IND_VOL
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) 
        {
            results->storeRequestResult(request, indVol); 
        }
    }
}

pair<double,double> EnergyCapFloor::price(bool isSmileOff) const
{
    static const string method = "EnergyCapFloor::price";

    double premium;
    double totalPremium = 0.0;

    // instrument parameters
    double actualPastAverage;
    int numberPastFixings = 0;
    
    double optionStrike = strike;
    string notionalType = energySwap->getNotionalType();
    EnergySwapScheduleBaseSP swapSchedule = energySwap->getSwapSchedule();
	
    int buySellAdjuster = isBuy ? 1 : -1;
    
    // market parameters
    EnergyFuturesCurveConstSP futuresCurve = energySwap->getFuturesCurve();
    
    DateTime today = energySwap->getValueDate();
    YieldCurveConstSP currencyCurve = energySwap->getYieldCurve();
    
    // pricing parameters
    int numFixings;
    double ratio;
    DateTime expiryDate;
    int numberFutures;
    double newOptionStrike;
    double discountFactor;
    double averageFuture;
    double averageVol;
    double yFrac;
    int i;

    double localNotional = notional;

    EnergyUnderlyerConstSP energyUnderlyer = swapSchedule->getEnergyUnderlyer();

    EnergyStreamScheduleFloatingConstSP streamSchedule = swapSchedule->getFloatingEnergyStreamSchedule();
    vector<EnergySwapCouponSP> paymentList = streamSchedule->getCoupons();

    int nCoupons = paymentList.size();

    // check if any coupons have finished fixing, but not paid out
    int firstFixingCoupon=0;
    DateTime paymentDate;

    for (i = 0; i < nCoupons; i++)
    {
        paymentDate = paymentList[i]->getCashflowDate();
        if (!today.isGreater(paymentDate))
        {
            discountFactor = currencyCurve->pv(today, paymentDate);

            numFixings = paymentList[i]->getNumFixings();
            double numPastFixings = paymentList[i]->getNumPastFixings();
            double numFutureFixings =  numFixings - numPastFixings;

            if (numFutureFixings == 0)
            {
                // today is between last fixing and payment dates
                double fixedPrice;
                firstFixingCoupon++;

                averageFuture = paymentList[i]->getActualPastAverage();// same as average?
            
                if (isCap)
                {
                    fixedPrice = Maths::max(averageFuture - optionStrike, 0.0);
                }
                else
                {
                    fixedPrice = Maths::max(optionStrike - averageFuture, 0.0);
                }

                fixedPrice *= discountFactor;

                if ( notionalType =="COUPON")
                {
                    localNotional = notional;
                }
                else if (notionalType == "FIXING")
                {
                    localNotional =notional * numFixings;
                }
                else if (notionalType == "DAILY")
                {    
                    DateTime startDate = paymentList[i]->getNotionalStartDate();
                    DateTime endDate = paymentList[i]->getNotionalEndDate(); 
                    localNotional = notional * ( endDate.daysDiff(startDate) + 1);
                }
                else if (notionalType == "TOTAL")
                {
                    localNotional = notional / nCoupons;                      
                }
                else
                {
                    throw ModelException(method, "Notional Type not implemented");
                }

                totalPremium += fixedPrice*localNotional*buySellAdjuster;
        
                // shouldn't we just return now???? coupons are in ascending order.
                return make_pair(totalPremium, 0.0);
            }
        }
        else
            firstFixingCoupon++; // this was not in fxdr code, which is a bug I think!!! Skipping coupons that have been paid.
    }

    
    // all of these coupons have >=1 future fixings and hence have volatility
    for (i = firstFixingCoupon; i < nCoupons; i++)
    {
        paymentDate = paymentList[i]->getCashflowDate();
        
        numFixings = paymentList[i]->getNumFixings();
        int numPastFixings = paymentList[i]->getNumPastFixings(); // this will never happen from swap coupon code
        int numFutureFixings = numFixings - numPastFixings;
        
        ratio = (/*static_cast*/double) numFutureFixings / (/*static_cast*/double) numFixings;

        expiryDate = paymentList[i]->getFixingDate(numFixings-1);        
    
        if (numPastFixings)
        {
            newOptionStrike = (numFixings*optionStrike - numPastFixings*paymentList[i]->getActualPastAverage()) / numFutureFixings;
        }
        else
            newOptionStrike = optionStrike;

        discountFactor = currencyCurve->pv(today, paymentDate);
        yFrac = expiryDate.daysDiff(today) / 365.0;
        
        string tempLabel;
        string prevLabel("");
        double tempRate;
		double tempVol;
        averageFuture = 0.0;
        averageVol = 0.0;
        DateTime futuresExpiry;
        DateTime futuresMaturity;

		vector<int> counts;
	    vector<double> futures;
	    vector<double> vols;
	    vector<double> t;
		int counter = 0;
		DateTime fixDate;


        for (int j = 0; j < numFixings; j++)
        {
			fixDate = paymentList[i]->getFixingDate(j);
			t.push_back(fixDate.daysDiff(today) / 365.0);
            tempLabel = paymentList[i]->getFixingLabel(j);
                
            if ( tempLabel != prevLabel)
            {
				prevLabel = tempLabel;
				
				if(j == 0)
					counter = 1;
				else
				{
					// saving previous values...
					counts.push_back(counter);
					futures.push_back(tempRate);
                    vols.push_back(tempVol);
					counter = 1;
				}
				

                //tempRate = paymentList[i]->getFixingRate(j); not this static rate in original swap
		        tempRate = futuresCurve->fixing(tempLabel);
				
                futuresExpiry = energyUnderlyer->optionExpiryDate(tempLabel);
                futuresMaturity = energyUnderlyer->expiryDate(tempLabel);
                
                // If contract has expired then need to look up maturity date rather than expiry date
                if ( !futuresExpiry.isGreater(today))
                    futuresExpiry = futuresMaturity;

                if ( isSmileOff )
                    tempVol = volData->getATMVol(futuresExpiry);
                else
                    tempVol = volData->getSmileVolByStrike(futuresExpiry, optionStrike);
            }
			else
				counter++;

            averageFuture += tempRate;
       
        }    

		counts.push_back(counter);
	    futures.push_back(tempRate);
        vols.push_back(tempVol);

        averageFuture /= numFixings;

		if ( counts.size() == 1)
            averageVol = oneDCalculateBasketVol( counts, futures, vols, t);
		else if (counts.size() == 2)
            averageVol = twoDCalculateBasketVol( counts, futures, vols, t);
		else
            averageVol = genericCalculateBasketVol( counts, futures, vols, t);

        //Note, Black takes variance
        premium = Black::price( isCap, averageFuture, 
                                newOptionStrike, 1, averageVol*averageVol*yFrac);
        premium *= discountFactor;

        if(notionalType == "COUPON")
        {
            localNotional = notional;
        }
        else if (notionalType == "FIXING")
        {
            localNotional = notional * numFixings;
        }
        else if (notionalType == "DAILY")
        {    // +1 because notional start is inclusive
            DateTime startDate = paymentList[i]->getNotionalStartDate();
            DateTime endDate = paymentList[i]->getNotionalEndDate(); 
            localNotional = notional * ( endDate.daysDiff(startDate) + 1);
        }
        else if ( notionalType == "TOTAL")
        {
            localNotional = notional / nCoupons;  
        }
        else
        {
            throw ModelException(method, "Notional Type not implemented");
        }

        totalPremium += ratio*premium*buySellAdjuster*localNotional;
    }    

    return make_pair(totalPremium, averageVol);

}

double EnergyCapFloor::oneDCalculateBasketVol(
	const vector<int>& weights,
	const vector<double>& rates,
	const vector<double>& vols,
	const vector<double>& t) const
{
	int i;

	int n = weights[0];
	double f = rates[0];
	double v = vols[0];

	double sum = 0.0;

	for (i=0; i<n; i++)
		sum += (2.0*(n-i)-1.0) * exp(v*v * t[i]);

	double m1 = f;
	double m2 = f*f*sum / (n*n);

	return sqrt(log(m2/(m1*m1))/t.back());
}

double EnergyCapFloor::twoDCalculateBasketVol(
	const vector<int>& weights,
	const vector<double>& rates,
	const vector<double>& vols,
	const vector<double>& t) const
{
	int i;

	int n1 = weights[0], n2 = weights[1];
	double f1 = rates[0], f2 = rates[1];
	double v1 = vols[0], v2 = vols[1];

	double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

	for (i=0; i<n1; i++)
		sum1 += (2.0*(n1-i)-1.0) * exp(v1*v1 * t[i]);

	for (i=n1; i<n1+n2; i++)
		sum2 += (2.0*(n1+n2-i)-1.0) * exp(v2*v2 * t[i]);

	for (i=0; i<n1; i++)
		sum3 += exp(v1*v2 * t[i]);

	double m1 = (n1*f1+n2*f2)/(n1+n2);
	double m2 = (f1*f1*sum1 + f2*f2*sum2 + 2.0*f1*f2*n2*sum3) / ((n1+n2)*(n1+n2));

	return sqrt(log(m2/(m1*m1))/t.back());
}


double EnergyCapFloor::genericCalculateBasketVol(
	const vector<int>& weights,
	const vector<double>& futures,
	const vector<double>& vols,
	const vector<double>& t) const
{
	double f1, f2, v1, v2, t1, t2;
	double sum1= 0.0, sum2=0.0;

	int i, j;
	int n=0;
	for (i=0; i<weights.size(); i++)
		n += weights[i];
	
	int i1 = 0;
	int s1 = weights[0];
	for (i=0; i<n; i++)
	{
		if (i>=s1)
		{
			i1++;
			s1 += weights[i1];
		}
		f1 = futures[i1];
		v1 = vols[i1];
		t1 = t[i];

		sum1 += f1;

		int i2 = 0;
		int s2 = weights[0];

		for (j=0; j<n; j++)
		{
			if (i>=s2)
			{
				i2++;
				s2 += weights[i2];
			}			
			f2 = futures[i2];
			v2 = vols[i2];
			t2 = t[j];

			sum2 += f1*f2*exp(v1*v2*Maths::min(t1,t2));
		}
	}

	return sqrt(log(sum2/(sum1*sum1))/t[n-1]);		
}



/** private class */
class EnergyCapFloorClosedForm : public ClosedFormEnergy::IProduct
{
public:

    EnergyCapFloorClosedForm(const EnergyCapFloor* cap): capFloor(cap){}

    void price(ClosedFormEnergy*   model,
               Control*        control, 
               CResults*       results) const;
    
private:

    pair<double,double> price(bool) const;
    

    const EnergyCapFloor*  capFloor;

};

/** Implementation of ClosedFormEnergy::IntoProduct interface */
ClosedFormEnergy::IProduct* EnergyCapFloor::createProduct(
                                     ClosedFormEnergy* model) const
{
    return new EnergyCapFloorClosedForm(this);
}


void EnergyCapFloorClosedForm::price(ClosedFormEnergy* model,
                               Control*       control, 
                               CResults*      results) const
{
    static const string method = "EnergyCapFloorClosedForm::price";
    
    pair<double,double> premium;

    try
    {
        // valueDate >= matDate is taken care of here
       
        premium = capFloor->price(model->isSmileOff());
      
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }

    results->storePrice(premium.first, "USD");

	if ( control->isPricing() ) 
        capFloor->addOutputRequests(control, results,
                               premium.first, premium.second);

}
    

class EnergyCapFloorHelper
{

public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyCapFloor, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedFormEnergy::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultCapFloor);

        FIELD(energySwap,        "Underlying Swap");
        FIELD(isCap,           "Is it a cap (true/false)");
        FIELD(volData,        "Energy Implied Vol Surf Wrapper");

        // Following three override on underlying instrument
        FIELD(strike,       "Contract Strike");
        FIELD_MAKE_OPTIONAL(strike);
        FIELD(notional,     "Notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(isBuy,        "Is it a Buy(true) and sell(false)");
        FIELD_MAKE_OPTIONAL(isBuy);
        
    }

    static IObject* defaultCapFloor()
    {
        return new EnergyCapFloor();
    }
};

CClassConstSP const EnergyCapFloor::TYPE = CClass::registerClassLoadMethod(
    "EnergyCapFloor", typeid(EnergyCapFloor), EnergyCapFloorHelper::load);

bool  EnergyCapFloorLoad() { return (EnergyCapFloor::TYPE != 0);   }


DRLIB_END_NAMESPACE

