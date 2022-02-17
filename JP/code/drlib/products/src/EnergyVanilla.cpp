//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyVanilla.cpp
//
//   Description : Energy Vanilla instrument
//
//   Author      : Sean Chen
//
//   Date        : Aug. 2, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/EnergyVanilla.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "edginc/Control.hpp"
#include "edginc/Black.hpp"


DRLIB_BEGIN_NAMESPACE


void EnergyVanilla::validatePop2Object()
{
    static const string method("EnergyVanilla::validatePop2Object");
}

void EnergyVanilla::GetMarket(const IModel*   model, 
                              const CMarketDataSP  market)
{
	static const string method("EnergyVanilla::GetMarket");

    market->GetReferenceDate(valueDate);

    energyUnderlyer.getData(model, market);
    futuresCurve.getData(model, market);

    yieldCurve.getData(model, market);

    expiryDate = energyUnderlyer->expiryDate(contractLabel, 1); // option expiry date
    contractDate = energyUnderlyer->expiryDate(contractLabel,0);  // future expiry date

    if (expiryDate>= contractDate)
	{
        throw ModelException(method, "Expiry date must always be before futures maturity");
    }

    // in FXDR, no payment date needed. Why here??????
    //if (expiryDate> paymentDate){
      //  throw ModelException(method, "Expiry date must always be before or on payment date");
    //}

    volData.getData(model, market);

    //instSettle->getMarket(model, market.get());
}

DateTime EnergyVanilla::getValueDate() const
{
    return valueDate;
}

/** when to stop tweaking */
DateTime EnergyVanilla::endDate(const Sensitivity* sensControl) const 
{
    //DateTime instEnd  = instSettle->settles(matDate, asset.get());
    //DateTime assetEnd = asset->settleDate(matDate);
    
    return contractDate;
}

/** when to begin tweaking */
DateTime EnergyVanilla::beginDate(const SensControl* sensControl) const 
{
    return valueDate;
}

string EnergyVanilla::discountYieldCurveName() const
{
	return yieldCurve.getName();
}

EnergyVanilla::EnergyVanilla(): CInstrument(TYPE), notional(1)
{    
}

EnergyVanilla::~EnergyVanilla()
{    
}

void EnergyVanilla::Validate()
{    
}

void EnergyVanilla::addOutputRequests(Control* control,
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


/** private class */
class EnergyVanillaClosedForm : public ClosedFormEnergy::IProduct
{
public:

    EnergyVanillaClosedForm(const EnergyVanilla* vanilla): vanilla(vanilla){}

    void price(ClosedFormEnergy*   model,
               Control*        control, 
               CResults*       results) const;
	
private:

	pair<double,double> priceEuropean(bool) const;
	pair<double,double> priceAmerican(bool) const;
	double priceVanillaTree(
            bool isEuropean,
            bool isCall,
            const DateTime& today,
            const DateTime& expiry,
            double fwdPrice, 
            double strike,
            double vol,
            long n, /* no_time_steps*/
            long* m) const;

	const EnergyVanilla*  vanilla;

};

/** Implementation of ClosedFormEnergy::IntoProduct interface */
ClosedFormEnergy::IProduct* EnergyVanilla::createProduct(
                                     ClosedFormEnergy* model) const
{
    return new EnergyVanillaClosedForm(this);
}


void EnergyVanillaClosedForm::price(ClosedFormEnergy* model,
                               Control*       control, 
                               CResults*      results) const
{
    static const string method = "EnergyVanillaClosedForm::price";
    
	pair<double,double> premium;

    try
	{

        // valueDate >= matDate is taken care of here
        if(vanilla->isEuropean)
        {
            premium = priceEuropean(model->isSmileOff());
        }
        else
        {
            premium = priceAmerican(model->isSmileOff());
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }

    results->storePrice(premium.first, vanilla->yieldCurve->getCcy());

	vanilla->addOutputRequests(control, results,
                               premium.first, premium.second);

}
    
pair<double,double> EnergyVanillaClosedForm::priceEuropean(bool isSmileOff) const
{
    static const string method = "EnergyVanilla::priceEuropean";

	DateTime expDate = vanilla->expiryDate;
    DateTime today = vanilla->getValueDate(); 

    double premium = 0.0, vol = 0.0;
    double fwdPrice = vanilla->futuresCurve->getFwdRateForLabel(vanilla->contractLabel);

    if (today >= expDate)
    {
        // Pricing dead instrument
     
        if (vanilla->isCall && fwdPrice > vanilla->strike)
        {
            premium = fwdPrice - vanilla->strike;
        }
        else if ( !vanilla->isCall && fwdPrice<vanilla->strike)
        {
            premium = vanilla->strike - fwdPrice;
        }
    }
    else
    {
        // interpolate the vol using our LN request

        if ( isSmileOff )
            vol = vanilla->volData->getATMVol(vanilla->expiryDate);
        else
            vol = vanilla->volData->getSmileVolByStrike(vanilla->expiryDate, vanilla->strike);

        //double discFactor = instSettle->pv(vanilla->getValueDate(),
         //                                    expiryDate, 
          //                                   vanilla->discount.get(), 
          //                                   vanilla->asset.get());//?????

		double yFrac = expDate.daysDiff(today) / 365.0;
        double df = vanilla->yieldCurve->pv(today, expDate);
	    //Note, Black takes variance
        premium = Black::price( vanilla->isCall, fwdPrice, 
                                        vanilla->strike, 1, vol*vol*yFrac);
        premium *= df;
    }

    // apply buy/sell adjustment
    if (!vanilla->isBuy)
    {
        premium *= -1;
    }
    
    return make_pair(premium *= vanilla->notional, vol);
}


pair<double,double> EnergyVanillaClosedForm::priceAmerican(bool isSmileOff) const
{
    static const string method = "EnergyVanilla::priceAmerican";

    // instrument parameters
    
    //DateTime maturityDate = vanilla->energyUnderlyer->expiryDate(vanilla->contractLabel);
	DateTime expDate = vanilla->expiryDate;
    DateTime today = vanilla->getValueDate();  

    double fwdPrice = vanilla->futuresCurve->getFwdRateForLabel(vanilla->contractLabel);
	double optionStrike = vanilla->strike;

    double premium = 0.0;
	double vol;

    if (today >= expDate)
    {
        // Pricing dead instrument

        if (vanilla->isCall && fwdPrice > vanilla->strike)
        {
            premium = fwdPrice - vanilla->strike;
        }
        else if ( !vanilla->isCall && fwdPrice<vanilla->strike)
        {
            premium = vanilla->strike - fwdPrice;
        }
        
    }
    else
    {

        // interpolate the vol using our LN request
        if ( isSmileOff)
            vol = vanilla->volData->getATMVol(vanilla->expiryDate);
        else
            vol = vanilla->volData->getSmileVolByStrike(vanilla->expiryDate, vanilla->strike);

        double yFrac = expDate.daysDiff(today) / 365.0;
        double df = vanilla->yieldCurve->pv(today, expDate);    

        long nTimeSteps=50, nNodes=(nTimeSteps+2)*(nTimeSteps+3)/2;
        long m=0, mDummy;

        double vol_shift = 0.001;
        double price1, price2, price3, price4, price5;

        price1 = priceVanillaTree(
            vanilla->isEuropean,
            vanilla->isCall,
            today,
            vanilla->expiryDate,
            fwdPrice, 
            optionStrike,
            vol,
            nTimeSteps, 
            &m);

        // reprice american with one extra node (oscillation)
        price2 = priceVanillaTree(
            vanilla->isEuropean,
            vanilla->isCall,
            today,
            vanilla->expiryDate,
            fwdPrice, 
            optionStrike,
            vol,
            nTimeSteps+1, 
            &mDummy);

        // price european on the tree
        price3 = priceVanillaTree(
            true,
            vanilla->isCall,
            today,
            vanilla->expiryDate,
            fwdPrice, 
            optionStrike,
            vol,
            nTimeSteps, 
            &mDummy);
    
        // price european with one extra node
        price4 = priceVanillaTree(
            true,
            vanilla->isCall,
            today,
            vanilla->expiryDate,
            fwdPrice, 
            optionStrike,
            vol,
            nTimeSteps+1, 
            &mDummy);

        // exact price of european to use as a control variate
 
        //double discFactor = instSettle->pv(vanilla->getValueDate(),
         //                                    expiryDate, 
          //                                   vanilla->discount.get(), 
          //                                   vanilla->asset.get());//?????

	    //Note, Black takes variance
        price5 = Black::price( vanilla->isCall, fwdPrice, 
                                        vanilla->strike, df, vol*vol*yFrac);
//        premium *= discFactor;
	      // interpolate the vol using our LN request

        double coeff = 1.0-(double)m/(double)nNodes;
        premium = (0.5*(price1+price2)+coeff*(price5-0.5*(price3+price4)));

	}

    // apply buy/sell adjustment
    if (!vanilla->isBuy)
    {
        premium *= -1;
    }
    
    return make_pair(premium *= vanilla->notional, vol);
}

double EnergyVanillaClosedForm::priceVanillaTree(
            bool isEuropean,
            bool isCall,
            const DateTime& today,
            const DateTime& expiry,
            double fwdPrice, 
            double strike,
            double vol,
            long n, /* no_time_steps*/
            long* m) const
{
    double df = vanilla->yieldCurve->pv(today, expiry);        
    double dt = expiry.daysDiff(today) / (365.0*n);
    double fwd_fwd_df = exp(log(df)/n);
    double dummy;

	DoubleArray spotArray(2*n+5);
    DoubleArray priceArray(n+3);

    if (dt>0.0)
    {

        double a = 1.0;
        double b_sq = a*a * (exp(vol*vol*dt) - 1.0);
        double u =( (a*a +b_sq + 1.0)+sqrt( (a*a+b_sq+1.0)*(a*a+b_sq+1.0)-(4.0*a*a) ) )/(2.0*a);
        double p = (a - 1.0/u) / (u - 1.0/u);

        *m = 0;
        
        double intermediatePrice = fwdPrice * exp(-(n+2)*log(u));
        spotArray[0] = intermediatePrice;
        long i;
        for (i=1; i<= 2*n+4; i++)
        {
            intermediatePrice *= u;
            spotArray[i] = intermediatePrice;
        }

        for (i=0; i<=n+2; i++)
        {
            double x = spotArray[2*i]-strike;
			if (!isCall) x *= -1;
            if (x>=0) *m += 1;

            priceArray[i] = Maths::max(x, 0.0);
        }
        
        for (i=n+1; 2<=i; i--)
        {
            for (long j=0; j<=i; j++)
            {                
                double x = (spotArray[n+2+2*j-i]-strike);
				if (!isCall) x *= -1;
                double y = fwd_fwd_df * (p*priceArray[j+1] + (1.0-p)*priceArray[j]);

                if (x>=y) *m += 1;
                
                if (isEuropean)
                    priceArray[j] = Maths::max(0.0,y);
                else
                    priceArray[j] = Maths::max(x,y);
            }

            if (i==4)
                dummy = priceArray[2];
        }
	}
    else
    {
        throw ModelException();      
    }

    return priceArray[1];
}


class EnergyVanillaHelper
{

public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyVanilla, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedFormEnergy::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultVanilla);

        FIELD(energyUnderlyer,        "Underlyer Wrapper");
        FIELD(futuresCurve,        "Futures Curve Wrapper");
        FIELD(contractLabel,        "Contract Label");
		FIELD(volData,        "Energy Implied Vol Surf Wrapper");
		FIELD(yieldCurve,        "Yield Curve Wrapper");
        FIELD(strike,        "Contract Strike");
        FIELD(isCall,           "Is it a call option(true/false)");
        FIELD(isEuropean,     "Is it an European option(true/false)");
    	FIELD(isBuy,           "Is is a Buy(true) and sell(false)");
	    

	    // optional with default rules
		FIELD(notional,         "Option notional");
		FIELD_MAKE_OPTIONAL(notional);

        FIELD(valueDate,        "valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);

	    FIELD(expiryDate,     "Option Expiry Date");
        FIELD_MAKE_OPTIONAL(expiryDate);

	    FIELD(paymentDate,     "Payment Date");
        FIELD_MAKE_OPTIONAL(paymentDate);

	    FIELD(contractDate,     "Future Maturity Date");
        FIELD_MAKE_TRANSIENT(contractDate);

       // FIELD(instSettle, "Instrument settlement at maturity");
      //  FIELD_MAKE_OPTIONAL(instSettle);
        
    }

    static IObject* defaultVanilla()
	{
        return new EnergyVanilla();
    }
};

CClassConstSP const EnergyVanilla::TYPE = CClass::registerClassLoadMethod(
    "EnergyVanilla", typeid(EnergyVanilla), EnergyVanillaHelper::load);

bool  EnergyVanillaLoad() { return (EnergyVanilla::TYPE != 0);   }


DRLIB_END_NAMESPACE

