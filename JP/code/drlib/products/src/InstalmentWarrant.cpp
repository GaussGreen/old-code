//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InstalmentWarrant.cpp
//
//   Description : Strike adjusted for dividends and reset rates
//
//   Date        : 7/28/2006
//
//	Author: Keith Law
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE
 
class InstalmentWarrant : public Generic1Factor,
							public CClosedFormLN::IIntoProduct
{
protected:
	const double			startStrike;
	const bool				isCall;
	const DateTimeArray		resetDates;
	const DoubleArray		resetRates;
	const double			initialRate;
	const double			spotAtMaturity;
	const bool				adjustStrike;
	const bool				getPVdivs;

public:

	static CClassConstSP const TYPE;

	InstalmentWarrant(): Generic1Factor(TYPE), startStrike(0.0), isCall(false), 
        initialRate(0.0), spotAtMaturity(0.0), adjustStrike(true), getPVdivs(false) {}

     // override base implementation if required
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    virtual void Validate();

    virtual void validatePop2Object();

	/** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
 
    virtual DateTime getValueDate() const;

	//copied from Generic1Factor
    void getSensStrikes(const OutputNameConstSP&         outputName,
                        const CVolRequest*               volRequest,
                        const SensitiveStrikeDescriptor& sensStrikeDesc,
                        const DoubleArraySP&             sensitiveStrikes);


    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                                  const IModel*      model){
    if (MonteCarlo::TYPE->isInstance(model)){
        const MonteCarlo& mc = dynamic_cast<const MonteCarlo&>(*model);
        return mc.getSensitiveStrikes(this, outputName);
    }
    throw ModelException("Generic1Factor::getSensitiveStrikes",
                         "Sensitive Strikes not supported");
	}
 
	bool sensShift(Theta* theta){
    // get the new date
    const DateTime& newDate = theta->rollDate(valueDate);
            
    /* If fwd start date falls between value date and theta date 
       have to set initial spot.*/
    if ( fwdStarting                           && 
         startDate.isGreaterOrEqual(valueDate) && 
         newDate.isGreaterOrEqual(startDate)   ) 
    {
        // returns the fwd price if shift is a theta fs
        initialSpot = asset->getThetaSpotOnDate(theta, startDate);
        
        // not fwd starting anymore 
        fwdStarting = false;
    }  
    
    // roll today 
    valueDate = newDate;
    
    return true; // continue to tweak components which implement Theta
	}
  
	void addOutputRequests(Control* control,
                           Results* results,
                           const double& indVol, 
						   const double& strike, 
						   const double& interest, 
						   const double& premium, 
						   const double& capital, 
						   const double& PVdivs,
                           const double& futStrike) const{

    // take care of additional outputs
    if ( control && control->isPricing() )
    {
        // IND_VOL
        OutputRequest* request = NULL;
        if (control->requestsOutput(OutputRequest::IW_IND_VOL, request)) 
        {
                     results->storeRequestResult(request, indVol); 
        }
       
		// strike
		if (control->requestsOutput(OutputRequest::IW_STRIKE, request)) 
        {
                     results->storeRequestResult(request, strike); 

        }

		// interest
		if (control->requestsOutput(OutputRequest::IW_INTEREST, request)) 
        {
                     results->storeRequestResult(request, interest);  

        }

		// put premium
		if (control->requestsOutput(OutputRequest::IW_PUT, request)) 
        {
                     results->storeRequestResult(request, premium); 

        }

		// capital
		if (control->requestsOutput(OutputRequest::IW_CAPITAL, request)) 
        {
                     results->storeRequestResult(request, capital);  

        }

		// PV dividends
		if (control->requestsOutput(OutputRequest::IW_PVDIVS, request)) 
        {
                     results->storeRequestResult(request, PVdivs); 

        }

        // strike at expiry, futStrike
		if (control->requestsOutput(OutputRequest::IW_FUTSTRIKE, request)) 
        {
                     results->storeRequestResult(request, futStrike); 

        }
    };
	
	} 

bool priceDeadInstrument(CControl* control, CResults* results, double adjStrike) const
{
    
    static string method = "InstrumentWarrant::priceDeadInstrument";

    double    value         = 0.0;
	double	put				= 0.0;
	double capital          = 0.0;

    DateTime matDate = resetDates.back();
	bool expired = valueDate >= matDate;
    
	if (!expired)
        return false; // not dead yet

    DateTime settlementDate = instSettle->settles(matDate, asset.get());

    if (valueDate >= settlementDate)
    {// settled already
        results->storePrice(0.0, discount->getCcy());
        addOutputRequests(control, results, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 
        return true;
    }


        put = GetIntrinsic(spotAtMaturity,
                                 adjStrike,
                                 isCall, 
                                 true /* isOption */);

		capital = spotAtMaturity - adjStrike;
		value = put + capital;

        // pv from settlement to today
        value *= discount->pv(valueDate, settlementDate);
      
    // store results
    results->storePrice(value, discount->getCcy());
    addOutputRequests(control, results, 0.0, adjStrike, 0.0, put, capital, 0.0, adjStrike); 

    return true;
}

private:
	friend class InstalmentWarrantHelper;
	friend class InstalmentWarrantClosedForm;
	
};


typedef smartPtr<InstalmentWarrant> instalmentWarrantSP;

class InstalmentWarrantHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(InstalmentWarrant, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultInstalmentWarrant);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        FIELD(startStrike,        "strike before any dividend or interest rate reset adjusted");
		FIELD(isCall, "true for Call, false for Put");
		FIELD(resetDates,        "Reset Dates");
		FIELD(resetRates, "Reset rates");
		FIELD(initialRate, "Initial Reset Rate");
		FIELD(spotAtMaturity, "Spot Price at Maturity, mandatory for expired instrument");
		FIELD(adjustStrike, "Y: Self-funding Instalment Warrant; N: Instalment Warrant");
		FIELD(getPVdivs, "Y: Instalment Warrant price will subtract Present Value of estimated dividends");

    }

    static IObject* defaultInstalmentWarrant(){
        return new InstalmentWarrant();
    }
};

CClassConstSP const InstalmentWarrant::TYPE = CClass::registerClassLoadMethod(
    "InstalmentWarrant", typeid(InstalmentWarrant), InstalmentWarrantHelper::load);

void InstalmentWarrant::GetMarket(const IModel*         model, 
                            const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    if (asset.usingCache())
    {// should always come in - just to cope with old regression convertion
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);
    }
    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
}


DateTime InstalmentWarrant::getValueDate() const
{
	return valueDate;
}


void InstalmentWarrant::validatePop2Object(){
    static const string method = "InstalmentWarrant::validatePop2Object";

    // validation against forward start
    if ( fwdStarting ) {
		throw ModelException(method, "Forward Starting feature not supported: Start date should not be greater than valueDate");
    }

    // validate reset dates are not empty
	if( resetDates.empty() )
        throw ModelException(method, "There must be at least One reset date");

   	// validate first reset date is greater than value date
   if( resetDates[0] < startDate )
        throw ModelException(method, "The first reset date must be greater than or equal to start date");

    // validate reset dates are in ascending order
	for(int i = 0; i< (resetDates.size()-1); i++)
	{
        if( resetDates[i] >= resetDates[i+1] )
        throw ModelException(method, "Reset dates ("+resetDates[i+1].toString()+") need to be in ascending order");
	}

 
   //validate that SpotAtStart is not zero or negative

    if ( !Maths::isPositive(initialSpot) )
    {
        throw ModelException(method,
                                     "Initial Spot value must be positive ");
    }
}
    
void InstalmentWarrant::Validate()
{
    static const string method = "InstalmentWarrant::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }
    
	AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);  
       
    // validate against struck until thoroughly tested and agreed methodology
    if (ccyTreatment == Asset::CCY_TREATMENT_STRUCK) {
        throw ModelException(method, "ccy struck type not supported.");
    }
  
    // validate spot at maturity is entered for dead instrument
   if((valueDate >= resetDates.back()) && Maths::isZero(spotAtMaturity))
   	   throw ModelException(method, "Trade matured, please enter spot at maturity");
 
}

class InstalmentWarrantClosedForm: public CClosedFormLN::IProduct{
private:
    const InstalmentWarrant*  instalmentWarrant; // a reference


public:
   InstalmentWarrantClosedForm(const InstalmentWarrant* instalmentWarrant): instalmentWarrant(instalmentWarrant){}


   // calculate strike adjusted with dividends and reset rates. 
   
   double calStrike(const DateTime valueDate) const
	{
	double strike = instalmentWarrant->startStrike;
    DateTimeArray resetDates = instalmentWarrant->resetDates;
    DoubleArray resetRates = instalmentWarrant->resetRates;
    DateTime startDate = instalmentWarrant->startDate;

    DateTime matDate = resetDates.back();
    DividendListSP divs = AssetUtil::getAllDivsBetweenDates(instalmentWarrant->asset.get(),
                                                        startDate,
                                                        matDate.rollDate(-1));
    DividendArray eqDiv = divs.get()->getArray();

    double initialRate = instalmentWarrant->initialRate;
    double firstRate = 0.0;
	double fracResetRate = 0.0;
	double intReduct = 0.0;
     
	int i = 0;
	int j = 0;

	// loop for all reset dates <= valueDate and exDiv dates <= valueDate
    // Don't need i<resetDates.size since the last reset date (matDate) will not go through the loop
    // Need j<eqDiv.size since if valueDate is matDate, the || condition will satisfy as eqDiv[j] may have
    // gone through the eqDiv range and becomes bad date which is < matDate, so need eqDiv.size() as control

	while ((j<eqDiv.size()&&(eqDiv[j].getExDate() <= valueDate))||(((resetDates[i] <= valueDate)&&(resetDates[i] < matDate))))
	{
        if((!!eqDiv.size())&& (j<eqDiv.size())){
            if(resetDates[i]<eqDiv[j].getExDate())
		    {
			    fracResetRate = resetRates[i]* (double)(resetDates[i+1].daysDiff(resetDates[i]))/365;

			    strike += fracResetRate * strike;
			    i++;
		    }

		    else if(resetDates[i]>eqDiv[j].getExDate())
		    {
                firstRate = (i==0?initialRate:resetRates[i-1]);

			    intReduct = firstRate * eqDiv[j].getDivAmount()* (double)(resetDates[i].daysDiff(eqDiv[j].getExDate()))/365;
			
			    strike = strike - eqDiv[j].getDivAmount() - intReduct;
			    j++;
		    }

		    //if(resetDates[i] == eqDiv[j].getExDate()), 'else' will not be reached since div DateTime is XCB, <EOD
            // put here for logic completeness.
            // the above 'if' and 'else if' can already handle the case resetDates[i]=eqDiv[j] since
            // new strike = old strike + old strike * resetRates[i] * yearFrac(resetDate[i],[i+1]) - exDiv - IntReduct
            // Since IntReduct = exDiv * resetRates[i] * yearFrac(resetDate[i],[i+1])
            // essentially, new strike = old strike - exDiv + (old strike - exDiv)*resetRates[i]*yearFrac(resetDates[i],[i+1])

		    else
		    {
                fracResetRate = resetRates[i]* (double)(resetDates[i+1].daysDiff(resetDates[i]))/365;

                intReduct = resetRates[i] * eqDiv[j].getDivAmount()* (double)(resetDates[i+1].daysDiff(eqDiv[j].getExDate()))/365; 

			    strike = strike + fracResetRate * strike - eqDiv[j].getDivAmount() - intReduct;

			    i++;
			    j++;
		    }
        }
        else { // take care of eqDiv.size = 0 or when j > eqDivSize
        
                fracResetRate = resetRates[i]* (double)(resetDates[i+1].daysDiff(resetDates[i]))/365;

			    strike += fracResetRate * strike;
			    i++;
        }
	}

	return strike;

	}

     // calculate expected future strike adjusted with dividends and reset rates. 
   
   double calFutureStrike(const DateTime valueDate) const
	{
	double strike = calStrike(valueDate);
  
    DateTimeArray resetDates = instalmentWarrant->resetDates;
    DoubleArray resetRates = instalmentWarrant->resetRates;
    DateTime startDate = instalmentWarrant->startDate;
    DoubleArray futureResetRates = resetRates;
    DateTime matDate = resetDates.back();

    DividendListSP divs = AssetUtil::getAllDivsBetweenDates(instalmentWarrant->asset.get(),
                                                        valueDate,
                                                        matDate.rollDate(-1));
    DividendArray eqDiv = divs.get()->getArray();
    
    double initialRate = instalmentWarrant->initialRate;
    double firstRate = 0.0;
	double fracResetRate = 0.0;
	double intReduct = 0.0;
	int k = 0;

    if (matDate>valueDate)
    {
        while (resetDates[k] <= valueDate) {k++;}
        int i = k;

        while (k < futureResetRates.size()-1)
        {
            futureResetRates[k] = instalmentWarrant->instSettle->pv(valueDate, resetDates[k],instalmentWarrant->discount.get(),instalmentWarrant->asset.get()) /
            instalmentWarrant->instSettle->pv(valueDate, resetDates[k+1],instalmentWarrant->discount.get(),instalmentWarrant->asset.get()) - 1 ;
            futureResetRates[k] = futureResetRates[k] * 365 / (double)(resetDates[k+1].daysDiff(resetDates[k])) + resetRates[k];
            k++;
        }
     
	    int j = 0;

	    // loop for all reset dates < matDate and exDiv dates < matDate
    
        while (j<eqDiv.size()||(resetDates[i] < matDate)|| (resetDates.size()==1 && j<eqDiv.size()))
	    {
            if((!!eqDiv.size())&& (j<eqDiv.size())){
                if(resetDates[i]<eqDiv[j].getExDate())
		        {
			        fracResetRate = futureResetRates[i]* (double)(resetDates[i+1].daysDiff(resetDates[i]))/365;

			        strike += fracResetRate * strike;
			        i++;
		        }

		        else if(resetDates[i]>eqDiv[j].getExDate())
		        {
                    firstRate = (i==0?initialRate:futureResetRates[i-1]);

			        intReduct = firstRate * eqDiv[j].getDivAmount()* (double)(resetDates[i].daysDiff(eqDiv[j].getExDate()))/365;
			
			        strike = strike - eqDiv[j].getDivAmount() - intReduct;
			        j++;
		        }

		     //if(resetDates[i] == eqDiv[j].getExDate()), 'else' will not be reached since div DateTime is XCB, <EOD
                // put here for logic completeness.
    
		        else
		        {
                    fracResetRate = futureResetRates[i]* (double)(resetDates[i+1].daysDiff(resetDates[i]))/365;

                    intReduct = futureResetRates[i] * eqDiv[j].getDivAmount()* (double)(resetDates[i+1].daysDiff(eqDiv[j].getExDate()))/365; 

			        strike = strike + fracResetRate * strike - eqDiv[j].getDivAmount() - intReduct;

			        i++;
			        j++;
		        }
            }
            else { // take care of eqDiv.size = 0 or when j>eqDiv.size()
            
                fracResetRate = futureResetRates[i]* (double)(resetDates[i+1].daysDiff(resetDates[i]))/365;

			    strike += fracResetRate * strike;
			    i++;
            }     
        
	    }
    }
	return strike;

	}
    // calculate interest
	double calInterest(const DateTime valueDate) const{

        bool adjustStrike = instalmentWarrant->adjustStrike;
        double startStrike = instalmentWarrant->startStrike;
        DateTimeArray resetDates = instalmentWarrant->resetDates;
        DoubleArray resetRates = instalmentWarrant->resetRates;
        DateTime startDate = instalmentWarrant->startDate;

        DateTime matDate = resetDates.back();     
        double initialRate = instalmentWarrant->initialRate;

		int n = 0;
		double interest = 0.0;
		DateTime untilDate;
		double strikeZero = startStrike;

		//find Next reset date
		while(resetDates[n] <= valueDate) n++;
 
		if(n==0) interest = initialRate * startStrike* (double)(resetDates[n].daysDiff(valueDate))/365;
        else if(n < resetDates.size())
		{
			// if resetDates[n-1] < valueDate, strikeZero is the strike of the last reset date, 
            // so can use calStrike above but replace function argument by last reset date (untilDate below)
			// if resetDates[n-1]==valueDate, strikeZero is the strike just before the value date
            // so can use calStrike above but replace function argument by 1 day before valueDate (untilDate below)

			untilDate = resetDates[n-1] < valueDate? resetDates[n-1].rollDate(-1):valueDate.rollDate(-1);
			
			strikeZero = adjustStrike?calStrike(untilDate):startStrike;				

			interest = resetRates[n-1]* strikeZero *(double)(resetDates[n].daysDiff(valueDate))/365;
		}	
        else interest = 0.0;

		return interest;
	}
		
    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results) const;

};

void InstalmentWarrantClosedForm::price(CClosedFormLN* model,
                               Control*       control, 
                               CResults*      results) const
{
    static const string method = "InstalmentWarrantClosedForm::price";
    try {
        double         discFactor;      
		double		   divDiscFactor;
        double         variance;        // variance between start date and mat.
        double         fwdAtStart;      // forward price at start date
        double         premium;         // put price 
        double         fwdPrice;        // forward at maturity
        double         ivol;            // indicative vol for put
		bool		    adjustStrike = instalmentWarrant->adjustStrike;	
		bool			getPVdivs = instalmentWarrant->getPVdivs;
		DateTimeArray	resetDates = instalmentWarrant->resetDates;
        DateTime        matDate; 
        double          adjStrike; // for capital calculation purpose and report purpose
        double          futStrike; // for put calculation
		DateTime		valueDate = instalmentWarrant->getValueDate();
		DividendArray	subDiv; // for PV divs calculation
		double			PVdivs = 0.0;
		double			interest;
		double			capital; // spot - adjusted strike
		double			SFI; //Total price: spot - strike + interest + put premium
   

		// find next reset date for Non-adjusted strike Instalment Warrant maturity purpose
		int n = 0;

        matDate = resetDates.back();

        if(!adjustStrike)
        {
            while(resetDates[n] <= valueDate) n++;
            matDate = resetDates[n];
        }

    	//matDate = adjustStrike?resetDates.back():resetDates[n];

		DateTime imntStartDate = instalmentWarrant->fwdStarting? 
            instalmentWarrant->startDate: instalmentWarrant->valueDate;
   
   
		// get estimated dividends from value date to next reset date for getPVdivs purpose
		
		 DividendListSP divs = AssetUtil::getAllDivsBetweenDates(instalmentWarrant->asset.get(),
                                                        valueDate,
                                                        matDate);

		subDiv = divs.get()->getArray();
		

		if(getPVdivs){
		for(int i=0; i<subDiv.size(); i++)
		{
		
			divDiscFactor = instalmentWarrant->instSettle->pv(valueDate,
                                             subDiv[i].getExDate(), 
                                             instalmentWarrant->discount.get(), 
                                             instalmentWarrant->asset.get());

			PVdivs += subDiv[i].getDivAmount() * divDiscFactor;
		}
		}


		//adjusting strike
		adjStrike = adjustStrike?calStrike(valueDate):instalmentWarrant->startStrike;

        // future strike
        futStrike = adjustStrike?calFutureStrike(valueDate):instalmentWarrant->startStrike;

		// valueDate >= matDate is taken care of here
        if(instalmentWarrant->priceDeadInstrument(control, results, adjStrike)){
            return; // dead instrument priced
        }
 
		// price live instrument
        fwdPrice = instalmentWarrant->asset->fwdValue(matDate);

		if (instalmentWarrant->fwdStarting)
        {
            fwdAtStart = instalmentWarrant->asset->fwdValue(instalmentWarrant->startDate);
        }


        // interest 
		interest = calInterest(valueDate);

		// Put Premium 
        // choose how to interpolate the vol - go for traditional route for now

		LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(
                                                   futStrike, 
                                                   imntStartDate, 
                                                   matDate,
                                                   instalmentWarrant->fwdStarting));
        
        // interrogate the model to see if neg fwd vars are allowed
        volRequest->allowNegativeFwdVar(model->negativeFwdVarAllowed());

        // interpolate the vol using our LN request
        CVolProcessedBSSP volBS(instalmentWarrant->asset->
                                getProcessedVol(volRequest.get()));
        // calculate the variance
        variance = volBS->CalcVar(imntStartDate, matDate);

        // calculate the indicative vol
        try {
            ivol = volBS->CalcVol(imntStartDate, matDate);
        }
        catch (exception& ) {
            ivol = 0.0;
        }

        // call Black model without discounting
        premium = Black::price(instalmentWarrant->isCall, fwdPrice, 
                               futStrike, 1.0, variance);
        // and discount
        discFactor = instalmentWarrant->instSettle->pv(instalmentWarrant->valueDate,
                                             matDate, 
                                             instalmentWarrant->discount.get(), 
                                             instalmentWarrant->asset.get());

        premium *= discFactor;

        double scalingFactor = InstrumentUtil::scalePremium(
                                        instalmentWarrant->oneContract,
                                        instalmentWarrant->fwdStarting,
                                        instalmentWarrant->notional,
                                        fwdAtStart,
                                        instalmentWarrant->initialSpot);

        premium *= scalingFactor;

		interest *= scalingFactor;

		capital = (instalmentWarrant->asset->fwdValue(valueDate) - adjStrike) * scalingFactor;

		
		PVdivs *= scalingFactor;
		
		
		//SFI
        SFI = capital + interest + premium - PVdivs;

        results->storePrice(SFI, instalmentWarrant->discount->getCcy());

        instalmentWarrant->addOutputRequests(control,
                          results,
                          ivol,
						  adjStrike,
						  interest,
						  premium,
						  capital,
						  PVdivs,
                          futStrike);

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* InstalmentWarrant::createProduct(
    CClosedFormLN* model) const{

		return new InstalmentWarrantClosedForm(this);
}
 
	bool InstalmentWarrantLoad()
{
   return (InstalmentWarrant::TYPE != 0);

}

DRLIB_END_NAMESPACE

