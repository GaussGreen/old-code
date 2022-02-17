//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDOIndexOption.cpp
//
//   Description : a pilot instrument for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CDOIndexOption.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Black.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/SCIDparameters.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

CDOIndexOption::~CDOIndexOption(){}

CDOIndexOption::CDOIndexOption(CClassConstSP clazz) : CInstrument(clazz){
}


string CDOIndexOption::discountYieldCurveName() const {
    return discount.getName();
}

/** Do some asset specific validation */
void CDOIndexOption::Validate() {
    static const string routine("CDOIndexOption::Validate");
	try 
    {
		if (coupon<=0) 
		    throw ModelException("coupon should be positive", routine);
		if (trancheMaturities->size()==0)
		    throw ModelException("No maturities given", routine);
        //verify that expiries are in increasing order
        for(int i=1; i < expiries.size() - 1; ++i)
           QLIB_VERIFY(expiries[i-1]<expiries[i],routine + ": expiries are not in increasing order");

        //verify that maturities are bigger than expiries
        for(int i=0; i < expiries.size(); ++i)
            for (int k=0; k<trancheMaturities->size(); k++)
                QLIB_VERIFY((*trancheMaturities)[k]->toDate(expiries[i]) > expiries[i], routine + ":an expiry is bigger than a maturity");

	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}

/** Get the asset , par curve and discount market data */
void CDOIndexOption::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);

    /*=========================================================================
     * GET THE DISCOUNT CURVE 
     *=======================================================================*/
    discount.getData(model, market);
}

DateTime CDOIndexOption::getValueDate() const {
    return valueDate;
}

double CDOIndexOption::payoff(double dl, double rd, double loss, double strike)
{
	if (flagLoss) return dl - strike*rd + loss;
		else return dl - strike*rd;
}

//returns the loss at expiry
double CDOIndexOption::getLossAndLegs(long expiryIndexInSim, DateTimeArray & simDates, DateTimeArray & maturities, DateTime& expiry, DoubleArray& RA, DoubleArray& DL, SCID * model) 
{
	SCIDparametersSP sCIDparamSP = model->getSCIDparameters();
    DoubleArray etlTimes(simDates.size());
    for (int i=0; i<etlTimes.size(); i++) etlTimes[i] = valueDate.yearFrac(simDates[i]);
    
    DoubleArray EL(etlTimes.size()), notionalEL(etlTimes.size()); 
  	double loss = sCIDparamSP->getFuturePortfolioEL(expiryIndexInSim, etlTimes, EL);
  	double notionalLoss = sCIDparamSP->getFuturePortfolioNotionalEL(expiryIndexInSim, etlTimes, notionalEL);
   
    for (int i=0; i < etlTimes.size(); ++i){
        EL[i] += loss;
        notionalEL[i] += notionalLoss;
    }
    
    sCIDparamSP->computeDefaultLeg(valueDate, discount.getSP(), simDates, EL, 1.0, expiry, maturities,DL);
    sCIDparamSP->computeRiskyAnnuityLeg(valueDate, discount.getSP(), simDates, notionalEL, 1.0, expiry, maturities,coupon, *dcc, RA);
    
    return loss;
}


/** 
 * This is where the pricer evaluates its MTM
 * 
 */

void CDOIndexOption::priceClosedForm(CResults* results, Control* control, SCID* model) {
    static const string method = "CDOIndexOption::priceClosedForm";

    try {
        if (! SCID::TYPE->isInstance(model))
            throw ModelException(method, "sCID instrument did not receive an SCID model");
        const string& ccy = discount->getCcy();

		SCIDparametersSP sCIDparamSP = model->getSCIDparameters();
		DateTimeArray tDate;
        IntArray nrFastPaths = model->getNrFastPaths();
		string MCmethod = model->getMCAlgorithm();
		int seed = model->getSeed();
        double loss;
        int nbPathsFull = model->getNbPathsFull();
        DoubleArray RA(trancheMaturities->size(),0.0), DL(trancheMaturities->size(),0.0);
        DoubleMatrix meanRA(trancheMaturities->size(),expiries.size()), meanDL(trancheMaturities->size(), expiries.size());
        DoubleMatrix forwardPrices(trancheMaturities->size(),expiries.size());
        //call and put values
        DoubleMatrixArray put(expiries.size()), call(expiries.size()), impliedVol(expiries.size()); 
        for (int i=0; i < expiries.size(); ++i){
            put[i] = CDoubleMatrixSP(new DoubleMatrix(strikes.size(), trancheMaturities->size()));
            call[i] = CDoubleMatrixSP(new DoubleMatrix(strikes.size(), trancheMaturities->size()));
            impliedVol[i] = CDoubleMatrixSP(new DoubleMatrix(strikes.size(), trancheMaturities->size()));
        }
        
        DoubleMatrixArray pathRA(expiries.size()), pathDL(expiries.size()); 
        for (int i=0; i < expiries.size(); ++i){
            pathRA[i] = CDoubleMatrixSP(new DoubleMatrix(trancheMaturities->size(), nbPathsFull));
            pathDL[i] = CDoubleMatrixSP(new DoubleMatrix(trancheMaturities->size(), nbPathsFull));
        }
        
        DoubleMatrix expiryLoss(expiries.size(), nbPathsFull);
        
        int expirySize = expiries.size();
        vector<long> expiryIndex(expirySize);
        vector<double> totalLoss(expirySize,0.0);
        vector<DateTimeArray> fastDates(expirySize);
        vector<DateTimeArray> maturities(expirySize);
        //computes the simulation dates for full and fast MC simulations
        for (int exp = 0; exp <  expirySize; ++exp) {
            // Compute Dates at which we are going to simulate the spread and losses
            if (exp == 0) 
                sCIDparamSP->push_backTimeLine(tDate, valueDate, expiries[exp], model->getFreqFullMC(), true);
            else
                sCIDparamSP->push_backTimeLine(tDate, tDate.back(), expiries[exp], model->getFreqFullMC(), false);
            expiryIndex[exp] = tDate.size()-1;
            
            maturities[exp].resize(trancheMaturities->size());
            for (int k=0; k<maturities[exp].size(); k++)
                maturities[exp][k] = (*trancheMaturities)[k]->toDate(expiries[exp]);
            // Compute Dates at which we are going to compute the TEL
            sCIDparamSP->push_backTimeLine(fastDates[exp], tDate.back(), maturities[exp].back(), model->getFreqFastMC(), true);
        }
        array<double> discTimes = sCIDparamSP->DateAsDouble(tDate);

        double timeSteps = model->getCFtimeSteps();
        //set parameters for Full MC simulation
        sCIDparamSP->setFullMC(seed+3, discTimes);
        //set specific parameters for fast MC simulation - Tranche options should add information about attachement points
        setFastMCParameters(sCIDparamSP);

        for (int mc=0; mc<nbPathsFull; mc++) {
            sCIDparamSP->FullMCSimulation();
            for (int exp = 0; exp <  expirySize; ++exp) {
                loss = discount->pv(expiries[exp])*getLossAndLegs(expiryIndex[exp], fastDates[exp], maturities[exp], expiries[exp], RA, DL, model);
                expiryLoss[exp][mc] = loss;
                totalLoss[exp] +=loss;
                for (int k=0; k<maturities[exp].size(); k++)
                {
                    meanDL[k][exp] +=DL[k];
                    meanRA[k][exp] +=RA[k];
                    (*pathRA[exp])[k][mc] = RA[k];
                    (*pathDL[exp])[k][mc] = DL[k];
                    for (int l=0; l<strikes.size(); l++)
                    {
                        (*call[exp])[l][k] += Maths::max(0.0,payoff(DL[k],RA[k], loss, strikes[l]));
                        (*put[exp])[l][k] += Maths::max(0.0,-payoff(DL[k],RA[k], loss, strikes[l]));
                    }
                }
            }
        }

        for (int exp = 0; exp <  expirySize; ++exp) {
            //computes implied volatility and updates prices for call and put 
            double variance, callVol, putVol, yearFrac = valueDate.yearFrac(expiries[exp]);
            for (int k=0; k<maturities[exp].size(); k++)
            {
		        if (flagLoss) forwardPrices[k][exp] = (totalLoss[exp] + meanDL[k][exp])/meanRA[k][exp];
			        else forwardPrices[k][exp] = meanDL[k][exp]/meanRA[k][exp];
                //computes the discount from value date to expiry
                for (int l=0; l<strikes.size(); l++)
                {				
                    // call version of impliedVariance that returns status
                    if (Black::impliedVariance(true, // option call or put
                                                forwardPrices[k][exp], // forward price
                                                strikes[l], // strike
                                                1,//discount->pv(expiries[exp]), // pv
                                                0.2 * 0.2 * yearFrac, // initial var guess
                                                (*call[exp])[l][k]/meanRA[k][exp], // option price
                                                1.0e-10 * yearFrac, // var accuracy
                                                variance)) // variance
                        callVol = sqrt(variance/yearFrac);
                    else
                        callVol = -1;
                    if (Black::impliedVariance(false, // option call or put
                                                forwardPrices[k][exp], // forward price
                                                strikes[l], // strike
                                                1,//discount->pv(expiries[exp]), // pv
                                                0.2 * 0.2 * yearFrac, // initial var guess
                                                (*put[exp])[l][k]/meanRA[k][exp], // option price
                                                1.0e-10 * yearFrac, // var accuracy
                                                variance)) // variance
                        putVol = sqrt(variance/yearFrac);
                    else
                        putVol = -1;
				    if(fabs(callVol-putVol)<0.001) 
                        (*impliedVol[exp])[l][k] = 0.5*(callVol+putVol);
                    else
                        (*impliedVol[exp])[l][k] = -1;
                    //scale call, and put to the number of full paths and compute volatility
                    /*(*call[exp])[l][k] *=discount->pv(expiries[exp])/(double)nbPathsFull; 
                    (*put[exp])[l][k] *=discount->pv(expiries[exp])/(double)nbPathsFull; 			
                    */
                    (*call[exp])[l][k] /=(double)nbPathsFull; 
                    (*put[exp])[l][k] /=(double)nbPathsFull; 			
                }
            }
        }

        //returns the result
		results->storePrice(0, ccy);
        //returns put results
        results->storeGreek(
            DoubleMatrixArraySP(new DoubleMatrixArray(put)),    
            "OptionResults", 
            OutputNameSP(new OutputName("put")));
        //returns call results
        results->storeGreek(
            DoubleMatrixArraySP(new DoubleMatrixArray(call)),    
            "OptionResults", 
            OutputNameSP(new OutputName("call")));
        //returns implied vols 
        results->storeGreek(
            DoubleMatrixArraySP(new DoubleMatrixArray(impliedVol)),    
            "OptionResults", 
            OutputNameSP(new OutputName("impliedVol")));
        //returns forward prices 
        results->storeGreek(
            CDoubleMatrixSP(new DoubleMatrix(forwardPrices)),    
            "OptionResults", 
            OutputNameSP(new OutputName("forwardPrices")));
        //returns path RAs 
        results->storeGreek(
            DoubleMatrixArraySP(new DoubleMatrixArray(pathRA)),    
            "OptionResults", 
            OutputNameSP(new OutputName("pathRA")));
        //returns path DLs 
        results->storeGreek(
            DoubleMatrixArraySP(new DoubleMatrixArray(pathDL)),    
            "OptionResults", 
            OutputNameSP(new OutputName("pathDL")));
        //returns path losses 
        results->storeGreek(
            CDoubleMatrixSP(new DoubleMatrix(expiryLoss)),    
            "OptionResults", 
            OutputNameSP(new OutputName("pathLoss")));
    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class CDOIndexOptionHelper {
public:
    static IObject* defaultCDOIndexOption();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CDOIndexOption, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(SCID::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultCDOIndexOption);

        FIELD(valueDate,                  "Valuation Date - default from market.");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(discount,                   "Discount curve");        
		FIELD(trancheMaturities,				 "Maturities");
		FIELD(expiries, 				 "Expiries");
		FIELD(strikes,   				 "Strikes");
		FIELD(dcc,						 "Day Count Convention");
		FIELD(flagLoss,					 "Shows if payoff takes into account the loss at the expiry");
		FIELD(coupon,					 "frequency of coupons (in Months)");
	}
};

IObject* CDOIndexOptionHelper::defaultCDOIndexOption() {
    return new CDOIndexOption();
}

CClassConstSP const CDOIndexOption::TYPE = 
    CClass::registerClassLoadMethod("CDOIndexOption", typeid(CDOIndexOption),CDOIndexOptionHelper::load);

/*=============================================================================
 * Pricing Models
 *===========================================================================*/
class CDOIndexOptionSCID : public SCID::IProduct {
private:
    const CDOIndexOption* instr; // a reference
    SCID* model;

public:
    CDOIndexOptionSCID(const CDOIndexOption* instr, SCID* model): instr(instr), model(model){}
    void price(SCID* model,
               Control*         control, 
               CResults*        results) const {
        const_cast<CDOIndexOption*>(instr)->priceClosedForm(results, control, model);
    }
};
    
/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCID::IProduct* CDOIndexOption::createProduct(SCID* model) const {
    return new CDOIndexOptionSCID(this, model);
}

/** Included in ProductsLib to force the linker to include this file */
bool CDOIndexOptionLoad() {
    return (CDOIndexOption::TYPE != 0);
}

DRLIB_END_NAMESPACE

