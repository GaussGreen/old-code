#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/FD1DMultiSolver.hpp"
#include "edginc/FD1DRegimes.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/ATMVolRequest.hpp"


DRLIB_BEGIN_NAMESPACE

////////////////////////// FD1DRegimes /////////////////////////
void FD1DRegimes::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD1DRegimes, clazz);
    SUPERCLASS(FD1DMulti);
    EMPTY_SHELL_METHOD(defaultFD1DRegimes);

	FIELD(stochFactor, "stochastic factors"); 
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(stochFactor);
}

//----------------------------------------------------------------
CClassConstSP const FD1DRegimes::TYPE = CClass::registerClassLoadMethod(
										"FD1DRegimes", typeid(FD1DRegimes), load);

//----------------------------------------------------------------
FD1DRegimes::FD1DRegimes( const CClassConstSP & type ) :
    FD1DMulti( type ),
    volType(VolSurface::TYPE->getName())
{
    changeOfVar = LOG_X;  //1: S; 2: log(S); 5: log(fwd)
}

//----------------------------------------------------------------
FD1DRegimes::~FD1DRegimes(){}

//----------------------------------------------------------------
void FD1DRegimes::validatePop2Object(){
    static const string method = "FD1DRegimes::validatePop2Object";
    try
    {
        FD1DMulti::validatePop2Object();

        if (isFwdInduction)
        {
            throw ModelException(method, "Fwd Induction hasn't implemented yet! ");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void FD1DRegimes::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
	FDModel::getMarket(market, instruments);

	string name = factors[0]->getName();
	stochFactor = MultiRegimeFactorSP::dynamicCast(market->GetData(name,MultiRegimeFactor::TYPE));

	stochFactor->validatePop2Object();
	nPdes = stochFactor->nbStates;
	iState = stochFactor->initState;
}

MarketDataFetcherSP FD1DRegimes::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(volType));
}

IModel::WantsRiskMapping FD1DRegimes::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** retrieving market data from MDF */    
void FD1DRegimes::retrieveFactor(){
    static const string method = "FD1DRegimes::retrieveFactor";
    try{
        FD1DMulti::retrieveFactor();

        if ((StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ||
            (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ){
                throw ModelException("FD1DRegimes::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }

		// get time metric
		ATMVolRequestSP		volRequest(new ATMVolRequest());
		CVolProcessedSP		vol(underlying->getProcessedVol(volRequest.get()));
        timeMetric = vol->GetTimeMetric();

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void FD1DRegimes::initModel()
{
	// temporary...
	nPdes = stochFactor->nbStates;

    //call the base;
    FD1DMulti::initModel();
}

void FD1DRegimes::finaliseModel(CControl*    control)
{
    //timeLine should be already be set up
    int iAsset;
    for (iAsset=0; iAsset < factors.size(); iAsset++) {
        dynamic_cast<const CAsset *>(
            factors[iAsset].get())->fwdValue(timeLine->StepDates, stepForward);
    }
       
    // call parent last
    FD1DMulti::finaliseModel( control );
}

//----------------------------------------------------------------
/** calculate the stuff which specifies the grid */
/**calculate FD grid and its underlying*/
/**alpha : truncation1D*/
void FD1DRegimes::setFdBounds(
    double alpha1,
    double & outLowB1, 
	double & outUpB1 ) const
{
    //get the vol from the underlying to set the FDBound         
    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];
	double volForBos = stochFactor->vol->expect(longest_T);
    double temp1 = alpha1* volForBos * sqrt(longest_T);
    double undSpot1 = underlying->fwdValue(timeLine->StepDates[0]);

    switch (changeOfVar) {
        case X:{
            //S
            outLowB1 = undSpot1 * exp(-temp1);
            outUpB1 = undSpot1* exp(temp1);
        }        
        break;
        case LOG_X: {
            //log(S)
            outLowB1 = log(undSpot1) - temp1;
            outUpB1 = log(undSpot1) + temp1;
        }
        break;
        case LOG_FWDX: {
            //log(fwd)
            outLowB1 = log(undSpot1) - temp1;
            outUpB1 = log(undSpot1) + temp1;
        }
        break;
    }
}

void FD1DRegimes::pdeCoeff(
    int step,
    DoubleArray & a,
    DoubleArray & c,
    DoubleArray & f,
    DoubleArray & g,
	DoubleArray & q,
	DoubleArray & jump) const
{
    if (step == timeLine->NumOfStep){
        throw ModelException("FD1DRegimes::pdeCoeff", 
            "pdeCoeff shouldn't be called at last step!.");
    }
            
	double dt_trading =timeLine->TradeYrFrac[step+1];
    double mu = log(stepForward[step+1]/stepForward[step])/dt_trading;

    double df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);
    double r = -log(df)/dt_trading;

	//fill all coeff of PDE
	int iPde, jPde, ix;
	double tmpSum;
	for (ix = 0; ix < dim; ix++){
		for (iPde = 0; iPde < nPdes; ++iPde){

			c[iPde * dim + ix] = 0.5 * stochFactor->volStates[iPde] * stochFactor->volStates[iPde];
			g[iPde * dim + ix] = 0.;

			tmpSum = 0.0;
			for (jPde = 0; jPde < nPdes; ++jPde){
				tmpSum += stochFactor->generator[jPde][iPde] * stochFactor->jumpsAtDrift[jPde][iPde];
				q[(iPde * nPdes + jPde) * dim + ix] = stochFactor->generator[jPde][iPde];
				jump[(iPde * nPdes + jPde) * dim + ix] = stochFactor->jumpsOnLogSpot[jPde][iPde];
			}

			if (stochFactor->rate.get()){
				a[iPde * dim + ix] = mu - r + stochFactor->rateStates[iPde] - c[iPde * dim + ix] - tmpSum; 
				f[iPde * dim + ix] = -stochFactor->rateStates[iPde] + stochFactor->generator[iPde][iPde];
			}else{
				a[iPde * dim + ix] = mu - c[iPde * dim + ix] - tmpSum; 
				f[iPde * dim + ix] = -r + stochFactor->generator[iPde][iPde];
			}
		}
	} 
}

/** calculate vol for current step for a set of spot levels
    returns number of vol calculated - one (flat for all node) or num */
int FD1DRegimes::GetStepVol(int step, vector<double>& vol, const double*, int, int){
	       
	if (vol.size() != 1)
    vol.resize(1);
    double tradeTime = timeLine->TradeTime[step];
	vol[0] = stochFactor->vol->expect(tradeTime);
	return 1;
}

//----------------------------------------------------------------
DRLIB_END_NAMESPACE
