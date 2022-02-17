
//
//   Group       : Credit QR
//
//   Filename    : SCIDRiskyZero.cpp
//
//   Description : a pilot instrument for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SCIDRiskyZero.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/SCIDparameters.hpp"

DRLIB_BEGIN_NAMESPACE

SCIDRiskyZero::~SCIDRiskyZero(){}

/* Default constructor */
SCIDRiskyZero::SCIDRiskyZero() : CInstrument(TYPE){
};

string SCIDRiskyZero::discountYieldCurveName() const {
    return discount.getName();
}

/** Do some asset specific validation */
void SCIDRiskyZero::Validate() 
{
    static const string routine("SCIDRiskyZero::Validate");
	try 
    {
/*		if (!((MCMethod=="CLOSED FORM")||(MCMethod=="FAST")||(MCMethod=="FULL1")||(MCMethod=="FULL2")||(MCMethod=="HYBRID")) )
		    throw ModelException("Unrecognized MCMethod", routine);
		if (timeSteps<1e-12)
		    throw ModelException("time steps too small", routine);
            */
        
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}

/** Get the asset , par curve and discount market data */
void SCIDRiskyZero::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);

    /*=========================================================================
     * GET THE DISCOUNT CURVE 
     *=======================================================================*/
    discount.getData(model, market);
}

DateTime SCIDRiskyZero::getValueDate() const {
    return valueDate;
}



void SCIDRiskyZero::GetParSpread(MaturityPeriodSP& freq, DoubleMatrix &modelPS, DoubleMatrix &marketPS, SCIDparametersSP &sCIDparamSP)
{
	modelPS.resize(maturities.size(),sCIDparamSP->getNbNames());
	marketPS.resize(maturities.size(),sCIDparamSP->getNbNames());
	DateTimeArray tDate;
	sCIDparamSP->push_backTimeLine(tDate, valueDate, maturities.back(), freq, true);

	DoubleMatrix modelSurvProba(tDate.size(), sCIDparamSP->getNbNames());
	modelSurvProba.fill(0);
	DoubleMatrix marketSurvProba = modelSurvProba;

	sCIDparamSP->NamesSurvProb(tDate, modelSurvProba);
	sCIDparamSP->MarketSurvProb(tDate, marketSurvProba);

	array<double> modelRiskyDiscount(tDate.size()), marketRiskyDiscount(tDate.size());
	for (int m=0; m<sCIDparamSP->getNbNames(); m++)
		{
			ICDSParSpreadsSP CDSm = sCIDparamSP->getCDSparSpread(m);
			for (int j=0; j<tDate.size(); j++) 
			{
				modelRiskyDiscount[j] = modelSurvProba[j][m];
				marketRiskyDiscount[j] = marketSurvProba[j][m];
			}
			EffectiveCurve modelCurve(valueDate, discount.getSP(), tDate, modelRiskyDiscount, EffectiveCurve::FLAT_FORWARD);
			EffectiveCurve marketCurve(valueDate, discount.getSP(), tDate, marketRiskyDiscount, EffectiveCurve::FLAT_FORWARD);
			for (int j=0; j<maturities.size(); j++)
			{
			    CashFlowArray cf = SwapTool::cashflows(valueDate, maturities[j], false, 1.0, 3, "M", CDSm->dayCountConv() );
				cf[cf.size()-1].amount -= 1.0; // ugly
				double modelPL = modelCurve.protectionPV(valueDate,valueDate,maturities[j],IDiscountCurveRisky::RECOVER_1);
				modelPL *= 1-sCIDparamSP->getRecovery(m);
				double marketPL = marketCurve.protectionPV(valueDate,valueDate,maturities[j],IDiscountCurveRisky::RECOVER_1);
				marketPL *= 1-sCIDparamSP->getRecovery(m);

				double modelRA = modelCurve.annuityPV(cf,valueDate,IDiscountCurveRisky::RECOVER_1);
				double marketRA = marketCurve.annuityPV(cf,valueDate,IDiscountCurveRisky::RECOVER_1);
				modelPS[j][m] = modelPL/modelRA;
				marketPS[j][m] = marketPL/marketRA;
			}
		}
}


void SCIDRiskyZero::priceClosedForm(CResults* results, Control* control, SCID* model) {
    static const string method = "SCIDRiskyZero::priceClosedForm";

    try {
        OutputRequest* request = 0;
        if (! SCID::TYPE->isInstance(model))
            throw ModelException(method, "sCID instrument did not receive a SCID model");
        const string& ccy = discount->getCcy();

		SCIDparametersSP sCIDparamSP = model->getSCIDparameters();
		DoubleMatrix survProba(maturities.size(), sCIDparamSP->getNbNames());
		survProba.fill(0);
		string MCMethod = model->getMCAlgorithm();
		int seed = model->getSeed();
		if (   (MCMethod=="FAST")||(MCMethod=="HYBRID")||(MCMethod=="FULL1")||(MCMethod=="FULL2")||(MCMethod=="MARKET")||(MCMethod=="CLOSED FORM") ) 
		{
			if (MCMethod=="MARKET")
				sCIDparamSP->MarketSurvProb(maturities, survProba);
			else
			if (MCMethod=="CLOSED FORM")
				sCIDparamSP->NamesSurvProb(maturities, survProba);
			else
			if (MCMethod=="FAST")
			{
				array<double> t = sCIDparamSP->DateAsDouble(maturities);
				double lastTime = t.back();
			
				survProba.fill(0.0);
				vector< DoubleMatrix > survProballWorlds(sCIDparamSP->getNbWorlds(), survProba);
				double timeSteps = model->getCFtimeSteps();
				sCIDparamSP->setFastMC(seed, timeSteps, int(lastTime/timeSteps), model->getNbPathsFastNoJump(), model->getNbPathsFastAtLeastOneJump());
				sCIDparamSP->ComputeSurvProbinAllWorld_fastMCtest(t,survProballWorlds);
				for (unsigned int i=0; i<survProballWorlds.size(); i++) 
					survProba.add(survProballWorlds[i],sCIDparamSP->getInitialWeights(i));
			}
			else if ((MCMethod=="FULL1")||(MCMethod=="FULL2"))
			{
				// Compute Dates at which we are going to simulate spreads and losses
				DateTimeArray tDate;
				sCIDparamSP->push_backTimeLine(tDate, valueDate, maturities.back(), model->getFreqFullMC(), true);
				DoubleArray maturitiesAsDouble = sCIDparamSP->DateAsDouble(maturities);
				DoubleArray defaultTime;
				vector<int> indexMaturities = DateTime::getCeilingProjection(maturities, tDate);

				sCIDparamSP->setFullMC(seed+3, tDate);
				int nbPathsFull = model->getNbPathsFull();
				for (int mc=0; mc<nbPathsFull; mc++)
				{
					sCIDparamSP->FullMCSimulation();
					if (MCMethod=="FULL1") sCIDparamSP->getSurvProba(0,indexMaturities,survProba,true);
					else
					{
						sCIDparamSP->getDefaultTime(defaultTime);
						for (size_t i=0; i<indexMaturities.size(); i++)
							for (int m=0; m<sCIDparamSP->getNbNames(); m++)
								if (maturitiesAsDouble[i]<defaultTime[m]) survProba[i][m]+=1.0;
					}
				}
				survProba.scale(1.0/nbPathsFull);
		    }
			else if (MCMethod=="HYBRID")
			{
			// Compute Dates at which we are going to simulate the spread and losses
				DateTimeArray tDate;
				sCIDparamSP->push_backTimeLine(tDate, valueDate, forwardDate, model->getFreqFullMC(), true);
				double forwardDateAsDouble = valueDate.yearFrac(forwardDate);
				DoubleArray maturitiesAsDouble = sCIDparamSP->DateAsDouble(maturities);
				DoubleArray defaultTime;
				DoubleMatrix survProbaTemp(survProba);

				sCIDparamSP->setFullMC(seed+3, tDate);
				int nbPathsFull = model->getNbPathsFull();
				for (int mc=0; mc<nbPathsFull; mc++)
				{
					sCIDparamSP->FullMCSimulation();
					sCIDparamSP->getDefaultTime(defaultTime);
					sCIDparamSP->getFutureSurvProba(tDate.size()-1, maturitiesAsDouble, survProbaTemp);
					for (int m=0; m<sCIDparamSP->getNbNames(); m++)
						if (defaultTime[m]>forwardDateAsDouble)
							for (int i=0; i<maturities.size(); i++)
								survProba[i][m] += survProbaTemp[i][m];
				}
				survProba.scale(1/double(nbPathsFull));
			}

			request = control->requestsOutput(OutputRequest::DEFAULT_PROBABILITY);
			survProba.scale(-1.0);
			survProba.scalarAdd(1.0);
    	    if (request) results->storeRequestResult(request, CDoubleMatrixSP (new DoubleMatrix(survProba)));
		}

		results->storePrice(0, ccy);
		request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
		MaturityPeriodSP freq = model->getFreqFastMC();
		DoubleMatrix spreadError, marketSP;
		GetParSpread(freq, spreadError, marketSP, sCIDparamSP);  // spreadError is the model par spread for now
		spreadError.add(marketSP,-1.0);  // and now really the error
		spreadError.scale(10000.0); // now in bps
        if (request) results->storeRequestResult(request, CDoubleMatrixSP (new DoubleMatrix(spreadError)));


    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class SCIDRiskyZeroHelper {
public:
    static IObject* defaultSCIDRiskyZero();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SCIDRiskyZero, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(SCID::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultSCIDRiskyZero);

        FIELD(valueDate,                  "Valuation Date - default from market.");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(discount,                   "Discount curve");        
		FIELD(maturities,				 "Maturities");
		FIELD(forwardDate,				 "Hybrid Date");
		FIELD(		 dcc,						 "Day Count Convention");
	}
};

IObject* SCIDRiskyZeroHelper::defaultSCIDRiskyZero() {
    return new SCIDRiskyZero();
}

CClassConstSP const SCIDRiskyZero::TYPE = 
    CClass::registerClassLoadMethod("SCIDRiskyZero", typeid(SCIDRiskyZero),SCIDRiskyZeroHelper::load);

/*=============================================================================
 * Pricing Models
 *===========================================================================*/
class SCIDRiskyZeroSCID : public SCID::IProduct {
private:
    const SCIDRiskyZero* instr; // a reference
    SCID* model;

public:
    SCIDRiskyZeroSCID(const SCIDRiskyZero* instr, SCID* model): instr(instr), model(model){}
    void price(SCID* model,
               Control*         control, 
               CResults*        results) const {
        const_cast<SCIDRiskyZero*>(instr)->priceClosedForm(results, control, model);
    }
};
    
/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCID::IProduct* SCIDRiskyZero::createProduct(SCID* model) const {
    return new SCIDRiskyZeroSCID(this, model);
}

/** Included in ProductsLib to force the linker to include this file */
bool SCIDRiskyZeroLoad() {
    return (SCIDRiskyZero::TYPE != 0);
}

DRLIB_END_NAMESPACE

