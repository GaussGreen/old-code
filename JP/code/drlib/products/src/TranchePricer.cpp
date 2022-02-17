//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : TranchePricer.cpp
//
//   Description : a pilot instrument for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/TranchePricer.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/JumpAffine.hpp"


DRLIB_BEGIN_NAMESPACE

TranchePricer::~TranchePricer(){}

/* Default constructor */
TranchePricer::TranchePricer() : CInstrument(TYPE){
};

string TranchePricer::discountYieldCurveName() const {
    return discount.getName();
}

/** Do some asset specific validation */
void TranchePricer::Validate() {
    static const string routine("TranchePricer::Validate");
	try 
    {
		if (kmin.size()!=kmax.size()) 
		    throw ModelException("not the same number of lower and upper attachments points", routine);
		if (coupon<=0) 
		    throw ModelException("coupon should be positive", routine);
		if (maturities.size()==0)
		    throw ModelException("No maturities given", routine);
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}

/** Get the asset , par curve and discount market data */
void TranchePricer::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);

    /*=========================================================================
     * GET THE DISCOUNT CURVE 
     *=======================================================================*/
    discount.getData(model, market);
}

DateTime TranchePricer::getValueDate() const {
    return valueDate;
}

/** 
 * This is where the pricer evaluates its MTM
 * 
 */

void TranchePricer::priceClosedForm(CResults* results, Control* control, SCID* model) {
    static const string method = "TranchePricer::priceClosedForm";

    try {
        OutputRequest* request = 0;
        if (! SCID::TYPE->isInstance(model))
            throw ModelException(method, "sCID instrument did not receive an SCID model");
        const string& ccy = discount->getCcy();

		SCIDparametersSP sCIDparamSP = model->getSCIDparameters();
		DateTimeArray tDate;
		DoubleMatrix ETLmean, ETLcorrelationWithPEL;
		string MCmethod = model->getMCAlgorithm();
		int seed = model->getSeed();
		if ((MCmethod!="FAST")&&(MCmethod!="HYBRID")&&(MCmethod!="FULL")) MCmethod= "FAST";

		if (  (MCmethod=="FAST")|| ( (MCmethod=="HYBRID")&&(forwardDate==valueDate)) )
		{
			// Compute Dates at which we are going to compute the TEL
			DateTime firstDate;
			if (forwardPricer) firstDate = forwardDate;
				       	  else firstDate = valueDate;
			sCIDparamSP->push_backTimeLine(tDate, firstDate, maturities.back(), model->getFreqFastMC(), true);
			array<double> t = sCIDparamSP->DateAsDouble(tDate);
		
			ETLmean.resize(tDate.size(),kmin.size());
			ETLmean.fill(0.0);

			vector< DoubleMatrix > TELallWorlds(sCIDparamSP->getNbWorlds(), ETLmean);

			// Compute Tranche Expected Losses at all times t
			double timeSteps = model->getCFtimeSteps();
			sCIDparamSP->setFastMC(seed, timeSteps, int(t.back()/timeSteps), model->getNbPathsFastNoJump(), model->getNbPathsFastAtLeastOneJump());
			sCIDparamSP->setConvolution(kmin, kmax);
			sCIDparamSP->ComputeTELinAllWorlds(t, TELallWorlds, model->getConvolutionNoJump(), model->getConvolutionAtLeastOneJump());
			for (unsigned int i=0; i<TELallWorlds.size(); i++) 
				ETLmean.add(TELallWorlds[i],sCIDparamSP->getInitialWeights(i));
		}
		else if (MCmethod=="FULL")
		{
			// Compute Dates at which we are going to simulate spreads and losses
			int indexForward;
			DateTime lastDate;
			if (forwardPricer) lastDate = forwardDate;
				else lastDate = maturities.back();
			sCIDparamSP->push_backTimeLine(tDate, valueDate, lastDate, model->getFreqFullMC(), true);
			if (forwardPricer)
			{
				indexForward = tDate.size()-1;
				sCIDparamSP->push_backTimeLine(tDate, tDate.back(), maturities.back(), model->getFreqFullMC(), false);
			}
			else indexForward = 0;

			DoubleArray loss(tDate.size());
			ETLmean.resize(tDate.size(), kmin.size());  		    ETLmean.fill(0.0);
			ETLcorrelationWithPEL.resize(tDate.size(), kmin.size());ETLcorrelationWithPEL.fill(0.0);
			DoubleMatrix ETLsquareMean(tDate.size(),kmin.size());   ETLsquareMean.fill(0.0);
			DoubleMatrix ETLtimePEL(tDate.size(),kmin.size());      ETLtimePEL.fill(0.0);
			DoubleArray PELmcMean(tDate.size(),0.0);
			DoubleArray PELsquareMcMean(tDate.size(),0.0);
			DoubleArray PELcfMean(tDate.size());                   
			sCIDparamSP->PortfolioEL(tDate,&PELcfMean[0]);
			if (forwardPricer)
			{
				for (int i=0; i<indexForward; ++i) PELcfMean[i] = 0.0;
				for (int i=indexForward+1; i<tDate.size(); ++i) PELcfMean[i] -= PELcfMean[indexForward];
			}


			sCIDparamSP->setFullMC(seed+3, tDate);

			int nbPathsFull = model->getNbPathsFull();
			for (int mc=0; mc<nbPathsFull; mc++)
			{
				sCIDparamSP->FullMCSimulation();
				sCIDparamSP->getPortfolioLoss(loss);
				for (int j=indexForward; j<tDate.size(); j++)
				{
					double portfolioLoss = loss[j]-loss[indexForward];
					for (int i=0; i<kmin.size(); i++)
					{
						double trancheLoss = PayoffTranche(portfolioLoss, kmin[i], kmax[i]);
						ETLmean[j][i] += trancheLoss;
						ETLsquareMean[j][i] += trancheLoss*trancheLoss;
						ETLtimePEL[j][i] += trancheLoss*portfolioLoss;
					}
					PELmcMean[j] += portfolioLoss;
					PELsquareMcMean[j] += portfolioLoss*portfolioLoss;
				}
			}
			ETLmean.scale(1.0/nbPathsFull);
			ETLsquareMean.scale(1.0/nbPathsFull);
			ETLtimePEL.scale(1.0/nbPathsFull);
			for (int j=indexForward; j<tDate.size(); j++) 
			{
				PELmcMean[j] /= double(nbPathsFull);
				PELsquareMcMean[j] /= double(nbPathsFull);
			}
			for (int j=indexForward; j<tDate.size(); j++)
			{
				double portfolioVariance = PELsquareMcMean[j] - Maths::square(PELmcMean[j]);
				for (int i=0; i<kmin.size(); i++)
				{
					ETLcorrelationWithPEL[j][i] = ETLtimePEL[j][i] - ETLmean[j][i]*PELmcMean[j]; // PELmcMean should maybe be replaced by PELcfMean[j]
					double trancheVariance = ETLsquareMean[j][i]-Maths::square(ETLmean[j][i]);
					double alpha = 0.0;
					if (trancheVariance>0)
					{
						alpha = ETLcorrelationWithPEL[j][i]/portfolioVariance;
						ETLcorrelationWithPEL[j][i] /= sqrt(trancheVariance*portfolioVariance);
					}
					//if (model->useControlVariates())
					//	ETLmean[j][i] -= alpha * (PELmcMean[j] - PELcfMean[j]);
				}
			}
	    }
		else if (MCmethod=="HYBRID")
		{
			// Compute Dates at which we are going to simulate the spread and losses
			sCIDparamSP->push_backTimeLine(tDate, valueDate, forwardDate, model->getFreqFullMC(), true);
			array<double> discTimes = sCIDparamSP->DateAsDouble(tDate);
			// Compute Dates at which we are going to compute the TEL
			sCIDparamSP->push_backTimeLine(tDate, tDate.back(), maturities.back(), model->getFreqFastMC(), false);
			DoubleArray etlTimes(tDate.size()-discTimes.size());
			for (int i=0; i<etlTimes.size(); i++) etlTimes[i] = valueDate.yearFrac(tDate[i+discTimes.size()]);

			sCIDparamSP->setFullMC(seed+3, discTimes);
			double timeSteps = model->getCFtimeSteps();
			sCIDparamSP->setFastMC(seed, timeSteps, int(etlTimes.back()/timeSteps), model->getNbPathsFastNoJump(), model->getNbPathsFastAtLeastOneJump());
			sCIDparamSP->setConvolution(kmin, kmax);

			ETLmean.resize( tDate.size(), kmin.size());
			ETLmean.fill(0);
			DoubleArray loss(discTimes.size());
			int nbPathsFull = model->getNbPathsFull();
			for (int mc=0; mc<nbPathsFull; mc++)
			{
				sCIDparamSP->FullMCSimulation();
				if (!forwardPricer)
				{
					sCIDparamSP->getPortfolioLoss(loss);
					for (int j=0; j<discTimes.size(); j++)
						for (int i=0; i<kmin.size(); i++)
							ETLmean[j][i] += PayoffTranche(loss[j], kmin[i], kmax[i]);
				}
				sCIDparamSP->getFutureETL(discTimes.size()-1, etlTimes, discTimes.size(), ETLmean, (!forwardPricer), model->getConvolutionNoJump(), model->getConvolutionAtLeastOneJump());
			}
			ETLmean.scale(1/double(nbPathsFull));
		}
		else throw ModelException("Unknown computation method", method);

	// Put these expected losses in Effective Curves
		DoubleMatrix PL(kmin.size(),maturities.size()), RA(kmin.size(),maturities.size()), TELoutput(kmin.size()+2, maturities.size());
		vector<CashFlowArray> cf(maturities.size());
		for (unsigned int i=0; i<cf.size(); i++)
		{
			cf[i] = SwapTool::cashflows(valueDate, maturities[i], false, 1.0, coupon, "M", &(*dcc));
			cf[i][cf[i].size()-1].amount -= 1.0; // ugly
		}

		DateTimeArray tDateWithToday;
		int index=0;
		if (tDate[0]!=valueDate) index=1;

		tDateWithToday.resize(tDate.size()+index);
		tDateWithToday[0] = valueDate;
		for (int i=0; i<tDate.size(); i++) tDateWithToday[i+index]=tDate[i];
		array<double> riskyDiscount(tDateWithToday.size());

		for (int i=0; i<kmin.size(); i++)
		{
			double div = 1.0 / ( kmax[i] - kmin[i] );
			riskyDiscount[0] = 1.0;
			for (int j=0; j<tDate.size(); j++) riskyDiscount[j+index] = Maths::max(1e-12,1 - ETLmean[j][i]*div);

			EffectiveCurve trancheCurve(valueDate, discount.getSP(), tDateWithToday, riskyDiscount, EffectiveCurve::FLAT_FORWARD);
			for (int j=0; j<maturities.size(); j++)
			{
				PL[i][j] = trancheCurve.protectionPV(valueDate,valueDate,maturities[j],IDiscountCurveRisky::RECOVER_1);
				RA[i][j] = trancheCurve.annuityPV(cf[j],valueDate,IDiscountCurveRisky::RECOVER_1);
				TELoutput[i][j] = 1 - trancheCurve.survivalProb(maturities[j]);
			}
		}
	    for (int j=0; j<maturities.size(); j++)
			TELoutput[kmin.size()][j] = sCIDparamSP->getMarketPortfolioEL(maturities[j]);
		sCIDparamSP->PortfolioEL(maturities, &TELoutput[kmin.size()+1][0]) ;

		results->storePrice(0, ccy);
		request = control->requestsOutput(OutputRequest::TRANCHE_CONTINGENT_LEG_PRICE);
        if (request) results->storeRequestResult(request, CDoubleMatrixSP (new DoubleMatrix(PL)));
		request = control->requestsOutput(OutputRequest::TRANCHE_FEE_LEG_PRICE);
        if (request) results->storeRequestResult(request, CDoubleMatrixSP (new DoubleMatrix(RA)));
		request = control->requestsOutput(OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE);
        if (request) results->storeRequestResult(request, CDoubleMatrixSP (new DoubleMatrix(TELoutput)));

		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(ETLcorrelationWithPEL)), 
						    "TranchePricer", 
							OutputNameSP(new OutputName("CorrelationTrancheIndex")));

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void TranchePricer::priceClosedForm(CResults* results, Control* control, SCIDtree* model) {
	static const string method = "TranchePricer::priceClosedForm";

	try {
		OutputRequest* request = 0;
		if (! SCIDtree::TYPE->isInstance(model))
			throw ModelException(method, "sCIDtree instrument did not receive an SCIDtree model");
		const string& ccy = discount->getCcy();

		SCIDtreeParametersSP sCIDtreeParamSP = model->getSCIDtreeParameters();
		DateTimeArray tDate;		
		sCIDtreeParamSP->push_backTimeLine(tDate, valueDate, maturities.back(), model->getFreqFastMC(), false);
		YieldCurveSP discountsp = discount.getSP();
		sCIDtreeParamSP->setFastMC(model->getSeed(),
								   model->getCFtimeSteps(), 
								   valueDate.yearFrac(tDate.back()),
								   model->getNbPathsFastNoJump(), 
								   model->getNbPathsFastAtLeastOneJump(),
								   model->getConvolutionNoJump(),
								   model->getConvolutionAtLeastOneJump(),
								   kmin,
								   kmax,
								   coupon,
								   maturities,
								   discountsp,
								   dcc);
		sCIDtreeParamSP->setTELtimesInTree(tDate);
	
// Calibration (or price computations)
		DoubleArray marketPEL(maturities.size());	
		for (int i=0; i<maturities.size(); i++)
			marketPEL[i] = sCIDtreeParamSP->getMarketPortfolioEL(maturities[i]);
		results->storeGreek(DoubleArraySP(new DoubleArray(marketPEL)), "TranchePricer", OutputNameSP(new OutputName("marketPEL")));
		string MCmethod = model->getMCAlgorithm();
		vector< SimpleTreeSP > leaves;

		sCIDtreeParamSP->InitializeTree();
		sCIDtreeParamSP->computeAllTELandPEL(false);
		sCIDtreeParamSP->computeAllLegsFromTEL();

		leaves.clear();
		sCIDtreeParamSP->getLeaves(leaves);

		DoubleMatrix TELmean;
		DoubleArray PELmean;
		DoubleMatrix DLmean,RAmean;
		sCIDtreeParamSP->readAverageLegs(leaves,DLmean,RAmean);
		sCIDtreeParamSP->readAveragePEL(leaves,PELmean);
		sCIDtreeParamSP->readAverageTEL(leaves,TELmean);

		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(TELmean)), "TranchePricer", OutputNameSP(new OutputName("meanTEL")));
		results->storeGreek(DoubleArraySP(new DoubleArray(PELmean)), "TranchePricer", OutputNameSP(new OutputName("meanPEL")));
		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(DLmean)), "TranchePricer", OutputNameSP(new OutputName("meanDL")));
		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(RAmean)), "TranchePricer", OutputNameSP(new OutputName("meanRA")));


	} catch (exception& e) {
		throw ModelException(&e, method);
	}
}



/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class TranchePricerHelper {
public:
    static IObject* defaultTranchePricer();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TranchePricer, clazz);
        SUPERCLASS(CInstrument);
		IMPLEMENTS(SCIDtree::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultTranchePricer);

        FIELD(valueDate,                  "Valuation Date - default from market.");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(discount,                   "Discount curve");        
		FIELD(kmin,						 "Lower attachment points");
		FIELD(kmax,						 "Upper attachment points");
		FIELD(maturities,				 "Maturities");
		FIELD(forwardDate,				 "Hybrid Date");
		FIELD(coupon,					 "frequency of coupons (in Months)");
		FIELD(dcc,						 "Day Count Convention");
		FIELD(forwardPricer,			 "Are we pricing a forward tranche (TRUE) or a standard one (FALSE)");

	}
};

IObject* TranchePricerHelper::defaultTranchePricer() {
    return new TranchePricer();
}

CClassConstSP const TranchePricer::TYPE = 
    CClass::registerClassLoadMethod("TranchePricer", typeid(TranchePricer),TranchePricerHelper::load);

/*=============================================================================
 * Pricing Models
 *===========================================================================*/
class TranchePricerSCID : public SCID::IProduct {
private:
    const TranchePricer* instr; // a reference
    SCID* model;
public:
    TranchePricerSCID(const TranchePricer* instr, SCID* model): instr(instr), model(model){}
    void price(SCID* model, Control* control, CResults* results) const 
	{
        const_cast<TranchePricer*>(instr)->priceClosedForm(results, control, model);
    }
};
    

class TranchePricerSCIDtree : public SCIDtree::IProduct {
private:
	const TranchePricer* instr; // a reference
	SCIDtree* model;
public:
	TranchePricerSCIDtree(const TranchePricer* instr, SCIDtree* model): instr(instr), model(model){}
	void price(SCIDtree* model, Control* control, CResults* results) const 
	{
		const_cast<TranchePricer*>(instr)->priceClosedForm(results, control, model);
	}
};


/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCID::IProduct* TranchePricer::createProduct(SCID* model) const {
    return new TranchePricerSCID(this, model);
}

/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCIDtree::IProduct* TranchePricer::createProduct(SCIDtree* model) const 
{
	return new TranchePricerSCIDtree(this, model);
}


/** Included in ProductsLib to force the linker to include this file */
bool TranchePricerLoad() {
    return (TranchePricer::TYPE != 0);
}

DRLIB_END_NAMESPACE

