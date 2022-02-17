
#include "edginc/config.hpp"
#include "edginc/TranchePricerCalibrator.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/Optimizer.hpp"


DRLIB_BEGIN_NAMESPACE

TranchePricerCalibrator::~TranchePricerCalibrator(){}

/* Default constructor */
TranchePricerCalibrator::TranchePricerCalibrator() : CInstrument(TYPE){
};

string TranchePricerCalibrator::discountYieldCurveName() const {
    return discount.getName();
}

/** Do some asset specific validation */
void TranchePricerCalibrator::Validate() {
    static const string routine("TranchePricerCalibrator::Validate");
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
void TranchePricerCalibrator::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
}

DateTime TranchePricerCalibrator::getValueDate() const {
    return valueDate;
}



void TranchePricerCalibrator::priceClosedForm(CResults* results, Control* control, SCIDtree* model) {
	static const string method = "TranchePricerCalibrator::priceClosedForm";

	try {
		OutputRequest* request = 0;
		if (! SCIDtree::TYPE->isInstance(model))
			throw ModelException(method, "sCID instrument did not receive an SCID model");
		const string& ccy = discount->getCcy();

		SCIDtreeParametersSP sCIDtreeParamSP = model->getSCIDtreeParameters();
		DateTimeArray tDate;		
		sCIDtreeParamSP->push_backTimeLine(tDate, valueDate, maturities.back(), model->getFreqFastMC(), false);

		YieldCurveSP discountSP = discount.getSP();
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
								   discountSP,
								   dcc);
		sCIDtreeParamSP->setTELtimesInTree(tDate);
	
// Calibration (or price computations)
		DoubleArray marketPEL(maturities.size());	
		for (int i=0; i<maturities.size(); i++)
			marketPEL[i] = sCIDtreeParamSP->getMarketPortfolioEL(maturities[i]);
		results->storeGreek(DoubleArraySP(new DoubleArray(marketPEL)), "TranchePricerCalibrator", OutputNameSP(new OutputName("marketPEL")));
		vector< SimpleTreeSP > leaves;

		DoubleArrayArray selectedPos;
		sCIDtreeParamSP->guessPositionAtAllMaturities(calibrateParSpread,marketPEL,worldDx, smoothness[0], selectedPos);
		results->storeGreek(DoubleArrayArraySP( new DoubleArrayArray(selectedPos)), "TranchePricerCalibrator", OutputNameSP(new OutputName("selectedPos")));
		sCIDtreeParamSP->a_t.clear();
		sCIDtreeParamSP->a_t.resize(maturities.size());
		for (int i=0; i<maturities.size(); ++i)
		{
			if (i==0) sCIDtreeParamSP->a_t[i].resize(1);
			else
			{
				leaves.clear();
				sCIDtreeParamSP->getLeaves(leaves,i);
				sCIDtreeParamSP->a_t[i].resize(leaves.size());
			}
			for (int j=0; j<sCIDtreeParamSP->a_t[i].size(); ++j)
				sCIDtreeParamSP->a_t[i][j] = selectedPos[i];
			sCIDtreeParamSP->CalibrateTreeToTranches(calibrateParSpread,marketPEL,smoothness[i],i+1,i+1);
		}


		leaves.clear();
		sCIDtreeParamSP->getLeaves(leaves);

		DoubleMatrix valuesAtTreeTime;
		sCIDtreeParamSP->readTree(valuesAtTreeTime);
		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(valuesAtTreeTime)), "TranchePricerCalibrator", OutputNameSP(new OutputName("worlds")));

		DoubleMatrix TELmean;
		DoubleArray PELmean;
		DoubleMatrix DLmean,RAmean;
		sCIDtreeParamSP->readAverageLegs(leaves,DLmean,RAmean);
		sCIDtreeParamSP->readAveragePEL(leaves,PELmean);
		sCIDtreeParamSP->readAverageTEL(leaves,TELmean);
		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(TELmean)), "TranchePricerCalibrator", OutputNameSP(new OutputName("meanTEL")));
		results->storeGreek(DoubleArraySP(new DoubleArray(PELmean)), "TranchePricerCalibrator", OutputNameSP(new OutputName("meanPEL")));
		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(DLmean)), "TranchePricerCalibrator", OutputNameSP(new OutputName("meanDL")));
		results->storeGreek(CDoubleMatrixSP(new DoubleMatrix(RAmean)), "TranchePricerCalibrator", OutputNameSP(new OutputName("meanRA")));

	} catch (exception& e) {
		throw ModelException(&e, method);
	}
}



class TranchePricerCalibratorHelper {
public:
    static IObject* defaultTranchePricerCalibrator();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TranchePricerCalibrator, clazz);
        SUPERCLASS(CInstrument);
		IMPLEMENTS(SCIDtree::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultTranchePricerCalibrator);

        FIELD(valueDate,                  "Valuation Date - default from market.");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(discount,                   "Discount curve");        
		FIELD(kmin,						 "Lower attachment points");
		FIELD(kmax,						 "Upper attachment points");
		FIELD(maturities,				 "Maturities");
		FIELD(coupon,					 "frequency of coupons (in Months)");
		FIELD(dcc,						 "Day Count Convention");
		FIELD(calibrateParSpread,        "Market Par Spread");
		FIELD(smoothness,				 "smoothness coefficient added in the Quad Prog");
		FIELD(worldDx,				     "worlds level defined every Dx");

	}
};

IObject* TranchePricerCalibratorHelper::defaultTranchePricerCalibrator() {
    return new TranchePricerCalibrator();
}

CClassConstSP const TranchePricerCalibrator::TYPE = 
    CClass::registerClassLoadMethod("TranchePricerCalibrator", typeid(TranchePricerCalibrator),TranchePricerCalibratorHelper::load);

/*=============================================================================
 * Pricing Models
 *===========================================================================*/

class TranchePricerCalibratorSCIDtree : public SCIDtree::IProduct {
private:
	const TranchePricerCalibrator* instr; // a reference
	SCIDtree* model;
public:
	TranchePricerCalibratorSCIDtree(const TranchePricerCalibrator* instr, SCIDtree* model): instr(instr), model(model){}
	void price(SCIDtree* model, Control* control, CResults* results) const 
	{
		const_cast<TranchePricerCalibrator*>(instr)->priceClosedForm(results, control, model);
	}
};

/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCIDtree::IProduct* TranchePricerCalibrator::createProduct(SCIDtree* model) const 
{
	return new TranchePricerCalibratorSCIDtree(this, model);
}


/** Included in ProductsLib to force the linker to include this file */
bool TranchePricerCalibratorLoad() {
    return (TranchePricerCalibrator::TYPE != 0);
}

DRLIB_END_NAMESPACE

