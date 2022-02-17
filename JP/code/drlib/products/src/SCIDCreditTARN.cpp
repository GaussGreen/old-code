//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : SCIDCreditTARN.cpp
//
//   Description : a pilot instrument for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SCIDCreditTARN.hpp"
//#include "edginc/EffectiveCurve.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/SCIDparameters.hpp"

DRLIB_BEGIN_NAMESPACE

SCIDCreditTARN::~SCIDCreditTARN(){}

/* Default constructor */
SCIDCreditTARN::SCIDCreditTARN() : CInstrument(TYPE){
};

/** Do some asset specific validation */
void SCIDCreditTARN::Validate() {
    static const string routine("SCIDCreditTARN::Validate");

	try 
    {
		if (TriggerDates.size()!=Triggers.size())
		    throw ModelException("Trigger Dates not of the same size as Triggers", routine);
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}

/** Get the asset , par curve and discount market data */
void SCIDCreditTARN::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);

    /*=========================================================================
     * GET THE DISCOUNT CURVE 
     *=======================================================================*/
    discount.getData(model, market);
}

/** 
 * This is where the pricer evaluates its MTM
 * 
 */

void SCIDCreditTARN::priceClosedForm(CResults* results, Control* control, SCID* model) {
    static const string method = "SCIDCreditTARN::priceClosedForm";

    try {
        OutputRequest* request = 0;
        if (! SCID::TYPE->isInstance(model))
            throw ModelException(method, "sCID instrument did not receive an SCID model");
        const string& ccy = discount->getCcy();

		SCIDparametersSP sCIDparamSP = model->getSCIDparameters();
		DateTimeArray tDate(1, valueDate);
		DateTime firstDate = valueDate, lastDate;
		int startCount=0;
		while ((startCount<TriggerDates.size())&&(TriggerDates[startCount]<=valueDate)) startCount++;
		QLIB_VERIFY(startCount<TriggerDates.size()-1,"no TriggerDates after valueDate");

		DateTimeArray goodTrigDate;
		for (int i=startCount; i<TriggerDates.size(); i++)
		{
			lastDate = TriggerDates[i];
			goodTrigDate.push_back(lastDate);
			sCIDparamSP->push_backTimeLine(tDate,firstDate, lastDate, model->getFreqFullMC(), false);
			firstDate = lastDate;
		}
		vector<int> indices = DateTime::getIndexes(tDate,goodTrigDate); 
		DoubleArray  redeem(TriggerDates.size()+1, 0.0);
		IntArray nbDefNames(tDate.size());

		int nbPaths = model->getNbPathsFull();

		sCIDparamSP->setFullMC(model->getSeed(), tDate);

		for (int mc=0; mc<nbPaths; mc++)
		{
			sCIDparamSP->FullMCSimulation();
			sCIDparamSP->getNbDefaultedNames(nbDefNames);
			bool triggered = false;
			for (int i=0; (i<goodTrigDate.size())&&(!triggered); i++)
			{
				if (nbDefaultedNames + nbDefNames[indices[i]] <= Triggers[i+startCount])
				{
					redeem[i]++;
					triggered=true;
				}
			}
			if (!triggered) redeem[redeem.size()-1]++;
		}
		for (int i=0; i<redeem.size(); i++) redeem[i] /= nbPaths;


		results->storePrice(0, ccy);
		request = control->requestsOutput(OutputRequest::REDEEMING_PROBABILITY);
        if (request) results->storeRequestResult(request, DoubleArraySP (new DoubleArray(redeem)));

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class SCIDCreditTARNHelper {
public:
    static IObject* defaultSCIDCreditTARN();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SCIDCreditTARN, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(SCID::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultSCIDCreditTARN);

	FIELD(valueDate,					 "Date of the Trade");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(discount,                   "Discount curve");        
	FIELD(nbDefaultedNames,			 "Number of Defaulted Names by Trade Date");
	FIELD(Triggers,					 "Triggers Level");
	FIELD(TriggerDates,				 "Trigger Dates");

	}
};

IObject* SCIDCreditTARNHelper::defaultSCIDCreditTARN() {
    return new SCIDCreditTARN();
}

CClassConstSP const SCIDCreditTARN::TYPE = 
    CClass::registerClassLoadMethod("SCIDCreditTARN", typeid(SCIDCreditTARN),SCIDCreditTARNHelper::load);

/*=============================================================================
 * Pricing Models
 *===========================================================================*/
class SCIDCreditTARNSCID : public SCID::IProduct {
private:
    const SCIDCreditTARN* instr; // a reference
    SCID* model;

public:
    SCIDCreditTARNSCID(const SCIDCreditTARN* instr, SCID* model): instr(instr), model(model){}
    void price(SCID* model,
               Control*         control, 
               CResults*        results) const {
        const_cast<SCIDCreditTARN*>(instr)->priceClosedForm(results, control, model);
    }
};
    
/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCID::IProduct* SCIDCreditTARN::createProduct(SCID* model) const {
    return new SCIDCreditTARNSCID(this, model);
}

/** Included in ProductsLib to force the linker to include this file */
bool SCIDCreditTARNLoad() {
    return (SCIDCreditTARN::TYPE != 0);
}

DRLIB_END_NAMESPACE

