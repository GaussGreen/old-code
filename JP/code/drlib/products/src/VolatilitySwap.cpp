//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolatilitySwap.cpp
//
//   Description : Volatility Swap 
//
//   Author      : Bruno O Melka
//
//   Date        : 04 Apr 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolVarSwap.hpp"

DRLIB_BEGIN_NAMESPACE

//**********************************************************************//
//****************   INSTRUMENT CLASS AND METHODS    *******************//
//**********************************************************************//

class VolatilitySwap:	public VarianceSwap,
						virtual public VanVSModel::IIntoProduct {

public:
	friend class VolSwapProd;
    static CClassConstSP const TYPE;
	static void load(CClassSP& clazz);
	static IObject* defaultVolatilitySwap();

	VolatilitySwap() :VarianceSwap(TYPE) {
		payoffType = "FORWARD";
		cap = 1.0;
		noDivAdj = true;
		useMatytsinOnly = false;
	}

	/** retrieve market data */
    virtual void GetMarket(const IModel*       model, 
							const CMarketDataSP market);



    /** Implementation of FourierEngine::IntoProduct interface */
    virtual FourierProduct* createProduct(const FourierEngine* model) const;

	/** Implementation of VanVSModel::IntoProduct interface */
    virtual VanVSModel::IProduct* createProduct(const VanVSModel* model) const;

protected:
	CAssetWrapper   adjAsset;			// used with capModel of VanVSModel
	bool			useMatytsinOnly;	// uses only Matytsin's Formula in Fourier for pricing			

};

typedef smartPtr<VolatilitySwap> VolatilitySwapSP;

// retrieve market data
void VolatilitySwap::GetMarket(const IModel*   model, 
								const CMarketDataSP    market) {
	static const string method = "VolatilitySwap::GetMarket";	
    try {
		const IModel* varSwapModel = 0;
		adjAsset = CAssetWrapper(asset.getName());
		if (VanVSModel::TYPE->isInstance(model)) {
			VanVSModel* vsModel = dynamic_cast<VanVSModel*>(const_cast<IModel*>(model));
			CAsset::getAssetMarketData(vsModel->capModel.get(), market.get(), ccyTreatment, 
									   discount, adjAsset);
			varSwapModel = vsModel->varSwapModel.get();
		}
		else {
			varSwapModel = model;
		}
		CAsset::getAssetMarketData(varSwapModel, market.get(), ccyTreatment, 
								   discount, asset);
		market->GetReferenceDate(valueDate);
		discount.getData(varSwapModel, market);
		instSettle->getMarket(varSwapModel, market.get());
	}
	catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Implementation of FourierEngine::IntoProduct interface */
FourierProduct* VolatilitySwap::createProduct(const FourierEngine* model) const {
    static const string routine("VolatilitySwap::createProduct");

    DateTime maturity = samples[samples.size()-1].date;
               
    return new VarOptFP(this, true, useMatytsinOnly, maturity);
}

void VolatilitySwap::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolatilitySwap, clazz);
    SUPERCLASS(VarianceSwap);
    IMPLEMENTS(VanVSModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultVolatilitySwap);
	FIELD(adjAsset,"Underlying of adjustment");
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(adjAsset);
	FIELD(useMatytsinOnly, "uses only Matytsin's formula for pricing");
	FIELD_MAKE_OPTIONAL(useMatytsinOnly);
}

IObject* VolatilitySwap::defaultVolatilitySwap(){
    return new VolatilitySwap();
}

// Type registration
CClassConstSP const VolatilitySwap::TYPE = CClass::registerClassLoadMethod(
    "VolatilitySwap", typeid(VolatilitySwap), VolatilitySwap::load);


//**********************************************************************//
//****************   PRODUCT CLASS AND METHODS    **********************//
//**********************************************************************//

class VolSwapProd: virtual public VanVSModel::IProduct {
private:
    
    const VolatilitySwap* inst; // a reference

public:
    
    VolSwapProd(const VolatilitySwap* instrument): inst(instrument){}
    
    // This is the method responsible for pricing the VarSwap and the cap and aggregating
    // the results
    void price(VanVSModel*       model,
               Control*          control, 
               CResults*         results) const;
};


// This is the method responsible for pricing the VolSwap and the adjustment and aggregating
// the results
void VolSwapProd::price(VanVSModel*       model,               
						  Control*          control, 
                          CResults*         results) const {
    static const string method = "VanVSModel::price";
    try {

		// Create a copy of the deal with no strike and insuring no mean and no scale.
		VolatilitySwapSP vswap(copy(inst));
		vswap->strikeVol = 0.0;
		vswap->notional = 0.01;
		vswap->subtractMeanVol = false;
		vswap->dontScaleByStrike = true;

        // Price the var swap
        model->varSwapModel->Price(vswap.get(), control, results);

        // Price the convexity adjustment
        CControlSP ctrl(copy(control));
        ctrl->reset();  // want this as if nothing has happened
        CResultsSP adjResults(new Results);
		vswap->asset = inst->adjAsset;
        model->capModel->Price(vswap.get(), ctrl.get(), adjResults.get());      
        
        double value = adjResults->retrievePrice() - inst->strikeVol;
		if (!inst->useMatytsinOnly) {
			value += sqrt(results->retrievePrice());
		}
		value *= 100.0 * inst->notional;
        results->storePrice(value, results->getCcyName());

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
};

VanVSModel::IProduct* VolatilitySwap::createProduct(const VanVSModel* model) const {
    return new VolSwapProd(this);
}

bool VolatilitySwapLoad()
{
    return (VolatilitySwap::TYPE != 0);
}


DRLIB_END_NAMESPACE
