//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ImpliedLossAddIn.cpp
//
//   Description :	Add in that returns par spreads, annuities or losses generated 
//					by ImpliedLossModel
//                 
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ImpliedLossAddIn.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

void ImpliedLossAddIn::load(CClassSP& clazz)
{
        clazz->setPublic();
        REGISTER(ImpliedLossAddIn, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultImpliedLossAddIn);
        
		FIELD(market, "market");
		
		FIELD(cdoQuotes, "CDO tranche quotes");

		FIELD(indexSpreads, "Index swap spreads");

		FIELD(model, "Interpolation model");

		FIELD(expiries,"Output expiries");

		FIELD(lowStrikes, "Output low strikes");

		FIELD(highStrikes, "Output high strikes");

		FIELD(outputTypes, "outputType = 'LOSS', 'SPREAD', 'UPFRONT',  'ANNUITY', 'RISKY_ZERO' or 'CONT_LEG'");

		FIELD(coupons,"coupons when outputType = UPFRONT [default = 500bps]");
		FIELD_MAKE_OPTIONAL(coupons);
		
        
        Addin::registerObjectMethod("IMPLIED_LOSS_PRICER",
                                    Addin::RISK,
                                   "Returns results from ImpliedLossModel",
                                    false,
                                    Addin::expandMulti,
                                    &ImpliedLossAddIn::run);
}



IObject* ImpliedLossAddIn::defaultImpliedLossAddIn()
{
    return new ImpliedLossAddIn();
}

ImpliedLossAddIn::ImpliedLossAddIn(): CObject(TYPE), coupons(0)
 {}

CClassConstSP const ImpliedLossAddIn::TYPE = CClass::registerClassLoadMethod(
    "ImpliedLossAddIn", typeid(ImpliedLossAddIn), load);

// * for class loading (avoid having header file) */
bool ImpliedLossAddInLoad() {
    return (ImpliedLossAddIn::TYPE != 0);
}




IObjectSP ImpliedLossAddIn::run()
{
    static const string method = "ImpliedLossAddIn::run";

    try
    {
		// Dummy model (note that ImpliedLossModel does not derive off CModel yet)
		NonPricingModel* dummy = new NonPricingModel();

		// get market data -----------------------------------------
		indexSpreads.getData(dummy, market.get());
		cdoQuotes.getData(dummy, market.get());
		
		// cascade to cdoQuotes
		cdoQuotes->getMarket(dummy, market.get());

		delete dummy;

        // Require a model for calculating fees
        // so just provide the default closed form
        ClosedFormForwardRatePricerSP cfPricer =
            ClosedFormForwardRatePricerSP(
                new ClosedFormForwardRatePricer());

		if(coupons.get() == NULL)
		{
			// return output values ---------------------------------
			return model->getValues(
				*cdoQuotes.getSP(),
				*indexSpreads.getSP(),
				*expiries,
				*lowStrikes,
				*highStrikes,
				*outputTypes,
				cfPricer);
		}
		else
		{
			return model->getValues(
				*cdoQuotes.getSP(),
				*indexSpreads.getSP(),
				*expiries,
				*lowStrikes,
				*highStrikes,
				*outputTypes,
				*coupons,
				cfPricer);
		}


    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}
 




DRLIB_END_NAMESPACE

