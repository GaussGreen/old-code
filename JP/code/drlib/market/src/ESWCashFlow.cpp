//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWCashFlow.hpp
//
//   Description   cash flow class for ESW
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ESWCashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

/////////////////////////////////////////////
////////////////// ESWCashFlow /////////
/////////////////////////////////////////////
// helpers
void ESWCashFlow::load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(ESWCashFlow, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultESWCashFlow);
		FIELD(LegName, "Name of CashFlow leg");
		FIELD_MAKE_OPTIONAL(LegName);
        FIELD(SpotFX,  "FX spot for ESW cash flow leg");
        FIELD(PayDates,  "Pay date array");
        FIELD(Amounts,  "Pay amount array");
        FIELD(discount,  "Discount ccy");
}


CClassConstSP const ESWCashFlow::TYPE = CClass::registerClassLoadMethod(
    "ESWCashFlow", typeid(ESWCashFlow), load);

// constructor
ESWCashFlow::ESWCashFlow(): CObject(TYPE), LegName("CF") {}

// validation
void ESWCashFlow::Validate()
{
	static string method = "ESWCashFlow::Validate";
	
	if (PayDates.size() != Amounts.size())
		throw ModelException(method, "date list and amount list are different sizes");

	if (SpotFX <= 0.0)
		throw ModelException(method, "invalid SpotFX input for cash flow leg");
}

// cash flow pv'ed inclusive of start and end dates
double ESWCashFlow::priceCF(const DateTime& valDate, const DateTime& endDate, 
                            const DateTime& callSettleDate, bool isCallable,
							const Control* control, CResults* results,
                            OutputRequestUtil::KnownCashFlows& knownCFs, 
                            DateTimeArray& payment_dates) const
{	
    int i;
    double result = 0.0;
    CashFlowArraySP cfCF(new CashFlowArray); // just for output printing

    for (i=0; i<PayDates.size(); i++)
    {// include future cf from value date exclusive and call date inclusive
        if (PayDates[i] > valDate && (!isCallable || PayDates[i] <=callSettleDate))
        {
            result += Amounts[i]*discount->pv(PayDates[i]);
            // prepare printing if needed
            if ((control)&&(control->isPricing()))
            {
				cfCF->push_back(CashFlow(PayDates[i], Amounts[i]));
            }

            knownCFs.addKnownCashFlow(discount->getCcy(), CashFlow(PayDates[i], Amounts[i]*SpotFX));
            
        }

        // add payment dates
        payment_dates.push_back(PayDates[i]);
    }
    // simply multiply by FXSpot input
    result *= SpotFX;
    // print if needed
    if ((control)&&(control->isPricing()) && cfCF->size()>0)
    {
		OutputNameConstSP cfOutput(new OutputName(LegName));
		results->storeGreek(cfCF, Results::DEBUG_PACKET, cfOutput);
    }
    return result;
}

DRLIB_END_NAMESPACE
