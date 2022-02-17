//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : PadCDS.cpp
//
//   Description : Basic CDS implementation for pricing par instruments
//
//   Author      : Gordon Stephens
//
//   Date        : 19 April 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ParCDS.hpp"
#include "edginc/Results.hpp"
#include "edginc/CDSHelper.hpp"

DRLIB_BEGIN_NAMESPACE

//--------
// IObject
//--------

CClassConstSP const ParCDS::TYPE = CClass::registerClassLoadMethod("ParCDS", typeid(ParCDS), load);

//class load information
void ParCDS::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ParCDS, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(ClosedForm::IIntoProduct);
    IMPLEMENTS(LastSensDate);
    EMPTY_SHELL_METHOD(defaultParCDS);
    FIELD(valueDate, "valuation date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(discount, "identifies the discount curve");
    FIELD_MAKE_OPTIONAL(discount);
    FIELD(effectiveDate, "the start date of the CDS");
    FIELD(benchmark, "the par instruments maturity");
    FIELD(cdsParSpreads, "identifies the par curve");
    FIELD(spread, "the fair fee spread for the CDS");
    FIELD(frequency, "the frequency of fee payments");
}

IObject* ParCDS::defaultParCDS()
{
    return new ParCDS();
}

void ParCDS::validatePop2Object()
{
    static const string method = "ParCDS::validatePop2Object";

    try
    {
        if (!cdsParSpreads)
        {
            throw ModelException(method, "par spreads not set");
        }

        if (!discount)
        {
            throw ModelException(method, "discount curve not set");
        }

        if (!benchmark)
        {
            throw ModelException(method, "benchmark not set");
        }

    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

//------------
// CInstrument
//------------

void ParCDS::Validate()
{
    //empty
}

/** what's today ? */
DateTime ParCDS::getValueDate() const
{
    return valueDate;
}

// copy market data relevant to the instrument
void ParCDS::GetMarket(const IModel* model, const CMarketDataSP market)
{
    static const string method = "ParCDS::GetMarket";

    try
    {
        market->GetReferenceDate(valueDate);
        discount.getData(model, market);
        cdsParSpreads.getData(model, market);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


/** Returns the name of the instrument's discount currency */
string ParCDS::discountYieldCurveName() const {
    return discount.getName();
}


//-------------------------
// ClosedForm::IIntoProduct
//-------------------------

/** private class */
class ParCDSClosedForm: public ClosedForm::IProduct{
private:
    const ParCDS* cds; // a reference

public:
    ParCDSClosedForm(const ParCDS* cds): cds(cds){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        cds->price(control, results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* ParCDS::createProduct(
    ClosedForm* model) const{
    return new ParCDSClosedForm(this);
}

//-------------
// LastSensDate
//-------------

/** when to stop tweaking */
DateTime ParCDS::endDate(const Sensitivity* sensControl) const
{
    return benchmark->toDate(valueDate);
}

//-------
// ParCDS
//-------

//default constructor
ParCDS::ParCDS(): CInstrument(TYPE) {}

ParCDS::ParCDS(const string&  discountName,
               const string&  cdsParSpreadsName,
               DateTime       effectiveDate,
               ExpirySP       benchmark,
               double         spread,
               int            frequency):
    CInstrument(TYPE),
    discount(discountName),
    effectiveDate(effectiveDate),
    benchmark(benchmark),
    cdsParSpreads(cdsParSpreadsName),
    spread(spread),
    frequency(frequency)
{
    //empty
}

//closed form pricing
void ParCDS::price(Control* control, CResults* results) const
{
    static const string method = "ParCDS::GetMarket";

    try
    {
        DateTime maturity = benchmark->toDate(valueDate);

        //build fee payments
        CashFlowArraySP feePayments = CDSHelper::calculateFeePayments( valueDate,
                                                       maturity,
                                                       frequency,
                                                       cdsParSpreads->dayCountConv(),
                                                       1, //notional
                                                       spread);

        double value = 0.0;
    
        //get the default rates
        DefaultRatesSP psDefRates = cdsParSpreads->defaultRates();
        
        //price the par CDS
        CDSHelper::CredDefSwapFromDefRatesObject(
            valueDate,
            effectiveDate,
            feePayments,
            1., //notional
            cdsParSpreads->getRecovery(),
            false, //pay accrual fee
            effectiveDate,
            maturity,
            discount.getSP(),
            psDefRates.get(),
            cdsParSpreads->dayCountConv(),
            false, /* isE2C */
            &value);
        
        results->storePrice(value, discount->getCcy());
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE

