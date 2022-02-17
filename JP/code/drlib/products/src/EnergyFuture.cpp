//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFuture.cpp
//
//   Description : Energy Fututes. Based on drcommodityfuture.h and
//                 drcommodityfuture.cpp in FXLIB.
//
//   Author      : Sean Chen
//
//   Date        : April 18, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenEnergyFuturePrice.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/EnergyContractLabel.hpp"
#include "edginc/EnergyFuture.hpp"



DRLIB_BEGIN_NAMESPACE

EnergyFuture::EnergyFuture(): CInstrument(TYPE), numContracts(1)
{    
}

EnergyFuture::~EnergyFuture()
{
}


string EnergyFuture::getName() const
{
    return name;
}

string EnergyFuture::discountYieldCurveName() const
{
    if (energyUnderlyerWrapper.getMO().get())
        return getEnergyUnderlyer()->getPricingCurrency();
    else
        return "";
}

double EnergyFuture::getNumContracts() const
{
    return numContracts;
}

string EnergyFuture::getBenchmark() const
{
    return benchmark;
}

DateTime EnergyFuture::getValueDate() const
{
    return valueDate;
}

double EnergyFuture::getPrice() const
{
    return price;
}

ExpirySP EnergyFuture::getExpiry(const SensControl* sensControl) const
{
    return ExpirySP(new EnergyContractLabel(getBenchmark()));
}
EnergyFuturesCurveConstSP EnergyFuture::getEnergyFuturesCurve() const
{
    return curveWrapper.getSP();
}

EnergyUnderlyerConstSP EnergyFuture::getEnergyUnderlyer() const
{
    return energyUnderlyerWrapper.getSP();
}

DateTime EnergyFuture::getValuedDate(const SensControl* sensControl) const
{
    return valueDate;
}


DateTime EnergyFuture::endDate(const Sensitivity* sensControl) const
{
    return maturityDate;
}

DateTime EnergyFuture::beginDate(const SensControl* sensControl) const
{
    return valueDate;
}

void EnergyFuture::GetMarket(const IModel* theModel, const CMarketDataSP theSP)
{
    theSP->GetReferenceDate(valueDate);
    energyUnderlyerWrapper.getData(theModel,theSP);
    curveWrapper.getData(theModel,theSP);
    maturityDate = energyUnderlyerWrapper->expiryDate(benchmark);
}

void EnergyFuture::validatePop2Object()
{
}

// Helper class
class EnergyFutureHelper
{

public: 

    static IObject* defaultEnergyFuture()
    {
        return new EnergyFuture();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic();
        REGISTER(EnergyFuture, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedFormEnergy::IIntoProduct);
        IMPLEMENTS(MonteCarlo::IIntoProduct);
        IMPLEMENTS(FirstSensDate);
		IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultEnergyFuture);
        FIELD(name,        "Name");
		FIELD_MAKE_OPTIONAL(name);
        FIELD(price,     "Contract Price");
        FIELD(benchmark,     "Contract Benchmark Label");
        FIELD(numContracts,     "Number of Contracts");
		FIELD_MAKE_OPTIONAL(numContracts);
        FIELD(description,      "Energy Future Description");
		FIELD_MAKE_OPTIONAL(description);
        FIELD(valueDate,     "Value Date");
        FIELD(energyUnderlyerWrapper,  "Energy Underlyer Wrapper");
        FIELD(curveWrapper,  "Energy Futures Curve Wrapper");
        FIELD(maturityDate,     "Contract Maturity");
        FIELD_MAKE_TRANSIENT(maturityDate);
    }
};


CClassConstSP const EnergyFuture::TYPE = CClass::registerClassLoadMethod(
   "EnergyFuture", typeid(EnergyFuture), EnergyFutureHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(EnergyFutureWrapper);

bool  EnergyFutureLoad() { return (EnergyFuture::TYPE != 0);   }

// Pricing class internally used.

class EnergyFutureClosedFormProd: public ClosedFormEnergy::IProduct
{
private:
    const EnergyFuture*  futureP; // a reference

public:
    EnergyFutureClosedFormProd(const EnergyFuture* theFuture): futureP(theFuture){}

    void price(ClosedFormEnergy* model,Control* control,CResults* results) const;
};

ClosedFormEnergy::IProduct* EnergyFuture::createProduct(ClosedFormEnergy* model) const
{
        return new EnergyFutureClosedFormProd(this);
}

void EnergyFutureClosedFormProd::price(ClosedFormEnergy* model,Control* control, 
                                   CResults* results) const
{
        static const string method = "EnergyFutureClosedFormProd::price";

        try 
        {    
            DateTime theValueDate = futureP->getValueDate();
        
            DateTime theMatDate = futureP->getEnergyUnderlyer()->expiryDate(futureP->getBenchmark());
            double mtmValue;

            if (theValueDate >= theMatDate ) 
            {
                mtmValue = 0.0;
            }
            else
            {
                double marketPrice = futureP->getEnergyFuturesCurve()->fixing(theMatDate);
                // mtm value > 0 if market price on curve > price on contract
                mtmValue = ( marketPrice - futureP->getPrice() )
                                * futureP->getNumContracts() 
                                * futureP->getEnergyUnderlyer()->getNotional();
            }
            results->storePrice(mtmValue,futureP->getEnergyUnderlyer()->getPricingCurrency());

        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
};


// view class for Monte Carlo state variables framework
class EnergyFutureMC : public MCProductClient 
{

private:
    const EnergyFuture*                  instrument;
    SVGenEnergyFuturePrice::IStateVarSP          jsv;
    SVGenEnergyFuturePriceSP                     jsvGen;

    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) 
	{
         jsv=jsvGen->getEnergyFuturePriceSV(newPathGen);
    }

public:
    EnergyFutureMC(const EnergyFuture* inst, const SimSeriesSP&   simSeries,
                   InstrumentSettlementSP&   instSettle) :
    MCProductClient(IMultiMarketFactors::asMulti(inst->energyUnderlyerWrapper.getSP()).get(),
                      inst->getValueDate(),
                      ((inst->energyUnderlyerWrapper).getSP())->getYieldCurve(),
                      IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)),
                      simSeries,
                      IPastValuesSP(IPastValues::Util::makeTrivial(inst->valueDate, 0.0)),
                      instSettle.get(),
                      inst->getValueDate()), instrument(inst)
     {
          DateTimeArray dates=simSeries->getAllDates();
          DateTime      maturity=inst->getEnergyUnderlyer()->expiryDate(inst->getBenchmark());
          DateTime::MonthDayYear mdy = maturity.toMDY();
          EnergyContractLabelArray labels;
          EnergyContractLabelSP eclptr(new EnergyContractLabel(mdy.month,mdy.year));
          labels.push_back(eclptr);
          SVGenEnergyFuturePriceSP genPtr(new SVGenEnergyFuturePrice(dates,labels));
          jsvGen=genPtr;
     }

     virtual void collectStateVars(IStateVariableCollectorSP svCollector) const 
	 {
           svCollector->append(jsvGen.get());
     }

	 // S. Chen Prices changed to IMCPrices
     void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices) 
	 {
           double spot=jsv->getSpot();
           prices.add(spot);
     }

};

IMCProduct* EnergyFuture::createProduct(const MonteCarlo* model) const  
{
    SimSeriesSP simSeries(new SimSeries(1));
    DateTimeArray dates;
    DateTime      maturity=this->getEnergyUnderlyer()->expiryDate(this->getBenchmark());
    dates.push_back(this->getValueDate());
    dates.push_back(maturity);
    simSeries->addDates(dates);
    InstrumentSettlementSP instSettle(new CashSettleDate(maturity));
    return new EnergyFutureMC(this, simSeries, instSettle);

}


DRLIB_END_NAMESPACE

