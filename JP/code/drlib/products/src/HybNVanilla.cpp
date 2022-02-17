//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HybNVanilla.cpp
//
//   Description : An European option with notional depending on rates performance.
//
//   Author      : Lei Fang
//
//   Date        : 28 Nov 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/FunctionWrapper.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VegaMatrixLite.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/IndexSpecIR.hpp"

DRLIB_BEGIN_NAMESPACE

class HybNVanilla : public Generic1Factor,
                    virtual public FDModel::IIntoProduct,
                    virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz)
    {
        clazz->setPublic();
        REGISTER(HybNVanilla, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultHybNVanilla);
        FIELD(ir, "interest rate leg");
        FIELD(rates, "rates");
        FIELD(leverageFactors, "leverage factors");
        FIELD(interps, "interpolation styles between rates");
        FIELD(rateAtMaturity, "swap rate at maturity");
        FIELD_MAKE_OPTIONAL(rateAtMaturity);
        FIELD(isCall, "is it a call option");
        FIELD(exerciseSchedule, "exercise schedule");
        FIELD(spotAtMaturity, "spot at maturity");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
        FIELD(resetDates, "maturity");
        FIELD_MAKE_TRANSIENT(resetDates);
    }

    static IObject* defaultHybNVanilla() { return new HybNVanilla(); }

    virtual void GetMarket(const IModel* model, const CMarketDataSP market);

    virtual void Validate();
 
    virtual void validatePop2Object();

	virtual bool sensShift(Theta* shift);

    virtual DateTime getValueDate() const;

    virtual FDProductSP createProduct(FDModel* model) const;

    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;

private:
    HybNVanilla() : Generic1Factor(TYPE), rateAtMaturity(0.0), spotAtMaturity(0.0) {}

    // registered fields
    IndexSpecIRSP ir;
    DoubleArray rates;
    DoubleArray leverageFactors;
    StringArray interps;
    double rateAtMaturity;
    bool isCall;
    ScheduleSP exerciseSchedule;
    double spotAtMaturity;
    // transient fields
    DateTimeArray resetDates;

    friend class HybNVanillaFDProd;
};

CClassConstSP const HybNVanilla::TYPE = CClass::registerClassLoadMethod(
    "HybNVanilla", typeid(HybNVanilla), HybNVanilla::load);

// product class
class HybNVanillaFDProd : public FDProduct
{
public:
    HybNVanillaFDProd(const HybNVanilla* inst, FDModel* m) :
        FDProduct(m), inst(inst)
    {
        eqIndex = model->createProduct( IProdCreatorSP( new
            IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

        irIndex = model->createProduct(inst->ir);
		irIndex->addModelResetDates(inst->ir->resetDates, inst->ir->resetDates);
    }

    virtual string getCcyTreatment() const { return inst->ccyTreatment; }
 
    virtual void init(Control* control) const;

	virtual void initProd() {}

    virtual void update(int& step, FDProduct::UpdateType type);

    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results);

    static const string INTERP_LINEAR;
    static const string INTERP_STAIRS;

private:
    const HybNVanilla* inst;

    TreeSliceSP value;
    FDProductSP irIndex;
    FDProductSP eqIndex;
};

const string HybNVanillaFDProd::INTERP_LINEAR = "L";
const string HybNVanillaFDProd::INTERP_STAIRS = "S";

namespace
{
class Interpolate : public SliceMarker< Interpolate >
{
public:
    Interpolate(
        const TreeSlice & rate,
        const DoubleArray & rates,
        const DoubleArray & leverageFactors,
        const StringArray & interps )
        :
        rate( rate ),
        rates( rates ),
        leverageFactors( leverageFactors ),
        interps( interps )
    {}

    // TreeSlice "expression template" primitives
    static const int sliceCount = 1;
    template< typename S >
    const S** listSlices(const S** list) const
    {
        return rate.listSlices( list );
    }
    double calc() const
    {
        return apply( rate.calc() );
    }
    void printDebug(char *s) const
    {
        strcat(s, "(Interpolate)");
    }

private:
    double apply( double rate ) const
    {
        // assume rates is an increasing array, bisection search
        int numRates = rates.size();
        if ( rate >= rates[ numRates - 1 ] ) return leverageFactors[ numRates - 1 ];

        int low = 0, high = numRates - 1, mid = 0;
        while ( high - low > 1 )
        {
            mid = (high + low) / 2;
            if ( rate < rates[mid] )
            {
                high = mid;
            } else {
				low = mid;
            }
        }        
        if ( rate < rates[mid]) --mid;

        // do interpolation bwtween mid and mid+1
        if ( interps[mid] == HybNVanillaFDProd::INTERP_STAIRS )
        {
            return leverageFactors[mid];
        } else if ( interps[mid] == HybNVanillaFDProd::INTERP_LINEAR ) {
            return ( ( rates[mid + 1] - rate ) * leverageFactors[mid] + 
                ( rate - rates[mid] ) * leverageFactors[mid + 1] ) /
                ( rates[mid+1] - rates[mid] );
		} else {
			throw ModelException("HybNVanilla::Interpolate ", "unrecognized interpolation method");
		}
    }

    const TreeSlice & rate;
    const DoubleArray & rates;
    const DoubleArray & leverageFactors;
    const StringArray & interps;
};
}

void HybNVanillaFDProd::init(Control* control) const
{
    DateTimeArray critDates = inst->exerciseSchedule->getDates();
    model->addCritDates( critDates );

    DateTimeArray segDates;
    segDates.resize(2);
    segDates[0] = inst->valueDate;
    segDates[1] = inst->exerciseSchedule->lastDate();
    IntArray density(1,1);
    model->initSegments( segDates, density );
}

void HybNVanillaFDProd::update(int& step, FDProduct::UpdateType type)
{
    DateTime currDate = model->getDate( step );
    int i;
    for (i = 0; i < inst->ir->resetDates.size(); ++i) {
        if (inst->ir->resetDates[i] == currDate) {
            break;
        }
    }
	
	if (i == inst->ir->resetDates.size() - 1)
	{
		startDEV( value = model->createSlice( inst->discount->getName() ) );
		value->name = "HybNVanilla";
	}

    if (i < inst->ir->resetDates.size()) {
		const TreeSlice & spot = eqIndex->getValue(step);
        const TreeSlice & rate = irIndex->getValue(step);

        if (inst->isCall) {
            *value = smax(spot - inst->exerciseSchedule->lastValue(), 0.) *
                Interpolate( rate, inst->rates, inst->leverageFactors, inst->interps );
        } else {
            *value = smax(inst->exerciseSchedule->lastValue() - spot, 0.) *
                Interpolate( rate, inst->rates, inst->leverageFactors, inst->interps );
		}
    }
}

void HybNVanillaFDProd::recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
{
    double price = model->getPrice0( *value );
    results->storePrice( price, disc->getCcy() );
}

void HybNVanilla::GetMarket(const IModel* model, const CMarketDataSP market)
{
    market->GetReferenceDate(valueDate);

    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discount, asset);
 
    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
    
    ir->setup(model, market.get());
    if (ir->factor.getName().empty()) ir->factor = YieldCurveWrapper(asset->getYCName());
    ir->factor.getData(model, market);
    ir->addResetDates(resetDates);
}

void HybNVanilla::Validate()
{
    static const string method = "HybNVanilla::Validate";
    
	DateTime matDate = exerciseSchedule->lastDate();
	DateTime settlementDate = instSettle->settles(matDate, asset.get());
	if ( (valueDate >= matDate) && (valueDate < settlementDate) ) 
	{
		if ( (0.0 == spotAtMaturity) || (0.0 == rateAtMaturity) )
		{
			throw ModelException(method, "Spot and swap rate at maturity must be provided");
		}
	}
}

void HybNVanilla::validatePop2Object()
{
    static const string method = "HybNVanilla::validatePop2Object";
    
	// set reset dates to be maturity date
    resetDates = exerciseSchedule->getDates();

    // can't get exercise schedule from Market - fail if it is NULL
    if( !exerciseSchedule )
    {
        throw ModelException(method, "Exercise schedule is NULL");
    }

    // can't get instrument settlement from Market - fail if it is NULL
    if( !instSettle )
    {
        throw ModelException(method, "Instrument settlement is NULL");
    }

    // check that we have at least one entry in the exercise schedule
    int numDates = exerciseSchedule->length();
    if ( numDates < 1)
    {
        throw ModelException(method, "Exercise schedule is empty");
    }
 
	int rsize = rates.size(), lsize = leverageFactors.size(), isize = interps.size();
	if ( (rsize != lsize) || (rsize != isize) ) 
	{
		throw ModelException(method, "All array sizes should be equal");
	} else if (rsize < 1) {
		throw ModelException(method, "Array size should be at least 1");
	}

	if (rates[0] != 0.0) 
	{
		throw ModelException(method, "The first element of rates array should be 0");
	}

	for ( int i = 0; i < rsize - 1; ++i )
	{
		if ((rates[i] >= rates[i+1]))
		{
			throw ModelException(method, "The rates array must be an increaing array");
		}
	}

	for ( int i = 0; i < isize; ++i )
	{
		if ((interps[i] != HybNVanillaFDProd::INTERP_LINEAR) &&
			(interps[i] != HybNVanillaFDProd::INTERP_STAIRS))
		{
			throw ModelException(method, "The interpolation methods can only be linear or stairs");
		}
		if (i == isize - 1)
		{
			if (interps[i] != HybNVanillaFDProd::INTERP_STAIRS)
			{
				throw ModelException(method, "The last interpolation method must be S");
			}
		}
	}
}

bool HybNVanilla::sensShift(Theta* shift)
{
	DateTime newDate = shift->rollDate(valueDate);
	DateTime matDate = exerciseSchedule->lastDate();
	
	if ( ( newDate >= matDate && valueDate < matDate ) ||
		( valueDate == matDate && Maths::isZero(spotAtMaturity) ) )
		spotAtMaturity = asset->getThetaSpotOnDate(shift, matDate);

	valueDate = newDate;
	return true;
}

DateTime HybNVanilla::getValueDate() const
{
    return valueDate;
}

bool HybNVanilla::priceDeadInstrument(CControl* control, CResults* results) const
{
    static const string method = "HybNVanilla::priceDeadInstrument";
    
    DateTime matDate = exerciseSchedule->lastDate();
    bool expired = (valueDate >= matDate);
    if (!expired)
        return false;

    DateTime settlementDate = instSettle->settles(matDate, asset.get());
    if (valueDate >= settlementDate) 
	{
        results->storePrice(0.0, discount->getCcy());
        return true;
	}

    double eqStrike = exerciseSchedule->lastValue();
    double value = GetIntrinsic(spotAtMaturity, eqStrike, isCall, true);

	// cannot use Interpolate, copy Interpolate::apply here
    double notional = 0.0;
    int numRates = rates.size();
    if ( rateAtMaturity >= rates[ numRates - 1 ] ) {
		notional = leverageFactors[ numRates - 1 ];
	} else {
		int low = 0, high = numRates - 1, mid = 0;
		while ( high - low > 1 ) {
			mid = (high + low) / 2;
			if ( rateAtMaturity < rates[mid] ) {
				high = mid;
			} else {
				low = mid;
			}
		}        
		if ( rateAtMaturity < rates[mid]) --mid;
		if ( interps[mid] == HybNVanillaFDProd::INTERP_STAIRS ) {
			notional = leverageFactors[mid];
		} else if ( interps[mid] == HybNVanillaFDProd::INTERP_LINEAR ) {
			notional = ( ( rates[mid + 1] - rateAtMaturity ) * leverageFactors[mid] + 
                ( rateAtMaturity - rates[mid] ) * leverageFactors[mid + 1] ) /
                ( rates[mid+1] - rates[mid] );
		}
	}

    value *= notional * discount->pv(valueDate, settlementDate);
	results->storePrice(value, discount->getCcy());
    return true;
}

DateTime HybNVanilla::endDate(const Sensitivity* sensControl) const
{
    DateTime matDate = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

FDProductSP HybNVanilla::createProduct(FDModel* model) const
{
    return FDProductSP( new HybNVanillaFDProd(this, model) );
}

// for class loading
bool HybNVanillaLoad()
{
	return (HybNVanilla::TYPE != 0);
}

DRLIB_END_NAMESPACE
