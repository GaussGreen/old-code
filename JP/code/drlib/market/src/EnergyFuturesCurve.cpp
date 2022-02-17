//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFuturesCurve.cpp
//
//   Description : Energy Fututes Curve. Based on drstdcc.h and
//                 drstdcc.cpp in FXLIB.
//
//   Author      : Sean Chen
//
//   Date        : April 18, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

EnergyFuturesCurve::EnergyFuturesCurve(): EnergyCurve(TYPE)
{
}

EnergyFuturesCurve::~EnergyFuturesCurve()
{
}

void EnergyFuturesCurve::validatePop2Object()
{
    static const string method = "EnergyFuturesCurve::validatePop2Object()";
    if (expiryLabels.empty() && futureMaturityDates.empty())
        throw ModelException(method,
            "You must provide either the expiry labels or the future "
            "maturity dates.");

    if (!expiryLabels.empty() && !futureMaturityDates.empty())
        throw ModelException(method, 
            "You cannot provide the expiry labels and the future maturity "
            "dates at the same time.  Pls provide either one of them.");

    if (!expiryLabels.empty() && expiryLabels.getLength() != rates.getLength())
        throw ModelException(method,
            "Number of expiry labels does not match the number of rates");

    if (!futureMaturityDates.empty() && futureMaturityDates.getLength() != rates.getLength())
        throw ModelException(method, 
            "Number of future maturity dates does not match the number of rates");

    /*
     *  if expiry labels instead of future maturity dates are specified,
     *  convert expiry labels to future maturity dates
     *  - yet it is hard to do here since the energyUnderlyer is needed
     */
}

const DateTimeArray& EnergyFuturesCurve::getFutureMaturityDates() const
{
    if (futureMaturityDates.empty())
    {
        futureMaturityDates.resize(expiryLabels.getLength());
        EnergyUnderlyerConstSP thisEnergyUnderlyer = getEnergyUnderlyer();
        for (int i = 0; i < expiryLabels.getLength(); ++i)
        {
            futureMaturityDates[i] = thisEnergyUnderlyer->expiryDate(expiryLabels[i]);
        }
    }
    return futureMaturityDates;
}

void EnergyFuturesCurve::buildLinearInterpolant()
{
    try
    {
        int numContracts = !expiryLabels.empty() ? expiryLabels.getLength() : futureMaturityDates.getLength();
        if (numContracts == 0)
            throw ModelException("EnergyFuturesCurve::buildLinearInterpolant()", "No points defined");
    
        expiries = ExpiryArraySP(new ExpiryArray(0));
        const DateTimeArray& futdates = getFutureMaturityDates();

        // expiry dates have to be doubles for interpolator
        DateTime tempDate;
        if (!expiryLabels.empty()) { // use expiry labels
            for (int i = 0; i < numContracts; ++i) {
                tempDate = energyUnderlyerWrapper->expiryDate(expiryLabels[i]);
                expirysInDouble.push_back(tempDate.getDate());
                expiries->push_back(ExpirySP(new EnergyContractLabel(expiryLabels[i])));
            }
        }
        else { // use future maturity dates:
            for (int i = 0; i < numContracts; ++i) {
                expirysInDouble.push_back(futureMaturityDates[i].getDate());
                // empty expiry arrays in this case
            }

        }

        // building an interpolant
        LinearInterpolatorSP aInterpolatorSP = 
                 LinearInterpolatorSP( new LinearInterpolator() );
        // new() is called implicitly, return LinearInterpolantConstSP(new LinearInterpolant(xdata, fdata));
        interpolantSP = Interpolator::InterpolantSP(const_cast<Interpolator::Interpolant*>(
                            aInterpolatorSP->computeInterp(expirysInDouble, rates).get()));
    
    }
    catch(exception& e)
    {
        throw ModelException(e); 
    }
}


double EnergyFuturesCurve::interpolate(const DateTime& theContractDate) const
{
      return interpolantSP->value(theContractDate.getDate());
}

double EnergyFuturesCurve::fixing(const DateTime& theContractDate) const
{
      return interpolantSP->value(theContractDate.getDate());
}

double EnergyFuturesCurve::interpolate(const string& theContractLabel) const
{
      return interpolantSP->value(energyUnderlyerWrapper.getSP()->expiryDate(theContractLabel).getDate());
}

double EnergyFuturesCurve::fixing(const string& theContractLabel) const
{
      return interpolantSP->value(energyUnderlyerWrapper.getSP()->expiryDate(theContractLabel).getDate());
}

//int EnergyFuturesCurve::getSize() const
//{
//    return expiryLabels.size();
//}
//
//string EnergyFuturesCurve::getExpiryLabel(int i) const
//{
//    return expiryLabels[i];
//}

double EnergyFuturesCurve::getFwdRate(int i) const
{
    return rates[i];
}

double EnergyFuturesCurve::getFwdRateForLabel(const string& label) const
{
	for (int i=0;i<expiryLabels.size(); i++)
	{
		if(label == expiryLabels[i])
			return rates[i];
	}
	return -999.0;
}

void EnergyFuturesCurve::getMarket(const IModel* model, const MarketData* market)
{
    energyUnderlyerWrapper.getData(model, market);
    buildLinearInterpolant();
	if (!energyInstVolBaseWrapper.isEmpty()) // energyInstVolBaseWrapper is optional
        energyInstVolBaseWrapper.getData(model, market);
}

/** Returns name identifying the curve for the tweak */
string EnergyFuturesCurve::sensName(const EnergyFuturesCurveParallel*) const
{
    return EnergyCurve::getName(); 
}

/** Returns name identifying the curve for the tweak */
string EnergyFuturesCurve::sensName(const EnergyFuturesCurvePointwise*) const
{
    return EnergyCurve::getName(); 
}

string EnergyFuturesCurve::getName() const
{
    return name;
}

ExpiryWindowArrayConstSP EnergyFuturesCurve::sensQualifiers(const EnergyFuturesCurvePointwise*) const
{
    if (expiries->empty())
        throw ModelException("EnergyFuturesCurve::sensQualifiers(pointwise)",
            "Current implementation does not support pointwise tweak when "
            "future maturity dates are provided.  Please provide expiry "
            "labels instead.");

    return ExpiryWindow::series(expiries); 
}

// Shift method for EnergyDelta tweaks which is different from any other types of tweaking. Within one
// EnergyDelta tweak, there are two shocks and for each one of the shock, there are two shifts at each 
// point of tweaked curve. To fit in the EDR sensShift framework, we need do a lot more work...
/** Shocks the object using given shock */
TweakOutcome EnergyFuturesCurve::sensShift(const PropertyTweak<EnergyFuturesCurveParallel>& tweak)
{
    static const string method = "EnergyFuturesCurve::sensShift(parallel)";
    
    double shockSize = tweak.coefficient;
    
    if(weights.empty())
	    throw ModelException( method, "No weights for Energy Delta tweaking");
	
   
    double incRate;
    
    for (int i = 0; i < rates.size(); i++)
    {
        incRate = weights[i]*shockSize*rates[i];
        rates[i] += incRate;
    }
              
    LinearInterpolatorSP aInterpolatorSP = 
         LinearInterpolatorSP( new LinearInterpolator() );

    // delete old, replace with new
    interpolantSP.reset(const_cast<Interpolator::Interpolant*>(
                 (aInterpolatorSP->computeInterp(expirysInDouble, rates)).get()));

    return TweakOutcome(tweak.coefficient, false);
}

TweakOutcome EnergyFuturesCurve::sensShift(const PropertyTweak<EnergyFuturesCurvePointwise>& tweak)
{
    static const string method = "EnergyFuturesCurve::sensShift(pointwise)";

    if (expiries->empty())
        throw ModelException("EnergyFuturesCurve::sensQualifiers(pointwise)",
        "Current implementation does not support pointwise tweak when "
        "future maturity dates are provided.  Please provide expiry "
        "labels instead.");

    try
    {
		int i;
		double incRate;

        if (!Maths::isZero(tweak.coefficient))
		{
            i = tweak.qualifier->expiry->search(expiries.get());
			incRate = tweak.coefficient * rates[i];
			rates[i] += incRate;
            interpolantSP->editY(i, rates[i]); // making sure the interpolant is also updated
		}
  
		return TweakOutcome(incRate, false);
	}
    catch (exception& e) 
    {           
         throw ModelException(e, method,
                             "EnergyFuturesCurvePointwise failed for " +
                             getName());
    }     
}


EnergyUnderlyerConstSP EnergyFuturesCurve::getEnergyUnderlyer() const
{
    return energyUnderlyerWrapper.getSP();
}

EnergyInstVolBaseConstSP EnergyFuturesCurve::getEnergyInstVolBase() const
{
    if (!energyInstVolBaseWrapper.isEmpty()) {
	    return energyInstVolBaseWrapper.getSP();
    }
    else {
        throw ModelException("EnergyFuturesCurve::getEnergyInstVolBase", 
            "No energyInstVolBaseWrapper is specified for this energy future curve!");
        return EnergyInstVolBaseConstSP(NULL);
    }
}

class EnergyFuturesCurveHelper
{

public: 

    static IObject* defaultEnergyFuturesCurve()
    {
        return new EnergyFuturesCurve();
    }

    /** Invoked when Class is 'loaded' */

    static void load(CClassSP& clazz)
    {
        clazz->setPublic();
        REGISTER(EnergyFuturesCurve, clazz);
        SUPERCLASS(EnergyCurve);
        IMPLEMENTS(IGetMarket);
		IMPLEMENTS(IMarketFactor);
        IMPLEMENTS(ITweakableWithRespectTo<EnergyFuturesCurveParallel>);
        IMPLEMENTS(ITweakableWithRespectTo<EnergyFuturesCurvePointwise>);

        EMPTY_SHELL_METHOD(defaultEnergyFuturesCurve);
        FIELD(name,      "Energy Future Curve name");
        FIELD(energyUnderlyerWrapper,  "Energy Underlyer Wrapper");
		FIELD(energyInstVolBaseWrapper, "Energy Inst Vol Base Wrapper");
        FIELD_MAKE_OPTIONAL(energyInstVolBaseWrapper); // only use by diffusion engine
        FIELD(futureMaturityDates, "Future maturity dates.  If you do not " 
            "have future maturity dates, you can provide expiry labels as "
            "an alternative.");
        FIELD_MAKE_OPTIONAL(futureMaturityDates);
        FIELD(expiryLabels, "Contract expiry labels.  If you do not have "
            "expiry labels, you can provide future maturity dates as an "
            "alternative.");
        FIELD_MAKE_OPTIONAL(expiryLabels);
        FIELD(rates,     "Contract Rates");
		FIELD(weights,   "Weights");
		FIELD_MAKE_OPTIONAL(weights); //But mandatory if doing Delta tweaking 
        FIELD(expirysInDouble, "Expiries as double");
        FIELD_MAKE_TRANSIENT(expirysInDouble);
        FIELD(interpolantSP, "Interpolater handle");
        FIELD_MAKE_TRANSIENT(interpolantSP);
        FIELD(expiries, "Expiry type for expiryLabels");
        FIELD_MAKE_TRANSIENT(expiries);
     

        Addin::registerConstructor("ENERGY_FUTURES_CURVE",
                         Addin::MARKET,
                         "Create Energy Future Curve",
                         EnergyFuturesCurve::TYPE);

    }
};

CClassConstSP const EnergyFuturesCurve::TYPE = CClass::registerClassLoadMethod(
   "EnergyFuturesCurve", typeid(EnergyFuturesCurve), EnergyFuturesCurveHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(EnergyFuturesCurveWrapper);

DEFINE_TEMPLATE_TYPE(EnergyFuturesCurveArray);

// * for class loading (avoid having header file) */
bool EnergyFuturesCurveLoad() {
	return (EnergyFuturesCurve::TYPE != 0);
}


DRLIB_END_NAMESPACE
