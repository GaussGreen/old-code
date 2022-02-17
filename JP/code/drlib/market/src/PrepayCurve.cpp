//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : PrepayCurve.cpp
//
//   Description : A prepayment curve for ABS
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/PrepayCurve.hpp"
#include "edginc/Format.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Equals.hpp"
#include "edginc/RollingSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

PrepayCurve::PrepayCurve(CClassConstSP clazz) : DecretionCurve(TYPE)
{
    isFlatBalances = true;
    prepaySettle = RollingSettlementSP(new RollingSettlement());
}


PrepayCurve::~PrepayCurve() {}


void PrepayCurve::validatePop2Object() 
{
    static const string method ="PrepayCurve::validatePop2Object";
    try{
        if(!dates.get() || dates->empty()) {
            throw ModelException(method, name + " prepayment dates are null or empty");
        }

        if(!rates.get() || rates->empty()) {
            throw ModelException(method, name + " prepayment balances are null or empty");
        }

        if(!factorDates.get() || factorDates->empty()) {
            throw ModelException(method, name + " factor dates are null or empty");
        }

        if(!factors.get() || factors->empty()) {
            throw ModelException(method, name + " factor history is null or empty");
        }

        /** ensure dates are increasing */
        DateTime::ensureStrictlyIncreasing(*dates, "prepayment dates", false);
        DateTime::ensureStrictlyIncreasing(*factorDates, "factor dates", false);

        /** factors */
        int sizeDates, sizeRates;
        sizeDates = factorDates->size();
        sizeRates = factors->size();

        if(sizeDates != sizeRates) {
            throw ModelException(method, name + " has a different number of "
                                 "factorDates (" + Format::toString(sizeDates) + ") than"
                                 " factors (" + Format::toString(sizeRates) + ")" );
        }

        /** balances */
        balances = DoubleArraySP(new DoubleArray(*rates));

        sizeDates = dates->size();
        sizeRates = rates->size();

        if(sizeDates != sizeRates) {
            throw ModelException(method, name + " has a different number of "
                                 "dates (" + Format::toString(sizeDates) + ") than"
                                 " balances (" + Format::toString(sizeRates) + ")" );
        }

        for(int i = 0; i < sizeDates; i++)
        {
            if((*balances)[i] < DBL_EPSILON) 
            {
                (*balances)[i] = DBL_EPSILON;
            }
        }

        if(!isStepBalances()) {
            /** assume exponential distribution between balances, i.e.
             ** b2 = b1 * exp(-speed * dt)
             */
            speeds = DoubleArraySP(new DoubleArray(sizeDates));
            speeds->back() = 0;
            double b1 = (*balances)[0];
            DateTime d1 = (*dates)[0];
            for(int i = 0; i < sizeDates - 1; i++) {
                double b2 = (*balances)[i+1];
                DateTime d2 = (*dates)[i+1];
                (*speeds)[i] = log(b1 / b2) / Actual365F::yearFraction(d1, d2);
                b1 = b2;
                d1 = d2;
            }
        }

    }
    catch(exception& e) {
        throw ModelException(e, method, "Validation failed for " + name);
    }
}

//------------------------------------------
// IDecretionCurve methods
//------------------------------------------

/** pv, which actually is balance here, can be relative to a start date */
double PrepayCurve::pv(const DateTime& startDate,
                          const DateTime& endDate) const
{
    static const string method = "PrepayCurve::pv";
    try {
        double p1, p2, p;
        p1 = pv(startDate);
        p2 = pv(endDate);
        p = p2 / p1;
        return p;
    }
    catch(exception &e) {
        throw ModelException(e, method);
    }
}
    
/** pv (balance) relative to initial balance */
double PrepayCurve::pv(const DateTime& endDate) const
{
    static const string method = "PrepayCurve::pv";
    try {
        double p1, p2, p;

        /** startDate is the initial date, ie. dates[0] */
        p1 = balances->front();
        if(Maths::isZero(p1)) {
            throw ModelException("Initial Balance is Zero", method);
        }

        if(endDate <= dates->front() ) {
            p2 = balances->front();
        }
        else if( endDate > dates->back() ) {
            p2 = balances->back();
        }
        else{
            int idx = endDate.findUpper(*dates) - 1;
            p2 = (*balances)[idx];
            if(!isStepBalances()) {
                DateTime d = (*dates)[idx];
                p2 *= exp(-(*speeds)[idx] * Actual365F::yearFraction(d, endDate));
            }
        }

        p = p2 / p1;
        return p;
    }
    catch(exception &e) {
        throw ModelException(e, method);
    }
}

/** return decretion speed on a date */
double PrepayCurve::getDecretionSpeed(const DateTime& date) const
{
    static const string method = "PrepayCurve::getDecretionSpeed";
    try {
        double speed;
        
        if( date <= dates->front() || date >= dates->back() ) {
            speed = 0;
        } else {
            int idx = date.findLower(*dates);
            speed = (*speeds)[idx];
        }

        return speed;
    }
    catch(exception &e) {
        throw ModelException(e, method);
    }
}

/**/
DateTimeArraySP PrepayCurve::getStepDates() const
{
    if (stepDates.get())
        return stepDates;
    if (dates->size() < 1)
        throw ModelException("empty dates array");
    
    stepDates = DateTimeArraySP(new DateTimeArray());

    int i = 0;
    double r = (*rates)[i];
    while (i < dates->size()) {
        if (r == (*rates)[i]) {
            ++i;
            continue;
        }
        stepDates->push_back((*dates)[i]);
        r = (*rates)[i];
        ++i;
    }
    stepDates->push_back((*dates)[i-1]);
    return stepDates;
}

/* for a given structure like 
   |  a |  b |  c |  d |  e |  f |
   |  1 |  1 | .8 | .6 | .6 | .6 |
   the getStepDates method will return c,d,f
   whereas the useful dates are b,c,d
   since they define the start->end where the steps occur
   so this method is designed to return the latter
   However, since the indexing into these dates for the pv
   factor is offset by -1, we need to +1 on the date pair
*/
DateTimeArraySP PrepayCurve::getAlternativeStepDates() const
{
    if (dates->size() < 1)
        throw ModelException("empty dates array");
    
    DateTimeArraySP altStepDates = DateTimeArraySP(new DateTimeArray());

    int i = 0;
    double r = (*rates)[i];
    while (++i < dates->size()-1)
    {
        if ((*rates)[i] != r)
        {
            //record the two dates
            altStepDates->push_back((*dates)[i]); //start
            altStepDates->push_back((*dates)[i+1]); //end
            r = (*rates)[i];
        }
    }

    //may have repeating dates, so get rid of them now, but consider time of day
    DateTime::removeDuplicates(*(altStepDates.get()),false);

    return altStepDates;
}

/** return if balances are stepwise or continuous */ 
bool PrepayCurve::isStepBalances() const
{
    return isFlatBalances;
}

/** Returns a reference to the value date - allows default implementations
 *  to work */
double PrepayCurve::getFactor(const DateTime& date) const
{
    static const string method = "PrepayCurve::getFactor";
    int idx = date.findLower(*factorDates);
    double f;
    if(idx < 0) 
    {
        f = factors->front();
    } else {
        f = (*factors)[idx];
    }
    return f;
}

//------------------------------------------
// Tweak methods
//------------------------------------------
/** Returns a name */
string PrepayCurve::sensName(const PrepayParallelTP*) const
{
    return name;
}

TweakOutcome PrepayCurve::sensShift(const PropertyTweak<PrepayParallelTP>& tweak)
{
    static const string method = "PrepayCurve::sensShift";
    
    try{
        if(!rates.get())
        {
            throw ModelException(method,
                                 "parallel tweak none-applicable. Check inputs. PrepayParallel failed for " +
                                 getName());

        } else {
            if (!Maths::isZero(tweak.coefficient)) {
                DateTimeArraySP datesS = DateTimeArraySP(new DateTimeArray(dates->getLength()));
                DoubleArraySP   ratesS = DoubleArraySP(new DoubleArray(rates->getLength()));
                int daysOffset = tweak.coefficient * 365;
                int len = dates->getLength();

                (*datesS)[0] = (*dates)[0].rollDate(daysOffset);
                for(int i = 0; i < len; i++)
                {
                    (*datesS)[i] = (*dates)[i].rollDate(daysOffset);
                    (*ratesS)[i] = (*rates)[i];
                }
                
                //fill in new rates
                for (int i = 0; i < dates->getLength(); i++) {
                    //for convenience
                    DateTime& datei = (*dates)[i];
                    double&   ratei = (*rates)[i];

                    if(datei <= datesS->front()) 
                    {
                        ratei = ratesS->front();

                    }
                    else if (datei >= datesS->back()) 
                    {
                        ratei = ratesS->back();
                    }
                    else {
                        int idx = datei.findUpper(*datesS) - 1;
                        double dr = (*ratesS)[idx+1] - (*ratesS)[idx];
                        DateTime d1 = (*datesS)[idx];
                        DateTime d2 = (*datesS)[idx+1];
                        double dx1, dx2;
                        dx1 = d1.yearFrac(datei);
                        dx2 = d1.yearFrac(d2);
                        ratei = (*ratesS)[idx] + dx1 / dx2 * dr;
                    }
                    (*balances)[i] = (ratei < DBL_EPSILON) ? DBL_EPSILON : ratei;
                }
                (*rates)[0]    = 1;
                (*balances)[0] = 1;
            }
        }
            
        // none of our components has a rho type sensitivity
        return TweakOutcome(tweak.coefficient, false); 
    } catch (exception& e){
        throw ModelException(e, method,
                             "PrepayParallel failed for " + getName());
    }
}

SettlementConstSP PrepayCurve::getSettlement() const
{
    return prepaySettle;
}

//---------------------------------------------------------
// MarketObject methods
//---------------------------------------------------------
void PrepayCurve::getMarket(const IModel* model, const MarketData *market)
{
    market->GetReferenceDate(today);
    if(prepaySettle.get())
    {
        prepaySettle->getMarket(model, market);
    }
}

int PrepayCurve::hashCode() const 
{
    int hcode = (size_t)getClass();
    hcode ^= dates->hashCode();
    hcode ^= rates->hashCode();
    return hcode;
}

//** Compare if equal */
bool PrepayCurve::equalTo(const IObject* curve) const
{
    try
    {
        if (this == curve)    // obvious first test
            return true;

        if (!curve || getClass() != curve->getClass())
            return false;
        
        const PrepayCurve* ic2 = STATIC_CAST(PrepayCurve, curve);

        if (!NumericArrayEquals<DoubleArray>(rates, ic2->rates))
            return false;

        if (!NumericArrayEquals<DoubleArray>(balances, ic2->balances))
            return false;

        return true;
    }
    catch (exception& e)
    {
        throw ModelException(e, "PrepayCurve::equalTo");
    }
}

/** clone */
IObject *PrepayCurve::clone() const
{
    PrepayCurve *copy = DYNAMIC_CAST(PrepayCurve, MarketObject::clone());

    copy->today = today;
    copy->isFlatBalances = isFlatBalances;
    copy->factorDates = DateTimeArraySP(new DateTimeArray(*factorDates));
    copy->factors = DoubleArraySP(new DoubleArray(*factors));
    copy->dates = DateTimeArraySP(new DateTimeArray(*dates));
    copy->rates = DoubleArraySP(new DoubleArray(*rates));
    copy->balances = DoubleArraySP(new DoubleArray(*balances));
    if(speeds.get()) 
    {
        copy->speeds = DoubleArraySP(new DoubleArray(*speeds));
    }
    if(prepaySettle.get())
    {
        copy->prepaySettle = SettlementSP(prepaySettle.clone());
    }

    return copy;
}


//---------------------------------------------------------
// Invoked when Class is 'loaded'
//---------------------------------------------------------
void PrepayCurve::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(PrepayCurve, clazz);
    SUPERCLASS(DecretionCurve);
	IMPLEMENTS(ITweakableWithRespectTo<PrepayParallelTP>);
    EMPTY_SHELL_METHOD(defaultPrepayCurve);
    FIELD(today, "today");
    FIELD_MAKE_TRANSIENT(today);

    FIELD(isFlatBalances, "is the balances step-wise?");
    FIELD_MAKE_OPTIONAL(isFlatBalances);
    FIELD(factorDates, "factor dates");
    FIELD(factors, "corresponding factor history");
    FIELD(dates, "dates");
    FIELD(rates, "decretion rates, balances for prepay, loss rate for loss");
    FIELD(speeds, "decretion speeds");
    FIELD_MAKE_TRANSIENT(speeds);
    FIELD(balances, "balances");
    FIELD_MAKE_TRANSIENT(balances);
    FIELD(prepaySettle, "delay of prepay in business days. default=0");
    FIELD_MAKE_OPTIONAL(prepaySettle);
}

IObject* PrepayCurve::defaultPrepayCurve(){
    return new PrepayCurve();
}

//---------------------------------------------------------
// Static variables
//---------------------------------------------------------

CClassConstSP const PrepayCurve::TYPE = CClass::registerClassLoadMethod(
    "PrepayCurve", typeid(PrepayCurve), load);

DEFINE_TEMPLATE_TYPE(PrepayCurveWrapper);



DRLIB_END_NAMESPACE
