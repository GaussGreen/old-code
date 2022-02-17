//----------------------------------------------------------------------------
//
//   Group       : Conertible Research
//
//   Filename    : ParSpreadCurve.cpp
//
//   Description : a par spread curve
//
//   Author      : Andre Segger
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_PARSPREADCURVE_CPP
#include "edginc/MaturityPeriod.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/ParSpreadCurve.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadCurve::ParSpreadCurve() : MarketObject(TYPE), upfronts(0) {};


ParSpreadCurve::ParSpreadCurve(const string&    parName,
                               const double&    parSpread,
                               const string&    maturity) : MarketObject(TYPE), name(parName), upfronts(0)
{
    // expiry array
    MaturityPeriodSP mat( new MaturityPeriod(maturity));
    expiries = ExpiryArraySP( new ExpiryArray(1));
    (*expiries)[0] = mat;

    // maturity period array
    spreads = DoubleArray(1);
    spreads[0] = parSpread;
}

ParSpreadCurve::ParSpreadCurve(const string&    parName,
                               const double&    parSpread,
                               const DateTime&  maturity) : MarketObject(TYPE), name(parName), upfronts(0)
{
    // expiry array
    BenchmarkDateSP mat( new BenchmarkDate(maturity));
    expiries = ExpiryArraySP( new ExpiryArray(1));
    (*expiries)[0] = mat;

    // maturity period array
    spreads = DoubleArray(1);
    spreads[0] = parSpread;
}

ParSpreadCurve::ParSpreadCurve(const string& name,
                               ExpiryArraySP expiries,
                               const CDoubleArray& spreads) 
    : MarketObject(TYPE), name(name), expiries(expiries), spreads(spreads), upfronts(0){}                 


//constructors with upfront fee
ParSpreadCurve::ParSpreadCurve(const string&    parName,
                               const double&    parSpread,
                               const double&    upfront,
                               const string&    maturity) : MarketObject(TYPE), name(parName)
{
    // expiry array
    MaturityPeriodSP mat( new MaturityPeriod(maturity));
    expiries = ExpiryArraySP( new ExpiryArray(1));
    (*expiries)[0] = mat;

    // maturity period array
    spreads = DoubleArray(1);
    spreads[0] = parSpread;

    // upfronts
    upfronts = CDoubleArraySP(new CDoubleArray(1));
    upfronts->push_back(upfront);
}

ParSpreadCurve::ParSpreadCurve(
    const string&       name,
    ExpiryArraySP       expiries,
    const CDoubleArray& spreads,
    CDoubleArraySP      upfronts) :
    MarketObject(TYPE), name(name), expiries(expiries), spreads(spreads), 
    upfronts(upfronts)
{}                 


static IObject* defaultParSpreadCurve(){
    return new ParSpreadCurve();
}


/** Get today and the par spreads curve from the market data cache */
void ParSpreadCurve::getMarket(const IModel* model, const MarketData *market){
    market->GetReferenceDate(today);
}


void ParSpreadCurve::validatePop2Object() {
    static const string method("ParSpreadCurve::validatePop2Object");
    try {
        if (!expiries.get() || expiries->empty()) {
            throw ModelException(method, name + " expiries are null or empty");
        }

        if (expiries->size() != spreads.size()) {
            throw ModelException(method, name + " has a different number of "
                                 "expiries (" +
                                 Format::toString(expiries->size()) +
                                 ") than spreads (" +
                                 Format::toString(spreads.size()) + ")");
        }
        

        // has upfronts input
        if(upfronts.get()) {
            int s = expiries->size();
            if(upfronts->size() != s) 
            {
                throw ModelException(method, name + " has different number of "
                                     "expiries (" + Format::toString(s) + 
                                     ") than upfronts (" +
                                     Format::toString(upfronts->size()) + ")");
            }
        }
        
        for (int i = 0; i < expiries->size(); i++) {
            if (Maths::isNegative(spreads[i])) {
                throw ModelException(method, "par spread at " + 
                                     (*expiries)[i]->toString() + 
                                     " is negative (" + 
                                     Format::toString(spreads[i]) + ")");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method, "validation failed for " + name);
    }
}

string ParSpreadCurve::getName() const {
    return name;
}

// same interpolation used as for CreditSpreads at present
double ParSpreadCurve::getCurrentSpread(const DateTime& valueDate,
                                        const DateTime& maturityDate,
                                        const BadDayConvention* bdc,
                                        const Holiday* hols) const {
    static const string method = "ParSpreadCurve::getCurrentSpread";

    double currentSpread = 0.0;

    if (maturityDate <= bdc->adjust((*expiries)[0]->toDate(valueDate), hols)) {
        currentSpread = spreads[0];
    } else if (maturityDate >= bdc->adjust(expiries->back()->toDate(valueDate), hols)) {
        currentSpread = spreads.back();
    } else {
        DateTime previousDate = bdc->adjust((*expiries)[0]->toDate(valueDate), hols);
        DateTime currentDate;
        for (int i=1; i<expiries->size(); ++i) {
            currentDate = bdc->adjust((*expiries)[i]->toDate(valueDate), hols);
            if (maturityDate < currentDate ) {
                /* do a linear interp. Not the most sophisticated way of calculatin rates, but sufficient
                   for out purposes. */
                currentSpread = spreads[i-1] +
                                double((maturityDate.getDate() - previousDate.getDate())) /
                                double((currentDate.getDate() - previousDate.getDate())) *
                                (spreads[i] - spreads[i-1]);
                break;
            }
            previousDate = currentDate;
        }
    }
    return currentSpread;
}

// Bucketwise shift - modifies the object
// Outside the usual sensitivity framework
void ParSpreadCurve::bucketShift(const DoubleArray& shifts)
{
    static const string method = "ParSpreadCurve::bucketShift";
    //validate that the shifts are the correct size
    if (shifts.size() != spreads.size())
    {
        throw ModelException(method,
                             "Incorrect shifts array size for " + getName());
    }

    //apply the shifts
    for (int i=0; i<shifts.size(); i++)
    {
        spreads[i] += shifts[i];
    }
}

/** Returns name identifying CDS Par Curve for PAR_SPREAD_RHO_POINTWISE */
string ParSpreadCurve::sensName(const ParSpreadPointwise*) const{
    return name;
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
TweakOutcome ParSpreadCurve::sensShift(const PropertyTweak<ParSpreadPointwise>& tweak) {
    static const string method = "ParSpreadCurve::sensShift";
    try{
        if (!Maths::isZero(tweak.coefficient)) {
            int i = tweak.qualifier->expiry->search(expiries.get());
            spreads[i] += tweak.coefficient;
        }

        return TweakOutcome(tweak.coefficient, false); // none of our components has a rho type sensitivity
    } catch (exception& e){
        throw ModelException(e, method,
                             "ParSpreadPointwise failed for " +
                             getName());
    }
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  par curve */
ExpiryWindowArrayConstSP ParSpreadCurve::sensQualifiers(const ParSpreadPointwise*) const {
    return ExpiryWindow::series(expiries);
}

//---------------  ParSpreadUpfronts tweak --------------------------------------------//
/** Returns name identifying CDS Par Curve for PAR_SPREAD_UPFRONT_POINTWISE */
string ParSpreadCurve::sensName(const ParSpreadUpfronts*) const
{
    return name;
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
TweakOutcome ParSpreadCurve::sensShift(const PropertyTweak<ParSpreadUpfronts>& tweak) 
{
    static const string method = "ParSpreadCurve::sensShift";
    
    try{
        if(!upfronts.get())
        {
            upfronts = DoubleArraySP(new DoubleArray(expiries->size()));
            for(int i = 0; i < upfronts->size(); i++)
            {
                (*upfronts)[i] = 0;
            }
        }
        
        if (!Maths::isZero(tweak.coefficient)) {
            int i = tweak.qualifier->expiry->search(expiries.get());
            (*upfronts)[i] += tweak.coefficient;
        }

        // none of our components has a rho type sensitivity
        return TweakOutcome(tweak.coefficient, false); 
    } catch (exception& e){
        throw ModelException(e, method,
                             "ParSpreadUpfronts failed for " +
                             getName());
    }
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  par curve */
ExpiryWindowArrayConstSP ParSpreadCurve::sensQualifiers(const ParSpreadUpfronts*) const 
{
    return ExpiryWindow::series(expiries);
}

//--------------- ParSpreadUpfrontParallel tweak --------------------------------------------//
/** Returns name identifying CDS Par Curve for PAR_SPREAD_UPFRONT_PARALLEL */
string ParSpreadCurve::sensName(const ParSpreadUpfrontParallelTP*) const{
    return name;
}

/** Shifts the object using given shift */
TweakOutcome ParSpreadCurve::sensShift(const PropertyTweak<ParSpreadUpfrontParallelTP>& tweak)
{
   static const string method = "ParSpreadCurve::sensShift";
   try {
       if(!upfronts.get())
       {
           upfronts = DoubleArraySP(new DoubleArray(expiries->size()));
           for(int i = 0; i < upfronts->size(); i++)
           {
               (*upfronts)[i] = 0;
           }
       }

       if (!Maths::isZero(tweak.coefficient)){
           for (int i = 0; i < upfronts->size(); i++) {
               (*upfronts)[i] += tweak.coefficient;
          }
       }

       return TweakOutcome(tweak.coefficient, false); // none of our components has a rho type sensitivity
   }
   catch (exception &e) {
      throw ModelException(e, method,
                           "ParSpreadUpfrontParallel tweak failed for " + getName());
   }
}


//--------------- ParSpreadParallel tweak --------------------------------------------//
/** Returns name identifying CDS Par Curve for PAR_SPREAD_RHO_PARALLEL */
string ParSpreadCurve::sensName(const ParSpreadParallel*) const{
    return name;
}

/** Shifts the object using given shift */
TweakOutcome ParSpreadCurve::sensShift(const PropertyTweak<ParSpreadParallel>& tweak){
   static const string method = "ParSpreadCurve::sensShift";
   try {
       if (!Maths::isZero(tweak.coefficient)){
          for (int i = 0; i < spreads.getLength(); i++) {
              spreads[i] += tweak.coefficient;
          }
          // zc.useParallelCurve(shiftSize);   // switch to parallel curve
       }

       return TweakOutcome(tweak.coefficient, false); // none of our components has a rho type sensitivity
   }
   catch (exception &e) {
      throw ModelException(e, method,
                           "ParSpreadParallel tweak failed for " + getName());
   }
}


/** Returns name identifying CDS Par Curve */
string ParSpreadCurve::sensName(const ParSpreadParallelRelative*) const{
    return name;
}

/** Shifts the object using given shift */
TweakOutcome ParSpreadCurve::sensShift(
     const PropertyTweak<ParSpreadParallelRelative>& tweak) 
{
   static const string method = "ParSpreadCurve::sensShift";
   try {
       if (!Maths::isZero(tweak.coefficient)){
          for (int i = 0; i < spreads.getLength(); i++) {
              spreads[i] *= (1 + tweak.coefficient);
          }
       }

       return TweakOutcome(tweak.coefficient, false); // none of our components has a rho type sensitivity
   }
   catch (exception &e) {
      throw ModelException(e, 
                           method,
                           "ParSpreadParallelRelative tweak failed for " + 
                           getName());
   }
}



// Returns name identifying CDS Par Curve for PAR_SPREAD_LEVEL
string ParSpreadCurve::sensName(ParSpreadLevel* shift) const {
    return name;
}

bool ParSpreadCurve::sensShift(ParSpreadLevel* shift) {
    static const string method = "ParSpreadCurve::sensShift";
    try {
        double spreadLevel = shift->getShiftSize();

        for (int j = 0; j < spreads.size(); j++) {
            spreads[j] = spreadLevel;
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method,
                             "ParSpreadLevel scenario failed for " + getName());
    }
    return false;  // all done
}


// Flat override support
bool ParSpreadCurve::sensShift(QuasiContractualBaseCorrelation* shift) {
    static const string method = "ParSpreadCurve::sensShift";
    try {
        double spreadLevel = shift->getCreditSpreadsLevel();

        for (int j=0; j < spreads.size(); ++j) {
            spreads[j] = spreadLevel;
        }
    }
    catch (exception& e ) {
        throw ModelException(e, 
                             method,
                             "QuasiContractualBaseCorrelation greek failed "
                             "for " + getName());
    }
    return false;  // all done
}

// ParSpreadLevelAtGivenDate::IShift implementation
string ParSpreadCurve::sensName(ParSpreadLevelAtGivenDate* shift) const {
    return name;
}

// ParSpreadLevelAtGivenDate::IShift implementation
bool ParSpreadCurve::sensShift(ParSpreadLevelAtGivenDate* shift) {
    static const string method = "ParSpreadLevelAtGivenDate::sensShift";
    try {
        double levelDate = shift->getLevelDate().getDate();

        DoubleArray expiriesAsDouble(expiries->size());
        for (int i = 0; i < expiries->size(); i++) {
            expiriesAsDouble[i] = (*expiries)[i]->toDate(today).getDate();
		}
        
        // interpolate the spreads linearly
        LinearInterpolator interpolator;
        Interpolator::InterpolantSP spreadsInterpolator = 
            Interpolator::InterpolantSP::constCast(
                interpolator.computeInterp(expiriesAsDouble, spreads));
        
        double spreadLevel = spreadsInterpolator->value(levelDate);

        for (int j = 0; j < spreads.size(); j++) {
            spreads[j] = spreadLevel;
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method,
                             "ParSpreadLevel scenario failed for " + getName());
    }
    return false;  // all done
}

// Returns name identifying CDS Par Curve for PAR_SPREAD_PARALLEL_SHIFT 
string ParSpreadCurve::sensName(ParSpreadParallelShift* shift) const{
    return name;
}

/** Shifts the object using given shift */
bool ParSpreadCurve::sensShift(ParSpreadParallelShift* shift){
   static const string method = "ParSpreadCurve::sensShift";
   try {
       double shiftSize = shift->getShiftSize();
       if (!Maths::isZero(shiftSize)){
          for (int i = 0; i < spreads.getLength(); i++) {
              spreads[i] += shiftSize;
          }
       }

       return false; // all done
   }
   catch (exception &e) {
       throw ModelException(e, method,
                           "ParSpreadParallelShift failed for " + getName());
   }
}

// Returns name identifying CDS Par Curve for PAR_SPREAD_PROP_SHIFT
string ParSpreadCurve::sensName(ParSpreadPropShift* shift) const {
    return name;
}

// Shifts the object using given shift (see ParSpreadPropShift::Shift)
bool ParSpreadCurve::sensShift(ParSpreadPropShift* shift) {
    static const string method = "ParSpreadCurve::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            for (int j = 0; j < spreads.size(); j++) {
                spreads[j] *= (1.0 + shiftSize);
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method,
                             "ParSpreadPropShift scenario failed for " +
                             getName());
    }
    return false;  // all done
}


// Returns name identifying CDS Par Curve for additive or multiplicative weighted shift
string ParSpreadCurve::sensName(ParSpreadWeightedShift* shift) const {
    return name;
}

// Shifts the object using given shift
bool ParSpreadCurve::sensShift(ParSpreadWeightedShift* shift) {
    shift->shiftArray(expiries, &spreads, today);
    return false; // nothing else to tweak
}


// Returns name identifying CDS Par Curve for theta shift
string ParSpreadCurve::sensName(Theta* shift) const {
    return name;
}

bool ParSpreadCurve::sensShift(Theta* shift) {
    today = shift->rollDate(today);
    return false; // nothing else to tweak
}

/** get spreads*/
DoubleArrayConstSP ParSpreadCurve::getParSpreads() const {
    //copy & wrap them in a SP
    DoubleArrayConstSP spds = DoubleArrayConstSP(new DoubleArray(spreads));
    return spds;
}

DoubleArrayConstSP ParSpreadCurve::getUpfronts() const {
    return upfronts;
}

ExpiryArrayConstSP ParSpreadCurve::getExpiries() const {
    return expiries;
}

void ParSpreadCurve::acceptWrapperNameCollector(
    const ParSpreadCurve* parSpreadCurve,
    WrapperNameCollector* collector)
{
    collector->addName(parSpreadCurve->getName());
}

void ParSpreadCurve::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ParSpreadCurve, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(ITweakableWithRespectTo<ParSpreadParallel>);
    IMPLEMENTS(ITweakableWithRespectTo<ParSpreadParallelRelative>);
    IMPLEMENTS(ITweakableWithRespectTo<ParSpreadPointwise>);
    IMPLEMENTS(ITweakableWithRespectTo<ParSpreadUpfronts>);
    IMPLEMENTS(ITweakableWithRespectTo<ParSpreadUpfrontParallelTP>);
    IMPLEMENTS(QuasiContractualBaseCorrelation::IShift);
    IMPLEMENTS(ParSpreadLevel::IShift);
    IMPLEMENTS(ParSpreadLevelAtGivenDate::IShift);
    IMPLEMENTS(ParSpreadParallelShift::IShift);
    IMPLEMENTS(ParSpreadPropShift::IShift);
    IMPLEMENTS(ParSpreadWeightedShift::IShift);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(Theta::IShift);
    EMPTY_SHELL_METHOD(defaultParSpreadCurve);
    ClassSetAcceptMethod(ParSpreadCurve::acceptWrapperNameCollector);
    FIELD(name,    "name of the CDS par spread curve");
    FIELD(expiries,       "expiry dates of the spreads");
    FIELD(spreads, "CDS par spreads");
    FIELD(today, "today");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(upfronts, "upfront fees in decimal assuming national is 1.0");
    FIELD_MAKE_OPTIONAL(upfronts);

    Addin::registerConstructor("PAR_SPREAD_CURVE",
        Addin::MARKET,
        "Creates a par CDS spread curve",
        TYPE);

    // Register "spreads" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "spreads",
        new Range(OpenBoundary(0.0),  OpenBoundary(10.0)));


}

CClassConstSP const ParSpreadCurve::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadCurve", typeid(ParSpreadCurve), load);

DEFINE_TEMPLATE_TYPE(ParSpreadCurveWrapper);

DRLIB_END_NAMESPACE
