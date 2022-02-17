//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : NumericalIntegrationLN.cpp
//
//   Description : Numerical Integration Algorithm 
//
//   Author      : Andrew J Swain
//
//   Date        : 5 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/NumericalIntegrationLN.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/ImpliedSample.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"


DRLIB_BEGIN_NAMESPACE

static const int HI_MAX_STD_DEV = 6;
static const int LO_MAX_STD_DEV = 12;  // tend to need more downside than upside
static const int FINE_STD_DEV   = 20;

void NumericalIntegrationLN::validatePop2Object() {
    static const string method = "NumericalIntegrationLN::validatePop2Object";
    try {
        if (steps < 3) {
            throw ModelException(method, "need at least 3 steps- given ("+
                                 Format::toString(steps) + ")");
        }

        if (steps%2 == 0) {
            throw ModelException(method, 
                                 "require odd number of steps - given ("+
                                 Format::toString(steps) + ")");
        }
   
        // cutoffs must be 0 <= lo < hi <= 1
        if (Maths::isNegative(loCutOff) || hiCutOff > 1.0 ||
            loCutOff >= hiCutOff) {
            throw ModelException(method, 
                                 "low cut off (" + Format::toString(loCutOff)+ 
                                 ") and high cut off (" + 
                                 Format::toString(hiCutOff) + ") must be in "
                                 "range  0 <= lo < hi <= 1");
        }          
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Constructor takes type of vol to use */
NumericalIntegrationLN::NumericalIntegrationLN(const string& volType):
    CModelLN(TYPE, volType), hiCutOff(0.0), loCutOff(0.0) {}


/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP NumericalIntegrationLN::createMDF() const{
    return MarketDataFetcherSP(new MarketDataFetcherLNSpline(getVolType()));
}


/** calculate single price and store result in results */
void NumericalIntegrationLN::Price(CInstrument*  instrument, 
                                   CControl*     control, 
                                   CResults*     results){
    static const string method = "NumericalIntegrationLN::Price";
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support NumericalIntegrationLN::IntoProduct");
    }
    IProduct* product = 0;
    try {
        if (instrument->priceDeadInstrument(control, results)) {
            return; // done for a dead instrument
        }
        
        // cast to NumericalIntegrationLN::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        product = intoProd.createProduct(this);

        product->price(this, control, results);
    } 
    catch (exception& e) {
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}

/** integrate payoff using Simpson's rule */
double NumericalIntegrationLN::integrate(
    NumericalIntegrationLN::IProduct* product) const {
    static const string method = "NumericalIntegrationLN::integrate";
    try {
        PDFCalculatorSP pdf(product->pdfCalculator());

        // get range of integration
        double          loBound = lowBound(product);
        double          hiBound = highBound(product);

        double          stepSize = (hiBound-loBound)/steps;
        double          level = loBound;

        DoubleArray     levels(steps);
        DoubleArray     density(steps);
        DoubleArray     payoff(steps);

        int i;
        for (i = 0; i < steps; i++) {
            levels[i] = level;
            level += stepSize;
        }

        // get distribution
        DateTime time = product->time();

        try {
            pdf->localDensity(levels, time, density);
        }
        catch (PDFCalculatorException& e){
            throw PDFCalculatorException(e, method);
        }

        // and evaluate
        double integral = 0.0;
        for (i = 0; i < steps; i++) {
            payoff[i] = density[i] * product->payoff(levels[i]);
        }

        // use Simpsons rule to get value
        for (i = 0; i < steps - 2; i+=2) {
            integral += payoff[i] + 4.0* payoff[i+1] + payoff[i+2];
        }

        integral *= stepSize/3.0;
        return integral;
    } 
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



double NumericalIntegrationLN::highBound(IProduct* product) const {
    static const string method = "NumericalIntegrationLN::highBound";
    try {
        DateTime time   = product->time();
        double   fwd    = product->centre(time); 
        double   var    = product->variance(fwd, time);
        double   target = fabs(N1Inverse(hiCutOff));

        PDFCalculatorSP pdf(product->pdfCalculator());
        DateTimeArray   datesTo(1, time);
        DoubleArray     sqrtVar(1, sqrt(var));
        DoubleArray     fwds(1, fwd);

        // cheat Bracket - not interested in fwd starting
        // so don't need today or spot - also only works for +ve std dev
        Brackets        bracket(target,
                                pdf.get(),
                                time,
                                datesTo,
                                fwds,
                                fwd,
                                sqrtVar,
                                HI_MAX_STD_DEV,
                                FINE_STD_DEV,
                                false);

        double hi = bracket.upStrikes[0];

        return hi;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double NumericalIntegrationLN::lowBound(IProduct* product) const {
    static const string method = "NumericalIntegrationLN::lowBound";
    try {
        DateTime time   = product->time();
        double   fwd    = product->centre(time);
        double   var    = product->variance(fwd, time);
        double   target = fabs(N1Inverse(loCutOff));

        PDFCalculatorSP pdf(product->pdfCalculator());
        DateTimeArray   datesTo(1, time);
        DoubleArray     sqrtVar(1, sqrt(var));
        DoubleArray     fwds(1, fwd);

        // cheat Bracket - not interested in fwd starting
        // so don't need today or spot - also only works for +ve std dev
        Brackets        bracket(target,
                                pdf.get(),
                                time,
                                datesTo,
                                fwds,
                                fwd,
                                sqrtVar,
                                LO_MAX_STD_DEV,
                                FINE_STD_DEV,
                                false);

        double lo = bracket.downStrikes[0];

        return lo;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void NumericalIntegrationLN::sensRestore(Theta* shift) {
    isTimeShift = false;
}

/** Does special things for Theta-type tweaks */
bool NumericalIntegrationLN::sensShift(Theta* shift) {

    DateTime dummyDate(0,0);
    DateTime newDate = shift->rollDate(dummyDate);

    if (!newDate.equals(dummyDate)){
        // i.e. this is really is theta or similar, not the roll to now
        isTimeShift = true;
    }
    return true;
}

bool NumericalIntegrationLN::doingTimeShift() const {
    return isTimeShift;
}

NumericalIntegrationLN::NumericalIntegrationLN():CModelLN(TYPE), isTimeShift(false) {};

class NumericalIntegrationLNHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(NumericalIntegrationLN, clazz);
        SUPERCLASS(CModelLN);
        IMPLEMENTS(Theta::RestorableShift);
        EMPTY_SHELL_METHOD(defaultNumericalIntegrationLN);
        FIELD(steps, "steps");
        FIELD(hiCutOff, "hiCutOff (probability)");
        FIELD(loCutOff, "loCutOff (probability)");
        FIELD(isTimeShift, "");
        FIELD_MAKE_TRANSIENT(isTimeShift);
        
    }

    // for NumericalIntegrationLN::IIntoProduct
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(NumericalIntegrationLN::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultNumericalIntegrationLN(){
        return new NumericalIntegrationLN();
    }
};

CClassConstSP const NumericalIntegrationLN::TYPE = 
					CClass::registerClassLoadMethod("NumericalIntegrationLN", 
					typeid(NumericalIntegrationLN), 
					NumericalIntegrationLNHelper::load);
bool  NumericalIntegrationLNLoad() {
    return (NumericalIntegrationLN::TYPE != 0);
}



CClassConstSP const NumericalIntegrationLN::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("NumericalIntegrationLN::IIntoProduct",
                                    typeid(NumericalIntegrationLN::IIntoProduct), 
                                    NumericalIntegrationLNHelper::loadIntoProduct);

DRLIB_END_NAMESPACE
