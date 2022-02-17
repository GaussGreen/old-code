//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LiquiditySpreadCurve.cpp
//
//   Description : a liquidity spread curve
//
//   Author      : André Segger
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_LIQUIDITYSPREADCURVE_CPP
#include "edginc/LiquiditySpreadCurve.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

void LiquiditySpreadCurve::validate()
{
    static const string method = "CDSParSpreads::validate";
    try
    {
        /*
        ExpiryArrayConstSP expiries = parSpreads->getExpiries();

        // validate that all of the maturities are on cycle with the swap frequency.
        // Use a not-so-arbitrary date as a base date to test from. The on-cycle
        // routine is pretty crappy and automatically fails for dates after the
        // 28th of the month.
        DateTime pseudoValueDate("01-JAN-2001", DateTime::START_OF_DAY);
        DateTimeArray tmpSwapDates;        
        int i;
        for (i=0; i< expiries->size(); i++) {
            
            tmpSwapDates = SwapTool::simpleSwapDates(
                pseudoValueDate,
                (*expiries)[i]->toDate(pseudoValueDate),
                parSwapFreq);
            
            if (SwapTool::onCycle(
                pseudoValueDate,     
                tmpSwapDates[tmpSwapDates.size()-1],   
                (tmpSwapDates.size()-1)*(4/parSwapFreq),
                "Q") == false)
            {
                throw ModelException(
                    method,
                    "Par swap maturities must be on cycle with the par swap frequency");
            }
        } 
        */
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

LiquiditySpreadCurve::LiquiditySpreadCurve() : MarketObject(TYPE) {};

LiquiditySpreadCurve::LiquiditySpreadCurve(const CClassConstSP& clazz): MarketObject(clazz) {
}

LiquiditySpreadCurve::LiquiditySpreadCurve(const string& credName,
                                           const ExpiryArray* expiries,
                                           const DoubleArray& rates): MarketObject(TYPE),
                                                                      name(credName),
                                                                      expiries(copy(expiries)),
                                                                      liquiditySpreads(rates)
{
    // empty
}

LiquiditySpreadCurve::LiquiditySpreadCurve(const string&        name,
                                             const ExpiryArray* expiries,
                                             const DoubleArray& rates,
                                             const int          parSwapFreq,
                                             const string&      parDCC,
                                             const double       parRecovery,
                                             const bool         parAccrueFee):MarketObject(TYPE),
                                                                              name(name),
                                                                              expiries(copy(expiries)),
                                                                              liquiditySpreads(rates),
                                                                              parSwapFreq(parSwapFreq),
                                                                              parDCC(parDCC),
                                                                              parRecovery(parRecovery),
                                                                              parAccrueFee(parAccrueFee)
{
    // empty
}



static IObject* defaultCreditSpreadCurve(){
    return new LiquiditySpreadCurve();
}

/** Returns name identifying yield curve for rho parallel */
string LiquiditySpreadCurve::sensName(LiquiditySpreadRhoParallel* shift) const{
    return getName();  // or is it getCcy() ? who knows ?
}

/** Shifts the object using given shift */
bool LiquiditySpreadCurve::sensShift(LiquiditySpreadRhoParallel* shift){
   static const string method = "LiquiditySpreadCurve::sensShift";
   try {
       double shiftSize = shift->getShiftSize();
       if (!Maths::isZero(shiftSize)){
          for (int i = 0; i < liquiditySpreads.getLength(); i++) {
             liquiditySpreads[i] += shiftSize;
          }
          // zc.useParallelCurve(shiftSize);   // switch to parallel curve
       }

       return false; // none of our components has a rho type sensitivity
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

/** Restores the object to its original form */
void LiquiditySpreadCurve::sensRestore(LiquiditySpreadRhoParallel* shift){
   double shiftSize = shift->getShiftSize();
   if (!Maths::isZero(shiftSize)){
      for (int i = 0; i < liquiditySpreads.getLength(); i++) {
         liquiditySpreads[i] -= shiftSize;
      }
      // zc.restore();
   }
}

/** Returns the name of the yield curve - used to determine whether 
    to tweak the object */
string LiquiditySpreadCurve::sensName(const LiquiditySpreadPointwise*) const {
    return getName();  
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryWindowArrayConstSP LiquiditySpreadCurve::sensQualifiers(const LiquiditySpreadPointwise*) const {
    return ExpiryWindow::series(expiries);
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
TweakOutcome LiquiditySpreadCurve::sensShift(const PropertyTweak<LiquiditySpreadPointwise>& tweak) {
    static const string method = "";
    try {
        if (!Maths::isZero(tweak.coefficient)) {
            int i = tweak.qualifier->expiry->search(expiries.get());
            liquiditySpreads[i] += tweak.coefficient;
        }

        return TweakOutcome(tweak.coefficient, false);
    }
    catch (exception& e) {
        throw ModelException(e, "LiquiditySpreadCurve::sensShift(LiquiditySpreadPointwise)");
    }
}

/** Restores the object to its original form */
void LiquiditySpreadCurve::sensRestore(const PropertyTweak<LiquiditySpreadPointwise>& tweak) {
    if (!Maths::isZero(tweak.coefficient)){
        int i = tweak.qualifier->expiry->search(expiries.get());
        liquiditySpreads[i] -= tweak.coefficient;
    }
}

ExpiryArrayConstSP LiquiditySpreadCurve::getExpiries()const
{
    return expiries;
}

DateTimeArray LiquiditySpreadCurve::getExpiryDates(const DateTime& today)const
{
    static const string method = "LiquiditySpreadCurve::getExpiryDates";
    try
    {
        // first turn expiries into absolute dates
        int numExpiries = expiries->size();
        DateTimeArray benchmarks(numExpiries);

        for (int i = 0; i < numExpiries; i++) {
            benchmarks[i] = (*expiries)[i]->toDate(today);
        }
        
        return benchmarks;
    }
    catch (exception &e) 
    {
        throw ModelException(e, method, 
                             "Failed to convert offsets to dates");
    }
}


void LiquiditySpreadCurve::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(LiquiditySpreadCurve, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultCreditSpreadCurve);
    FIELD(name, "name of the liquidity spread curve");
    FIELD(expiries, "liquidity spread curve dates");
    FIELD(liquiditySpreads, "liquidity spreads");
    FIELD(parSwapFreq       ,"par swap frequency");
    FIELD(parDCC            ,"day count convention");
    FIELD(parRecovery       ,"recovery");
    FIELD(parAccrueFee      ,"pay accrued on default");
}

CClassConstSP const LiquiditySpreadCurve::TYPE = CClass::registerClassLoadMethod(
    "LiquiditySpreadCurve", typeid(LiquiditySpreadCurve), load);

DEFINE_TEMPLATE_TYPE(LiquiditySpreadCurveWrapper);
   
DRLIB_END_NAMESPACE
