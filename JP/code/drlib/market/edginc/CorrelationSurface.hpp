//----------------------------------------------------------------------------
//
//   Group       : IR Derivatives Research
//
//   Filename    : CorrelationSurface.hpp 
//
//   Description : CorrelationSurface - term structure vs
//                 strike to represent correlation smile surface
//
//   Author      : Steve Marks
//
//   Date        : Nov 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CORRELATION_SURFACE_HPP
#define EDR_CORRELATION_SURFACE_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/CorrelationBase.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/VolSurface.hpp"

DRLIB_BEGIN_NAMESPACE

class CorrelationSurface;
DECLARE(CorrelationSurface);

class MARKET_DLL CorrelationSurface : public CorrelationBase
{
public:
    static CClassConstSP const TYPE;

    virtual ~CorrelationSurface() {};

    // CorrelationBase Interface 
    // ??? not sure what to do with this
    virtual void configureForSensitivities(CClassConstSP clazz1,
                                           CClassConstSP clazz2) {
        throw ModelException(__FUNCTION__, " function not implemented");
    }

    /** Returns true if this correlation is [really] sensitive to the
        supplied sensitivity */
    virtual bool isSensitiveTo(const IPerNameSensitivity* sens) const {
        throw ModelException(__FUNCTION__, " function not implemented");
    }

    // MarketObject interface
    virtual void getMarket(const IModel* model, const MarketData* market);
    virtual string getName() const { return name; }

    // local functions
    double getCorrelation(DateTime expiry, double strike) const;

protected:
    CorrelationSurface(CClassConstSP clazz) : CorrelationBase(clazz) {};
private:
    static IObject* defaultConstructor(void) { return new CorrelationSurface(TYPE); }
    static void load(CClassSP& clazz);

    // exported variables
    string             name;
    DoubleArray        strikes;
    CDoubleMatrixSP    correlation;
    ExpiryArraySP      expiries;
    DateTime           baseDate;

    string             asset1;
    string             asset2;
    
    // reuse existing/generic surface structure for interpolation etc.

};

DRLIB_END_NAMESPACE
#endif
