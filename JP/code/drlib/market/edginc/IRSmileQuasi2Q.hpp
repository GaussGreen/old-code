//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRSmileQuasi2Q.hpp
//
//   Description : SmileQuasi2Q - as Smile2Q but Q parameters can be surfaces.
//
//   Author      : Anwar E Sidat
//
//   Date        : 25-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IRSmileQuasi2Q_HPP
#define QLIB_IRSmileQuasi2Q_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRExoticParam.hpp"
#include "edginc/DoubleMatrix.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MarketData);
FORWARD_DECLARE(Surface);

/** Quasi2Q smile class.
 */
class MARKET_DLL IRSmileQuasi2Q : public IRExoticParam
{
public:
    static CClassConstSP const TYPE;

    IRSmileQuasi2Q();
    virtual ~IRSmileQuasi2Q();

    /** Get methods for input parameters. */
    virtual DoubleMatrix qLeft() const    { return QLeft; }
    virtual DoubleMatrix qRight() const   { return QRight; }
    virtual DoubleMatrix fwdShift() const { return FwdShift; }

    /** Returns name of model. */
    virtual string getName() const;

    /** overrides default */
    virtual void validatePop2Object();

    /** Initialise data. */
    void initialiseData(const MarketData* pMarket, const IModel* pModel = 0);

    /** Get expiries */
    const ExpiryArray& getExpiries() const { return Expiries; }

    /** Get tenors */
    const ExpiryArray& getTenors() const { return Tenors; }

    /** Get value off surface. */
    double getQLeft(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getQRight(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getFwdShift(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;

protected:

    IRSmileQuasi2Q(const CClassConstSP& clazz);
    IRSmileQuasi2Q(const IRSmileQuasi2Q& irv);
    IRSmileQuasi2Q& operator=(const IRSmileQuasi2Q& irv);

    //Fields
    string        Key;        // Handle name
    ExpiryArray   Expiries;
    ExpiryArray   Tenors;
    DoubleMatrix  QLeft;
    DoubleMatrix  QRight;
    DoubleMatrix  FwdShift;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IRSmileQuasi2Q(); }

    // Get value of smile parameter given expiry/maturity
    double getValue(
        const SurfaceSP&    surface,
        const DateTime&     baseDate,
        const DateTime&     expiryDate,
        const Expiry*       pMaturityTenor,
        const double        scalarValue)
        const;

    // Private members
    bool      bInitialised;
    SurfaceSP ptrSurfaceQLeft;
    SurfaceSP ptrSurfaceQRight;
    SurfaceSP ptrSurfaceFwdShift;
};

typedef smartConstPtr<IRSmileQuasi2Q> IRSmileQuasi2QConstSP;
typedef smartPtr<IRSmileQuasi2Q>      IRSmileQuasi2QSP;
typedef MarketWrapper<IRSmileQuasi2Q> IRSmileQuasi2QWrapper;

#ifndef QLIB_IRSmileQuasi2Q_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRSmileQuasi2Q>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRSmileQuasi2Q>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRSmileQuasi2Q>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRSmileQuasi2Q>);
#endif

// support for wrapper class
#ifndef QLIB_IRSmileQuasi2Q_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRSmileQuasi2Q>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRSmileQuasi2Q>);
#endif

DRLIB_END_NAMESPACE

#endif
