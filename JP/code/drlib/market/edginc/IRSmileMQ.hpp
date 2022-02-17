//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRSmileMQ.hpp
//
//   Description : SmileMQ - MultiQ Smile.
//
//   Author      : Anwar E Sidat
//
//   Date        : 05-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IRSmileMQ_HPP
#define QLIB_IRSmileMQ_HPP

#include "edginc/IRExoticParam.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Surface.hpp"

DRLIB_BEGIN_NAMESPACE

// To avoid redundant file includes.
FORWARD_DECLARE(Expiry);
FORWARD_DECLARE(Surface);

/** MQ smile class.
 */
class MARKET_DLL IRSmileMQ : public IRExoticParam
{
public:
    static CClassConstSP const TYPE;

    IRSmileMQ();
    virtual ~IRSmileMQ();

    /** Returns name of model. */
    virtual string getName() const;

    /** Overrides default */
    virtual void validatePop2Object();

    /** Initialise data. */
    void initialiseData(const MarketData* pMarket, const IModel* pModel = 0);

    /** Get expiries */
    const ExpiryArray& getExpiries() const { return Expiries; }

    /** Get tenors */
    const ExpiryArray& getTenors() const { return Tenors; }

    /** Get value off surface. */
    double getSkew(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getVolOfVol(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getBBRP(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getBBVP(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getDeltaLeft(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getDeltaRight(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getTauLeft(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getTauRight(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getLiquidityFlag(const DateTime& baseDate, const DateTime& expiryDate, const Expiry* pMaturityTenor) const;
    double getNormalCutoff() const { return NormalCutoff; }
    double getNCK() const { return NCK; }

    /** Get surfaces (these access functions may be changed to use enums later). */
    SurfaceSP getSurfaceSkew() const { return ptrSurfaceSkew; }
    SurfaceSP getSurfaceVolOfVol() const { return ptrSurfaceVolOfVol; }
    SurfaceSP getSurfaceBBRP() const { return ptrSurfaceBBRP; }
    SurfaceSP getSurfaceBBVP() const { return ptrSurfaceBBVP; }
    SurfaceSP getSurfaceDeltaLeft() const { return ptrSurfaceDeltaLeft; }
    SurfaceSP getSurfaceTauLeft() const { return ptrSurfaceTauLeft; }
    SurfaceSP getSurfaceDeltaRight() const { return ptrSurfaceDeltaRight; }
    SurfaceSP getSurfaceTauRight() const { return ptrSurfaceTauRight; }
    SurfaceSP getSurfaceLiquidityFlag() const { return ptrSurfaceLiquidityFlag; }

protected:

    IRSmileMQ(const CClassConstSP& clazz);
    IRSmileMQ(const IRSmileMQ& irv);
    IRSmileMQ& operator=(const IRSmileMQ& irv);

    //Fields
    string         Key;        // Handle name
    ExpiryArray    Expiries;
    ExpiryArray    Tenors;
    DoubleMatrix   Skew;
    DoubleMatrix   VolOfVol;
    DoubleMatrix   BBRP;
    DoubleMatrix   BBVP;
    DoubleMatrix   DeltaLeft;
    DoubleMatrix   TauLeft;
    DoubleMatrix   DeltaRight;
    DoubleMatrix   TauRight;
    DoubleMatrix   LiquidityFlag;
    double         NormalCutoff;
    double         NCK;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IRSmileMQ(); }
    
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
    SurfaceSP ptrSurfaceSkew;
    SurfaceSP ptrSurfaceVolOfVol;
    SurfaceSP ptrSurfaceBBRP;
    SurfaceSP ptrSurfaceBBVP;
    SurfaceSP ptrSurfaceDeltaLeft;
    SurfaceSP ptrSurfaceDeltaRight;
    SurfaceSP ptrSurfaceTauLeft;
    SurfaceSP ptrSurfaceTauRight;
    SurfaceSP ptrSurfaceLiquidityFlag;
};

typedef smartConstPtr<IRSmileMQ> IRSmileMQConstSP;
typedef smartPtr<IRSmileMQ>      IRSmileMQSP;
typedef MarketWrapper<IRSmileMQ> IRSmileMQWrapper;

#ifndef QLIB_IRSmileMQ_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRSmileMQ>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRSmileMQ>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRSmileMQ>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRSmileMQ>);
#endif

// support for wrapper class
#ifndef QLIB_IRSmileMQ_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRSmileMQ>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRSmileMQ>);
#endif

DRLIB_END_NAMESPACE

#endif
