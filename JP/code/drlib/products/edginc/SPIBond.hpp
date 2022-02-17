//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIBond.hpp
//
//   Description : Bond interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_BOND_HPP
#define EDR_SPI_BOND_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Theta.hpp"
#include "edginc/SPIUtil.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

// this is the actual interface of what a bond object 
// needs to do internally
class PRODUCTS_DLL IBondSPI {
public:
    virtual void init(const DateTime&      today,
                      const DateTime&      maturity,
                      const DateTimeArray* simDates) = 0;

    virtual void getYieldCurveData(const IModel*          model, 
                                   const CMarketDataSP    market) = 0;

    virtual void getBondPrices(const YieldCurve*    disc,
                               DoubleArray&         bondPrices,
                               const DateTimeArray* forDates = 0) = 0;

    virtual void getFutureBondPrices(vector<SVExpectedDiscFactorSP>& dfSVBond,
                                     int                             iFirstFutureStep,
                                     DoubleArray&                    bondPrices) = 0;

    virtual double getBondPriceToday(const DateTime& today) = 0;

    virtual DateTimeArray getEssentialDates() const = 0;

    // This may return 0 which means it is not (piecewise)linear
    // Otherwise the Schedule returned IS the bond floor level
    virtual const Schedule* getLinearBondFloor() const = 0;

    virtual void roll(const Theta::Util&   thetaUtil,
                            const YieldCurve*    disc) = 0;

    // Slight design flaw due to incomplete spec - the yield curve used for
    // closing calculation is also needed for Loan Cost calculation. Rather than
    // take another Yield Curve into LoanCost classes I decided to fudge by
    // taking it from the bond. It's also now used for bond floor today for SPIs
    // which need a BF history. This design is getting very creaky. Maybe need 
    // a fixing curve object held high up (preferbaly with a flag to kill it when rolling)
    virtual const YieldCurve* getYC() = 0;

    virtual void setYC(YieldCurveWrapper& newYC) = 0;

    virtual void setYC(YieldCurveSP newYC) = 0;

    virtual ~IBondSPI();
};
typedef refCountPtr<IBondSPI> IBondSPISP;

// this is the external interface for abstraction so that users can bolt in 
// any Bond type they want - note this includes the SPIBondWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface IBondSPI as soon as possible
class PRODUCTS_DLL IBondSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    virtual IBondSPISP getBondSPI() = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IBondSPIInterface> IBondSPIInterfaceSP;

// Now we have concrete implementations of the IBondSPI interface
class PRODUCTS_DLL SPIBondStd : public SPIInterfaceIMS,
                    virtual public IBondSPI,
                    virtual public IBondSPIInterface {
private:

    // fields
    YieldCurveWrapper    yc;        // PURELY USED to compute bond price today
    string               stubTypeString;
    DateTime             datedDate; // accrue start date for first coupon
    CashFlowArraySP      cashFlows;  // coupon amounts
    double               redemptionAmt;
    string               dayCountConvString;
    CashFlowArraySP      pastBondValues;
    double               fundingSpread; // Sz
    bool                 hasLinearBondFloor; // optional - allows override of bond floor 
    ScheduleSP           linearBondFloor;    // optional - to be piecewise linear

    // transient
    DoubleArray          bondPrices; // 
    DateTime             today;
    DateTime             maturityDate; // needs info about lags so assigned late
    const DateTimeArray* simDates;
    vector<int>          bondDatesIdx; // bondDates indexed within simDates

    // derived
    DayCountConventionSP accruedDCC;
    int                  stubType;

    static const int stubNone;
    static const int stubSwap;
    static const int stubBond;

public: 
    // get the internal interface
    virtual IBondSPISP getBondSPI();

    // methods
    virtual void validatePop2Object();

    virtual DateTimeArray getEssentialDates() const;

    // This may return 0 which means it is not (piecewise)linear
    // Otherwise the Schedule returned IS the bond floor level
    virtual const Schedule* getLinearBondFloor() const;

    virtual void init(const DateTime& today,
                      const DateTime& maturity,
                      const DateTimeArray* simDates);

    void getYieldCurveData(const IModel*          model, 
                           const CMarketDataSP    market);

    const YieldCurve* getYC();

    // used when rolling to turn off fixing curve
    void setYC(YieldCurveWrapper& newYC);

    // used when retrospectively building bond floor history
    void setYC(YieldCurveSP newYC);

    void getBondPrices(const YieldCurve*    disc,
                       DoubleArray&         bondPrices,
                       const DateTimeArray* forDates = 0);

    void getFutureBondPrices(vector<SVExpectedDiscFactorSP>& dfSVBond,
                             int                             iFirstFutureStep,
                             DoubleArray&                    bondPrices);

    void roll(const Theta::Util&   thetaUtil,
              const YieldCurve*    disc);

    // This uses the special yield curve specified in the SPIBond
    double getBondPriceToday(const DateTime& today);

    static CClassConstSP const TYPE;

    SPIBondStd(); // for reflection

private:
    SPIBondStd(const SPIBondStd& rhs); // not implemented
    SPIBondStd& operator=(const SPIBondStd& rhs); // not implemented

    static IObject* defaultSPIBondStd();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIBondStd> SPIBondStdSP;

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPIBondWrapper : public CObject,
                                    virtual public IBondSPIInterface {
public:
    static CClassConstSP const TYPE;

    string          SPIBondType;
    SPIBondStdSP    bondStd;

public:

    virtual IBondSPISP getBondSPI();

    // validation
    void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
     
    // for reflection
    SPIBondWrapper();

    static IObject* defaultSPIBondWrapper();
};
typedef smartPtr<SPIBondWrapper> SPIBondWrapperSP;

DRLIB_END_NAMESPACE

#endif
