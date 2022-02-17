//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIDynamicBasket.hpp
//
//   Description : Dynamic basket interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_DYNBASK_HPP
#define EDR_SPI_DYNBASK_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/SPIAlgorithm.hpp"
#include "edginc/SPIBond.hpp"
#include "edginc/SPICoupons.hpp"
#include "edginc/SPILoanCost.hpp"
#include "edginc/SPILockIn.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/DateBuilder.hpp"
#include "edginc/Holiday.hpp"

#include "edginc/Maths.hpp"
// HUGE_VAL ?
#define VERY_BIG       1e+15; 
#define SAMPLE_ONLY      0
#define REBALANCE        1
#define FEE_NOTIFICATION 2
#define FEE_PAYMENT      4
#define FINAL_REBAL      8
#define AVERAGE          16

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL SPIDynamicBasket : public CObject {
private:
    // Potential rebalance dates plus a few at the start & end for lag
    ISPIDateBuilderSP     rebalanceDates; //(either all the dates or a parametric form)
    DateTimeArraySP       theRebalanceDates; // transient - built from parametric form
    DateTimeArray         sparseRebalanceDates; // rebalanceDates less some which are skipped $unregistered
    int                   maxNumSkipDates; // performance enhancer, optional, default = 0
    bool                  skipping; // might have maxNumSkipDates and monthly rebal so we're not 
                                    // actually skipping. Should then avoid turning off lag and scaling gap risk
    // these define the interval when lag adj is turned off when skipping steps
    int                   iStepEarlyLagAdj; // when skipping no lag adj after this date
    int                   iStepLateLagAdj; // when skipping no lag adj before this date
    IntArray              daysBetweenSimDates;

protected:
    SPIDynamicBasket(CClassConstSP clazz);// for reflection

public:
    static CClassConstSP const TYPE;

    // fields

    // components of the dynamic basket
    IBondSPIInterfaceSP      bond;

    // associated fees
    IFeesSPIInterfaceSP      feesSPI;

    // Coupons may be paid out of the basket
    ICouponsSPIInterfaceSP   couponsSPI;

    // Crash/Bounds/Rebalance decision
    IAlgorithmSPIInterfaceSP algorithm;

    // Accrual of loan cost. 
    ILoanCostSPIInterfaceSP  loanCost;

    // Lag between trade and knowing level at which traded
    IntArray              numPubLagDays;     // [numEquities], publication lag
    IntArray              numExecLagDays;    // [numEquities], execution lag
    bool                  overrideInitialLag; // we can override the lag for day one only
    int                   maxExecLagDays;    // useful
    int                   maxPubLagDays;     // useful

    // costs incurred when doing a full rebalance - proportional to amount of equity bought/sold
    // deducted from basket
    bool                  hasRebalCost;
    double                pastRebalCostRate;
    double                futureRebalCostRate;
    bool                  overrideInitialRebalCost;
    bool                  overrideFinalRebalCost;

    double                initialBasketLevel;

    // the following became necessary when we added coupons and fees on notional
    bool                  hasBFHistory; // do we have history of our bond floor
    CashFlowArraySP       bondFloorHistory;
    int                   lastHistoricBFIndex;  
    DoubleArray           cachedDFQuot; // for BF calcs when rolling
    mutable double        bondFloorToday; // ditto  

    // validation
    void validatePop2Object();

    void validateBondFloor(const DateTime today);

    void validateDates(const Holiday*       ccyHols, 
                       const IMultiFactors* assets,
                       const ObservationSourceArray& sources,
                       bool           excludeAssetHols);

    // For performance try to remove potential rebalance dates so
    // long as any gap left is <= maxNumSkipDates. 
    // Constrain to have some unskipped at the start and at the end.
    // Any past dates are of course unskipped.
    const DateTimeArray& prepareCoarseRebalDates(const DateTimeArray& someEssentialDates,
                                                 const DateTime&      today);

    const DateTimeArray& getRebalanceDates() const;

    const DateTimeArray& getAllRebalanceDates() const;

    // Dealing with actual - not sparse - rebalance dates 
    const DateTime& getNextRebalanceDate(const DateTime& from) const;

    const DateTime& getLastRebalDate() const;

    // Date on which basket is finally known
    const DateTime& getFinalDate() const;

    bool doLagAdj(int iStep) const;

    bool doRebalCost(int iStep) const;

    int getPubLag(int iAsset,
                  int iStep);

    void scaleGapRisk(int     iStep,
                      double& gaprisk) const;

    // if we've rolled forward we've stored DF quotients between all rebalance dates
    // now recast these as a set of DFQuotients between the current sparse rebal dates
    const DoubleArray unpickDiscountFactorQuotients() const;

    const DoubleArray getDiscountFactorsQuotients(bool sparse, const YieldCurve* disc) const;

    // Update contents for a time shift
    void roll(const Theta::Util&   thetaUtil,
              const YieldCurve*    disc,
              double               dayCountBasis);

    // note there's also a protected default constructor for use with the derived class
    SPIDynamicBasket(); // for reflection

private:
    SPIDynamicBasket(const SPIDynamicBasket& rhs); // not implemented
    SPIDynamicBasket& operator=(const SPIDynamicBasket& rhs); // not implemented

    static IObject* defaultSPIDynamicBasket();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIDynamicBasket> SPIDynamicBasketSP;

/****************************************************************************/

// empty derived class which is just an SPIDynamicBasket
// this enables us to split the interface at the IMS level
// into SPI and Super SPI. We restrict access to the full set of
// params for the former at the Aladdin/IMS level
class PRODUCTS_DLL SPIDynamicBasketSuper : public SPIDynamicBasket {
public:
    static CClassConstSP const TYPE;
    SPIDynamicBasketSuper();// for reflection

private:
    SPIDynamicBasketSuper(const SPIDynamicBasketSuper& rhs); // not implemented
    SPIDynamicBasketSuper& operator=(const SPIDynamicBasketSuper& rhs); // not implemented

    static IObject* defaultSPIDynamicBasketSuper();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
