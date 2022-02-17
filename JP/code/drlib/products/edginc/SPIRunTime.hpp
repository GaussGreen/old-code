//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIRunTime.hpp
//
//   Description : Run time classes for SPI and Rainbow SPI
//                 Handles the algorithm in a central place
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_RT_HPP
#define EDR_SPI_RT_HPP

#include "edginc/config.hpp"
#include "edginc/SPIDynamicBasket.hpp"
#include "edginc/SPIBondFloor.hpp"
#include "edginc/SPILockIn.hpp"
#include "edginc/SPICoupons.hpp"
#include "edginc/SPIFees.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/StateVariableCollector.hpp"
#include "edginc/StateVariableGen.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

// fwd declaration
class SPIRunTime;
typedef refCountPtr<SPIRunTime> SPIRunTimeSP;

class SyntheticPortfolioInsurance;

// to separate the contingency from the the coupon itself
class PRODUCTS_DLL SPICouponsRT {
public:
    ICouponsSPISP   coupons;

    SPICouponsRT(ICouponsSPISP coupons);

    virtual double getCoupon(SPIRunTime* dynBask, 
							 double BF, int iStep);

    virtual double getCumSE();

    virtual double getCumTE();

    virtual double getCumB();

    virtual ~SPICouponsRT();

protected:
    // store the cum basket values    
    double B;
};
typedef refCountPtr<SPICouponsRT> SPICouponsRTSP;

// coupons contingent on basket level
class PRODUCTS_DLL SPICouponsBasketRT : public SPICouponsRT {
public:
    SPICouponsBasketRT(ICouponsSPISP coupons);

    virtual double getCoupon(SPIRunTime* dynBask, 
                             double BF, int iStep);

    virtual ~SPICouponsBasketRT();
};

// coupons contingent on TE level
class PRODUCTS_DLL SPICouponsExposureRT : public SPICouponsRT {
public:
    SPICouponsExposureRT(ICouponsSPISP coupons);

    virtual double getCoupon(SPIRunTime* dynBask, 
                             double BF, int iStep);

    virtual double getCumSE();

    virtual double getCumTE();

    virtual ~SPICouponsExposureRT();

private:
    // store the cum exposure values    
    double SE;
    double TE;
};

/*****************************************************************************/
// to separate the contingency from the the fee itself
class PRODUCTS_DLL SPIFeesRT {
public:

    SPIFeesRT(SPIRunTime* dynBasket);

    virtual ~SPIFeesRT();

    virtual double calculateFees(int iStep, double& feeToPay, double pBF, double pBL);

protected:
    virtual double getCtgtFee(int iStep, double pBF, double pBL);

    virtual double getCtgtPaidFee(int iStep);

    IFeesSPISP   fees;
    SPIRunTime* dynBask;
    int firstStep;
    int lastStep;
};
typedef refCountPtr<SPIFeesRT> SPIFeesRTSP;

// fees contingent on SE level
class PRODUCTS_DLL SPIFeesExposureRT : public SPIFeesRT {
public:

    SPIFeesExposureRT(SPIRunTime* dynBasket,
                        ISPIBondFloorSP bondFloor, ILockInSPI* lockIn);

    virtual ~SPIFeesExposureRT();
    
protected:
    virtual double getCtgtFee(int iStep, double pBF, double pBL);
    
    /******
    // note to save time we rely on getFee having been called first to set the threshold flag
    *******/
    virtual double getCtgtPaidFee(int iStep);
    
private:

    ISPIBondFloorSP bFloor; 
    ILockInSPI* lock;
    bool thresholdBreached;  
};

// With as much pre-processed as poss - for runtime use in payoff
// Possibly better to have an array of structures which hold all the feeFactor, etc for a given step?
class PRODUCTS_DLL SPIRunTime {
public:
    InstrumentSettlementConstSP  settlement;
    const SPIDynamicBasket*  dynBask;
    IAlgorithmSPISP          algo;
    IBondSPISP               bond;
    IFeesSPISP               fees;
    ICouponsSPISP            coupons;
    ILoanCostSPISP           loanCost;
    ICutoffSPISP             cutoff;
    ILockInSPI*              lockIn;
    ISPIBondFloorSP          bondFloor;
    int                      numAlgAssets;     // convenient
    const YieldCurve*        disc;       // discount curve - needed for vol interp pv'ing strike
    int           numSteps;      // convenient
    int           iFirstFutureStep; // ditto
    double        B;     // dynamic basket level
    double        BF;   // bond floor level
    double        BL;   // lock in level
    double        SC;   // sustainable crash
    double        UE;   // unbalanced exposure
    double        UC;   // unbalanced crash
    double        SE;   // sustainable exposure
    double        TE;   // target exposure
    double        ZSpecial; // value priced off the EURIBOR curve
    int           haveZSpecial;  // -1 means not, else >=0 and it's the index for the "special" step
    double        nZ;    // current allocation to bond (so set at previous step)
    DoubleArray   nE;    // [numEquities] current allocation to each equity (so set at previous step)

    DoubleArray   A;     // [numEquities] - purely for debug
    double        LC; // loan cost
    double        LB; // loan balance
    DoubleArray   RC;     // [numEquities] - purely for debug
    DoubleArray   RCFactor;     // [numSimDates] - to handle different rates in past/future
    int           cutoffStep;  // step at which cutoff occurs
    bool          isCutoff;
    double        sumB; // for averaging out of basket
    double        couponAmt;
    double        payoff;
    double        sumCoupons; // cumulative coupon
    double        accumulatedFeeToPay;
    double        paidFeeToPay; // for reporting
    double        currentFeeAmount;
    bool          isRebal;
    bool          debugOn;
    double        notional;
    double        matCashFlow;

    // Handling trading lag adjustment :-
    // The lag interval may cover several rebalancing events so need to keep up to #lag days values around (for
    // the case of highest freq rebal being one rebalance per day during the lag interval)
    // XXX NOT SURE which is more efficient. Having to refresh these structures each iter with the "SoFar" construct
    // OR to have arrays across all time steps which don't need to be refreshed but the indexing is "global" and
    // more values have to be stored?
    DoubleArray   E; // [numEquities] values at this step - offset for publication lag
    DoubleArray   Einit; // [numEquities] values at 0th step needed for vol interp
    DoubleMatrix  ELag; /* [numEquities][numLagDays], so [i][j] holds value of E for asset #i, where the j rotates
                            through 0, 1, ..., gx-1 being (iStep%gx) so can reliably retrieve values for up to the
                            previous gx days.
                            The numLagDays may be different for each equity so allocate max.
                            ELag is assigned to at each rebal event.
                        */
    DoubleMatrix  nEDiffLag;   /* [numAlgAssets][numLagDays]. so [i][j] holds value of nE-nE[-1] for asset #i, "j days ago" according
                                    to the indexing prescription given above.
                                    nEDiffLag is assigned to at each rebal event */
    BoolArray     isRebalLag;  /* [numLagDays]. so [j] holds value isRebal from "j days ago (0<=j<maxNumLagDays]" according
                                    to the indexing prescription given above.
                                    isRebalLag is assigned to each simulation date.  */

    IntArray              sampleFlags; // 0 - sample only, 1 - rebalance, 2 - fees, 4 - lock in, 8 - FINAL_REBAL
    int                   iStepLastRebal;     // iStep value for which sampleFlags is FINAL_REBAL
    int                   iStepFirstRebal;    // iStep before which we are simply gathering pub lag equity level info

	// run time objects
    SPICouponsRTSP couponsRT;
    SPIFeesRTSP feesRT;
    SPIPostCutoffRTSP postCutoff;

    DateTimeArray   gapRiskBucketDates; // derived from the Offsets relative to valueDate
    IntArray        gapRiskBucketIdx; // [total numsteps] = iBucket
    DoubleArray     gapRiskByBucket;  // [num buckets] summed into buckets for each iter

    SPIKnownFlowsSP          instKnownCashFlows;

    // fields for handling past
    double        BSoFar;
    double        BLSoFar;
    double        nZSoFar;
    DoubleArray   nESoFar;
    double        LCSoFar;
    double        LBSoFar;
    DoubleMatrix  ELagSoFar;
    DoubleMatrix  nEDiffLagSoFar;
    BoolArray     isRebalLagSoFar;
    DoubleArray   ASoFar;     // [numEquities] - purely for debug
    DoubleArray   RCSoFar;     // [numEquities] - purely for debug
    double        sumBSoFar; 
    double        sumCouponsSoFar; 
    bool          isCutoffSoFar;
    double        accumulatedFeeToPaySoFar;
    bool          terminatedEarly;
    int           terminalStep;
    double        terminalCoupon;
    double        terminalTotalCoupon;
    double        targetCoupon;
    double        bonusCoupon;
    double        terminalTotalRedemption;

    // start values from the "SoFar" 
    // I think I need a BSoFar to handle after mat lag-adj-only use of value()
    virtual void init(int  startIdx,
                      bool doingPast);

    // preserve values from past into "SoFar"
    void setSoFar();
                                           
    SPIRunTime(const SyntheticPortfolioInsurance* inst,
               ILockInSPI*                        lockIn,
               InstrumentSettlementConstSP        settlement,
               const DateTime*                    instPaymentDate,
               const DateTimeArray&               feePayDates);

    // single algorithm step - returns true if we terminate early at this step
    bool singleSPIStep(bool doingPast, int iStep, int endIdx, 
                       double& futurePayoff);

    bool getBondClosing(double& closing);

    // the i is step index
    double value(int i);

    // needed for gap risk only!
    double equityComponent();

    double rebalCost(int iStep, int iAsset);

    // this could be const if it wasn't for the debug array A[]!
    double lagAdj(int iStep, int iAsset);

    double unbalExposure(int i);

    double sustLoss(int iStep, double BF);

    void rebalanceAndPrepareForNextStep(int iStep, bool isRebal, 
                                        double crash, double exposure);

    // Vol interp = Cash(i)/( H(i) * Z(i) ) with H(i) = nE1(i)*E1(0)+nE2(i)*E2(0)
    // and Cash(i) = BF/nZ.Z/F/LC etc
    double volInterp(int               endIdx,      // endIdx of the past path gen
                        double            BFSoFar,
                        double            strike);

    /** I know these shouldn't be here, but ... */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){};
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const {};

    virtual double Z(int iStep) const;
    virtual double loanCostAccrueFactor(int iStep);
    virtual double feeTimeFactor(int iStep);
    virtual double feeDF(int iStep);
    virtual double couponDF(int iStep);
    virtual void refreshDiscountFactorsQuotients(DoubleArray& DFQuot) const;
private:
    SPIRunTime(const SPIRunTime& rhs); // not implemented
    SPIRunTime& operator=(const SPIRunTime& rhs); // not implemented

protected:
    DoubleArray   bondPrices;     // [numSimDates], bond prices
    DoubleArray   loanCostAccrueFactorArray; // [numSimDates]
    DoubleArray   feeTimeFactorArray;        // [numSimDates]

    // A shame to need 2 arrays like this, but coupons pay according to a different schedule to fees
    DoubleArray   feeDFArray;       // from (sample date+delay) to mat. Used for fees being paid
    DoubleArray   couponDFArray;    // from (coupon pay date) to mat. Used also for early redemption from coupon target being met
};

class SPIRunTimeSV : public SPIRunTime {
public:
    
    SPIRunTimeSV(const SyntheticPortfolioInsurance* inst,
                 ILockInSPI*                        lockIn,
                 InstrumentSettlementConstSP        settlement,
                 const DateTimeArray&               feePayDates);

    void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    void init(int  startIdx,
              bool doingPast);

    double Z(int iStep) const;
    double loanCostAccrueFactor(int iStep);
    double feeTimeFactor(int iStep);
    double feeDF(int iStep);
    double couponDF(int iStep);
    void refreshDiscountFactorsQuotients(DoubleArray& DFQuot) const;

private:
    SPIRunTimeSV(const SPIRunTimeSV& rhs); // not implemented
    SPIRunTimeSV& operator=(const SPIRunTimeSV& rhs); // not implemented

    // perhaps better to localise these in each appropriate class but 
    // we'll follow the use of SPIRunTime
    SVGenDiscFactorSP                 dfGenCoupon;
    SVDiscFactorSP                    dfSVCoupon;
    SVGenDiscFactorSP                 dfGenFee;
    SVDiscFactorSP                    dfSVFee;
    vector<SVGenExpectedDiscFactorSP> dfGenRebalDates;
    vector<SVExpectedDiscFactorSP>    dfSVRebalDates;
    vector<SVGenExpectedDiscFactorSP> dfGenBond;
    vector<SVExpectedDiscFactorSP>    dfSVBond;

    vector<int>               feeDatesMap; // to map from iStep to fee dates
};

DRLIB_END_NAMESPACE

#endif
