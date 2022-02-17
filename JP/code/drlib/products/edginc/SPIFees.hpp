//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIFees.hpp
//
//   Description : Fees interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_FEES_HPP
#define EDR_SPI_FEES_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/TypeConvert.hpp"
#include "edginc/RevertTypeConvert.hpp"
#include "edginc/SPIUtil.hpp"
#include "edginc/SPILoanCost.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

#define SPI_FEES_TYPE_STD   "Standard"
#define SPI_FEES_TYPE_KNP   "KeptAndPaid"
#define SPI_FEES_TYPE_PNL   "PerNotional"
#define SPI_FEES_TYPE_LIB   "Libor"

// this is the actual interface of what a spi fees object 
// needs to do internally
class PRODUCTS_DLL IFeesSPI {
public:

    virtual void init(int                  numAlgAssets,
                      const DateTimeArray& rebalDates,
                      int                  iStepFirstRebal,
                      double               Basis,
                      ILoanCostSPI*  loanCost,
                      const YieldCurve*    yc) = 0;

    virtual void refresh(int iFirstFutureStep) = 0;

    virtual double getFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const = 0;

    virtual double getPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const = 0;

    virtual double getContingentFeeAmount(int                iStep,
                                double             nZ,
                                double             Z,
                                const DoubleArray& nE,
                                const DoubleArray& E) const = 0;

    virtual double getContingentPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const = 0;

    virtual double getFeeAtMin(const DoubleArray& exposureUpperBound,
                               double             equityExposureMin) const = 0;

    virtual bool hasContingentFees() const = 0;

    virtual bool isThresholdBreached(double level) const = 0;

    virtual const DateTimeArray& getNotificationDates(const DateTimeArray& rebalDates) const = 0;

    virtual DateTimeArray getPaymentDates() const = 0;

    virtual const DateTime getFeePayDate(const DateTime& rebalDate) const = 0;

    virtual DateTimeArray getEssentialDates() const = 0;

    // default implementation does nothing
    virtual void crossValidate(double equityExposureMin) const {};

    // post getMarket validation (handles settlement of fees)
    virtual void Validate(const InstrumentSettlement* instSettle) = 0;

    virtual ~IFeesSPI();

protected:
    // helper for getting settlement holidays
    static HolidayConstSP getSettleHolidays(const InstrumentSettlement* instSettle);
};
DECLARE_REF_COUNT(IFeesSPI);

/*****************************************************************************/
// this is the external interface for abstraction so that users can bolt in 
// any fee type they want - note this includes the SPIFeesWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface IFeesSPI as soon as possible
class PRODUCTS_DLL IFeesSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    // get the internal interface
    virtual IFeesSPISP getFeesSPI() = 0;

    virtual const string feeType() = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IFeesSPIInterface> IFeesSPIInterfaceSP;

/*****************************************************************************/
// Fee based on basket value
class PRODUCTS_DLL SPIFeesStd : public SPIInterfaceIMS,
                   virtual public IFeesSPI,
                   virtual public IFeesSPIInterface {
private:
    double        bondFee;
    DoubleArray   equityFee;             // [numEquities]
    DoubleArray   feeTimeFactor;         // [numSteps] $unregistered
    double        contingentBondFee;
    DoubleArray   contingentEquityFee;   // [numEquities]
    double        contingentThreshold;
    bool          hasContingentFee;

    // The concept of "fee payment date" changed partway through (with the advent of kept & paid fees)
    // so for "std" fees it is now a meaningless concept.
    // paymentDates & isPayDaily are kept around so don't break IMS.
    bool          isPayDaily;
    DateTimeArray paymentDates;
    DateTimeArray realPaymentDates;  // since under the new definition there are none, this is empty :)

public:
    friend class SPIAlgorithmStd;

    static CClassConstSP const TYPE;
    SPIFeesStd(); // for reflection

    virtual IFeesSPISP getFeesSPI();

    virtual const string feeType();

    void init(int                  numAlgAssets,
              const DateTimeArray& rebalDates,
              int                  iStepFirstRebal,
              double               Basis,
              ILoanCostSPI*  loanCost,
               const YieldCurve*    yc);

    void refresh(int iFirstFutureStep);

    double getFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    double getContingentFeeAmount(int                iStep,
                                  double             nZ,
                                  double             Z,
                                  const DoubleArray& nE,
                                  const DoubleArray& E) const;


    virtual double getPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual double getContingentPaidFeeAmount(int                iStep,
                                              double             nZ,
                                              double             Z,
                                              const DoubleArray& nE,
                                              const DoubleArray& E) const;

    virtual bool hasContingentFees() const;

    virtual bool isThresholdBreached(double level) const;

    virtual double getFeeAtMin(const DoubleArray& exposureUpperBound,
                               double             equityExposureMin) const;

    virtual const DateTimeArray& getNotificationDates(const DateTimeArray& rebalDates) const;

    virtual DateTimeArray getPaymentDates() const;

    virtual const DateTime getFeePayDate(const DateTime& rebalDate) const;

    virtual DateTimeArray getEssentialDates() const;

    // post getMarket validation (handles settlement of fees)
    virtual void Validate(const InstrumentSettlement* instSettle);

private:
    SPIFeesStd(const SPIFeesStd& rhs); // not implemented
    SPIFeesStd& operator=(const SPIFeesStd& rhs); // not implemented

    static IObject* defaultSPIFeesStd();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIFeesStd> SPIFeesStdSP;

/*****************************************************************************/
// Fee based on basket value wit some kept and some paid
class PRODUCTS_DLL SPIFeesKeptAndPaid : public SPIInterfaceIMS,
                   virtual public IFeesSPI,
                   virtual public IFeesSPIInterface  {
private:
    double        keptBondFee;
    double        paidBondFee;
    DoubleArray   keptEquityFee;             // [numEquities]
    DoubleArray   paidEquityFee;             // [numEquities]
    DoubleArray   totalEquityFee;            // [numEquities], transient
    double        contingentKeptBondFee;
    double        contingentPaidBondFee;
    DoubleArray   contingentKeptEquityFee;             // [numEquities]
    DoubleArray   contingentPaidEquityFee;             // [numEquities]
    DoubleArray   totalContingentEquityFee;            // [numEquities], transient
    double        contingentThreshold;
    bool          hasContingentFee;
    DateTimeArray paidFeePaymentDates; // actually notification dates - historical so can't fix now
    DateTimeArray paidFeeSettlementDates; // the actual settlement dates
    CashSettlePeriodSP feeSettle; // the settlement object for the fees
    bool          hasSettlePeriod; // do fees have their own settlement?
    int           settlePeriod;      // will become a CashSettlePeriod with same hols as inst
    DoubleArray   feeTimeFactor;         // [numSteps] $unregistered

    // following for Credit SPI a fee which depends on basket level but which
    // isn't deducted from basket. It's just paid.
    // they are handled in getPaidFeeAmount so that they do not affect
    // the bond floor, vol interp or the basket level
    bool          hasPaidNotDeductedFee;
    DoubleArray   paidNotDeductedEquityFee;             // [numEquities]

public:
    friend class SPIAlgorithmStd;

    static CClassConstSP const TYPE;
    SPIFeesKeptAndPaid();// for reflection

    // validation
    void validatePop2Object();

    virtual IFeesSPISP getFeesSPI();

    virtual const string feeType();

    void init(int                  numAlgAssets,
              const DateTimeArray& rebalDates,
              int                  iStepFirstRebal,
              double               Basis,
              ILoanCostSPI*  loanCost,
              const YieldCurve*    yc);

    void refresh(int iFirstFutureStep);

    double getFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    double getContingentFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    virtual double getPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual double getContingentPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual bool isThresholdBreached(double level) const;

    virtual bool hasContingentFees() const;

    virtual double getFeeAtMin(const DoubleArray& exposureUpperBound,
                               double             equityExposureMin) const;

    virtual const DateTimeArray& getNotificationDates(const DateTimeArray& rebalDates) const;

    virtual DateTimeArray getPaymentDates() const;

    virtual const DateTime getFeePayDate(const DateTime& rebalDate) const;

    virtual DateTimeArray getEssentialDates() const;

    // post getMarket validation (handles settlement of fees)
    virtual void Validate(const InstrumentSettlement* instSettle);

private:
    SPIFeesKeptAndPaid(const SPIFeesKeptAndPaid& rhs); // not implemented
    SPIFeesKeptAndPaid& operator=(const SPIFeesKeptAndPaid& rhs); // not implemented

    static IObject* defaultSPIFeesKeptAndPaid();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIFeesKeptAndPaid> SPIFeesKeptAndPaidSP;

/***********************************************************************************/
// Fee rate is based on notional, not basket.
// So no split across equity/bond at least

class PRODUCTS_DLL SPIFeesPerNotional : public SPIInterfaceIMS,
                   virtual public IFeesSPI,
                   virtual public IFeesSPIInterface  {
private:
    double        keptFee;
    double        paidFee;
    double        contingentKeptFee;
    double        contingentPaidFee;
    double        contingentThreshold;
    bool          hasContingentFee;
    DateTimeArray paidFeePaymentDates; // actually notification dates - historical so can't fix now
    DateTimeArray paidFeeSettlementDates; // the actual settlement dates
    CashSettlePeriodSP feeSettle; // the settlement object for the fees
    bool          hasSettlePeriod; // do fees have their own settlement?
    int           settlePeriod;      // will become a CashSettlePeriod with same hols as inst
    DoubleArray   feeTimeFactor;         // [numSteps] $unregistered

    // following for Credit SPI a fee which depends on BASKET level but which
    // isn't deducted from basket. It's just paid.
    // they are handled in getPaidFeeAmount so that they do not affect
    // the bond floor, vol interp or the basket level
    bool          hasPaidNotDeductedFee;
    DoubleArray   paidNotDeductedEquityFee;             // [numEquities]

public:
    friend class SPIAlgorithmStd;

    static CClassConstSP const TYPE;
    SPIFeesPerNotional();// for reflection

    // validation
    void validatePop2Object();

    virtual IFeesSPISP getFeesSPI();

    virtual const string feeType();

    void init(int                  numAlgAssets,
              const DateTimeArray& rebalDates,
              int                  iStepFirstRebal,
              double               Basis,
              ILoanCostSPI*  loanCost,
              const YieldCurve*    yc);

    void refresh(int iFirstFutureStep);

    double getFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    // when we know it's fee on notional ... i.e. in bond floor calcs
    double getFeeAmount(int iStep) const;

    double getContingentFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    virtual double getPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual double getContingentPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual bool isThresholdBreached(double level) const;

    virtual bool hasContingentFees() const;

    virtual double getFeeAtMin(const DoubleArray& exposureUpperBound,
                               double             equityExposureMin) const;

    virtual const DateTimeArray& getNotificationDates(const DateTimeArray& rebalDates) const;

    virtual DateTimeArray getPaymentDates() const;

    virtual const DateTime getFeePayDate(const DateTime& rebalDate) const;

    virtual DateTimeArray getEssentialDates() const;

    // post getMarket validation (handles settlement of fees)
    virtual void Validate(const InstrumentSettlement* instSettle);

private:
    SPIFeesPerNotional(const SPIFeesPerNotional& rhs); // not implemented
    SPIFeesPerNotional& operator=(const SPIFeesPerNotional& rhs); // not implemented

    static IObject* defaultSPIFeesPerNotional();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIFeesPerNotional> SPIFeesPerNotionalSP;

/***********************************************************************************/
// Fee based on risky holdings with rate determined by Libor
// Fee rate is specified per equity as mult*Libor + spread
// currently disabled spread while we consider bond floor
class PRODUCTS_DLL SPIFeesLibor : public SPIInterfaceIMS,
                     virtual public IFeesSPI,
                   virtual public IFeesSPIInterface  {

private:
    DoubleArray     keptEquityFeeMult;     // [numEquities]
    DoubleArray     paidEquityFeeMult;     // [numEquities]
    DoubleArraySP   keptEquityFeeSpread;  // [numEquities]
    DoubleArraySP   paidEquityFeeSpread;  // [numEquities]
    DateTimeArray   paidFeePaymentDates; // actually notification dates - historical so can't fix now
    DateTimeArray paidFeeSettlementDates; // the actual settlement dates
    CashSettlePeriodSP feeSettle; // the settlement object for the fees
    bool            hasSettlePeriod; // do fees have their own settlement?
    int             settlePeriod;      // will become a CashSettlePeriod with same hols as inst
    DoubleMatrix    feeFactors;            // [numSteps] $unregistered
    DoubleMatrix    paidFeeFactors;        // [numSteps] $unregistered
    DoubleArray     timeFactors;
    ILoanCostSPI*   loanCost; // of course, this should be const...  $unregistered
    int             numAlgAssets;

public:
    friend class SPIAlgorithmStd;

    static CClassConstSP const TYPE;
    SPIFeesLibor();// for reflection

    // validation
    void validatePop2Object();

    virtual IFeesSPISP getFeesSPI();

    virtual const string feeType();

    void init(int                  numAlgAssets,
              const DateTimeArray& rebalDates,
              int                  iStepFirstRebal,
              double               Basis,
              ILoanCostSPI*  loanCost,
              const YieldCurve*    yc);

    void refresh(int iFirstFutureStep);

    double getFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    double getContingentFeeAmount(int                iStep,
                        double             nZ,
                        double             Z,
                        const DoubleArray& nE,
                        const DoubleArray& E) const;

    virtual double getPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual double getContingentPaidFeeAmount(int                iStep,
                                    double             nZ,
                                    double             Z,
                                    const DoubleArray& nE,
                                    const DoubleArray& E) const;

    virtual bool isThresholdBreached(double level) const;

    virtual bool hasContingentFees() const;

    virtual double getFeeAtMin(const DoubleArray& exposureUpperBound,
                               double             equityExposureMin) const;

    virtual const DateTimeArray& getNotificationDates(const DateTimeArray& rebalDates) const;

    virtual DateTimeArray getPaymentDates() const;

    virtual const DateTime getFeePayDate(const DateTime& rebalDate) const;

    virtual DateTimeArray getEssentialDates() const;

    // overriden version for Libor fees. Stops spread if min equity exposure
    virtual void crossValidate(double equityExposureMin) const;

    // post getMarket validation (handles settlement of fees)
    virtual void Validate(const InstrumentSettlement* instSettle);

private:
    SPIFeesLibor(const SPIFeesLibor& rhs); // not implemented
    SPIFeesLibor& operator=(const SPIFeesLibor& rhs); // not implemented

    static IObject* defaultSPIFeesLibor();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIFeesLibor> SPIFeesLiborSP;

/****************************************************************************/
// a wrapper for the SPIFeesWrapper so we can split the IMS interface into 2
// note we needed all this before abstraction in IMS
class PRODUCTS_DLL SPIFeesWrapperSuper : public CObject,
                            virtual public ITypeConvert {
public:
    static CClassConstSP const TYPE;

    string                spiFeesType;
    SPIFeesStdSP          feesStd;
    SPIFeesKeptAndPaidSP  feesKeptAndPaid;
    SPIFeesPerNotionalSP  feesPerNotional;
    SPIFeesLiborSP        feesLibor;

private:
    SPIFeesWrapperSuper();
    SPIFeesWrapperSuper(const SPIFeesWrapperSuper& rhs);
    SPIFeesWrapperSuper& operator=(const SPIFeesWrapperSuper& rhs);

public:
    /** create a proper SPIFeesWrapper */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultSPIFeesWrapperSuper();

    // constructor for use with the SPIFeesWrapperSuper class
    SPIFeesWrapperSuper(string type,
                        SPIFeesStdSP std,
                        SPIFeesKeptAndPaidSP  keptAndPaid,
                        SPIFeesPerNotionalSP  perNotional,
                        SPIFeesLiborSP        libor);
};

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPIFeesWrapper : public CObject,
                       virtual public IRevertTypeConvert,
                       virtual public IFeesSPIInterface  {
public:
    static CClassConstSP const TYPE;

    string                SPIFeesType;
    SPIFeesStdSP          feesStd;
    SPIFeesKeptAndPaidSP  feesKeptAndPaid;
    SPIFeesPerNotionalSP  feesPerNotional;
    SPIFeesLiborSP        feesLibor;

    virtual IFeesSPISP getFeesSPI();

    virtual const string feeType();

    // validation
    void validatePop2Object();

    /** revert the convert method of SPIFeesWrapperSuper
        needed for PYRAMID and IMS */
    virtual IObjectSP revert(const string& interfaceType) const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    // for reflection
    SPIFeesWrapper();

    // constructor for use with the SPIFeesWrapperSuper class
    SPIFeesWrapper(string type, SPIFeesStdSP std,
                    SPIFeesKeptAndPaidSP  keptAndPaid,
                    SPIFeesPerNotionalSP  perNotional,
                    SPIFeesLiborSP        libor);

    static IObject* defaultSPIFeesWrapper();
};
typedef smartPtr<SPIFeesWrapper> SPIFeesWrapperSP;

DRLIB_END_NAMESPACE

#endif

