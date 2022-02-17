//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CredDefSwaption.hpp
//
//   Description : Credit default swaption
//
//   Author      : André Segger
//
//   Date        : 06 June 2003
//
//
//----------------------------------------------------------------------------

#ifndef CREDDEFSWAPTION_HPP
#define CREDDEFSWAPTION_HPP
#include "edginc/Instrument.hpp"
#include "edginc/Generic1FactorCredit.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/ClosedFormFA.hpp"
#include "edginc/ClosedFormCDSPSandFA.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/DeltaToCredit.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/CreditSpreadVegaParallel.hpp"

DRLIB_BEGIN_NAMESPACE


/** Credit default swap instrument - fixed amounts at known future dates */

class PRODUCTS_DLL CredDefSwaption: public Generic1FactorCredit,
                       public ClosedFormCDSPS::IIntoProduct,
                       public ClosedFormFA::IIntoProduct,
                       public ClosedFormCDSPSandFA::IIntoProduct,
                       public LastSensDate,
                       public CreditSpreadRhoParallel::RestorableShift,
                       public CreditSpreadRhoPointwise::IRestorableShift,
                       public CreditSpreadVegaParallel::IShift,
                       public DeltaToCredit::IShift,
                       virtual public ObjectIteration::IOverride, 
                       virtual public  IGetMarket {
public:
    static CClassConstSP const TYPE;
    friend class CredDefSwaptionHelper;
    friend class CredDefSwaptionClosedFormCDSPS;
    friend class CredDefSwaptionClosedFormFA;
    friend class CredDefSwaptionClosedFormCDSPSandFA;

    virtual void Validate();

    /** used to override CREDIT_SPREAD_RHO on cdsParSpreads
        ObjectIteration::IOverride interface */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    /** instrument specific market data handling */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Implementation of ClosedFormCDSPS::IntoProduct interface */
    virtual ClosedFormCDSPS::IProduct* createProduct(ClosedFormCDSPS* model) const;

    /** Implementation of ClosedFormFA::IntoProduct interface */
    virtual ClosedFormFA::IProduct* createProduct(ClosedFormFA* model) const;

    /** Implementation of ClosedFormCDSPSandFA::IntoProduct interface */
    virtual ClosedFormCDSPSandFA::IProduct* createProduct(ClosedFormCDSPSandFA* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** SpreadRhoParallel support */
    virtual string sensName(CreditSpreadRhoParallel* shift) const;
    virtual bool sensShift(CreditSpreadRhoParallel* shift);
    virtual void sensRestore(CreditSpreadRhoParallel* shift);

    /** SpreadRhoPointwise support */
    virtual string sensName(CreditSpreadRhoPointwise* shift) const;
    virtual ExpiryArrayConstSP sensExpiries(CreditSpreadRhoPointwise* shift) const;
    virtual bool sensShift(CreditSpreadRhoPointwise* shift);
    virtual void sensRestore(CreditSpreadRhoPointwise* shift);

    /** credit spread vega support */
    virtual string sensName(CreditSpreadVegaParallel* shift) const;
    virtual bool sensShift(CreditSpreadVegaParallel* shift);

    /** Returns name identifying vol for vega parallel */
    virtual string sensName(DeltaToCredit* shift) const;

    /** Shifts the object using given shift */
    virtual bool sensShift(DeltaToCredit* shift);

private:
    CredDefSwaption();
    CredDefSwaption(const CredDefSwaption& rhs);
    CredDefSwaption& operator=(const CredDefSwaption& rhs);
    void priceParSpreads(CResults* results, Control* control) const;
    void priceFirmAsset(CResults* results, Control* control) const;
    void addRequests(Control* control, CResults* results, 
                     CashFlowArraySP cleanSpreadCurve, IObjectSP currentSpread,
                     const double bpValue, const double cdsValue,
                     const double fwdSpread)const;
    void buildCreditSpreadCurve()const;

    // Converts basis of default rate from continuous to annualised
    CashFlowArraySP annualiseDefaultRates(const CashFlowArray& defaultRates)const;

    // Inputs of the swaption
    DateTime                optionIssueDate;
    DateTime                optionMaturityDate;
    bool                    optionKOOnDefault;
    string                  optionType;
    double                  fwdSwapSpreadVol;
    bool                    isCall;

    // additional information for skewed distribution
    bool                    doEdgworthExpansion;
    double                  skewness;
    double                  kurtosis;

    // Inputs of the underlying CDS
    DateTime                swapEffectiveDate;
    DateTime                swapMaturityDate;
    int                     frequency;
    double                  swapStrike;
    bool                    useSwapRecovery;
    mutable double          swapRecovery;
    bool                    payAccruedFee;
    string                  dcc;
    string                  bdc;
    // mutable CashFlowArraySP feePayments;

    bool                    addProtectionForPut;

    // Transient 
    mutable CreditSpreadCurveSP  creditSpreadCurve;      // For calculating CREDIT_SPREAD_RHO'S
    mutable CreditSpreadCurveSP  csRestore;              // Backup - cannot restore in a simple way
    mutable DayCountConventionSP swpAccrualDCC;          // used if payAccruedFee = true 
    mutable BadDayConventionSP   BDC;                    // bad day conv for BM adjustment 
    bool                         priceViaParCurves;      // only false for CREDIT_SPREAD_RHO'S
    bool                         createE2Csensitivities; // whether to create debt/equity sensitivities
    mutable bool                 e2cBasePriceCalculated;
    mutable double               e2cBasePrice;
    mutable CVolBaseWrapper      vol;
};

typedef smartConstPtr<CredDefSwaption>  CredDefSwaptionConstSP;
typedef smartPtr<CredDefSwaption>       CredDefSwaptionSP;


DRLIB_END_NAMESPACE
#endif




