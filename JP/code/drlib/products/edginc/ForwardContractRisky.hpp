//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardContractRisky.hpp
//
//   Description   Forward Contract valued with a risky discount curve
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FORWARD_CONTRACT_RISKY_HPP
#define EDG_FORWARD_CONTRACT_RISKY_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/LastSensDate.hpp"

// added for risky curve
#include "edginc/Generic1FactorCredit.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CredDefSwap.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/ClosedFormFA.hpp"
#include "edginc/ClosedFormCDSPSandFA.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CForwardContractRisky : 
    public CClosedFormLN::IIntoProduct,
    public LastSensDate,
    public Generic1FactorCredit
                                                                                        
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultForwardContract(){
        return new CForwardContractRisky();
    }

    // override base implementation if required
    virtual void GetMarket(const IModel*, const CMarketDataSP);
    
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

        
private:

    // for product to access instrument data
    friend class CForwardContractRiskyClosedFormProd;
    friend class ForwardContractNumerical;
    friend class CForwardContractRiskyHelper;

    CForwardContractRisky();

    DateTime                startDate;
    InstrumentSettlementSP  premiumSettle; // $unregistered
    bool                    fwdStarting;
    double                  initialSpot;
    double                  spotAtMaturity;
    ScheduleSP              exerciseSchedule;

    // old dividend re-invest (zeroing dividend) flag is back
    bool                    divReinvest; // $unregistered

    // instrument specs
    bool                    zeroDivs; // if true, ignore dividends
    bool                    zeroBorrow; // if true, ignore borrow cost

};

typedef smartConstPtr<CForwardContractRisky> ForwardContractRiskyConstSP;
typedef smartPtr<CForwardContractRisky> ForwardContractRiskySP;

DRLIB_END_NAMESPACE
#endif
