//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PVDiv.hpp
//
//   Description : value the synthetic dividend (or PVDiv) instrument
//
//
//----------------------------------------------------------------------------

#ifndef EDG_PV_DIV_HPP
#define EDG_PV_DIV_HPP

#include "edginc/Class.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/PVDivCreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CPVDiv : public CInstrument, 
               public virtual CClosedFormLN::IIntoProduct,
               public virtual Theta::IShift, 
               public MuParallel::IShift,
               public MuSpecial::IShift, 
               public MuPointwise::IShift,
               public virtual LastSensDate,
               public virtual CreditSupport::Interface
{
public:
    static CClassConstSP const TYPE;

    // override base implementation if required
    virtual void GetMarket(const IModel*, const CMarketDataSP);
    
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Returns rolls value date and sets samples for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    /** Returns name identifying this object for MU_PARALLEL */
    virtual string sensName(MuParallel* shift) const;
    
    /** Shifts the object using given shift. This is a wrapper for the
        DividendList MU_PARALLEL shift method */
    virtual bool sensShift(MuParallel* shift);

    /** Returns name identifying this object for MU_POINTWISE */
    virtual string sensName(MuPointwise* shift) const;
    
    /** Shifts the object using given shift. This is a wrapper for the
        DividendList MU_POINTWISE shift method */
    virtual bool sensShift(MuPointwise* shift);

        /** Returns name identifying this object for MU_S */
    virtual string sensName(MuSpecial* shift) const;

    /** Shifts the object using given shift. This is a wrapper for the
        DividendList MU_S shift method */
    virtual bool sensShift(MuSpecial* shift);

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual DateTime getValueDate() const;
        
    bool isCcyStruck() const;

    // returns true if PVDiv is fixed notional, invest on wgted avg with 1 avg date, which is
    // the edg definition of fwd starting.
    bool isLegacyFwdStarting() const;

    DateTime    getCumPayDate() const;

    virtual CreditSupportSP createCreditSupport(CMarketDataSP market);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultPVDiv();

    // for product to access instrument data
    friend class CPVDivClosedFormProd;
    // for credit exposure support
    friend class PVDivCreditSupport;

    CPVDiv();

    DateTime                 valueDate;
    InstrumentSettlementSP   instSettle;
    bool                     isFixedNotional; // true= num of share at each schedule date = N(i)/S(i), else = N(i)
    bool                     isCumulative; 
    SampleListSP             contractSched;
    SampleListSP             averageSched;   // only necessary for notional, type 2.
    CAssetWrapper            equityAsset;
    string                   ccyTreatment;
    YieldCurveWrapper        discount;      
    string                   averagingType; // either "InvestOnContractDates" (Type 1)
                                            // or "InvestOnWeightedAverage" (Type 2)
    bool                     useSynthetic; // to override divs in the underlying
    DividendListSP           synthDivs;    // overriding divs if useSynthetic flag true

    mutable DateTime         lastPayDate; /* transient field - populated in
                                              pricing call */
    class DivAdjuster;
    friend class DivAdjuster;
    void adjustDividend(double          sizeAdjust,
                        const DateTime& cumPayDate,
                        Dividend&       div) const;

};

typedef smartConstPtr<CPVDiv> PVDivConstSP;
typedef smartPtr<CPVDiv> PVDivSP;


DRLIB_END_NAMESPACE
#endif
