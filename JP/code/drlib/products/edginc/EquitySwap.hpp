//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquitySwap.hpp
//
//   Description   value a stock given a quoted price
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EQUITY_SWAP_HPP
#define EDG_EQUITY_SWAP_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ESWEquity.hpp"
#include "edginc/ESWDividend.hpp"
#include "edginc/ESWLibor.hpp"
#include "edginc/ESWCashFlow.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/MuPointwise.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/Theta.hpp"
#include "edginc/EquitySwapCreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CEquitySwap : public CInstrument, 
                    public CClosedFormLN::IIntoProduct,
                    public Theta::IShift,
                    public MuParallel::IShift,
                    public MuSpecial::IShift, 
                    public MuPointwise::IShift,
                    public CreditSupport::Interface,
                    public LastSensDate
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultEquitySwap(){
        return new CEquitySwap();
    }

    // override base implementation if required
    virtual void GetMarket(const IModel*, const CMarketDataSP);
    
    virtual void Validate();

    // initialisation
    virtual void Init();

    virtual void validatePop2Object();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** returns current value date */
    virtual DateTime getValueDate() const;

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

  
    bool sensShift(Theta* shift);

    /** when to stop tweaking */    
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual CreditSupportSP createCreditSupport(CMarketDataSP market){
        return CreditSupportSP(new EquitySwapCreditSupport(this, market));}

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    CEquitySwap();

    // not implemented
    CEquitySwap(const CEquitySwap& rhs);
    CEquitySwap& operator=(const CEquitySwap& rhs);

    // for product to access instrument data
    friend class CEquitySwapClosedFormProd;
    friend class EquitySwapCreditSupport;
    	
    DateTime                valueDate;
    
    bool					isFixedNotional;
    bool					excludeCashFlowsToday;
    bool					isCallable;
    bool                    validateSamples;
    string					swapType;  // STANDARD: standard swap 
                                       // BRZ: Brazil swap              
    mutable DateTime        callDate;
	
    // objects for easy i/o
    ESWEquitySP				eqLeg;
	ESWLiborSP				liborLeg;
	ESWDividendSP			divLeg;
	ESWCashFlowSP     		cashflowLeg;

	CAssetWrapper           asset;
    YieldCurveWrapper       discount;

	// sets all future sampleLevels to 0.
	static void zeroFutureSamples(
							 const DateTimeArray&  sampleDates,
							 DoubleArray&    sampleLevels,
							 const DateTime& valueDate,
							 int offset);

	void setLibAvType();
    void validateHistoricalSamples();

};

typedef smartPtr<CEquitySwap> CEquitySwapSP;

DRLIB_END_NAMESPACE
#endif
