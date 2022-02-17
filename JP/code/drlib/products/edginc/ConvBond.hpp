//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ConvBond.hpp
//
//   Description : Convertible Bond Instrument
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 27, 2001
//
//
//----------------------------------------------------------------------------

#ifndef CONVBOND_HPP
#define CONVBOND_HPP
#include "edginc/Instrument.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Bond.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/Asset.hpp"
#include "edginc/FD1F.hpp"
#include "edginc/FD1FGeneric.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/ImpliedYTM.hpp"
#include "edginc/ImpliedYTP.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/IRiskyPricer.hpp"
#include "edginc/ResetSchedule.hpp"
#include "edginc/InstrumentAsAsset.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/AdjCreditSpreadRhoParallel.hpp"
#include "edginc/AdjCreditSpreadRhoPointwise.hpp"


DRLIB_BEGIN_NAMESPACE


//////////////////////////////////////////////////////////////////////
//                                                                  //
//          Contingent Conversion : declaration                     //
//                                                                  //
// Handle contingent conversion                                     //
// conversion occurs if any of 3 condition is satisfied             //
//  a) periodic sample, stock price M out of N days above           //
//      trigger * conversionPrice                                   //
//      => convertible in the next period                           //
//  b) after above periods, stock ever hit trigger*conversionPrice  //
//      => convertible thereafter                                   //
//  c) between given start/end dates, average cvb price less than   //
//      trigger * avg parity in any rolling certain # days and		//
//		stock price on the day before conversion is within barrier	//
//      => convertible the day after the rolling window             //
//                                                                  //
//////////////////////////////////////////////////////////////////////

class ConvBond;

class PRODUCTS_DLL ContConversion : public CObject
{
public:
    static CClassConstSP const TYPE;

    /** validation */
    virtual void validatePop2Object();

    // some processing functions
    virtual void getMarket(const IModel* model, const MarketData* market, const ConvBond *cvb);

    // functions used in pricing by CVB product
    bool isTrivial() const;
    bool isHistEnabled(const DateTime &date) const;
    bool isTrigActiveAtMat() const { return periodTrigActiveAtMat; }
    // set the vol 
    void setVol(CVolProcessedBSConstSP vols);
    // insert period dates and floor start/end dates
    virtual void insertCritDates(DateTimeArray &critDates) const;

    // adjust price by trigger and floor
    void adjustPrices(const DateTime &date, 
                      int nStockStep, 
                      const double *s,  // stock grid
                      double *priceNCV, // non-convertible prices
                      double *priceCV,  // convertible prices
                      double callLevel,
                      const ConvBond *cvb) const;

    // function called by CVB
    bool rollDate(Theta* shift, const ConvBond *cvb);

private:

    void validateSample(const DateTimeArray *dates, DoubleArray *samples, double spot=0);
    void preprocess(const ConvBond *cvb);

	// simple util to see if a feature is non-trivial
	bool hasPeriod() const;
	bool hasPostPeriod() const;
	bool hasFloor() const;

    // relate to historical sample, trigger etc
    DoubleArray getHistQuotes(const DateTimeArray *dates, const DoubleArray *samples, DateTime end, int nDays) const;
    int getHistTrigCount(const DateTimeArray *dates, const DoubleArray *samples, const DateTime &end, int nDays, double trigger, int *nHistDays=0) const;
    double getHistSum(const DateTimeArray *dates, const DoubleArray *samples, const DateTime &end, int nDays, int *nHistDays=0) const;

    // trigger adj
    static double getAdjTrigger(int m, int n,  // m out of past n days
                                double trigger, // trigger level in $ terms
                                double fwdRate, // annual fwd rate
                                double vol);    // % vol for n days

    // find CVB info
    double getConversionPrice(const DateTime &date, const ConvBond *cvb) const;

private:
    friend class ContConversionHelper;

    ContConversion();

    string          badDayConvString;

    // periodic sampling to see if stock price M out of N days above trigger
    DateTime        firstPeriodStart;
    int             frequency;
    int             numPeriod;
    int             periodDaysM;
    int             periodDaysN;
    double          periodTrigger;
    string          periodTrigAdjType;
    
    // after all above periods, knock in trigger for conversion
    double          postPeriodTrigger;
    bool            postPeriodEnabled;
    bool            periodTrigActiveAtMat;

    // floor on cvb. if average cvb price over any rolling #floorDaysM days
    // is less than trigger * avg parity, and if the equity price on the 
    // conversion date is between up/down barriers, can convert in the ensuing 
    // window of floorDaysN days
    int             floorDaysM;
    int             floorDaysN;
    double          floorTrigger;
    double          upEqBarrier4Flr;
    double          dnEqBarrier4Flr;

    // historical record of stock price for periodic sample and for floor
    DateTimeArraySP eqSampleDates;  
    DoubleArraySP   eqSampleValues;

    // historical record of stock price for periodic sample and for floor
    DateTimeArraySP cvSampleDates;
    DoubleArraySP   cvSampleValues;

    ///////////////  INTERNAL FIELD //////////////////////////////////
    DateTime        valueDate;
    DateTime        bondMatDate;
    HolidayConstSP  hols;
    BadDayConventionSP   badDayConv;

    bool            curPeriodEnabled;   // flag to indicate current period is convertible
                                        // based on historical sample 
    bool            floorEnabled;       // historical triggering of floor provision. applicable to valueDate only         
    int             nHistPeriod;        // nb of periods with fully historical sampling
    DateTimeArraySP periodStarts;
    DoubleArraySP   periodTriggers;

    bool            ignoreFlrHist;      // ignore eq/cv history for floor provision $unregistered

    CVolProcessedBSConstSP vols;                // ATM vol for trigger adjustment $unregistered
};

typedef smartConstPtr<ContConversion> ContConversionConstSP;
typedef smartPtr<ContConversion> ContConversionSP;


// simple schedule to handle embedded warrant

class PRODUCTS_DLL AddOnConvRatio : public CObject {
public:
    static CClassConstSP const TYPE;

    /** validation */
    virtual void validatePop2Object();

    int nbSchedule() const;

    // obtain next trigger level strictly > spot. 
    // if no more trigger, return spotMax
    double getNextTrigger(double spot, double spotMax);

    // here ratio and cash are both input and output
    void adjustConvRatio(const DateTime& convDate, 
                        bool    isSpot,
                        double  spotOrLevel,    // level is value of sum of converted shares + cash, ie net conv value
                        double *ratio, 
                        double *cash) const;    // notice that cash is negative

	void insertCritDates(DateTimeArray &critDates) const;

private:
    friend class AddOnConvRatioHelper;

    AddOnConvRatio();

    DoubleArraySP   triggers;
	ScheduleSP		schedule1;
	ScheduleSP		schedule2;
	ScheduleSP		schedule3;

    // alternative simple interface for the schedules
    // extra shares apply between start/end dates, inclusive
    bool            interfaceEZ;
    DateTime        startDate;
    DateTime        endDate;
    DoubleArraySP   shares;
};

typedef smartConstPtr<AddOnConvRatio> AddOnConvRatioConstSP;
typedef smartPtr<AddOnConvRatio> AddOnConvRatioSP;

// class to function as data holder for putToStock feature. 
// with very basic function calls
class PRODUCTS_DLL PutToStock : public CObject {
public:
    static CClassConstSP const TYPE;

    /** validation */
    virtual void validatePop2Object();

public: // no need to hide
    PutToStock();

    static void load(CClassSP& clazz);
    static IObject* defaultPutToStock();

    bool          putToStockNotified;
    DateTime      putToStockPutDate;
    int           putToStockSampleDays;
    double        putToStockDiscount;

    // historical equity sample for putToStock if already notified
    DateTimeArraySP eqSampleDates;  
    DoubleArraySP   eqSampleValues;

    // transient field
    SampleListSP  putToStockSamples;
};

typedef smartConstPtr<PutToStock> PutToStockConstSP;
typedef smartPtr<PutToStock> PutToStockSP;

/** Convertible Bond Instrument */

class PRODUCTS_DLL ConvBond: public CInstrument, 
                public FD1F::IIntoProduct,
                public FDModel::IIntoProduct,
                public FD1FGeneric::IIntoProduct,
                public LastSensDate,
                public Theta::Shift,
                public ImpliedYTM::IFaceYTM,
                public ImpliedYTP::IFaceYTP,
                public ISensitiveStrikes,
                public IScaleOutputs,
                public IRiskyPricer,
                virtual public AdjCreditSpreadRhoParallel::IRestorableShift,
                virtual public AdjCreditSpreadRhoPointwise::IRestorableShift,
                virtual public IInstrumentAsAsset{
public:
    static CClassConstSP const TYPE;
    friend class ConvBondHelper;
    friend class ConvBondClosedForm;
    friend class OptOnConvBond;
    friend class TracerProd;
    friend class ConvBondAddin;
    friend class DECSVolAddin;
	friend class ConvBond1fProdDDE;

    static ConvBond* createDECS(YieldCurveWrapper         discount,
                                CreditCurveWrapper        creditSpreads,
                                CAssetWrapper             asset,
                                const DateTime&           valueDate,
                                const string&             ccyTreatment,
                                const double&             faceValue,
                                MaturityPeriodSP          bondMaturity,
                                HolidayWrapper            holidays,
                                int                       bondFrequency,
                                const string&             dayCountConvString,
                                const double&             coupon,
                                const double&             redemption,
                                const double&             spotAtStart,
                                const double&             downsideProtection,
                                const double&             decsPremium,
                                const bool                divPassThrough,
                                const DateTime&           knockInDate);

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate(); // real validation

    virtual void validatePop2Object(); // initialize a few params

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual bool sensShift(Theta* shift);

    /** Returns name identifying credit curve for adjusted credit rho parallel */
    string sensName(AdjCreditSpreadRhoParallel* shift) const;
    /** Sets riskyGrowth flag */
    bool sensShift(AdjCreditSpreadRhoParallel* shift);
    /** Restores the object to its original form */
    void sensRestore(AdjCreditSpreadRhoParallel* shift);

    /** Returns name identifying credit curve for adjusted credit rho pointwise */
    string sensName(AdjCreditSpreadRhoPointwise* shift) const;
    /** Sets riskyGrowth flag */
    bool sensShift(AdjCreditSpreadRhoPointwise* shift);
    /** Restores the object to its original form */
    void sensRestore(AdjCreditSpreadRhoPointwise* shift);
    ExpiryArrayConstSP sensExpiries(AdjCreditSpreadRhoPointwise* shift) const;

    double getStrike() const;
       
    double getMakeWholePayment( const DateTime & callDate ) const;

    void getCallLevel(const DateTime& callDate,
                      double   fwdVol,
                      bool    *isCallable,
                      bool    *isAdjustedForAccrued,
                      bool    *isHardCall,
                      double  *level) const; 

    void getCallRedemptionLevel(const DateTime& callDate,
                                bool    *isCallable,
                                bool    *isAdjustedForAccrued,
                                double  *level) const;
        
    void getPutLevel(const DateTime& putDate,
                     bool    *isPutable,
                     double  *level) const;

    void getSoftPutLevel(const DateTime& putDate,
                         bool    *isPutable,
                         double  *trigger,
                         double  *level) const;

    void getFirstPutInfo(bool *canCalcYield, DateTime *firstPutDate, double *firstPutLevel) const;

    // given put date, generate sample list for put to stock. all dates are business days
    // all dates are EOD (ok since this is only called if putNotified=false).
    // do NOT validate that all dates are past value date
    SampleListSP getPutToStockSample(const DateTime putDate) const;

    // return the multiplicative putToStock factor to apply to put level
    double calcPutToStockFactor(const DateTime putDate, const SampleList *samples, const CVolRequest *volRequest) const;
    
    double alreadyPutToStockValue(Control* control, Results* results) const;

    void getConversionInfo(const DateTime& convDate,
                           bool    *isConvertible,
                           double  *ratio,
                           double  *cash,
                           double  *parityFloor = 0) const;

	// interface for getting stock dependent conversion ratio. recommended interface over the older one above
    void getConversionInfo(const DateTime& convDate,
                           double  spotOrLevel,
                           bool    *isConvertible,
                           double  *ratio,
                           double  *cash,
                           bool    isSpot=true,
                           double  *parityFloor = 0) const;
	
	// return the strike such that the converted value (incl cash) match the input level
	double getImpliedStrike(const DateTime& convDate, 
                            double  level) const;

	bool hasAddOnConvRatios() const;

	bool hasContConversion() const;

    double getAccruedAtDate(const DateTime& aiDate) const;

    double getBondAccruedAtDate(const DateTime& aiDate) const;

    CashFlowArraySP getCashFlows(const DateTime& startDate) const;
    
    CashFlowArraySP getCoupons(const DateTime& startDate) const;

    void   getMarket(const IModel* model, const MarketData* market);

    bool canPriceClosedForm() const;

    void priceClosedForm(CControl* control, CResults* results, 
                         double& fairValue, double& bondValue) const;

    /** Implementation of CFDGridPass::IntoProduct interface */
    virtual FD1F::IProduct* createProduct(FD1F* model) const;

    /** Implementation of CFDGridPass::IntoProduct interface */
    virtual FD1FGeneric::IProduct* createProduct(FD1FGeneric* model) const;

    /** create a fd payoff tree - new for all fd/tree state variable interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    void recordOutputRequests(Control* control, Results* results, 
                              double fairValue, double bondFloor, bool isAnOCB) const;
    
    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual double yieldToMaturity(double price, bool clean) const;

    virtual double yieldToFirstPut(double price, bool clean) const;

    virtual bool   hasPut() const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    /** returns the coupon stream of the bond from which the fixed leg of the asset swap bond floor is calculated */
    CashFlowArraySP getAssetSwapLeg(const DateTime& swapMaturity, bool payAccruedAtMat) const;


    ConvBond(YieldCurveWrapper         discount,
             CreditCurveWrapper        creditSpreads,
             CAssetWrapper             asset,
             const DateTime&           valueDate,
             const string&             ccyTreatment,
             const double&             faceValue,
             MaturityPeriodSP          bondMaturity,
             HolidayWrapper            holidays,
             int                       bondFrequency,
             const string&             dayCountConvString,
             MaturityPeriodSP          callStartPeriod,
             const double&             callTriggerLevel,
             const double&             coupon,
             const double&             redemption,
             const double&             putLevel,
             MaturityPeriodSP          putStartPeriod,
             const double&             spotAtStart,
             const double&             conversionPremium,
             const DateTime&           knockInDate);

    void scaleOutputs(CControlSP control, ResultsSP unscaledResults);

    double alreadyCalledBlackApprox(Control* control, Results* results) const;
    
    /** Returns the bonds maturity date */
    DateTime getBondMaturityDate() const;

    // IRiskyPricer methods
    void setRisky(bool flag = true);
    bool isRisky() const { return riskyGrowth; }
    DateTime getEffMaturityDate() const { return effEndDate; }
    void setEffMaturityDate(const DateTime & maturityDate) {effEndDate = maturityDate; }

    DateTimeArraySP getPaymentDates() const;

    double getAssetSwapSpread(const double& bondPrice) const;

    double getZSpread(const double& bondPrice) const;

    static double calculateAssetSwapSpread(BondSP               bond,
                                           YieldCurveConstSP    yieldCurve,
                                           DateTime             valueDate,
                                           DateTime             workoutDate,
                                           double               workoutLevel,
                                           double               bondPrice,
                                           int                  frequency,
                                           bool                 stubAtFront,
					   DayCountConventionSP swapDCC);

    static double calculateZSpread(BondSP                bond,
				   YieldCurveConstSP     yieldCurve,
				   DateTime              valueDate,
				   DateTime              workoutDate,
				   double                workoutLevel,
				   double                bondPrice,
				   int                   frequency,
				   bool                  stubAtFront,
				   DayCountConventionSP  swapDCC);

    YieldCurveSP getRiskyCurve() const;

    YieldCurveSP getRiskyCurve(const DateTime* maturityDate) const;

    YieldCurveSP getRiskyCurve(const DateTime* maturityDate, bool useJointDdefProb) const;

    /** IInstrumentAsAsset interface
        Returns a date after which the instrument can no longer be used as an asset  */
    virtual DateTime maturityDate() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the 'coupons' or payments that this instrument will make during
        its lifetime. This can include historic payments. */
    virtual CashFlowArraySP getCoupons() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the dates on which the instrument has to be held in order to
        hold the right to the corresponding coupon as returned by getCoupons() */
    virtual DateTimeArraySP getExCouponDates() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the yield curve used for discounting */
    virtual YieldCurveConstSP getDiscount() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the accured interest (if any) to date */
    virtual double getAccrued() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    // for product to access instrument data
    friend class ConvBondTree1fProd;
    friend class ConvBondFDProd;
    friend class ConvBond1fProd;
	friend class ContConversion;

private:
    ConvBond();
    ConvBond(const ConvBond& rhs);
    ConvBond& operator=(const ConvBond& rhs);

    bool InterpSchedule(const DateTimeArray &dates, const DoubleArray &amounts,
        const string &interpType, bool matAdjust, DateTime interpDate, double *level) const; 

    // Gets asset out of FirmAsset if necessary
    CAssetSP  getEquityAsset() ;
    
    // -- Inputs -- // 

    // -- Instrument -- //                         
    YieldCurveWrapper         discount;
    CreditCurveWrapper        creditSpreads;
    DateTime                  valueDate;

    CAssetWrapper             asset;
    string                    ccyTreatment;

    CAssetWrapper             assetToUse; // $unregistered
    
    // Calls
    bool          redemptionAdjustForAccrued;
    bool          triggerAdjustForAccrued;
    bool          notResetCallTrigger; /* new */
    double        callOptimality; // model???
    string        softCallAdjustType; // model???  
    int           softCallDaysM;                   
    int           softCallDaysN;                   
    int           callNotification;
    string        callTreatment; // model???  
    
    ScheduleSP    softCallSchedule;
    ScheduleSP    softCallTriggerSchedule;
    ScheduleSP    callSchedule;

    ScheduleSP    softPutSchedule;
    ScheduleSP    softPutTriggerSchedule;

    // Puts
    bool          putAdjustForAccrued;             
    ScheduleSP    putSchedule;
    int           putNotification;

    PutToStockSP  putToStock;
    
    // Converison
    ScheduleSP    convCashSchedule;
    ScheduleSP    conversionRatios;
    ScheduleSP    parityFloorSchedule;
    bool          convertIntoIssuerShares;         

	ContConversionSP	contConversion;

	// additional conversions to handle the embedded warrant
    AddOnConvRatioSP    addOnConvRatios;

    // CashFlows
    BondSP     bond;                           

    // Misc
    InstrumentSettlementSP instSettle;             
    InstrumentSettlementSP premiumSettle;          
    bool          getCouponAtMat;                  
    bool          convertOnCall;                   
    double        currentConversionAvg;            
    double        currentConversionWeight;         
    bool          isPreferred;                     
    bool          riskyGrowth; // model???         
    string        makeWholeType;                   
    DateTime      makeWholeDate;                   
    double        makeWholeAmount;                 
    bool          alreadyCalled;                   
    DateTime      dateCallNotifSent;               
    bool          schedsArePcts;
    bool          inDefault;
    double        recoveryPct;
    bool          payAccruedUponDefault;
    DateTime      riskFreeCouponEndDate;   
    DateTime      effEndDate; // $unregistered

    // joint default probability for exchangeable bonds
    // currently a prototype implementation
    bool          useJointDefaultCorrelation;
    double        jointDefaultCorrelation;
    double        stockCreditSpread;

    // MEDS
    string      payoffType;
    bool        DECS;                                   
    bool        PERCS;                                  
    double      initialConvRatio;                       
    double      minConvRatio;                           
    double      initialPrice;                           
    double      convPrice;      
    
    bool        decsHasCutoff;
    double      decsCutoffLevel;
    DateTime    triggerStartDate;
    double      cappedDECSTrigger;

    bool        riskFreeCoupons;

    bool        continuousReset;
    int         delayReset;
    bool        includePutSpread;
    
    // French conventions
    bool frenchDivTreatment;                       
    bool frenchExtendedConv;                       
    MaturityPeriodSP frenchExtConvInt;             


    // Dividend pass-through feature - determines whether the stock dividends are paid as coupons
    bool    dividendPassThrough;
    double  dividendPassThroughPct;

    // Dividend adjusted feature
    bool dividendAdjusted;

    // Reset information
    ResetScheduleSP  resetSchedule;
    ResetScheduleSP  adjustedResetSchedule; // $unregistered
    bool             maxWithParity;
    double           manualResetConvPricePct;


    // state dependent reset - if the stock level hits the trigger during the observation period,
    // the reset level is reset to the reset level
    bool     triggerReset;
    DateTime resetObservationStartDate;
    DateTime resetObservationEndDate;
    double   resetTrigger;
    double   resetConversionRatio;
    string   resetType;

    // asset recovery information
    bool    useAssetRecovery;

    // DECS accelerated maturity - if the stock is above the trigger level, the issuer can
    // call the bond and the payoff is pv(coupons) + DECS payoff
    bool     accelerateDECS;
    double   accelerateDECSTriggerLevel;
    DateTime accelerateDECSStartDate;
    DateTime accelerateDECSEndDate;

    // whether we scale price or not for non preferred bonds
    bool dontScaleNonPrefPriceBy100;

    // debug
    bool DEBUG_SOD_Only;

    // whether stock is in Escrow
    bool    stockInEscrow;

    // contigent conversion trigger
    bool    contingentConversion;
    bool    triggerActiveAtMat;
    double  contingentConversionTrigger;
};


typedef smartConstPtr<ConvBond> ConvBondConstSP;
typedef smartPtr<ConvBond> ConvBondSP;


DRLIB_END_NAMESPACE
#endif
