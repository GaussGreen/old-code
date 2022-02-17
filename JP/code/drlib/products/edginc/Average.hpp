//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Average.hpp
//
//   Description : Average instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 9 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_AVERAGE_HPP
#define EDR_AVERAGE_HPP
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/AverageCreditSupport.hpp"
#include "edginc/Schedule.hpp"

DRLIB_BEGIN_NAMESPACE
class CVolRequestLN;

/** Shell structure for spot/hybrid/ratio payoffs */
class PRODUCTS_DLL Average: public CInstrument,
               public virtual LastSensDate,
               public virtual ISensitiveStrikes,
               public virtual CreditSupport::Interface {
public:
    static CClassConstSP const TYPE;

    /** "constructor" for average spot */
    static Average* makeAvgSpot(
        bool                        isCall,
        const DateTime&             maturity,
        double                      strike,
        const SampleList*           avgOut,
        const InstrumentSettlement* instSettle,
        const InstrumentSettlement* premiumSettle,
        const CAsset*               asset,
        string                      ccyTreatment,
        const YieldCurve*           discount,
        const DateTime&             valueDate,
        bool                        fwdStarting,
        const DateTime&             startDate,
        bool                        oneContract,
        double                      notional,
        double                      initialSpot);

    /** "constructor" for average ratio */
    static Average* makeAvgRatio(
        bool                        isCall,
        const DateTime&             maturity,
        double                      strike,
        const SampleList*           avgOut,
        const InstrumentSettlement* instSettle,
        const InstrumentSettlement* premiumSettle,
        const CAsset*               asset,
        string                      ccyTreatment,
        const YieldCurve*           discount,
        const DateTime&             valueDate,
        const SampleList*           avgIn,
        double                      notional);
    
    /** Price a spread, lowStrike and highStrike are passed by reference so the
        calling function can have the strikes converted to absolute values
        for fwdStarting options without recalculating the asset fwd value */
    static double priceSpotSpread(const DateTime& valueDate,
                                  const DateTime& startDate,
                                  const DateTime& maturityDate,
                                  bool  fwdStarting,
                                  bool  isCall,
                                  bool  oneContract,
                                  double notional,
                                  double initialSpot,
                                  double& lowStrike,
                                  double& highStrike,
                                  const InstrumentSettlement* instSettle,
                                  const Asset* asset,
                                  const YieldCurve* discount,
                                  const SampleList* avgOut);
    
    /* determine how to interpolate on a vol surface for a non fwd starting
     * average price option. Use effective strike (strike adjusted by 
     * historic samples) scaled by 1/(total future weights) for interpolation.
     */
    static CVolRequestLN* volInterpSpot(const DateTime& valueDate,
                                        const DateTime& startDate,
                                        const DateTime& maturity,
                                        bool  fwdStarting,
                                        double strike,
                                        const SampleList* avgOut);
        
    /* retrieve market data needed by Average - just valueDate, asset and
       discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** credit support */
    virtual CreditSupportSP createCreditSupport(CMarketDataSP market);

    /** returns ave-in sample list. returns 0 if does not exist.
        default returns 0;*/
    virtual SampleListSP getAveIn() const;

    /** expose below method for 3 hidden child classes */
    static DoubleArraySP getSensitiveStrikes(Average* inst,
                                            OutputNameConstSP outputName,
                                            const IModel*      model);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** input data validation */
    void validatePop2Object()
    {
        static const string method("Average::validatePop2Object");

        // can't get instrument settlement from Market - fail if it is NULL
        if( !instSettle )
        {
            throw ModelException(method, "Instrument settlement is NULL");
        }
        // check if isAmerican flag is TRUE, the exercise schedule cannot be empty
        // check if isAmerican flag is FALSE, the exercise schedule and maturity cannot be both empty
        if( isAmerican == true ){
            if ( !exerciseSchedule ){
                throw ModelException(method, "The schedule for american exercise is not specified");
            }
            else{
                maturity = exerciseSchedule->lastDate();
                strike = exerciseSchedule->lastValue();
            }
        }
        if( isAmerican == false){
             if (!!exerciseSchedule && !(maturity.empty())){
                 if( maturity != exerciseSchedule->lastDate() )
                    throw ModelException(method,"Maturity has to be the same date as the last date on the exercise schedule whenever exercise schedule is specified");
             }
            if ( !exerciseSchedule && maturity.empty() ) {
                throw ModelException(method, "Please, specify either exercise schedule or maturity date");
            }
            if (!!exerciseSchedule && maturity.empty()){
                maturity = exerciseSchedule->lastDate();
                strike = exerciseSchedule->lastValue();
            }
            if (!exerciseSchedule && !maturity.empty()){
                exerciseSchedule = ScheduleSP(new Schedule( DateTimeArray(1,maturity),
                                                            DoubleArray(1,strike),
                                                            "N" ));
            }
        }
    }

private:
    friend class AverageSpotHelper;
    friend class AverageRatioHelper;
    friend class AverageHybridHelper;

    friend class AverageCreditSupport;

    Average(const Average& rhs);
    Average& operator=(const Average& rhs);
    static void load(CClassSP& clazz);

	virtual bool getFwdStartDate(DateTime &startDate) const
	{ return false; }

protected:
    //// constructor
    Average(CClassConstSP clazz);

    //// how we interpolate vol - common method for MC
    virtual CVolRequestLN* volInterp(double strike) const = 0;   

    class MC; // interacts with MC
    friend class MC; // needs access to fields below

    // fields common to three average types
    bool         isCall;
    DateTime     maturity;
    double       strike;

    SampleListSP      avgOut;
    ScheduleSP        exerciseSchedule;    // FD: exercise schedule
    bool              isAmerican;          // FD: american or no

    InstrumentSettlementSP  instSettle;
    InstrumentSettlementSP  premiumSettle;


    CAssetWrapper      asset;
    string             ccyTreatment;
    YieldCurveWrapper  discount;
    DateTime           valueDate;
};

typedef smartPtr<Average> AverageSP;

// --------------------------------------------------- //
//// data class as external interface ////
class PRODUCTS_DLL AvgOutPerfMaker: public CObject {
public:
    static CClassConstSP const TYPE;  
    AvgOutPerfMaker():CObject(TYPE) {}; 

    void validatePop2Object();

    DateTime        getFirstDate();
    DateTime        getLastDate();
    DateTimeArray   getAvgDates();

private:
    friend class AvgOutPerfMakerHelper;
    friend class AvgOutPerf;

     // data
    IntArray    callPut;
    DoubleArray participation;
    DoubleArray strike;

    CashFlowArray    avgSamples;         //Average Schedule and Ovserved Values

    // data, like Generic Performance
/*    string      genPerfType;      // C/P/F/S, or CS/PS, or FBD; or as above
    DoubleArray strikesPct;    // 1, 2 or 3 strikes
    double      participation;        
    */
};

typedef smartPtr<AvgOutPerfMaker> AvgOutPerfMakerSP;

//// internal average out class ////
class PRODUCTS_DLL AvgOutPerf{
public:
    AvgOutPerf():isAlreadySetInst(false) {}; 
    // constructor
    AvgOutPerf(AvgOutPerfMakerSP AvgOutPerfMaker);

    void setInstInfo(const InstrumentSettlementSP instSettle,
		     const InstrumentSettlementSP premiumSettle,
		     const CAsset*                asset,
		     string                       ccyTreatment,
		     const YieldCurveConstSP      discount,
		     const DateTime               valueDate);

//    roll();

    //void makePayoff(bool isCall, double strike);
    
    double AvgValue(double const strikeScale, 
                    double const refLevel, 
                    double const notl,
                    const DateTime& startDate) const;
    
    SampleListSP  avgOut;       // avarage out schedule and weight.  (Allow to access to SampleList)

    // access class
    IntArray    getCallPutArray();
    DoubleArray getParticipation();
    DoubleArray getStrikeArray();
    DateTimeArray   getAvgDates();
/*
    double		getAdjFwd(double const strikeScale,
                          double const refLevel,
                          double const notl,
                          const InstrumentSettlement* instSettle,
                          const InstrumentSettlement* premiumSettle,
                          const CAsset*               asset,
                          string                      ccyTreatment,
                          const YieldCurve*           discount,
                          const DateTime&             valueDate,
                          const SampleListSP          avgIn) const;
*/
private:
//    friend class AvgOutPerfHelper;
     // data
    IntArray    callPut;
    DoubleArray participation;
    DoubleArray strike;

    InstrumentSettlementSP instSettle;
    InstrumentSettlementSP premiumSettle;
    const CAsset*               asset;
    string                      ccyTreatment;
    YieldCurveConstSP   discount;
    DateTime             valueDate;
    bool                 isAlreadySetInst;  //flag whether the Inst info is set.
};

typedef refCountPtr<AvgOutPerf> AvgOutPerfSP;       // not smartPtr because not inheritate CObject

//// internal average out class ////
class PRODUCTS_DLL AvgIn: public CObject {
public:
    static CClassConstSP const TYPE;  
    // constructor
    AvgIn():CObject(TYPE) {}; 
    AvgIn(SampleListSP avgIn, const CAsset* asset);

    // return Fwd and ATM BS variance.
	void	getFwdAndBSVar(const DateTime valDate,int iSample,double &fwd,double &bsVar);

    // return expectation average in with known sStart and sEnd, using Brownian Bridge.
	double getMeanCondAvgIn(const DateTime valDate, const double sStart,const double sEnd,const double varEnd);
    

private:
    friend class AvgInHelper;
     // data
    SampleListSP    avgIn;
    CAssetWrapper   asset; // $unregistered
};

typedef smartPtr<AvgIn> AvgInSP;

/*
//// internal average out class ////
typedef smartPtr<AvgOutPerf> AvgOutPerfSP;
class PRODUCTS_DLL AvgOutPerf{
public:
    static CClassConstSP const TYPE;  
    AvgOutPerf():CObject(TYPE) {}; 
    // constructor
    static AvgOutPerfSP makeAvgOutPerf(AvgOutPerfMakerSP AvgOutPerfMaker);

    void validatePop2Object(){};
    //void makePayoff(bool isCall, double strike);
    
    double AvgValue(double const strikeScale, 
                   double const refLevel, 
                   double const notl,
                const InstrumentSettlement* instSettle,
                const InstrumentSettlement* premiumSettle,
                const CAsset*               asset,
                string                      ccyTreatment,
                const YieldCurve*           discount,
                const DateTime&             valueDate,
                bool                        fwdStarting,
                const DateTime&             startDate,
                bool                        oneContract) const;
    
    SampleListSP  avgOut;       // avarage out schedule and weight.  (Allow to access to SampleList)

private:
    friend class AvgOutPerfHelper;

    // data
    IntArray    callPut;
    DoubleArray participation;
    DoubleArray strike;
//    string      genPerfType;      // C/P/F/S, or CS/PS, or FBD; or as above
//    DoubleArray strikesPct;    // 1, 2 or 3 strikes
//    double      participation;            
//    
//    // data for internal use
//    DoubleArray strikes;        // strike level, to calc Fwd, it could be duplicate.
//    BoolArray   isCalls;         // call or not array along strikes.
//    DoubleArray longShort;      // Long or short array along strikes.
//};
*/

struct IAverageSpot
{
    virtual double calcAvgVariance( double strike ) const = 0;
    virtual void priceLN(
        Control* control,
        CResults* results,
        double strike,
        double avgVariance ) const = 0;
};

DRLIB_END_NAMESPACE
#endif
