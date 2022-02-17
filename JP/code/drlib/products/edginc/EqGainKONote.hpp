//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EqGainKONote.hpp
//
//   Description   Equity Gain KO Note.  (Equity Turbo without Callable, but KO) 
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EQGAIN_KO_NOTE_HPP
#define EDG_EQGAIN_KO_NOTE_HPP

#include "edginc/DblBarrier.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/SwapLegIntFace.hpp"

DRLIB_BEGIN_NAMESPACE

class CEqGainKONote : public CDblBarrier,
                          virtual public BarrierBreach::IEventHandler,
                          virtual public KnownCashflows::IEventHandler 
{
public:
    static CClassConstSP const TYPE;

    virtual void Validate();
 
    FDProductSP createProduct(FDModel* model) const;

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

    CashFlowArraySP getKnownCashFlow() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    bool priceDeadInstrument(CControl* control, CResults* results) const;
    
    bool sensShift(Theta* shift);

    // the following internal class, AvgOutPerf and AverageOut are
    // proto-type of AvgOutPerf and AvgOutPerfMaker in average.hpp.
    // It's same name, but totally independet to them.

    //// data class as external interface ////
    class AvgOutPerf: public CObject {
    public:
        static CClassConstSP const TYPE;  

         // data, like Generic Performance
        string      genPerfType;      // C/P/F/S, or CS/PS, or FBD; or as above
        DoubleArray strikesPct;    // 1, 2 or 3 strikes
        double      participation;        
        CashFlowArray    avgSamples;         //Average Schedule and Ovserved Values

        AvgOutPerf():CObject(TYPE) {}; 
   
        void validatePop2Object();

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz){
            clazz->setPublic(); // make visible to EAS/spreadsheet
            clazz->setDescription("Average Out Peformance");
            REGISTER(CEqGainKONote::AvgOutPerf, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultAvgOutPerf);
            FIELD(genPerfType,"Generic Peformance. F,C,P,S,CS,PS,BDF");
            FIELD(strikesPct, "strike array");
            FIELD(participation, "participation")
            FIELD(avgSamples, "Avarage Sampling Schedule and observed level (if in past)");
        }
        
        static IObject* defaultAvgOutPerf(){
            return new CEqGainKONote::AvgOutPerf();
        }

    };

    typedef smartPtr<AvgOutPerf> AvgOutPerfSP;

    //// internal average out class ////
    class AverageOut{
    public:
        string      genPerfType;      // C/P/F/S, or CS/PS, or FBD; or as above
        DoubleArray strikesPct;    // 1, 2 or 3 strikes
        double      participation;        
        
        SampleListSP  avgOut;       // avarage out schedule and weight.
        
        // data for internal use
        DoubleArray strikes;        // strike level, to calc Fwd, it could be duplicate.
        BoolArray   isCalls;         // call or not array along strikes.
        DoubleArray longShort;      // Long or short array along strikes.
        
        // constructor
        AverageOut(AvgOutPerfSP avgOutPerf);

        void makePayoff(bool isCall, double strike, bool isLong){
            strikes.push_back(strike);
            isCalls.push_back(isCall);
            longShort.push_back(isLong? 1.0:-1.0);
        }
    };

    typedef refCountPtr<AverageOut> AverageOutSP;
    
    // It takes AvgOutPerf instance, not use the own avgOutPerf, which is public data.
    // It's safer to make a copy and modify the internal class.
    double avgValue(const AverageOutSP aop,
                    double const spot, 
                    double const SpotRef, 
                    const DateTime& stepDate) const;

    double avgValue(const AverageOutSP aop,
                    double const strikeScale, 
                    double const refLevel,
                    double const notl,
                    const DateTime& when) const ;

    //------------------------------------------//    
    // for event handling
    //------------------------------------------//    
    // KnownCashflows::IEventHandler interface
    void getEvents(const KnownCashflows* flows,
                   IModel* model, 
                   const DateTime& eventDate,
                   EventResults* events) const ;

    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach,
                   IModel* model, 
                   const DateTime& eventDate,
                   EventResults* events) const; 


private:
    friend class CEqGainKONoteHelper;
    friend class CEqGainKONoteFDProd;
    
    double          SpotFixing;

    // STB member data
    LiborLegSP liborLeg;
    FixedLegSP  fixedLeg;

    ScheduleSP  ParticipationSchedule;  //Participation for Equity Return
    DoubleArray CapStrikeArray;      //Strike Level for Cap Level

    ScheduleSP  CallSchedule;  //Call Schedule with strike level.
    bool        isCallable;     //Flag for Callable.  Default = false

    DoubleArray AddStrikeArray1;      //Additional Optin Strike 1.
    DoubleArray AddStrikeArray2;      //Additional Optin Strike 2.
    DoubleArray PartStrike1;            //Participation for Additional Option 1
    DoubleArray PartStrike2;            //Participation for Additional Option 2
    string      PayType1;               //Payoff type for Additional Strike 1
    string      PayType2;               //Payoff type for Additional Strike 2

    SampleListSP  pastValues;          //Historic values (Do not use!!)
    CashFlowArray    samples;         //Historic values Interface(Use this).
    
    bool        isCalled;               // if it's already called, true.
    DateTime    callDate;               // if it's already called, needs the date.

    int     segDensity;             // optional segment density to set more grid before next CritDates.
    
    bool    SwitchPlainSwap;         //Price with PlainSwap or Not.  Default = True (with PlainSwap)
    bool    RemoveSettlementOption;  //default = false, true shows only coupon part.
    
    double  scalingStrike;           //scaling factor for strike levels for fwd starting case.


    bool          hasAvgOut;          // has AverageOut Option?
    AvgOutPerfSP  avgOutPerf;         // Option at Maturity with Average Out.

protected:
    static void load(CClassSP& clazz);
    CEqGainKONote();
    CEqGainKONote(CClassConstSP clazz);

};


DRLIB_END_NAMESPACE
#endif
