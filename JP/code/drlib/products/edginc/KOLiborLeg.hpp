//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KOLiborLeg.cpp
//
//   Description   class for KO Libor Leg 
//
//
//   $Log: KOLiborLeg.cpp,v $
//----------------------------------------------------------------------------
#ifndef EDR_KO_LIBOR_LEG_HPP
#define EDR_KO_LIBOR_LEG_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/KOStubRule.hpp"

DRLIB_BEGIN_NAMESPACE
class KOLiborLeg;

class PRODUCTS_DLL KOLiborLegMaker: public CObject {
public:
    static CClassConstSP const TYPE;  
    friend class KOLiborLeg;

    KOLiborLegMaker(); 
    KOLiborLegMaker(const LiborLegSP floater, const KOStubRuleSP koRule);

    void validatePop2Object() ;

    //make a new copy with  ko settlement rule
    KOLiborLeg* makeKOLiborLeg(const InstrumentSettlementSP settle,
                               const DateTime valDate) const;

    // retrun payment dates in schedule
	DateTimeArray  getPayDates() const;

    // get Market for LigborLeg
    void getMarket(const IModel*     model, 
                   const MarketData* market,
                   const YieldCurveWrapper discount);

    // tweak theata and feed fixing level....
    void setFixingforThetaShift(const DateTime& valueDate, 
                              const YieldCurve* discount,
                              const DateTime& rollDate);

private:
    // Invoked when Class is 'loaded' 
    static void load(CClassSP& clazz);    
    static IObject* defaultKOLiborLegMaker();

    LiborLegSP          floater;
    KOStubRuleSP   koRule;

};

DECLARE(KOLiborLegMaker)

class PRODUCTS_DLL KOLiborLeg{
public:
    KOLiborLeg(const KOLiborLegMaker* maker,
               const InstrumentSettlementSP settle,
               const DateTime valDate);

    // constructor without maker class.
    KOLiborLeg(const LiborLegSP flt,
               const KOStubRuleSP koRule,
               const InstrumentSettlementSP settle,
               const DateTime valDate,
               const YieldCurveWrapper discount);

    // constructor, different settle schedule from floater's schedule.
    KOLiborLeg(LiborLegSP flt,
               KOStubRuleSP koRule,
               InstrumentSettlementSP settle,
               const DateTimeArray& payDates,
               const DateTimeArray& accrueDates,
               DateTime valDate,
               const YieldCurve* discount);

    // no settle version....
    KOLiborLeg(const LiborLegSP flt,
               const KOStubRuleSP koStub,
               const DateTime valDate,
               const YieldCurve* discount);

    KOLiborLeg(const LiborLegSP flt,
               const KOSettleSP koSettle,
               const DateTime valDate,
               const YieldCurve* discount);

    // return the critical dates.
    DateTimeArray getCritDates() const;

    // return AccrueDates.
    DateTimeArray getAccrueDates() const;

    CashFlowArraySP  makeKnownCashFlows(const DateTime hitDate, const bool isAlreadyHit);

	CashFlowArraySP	 getKnownCashFlows();

    // return simple cash flow array
	CashFlowArrayConstSP	 getCashFlows(const DateTime valDate, const YieldCurve* yc);

    // return cash flow when KO occur on hitDate
	CashFlow	  getKOCashFlow(const DateTime hitDate) const;

    // return KO value as of hitDate
    double        getKOValue(const DateTime hitDate, const YieldCurve* yc);

    // return non-KO (plain) value as of hitDate
	double		  getPV(const DateTime hitDate,const YieldCurve* yc);

    // is this necessary???
	double			getKOPV(const DateTime&   hitDate);

    // tweak spreads.  isTweak = false means roll the spreads back.
    void    tweakSpreads(bool isTweak);

protected:
    LiborLegSP      floater;
    KOSettleSP      koSettle;    
    DateTime        valDate;

    CashFlowArraySP     knownCFL;          // knownCFL

public:
    DoubleArray		getPVsAlongKOTimeStep(const DateTimeArray&    koDates,
                                          const double            scaling);

	CashFlowArraySP	getKOCashFlowArray(const DateTime& hitDate,
                                        const YieldCurve* discount,
                                        const string koStubRule,
                                        const bool isReverse);

private:
	bool			getCashFlowOnKO(const DateTime hitDate,
                                    const string koStubRule,
                                    double* pvValue,
                                    CashFlowArray*  coupon);

	double			getKOPV(const DateTime&   hitDate,
                            const string koStubRule);

    double			getKOPV(const DateTime&   baseDate,
                            const DateTime&   hitDate,
                            const YieldCurve* discount,
                            const string koStubRule);

//	void			setKOSettle(const InstrumentSettlementSP settle,
//                                const DateTimeArray payDateArray,
//                                const DateTimeArray accrueDates,
//                                const bool isRollingSettle,
//                                const bool isAccrueUpToSettle);

};

typedef refCountPtr<KOLiborLeg> KOLiborLegSP;

class KOLiborLegSV : public KOLiborLeg{
public:
    // constructor without maker class.
    KOLiborLegSV(const LiborLeg::LiborLegSVSP flt,
                const KOSettleSP koSettle,
                const DateTime valDate,
               const YieldCurve* discount);

    // return the PV of Cashflow from today to hitDate, 
    // with consideration about the KO settlement treatment (i.e. B or S or N)
	double			getKOPV(const DateTime&   hitDate);

private:
    LiborLeg::LiborLegSVSP      floaterSV;
};

typedef refCountPtr<KOLiborLegSV> KOLiborLegSVSP;

DRLIB_END_NAMESPACE

#endif
