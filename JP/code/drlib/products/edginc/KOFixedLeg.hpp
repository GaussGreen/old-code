//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KOLiborLeg.cpp
//
//   Description   class for KO Fixed Leg 
//
//
//   $Log: KOLiborLeg.cpp,v $
//----------------------------------------------------------------------------
#ifndef EDR_KO_FIXED_LEG_HPP
#define EDR_KO_FIXED_LEG_HPP

#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/SwapLegIntFace.hpp"
#include "edginc/KOStubRule.hpp"

DRLIB_BEGIN_NAMESPACE

class KOFixedLeg;

class PRODUCTS_DLL KOFixedLegMaker: public CObject {
public:
    static CClassConstSP const TYPE;  
    friend class KOFixedLeg;

    KOFixedLegMaker(); 
    KOFixedLegMaker(const FixedLegSP fix, const KOStubRuleSP koRule);

    void validatePop2Object() ;

    //make a new copy with  ko settlement rule
    KOFixedLeg* makeKOFixedLeg(const InstrumentSettlementSP settle) const;


    // retrun payment dates in schedule
	DateTimeArray  getPayDates() const;

private:
    // Invoked when Class is 'loaded' 
    static void load(CClassSP& clazz);    
    static IObject* defaultKOFixedLegMaker();

    FixedLegSP fixedLeg;
    KOStubRuleSP   koRule;
    
};

typedef smartPtr<KOFixedLegMaker> KOFixedLegMakerSP;

class PRODUCTS_DLL KOFixedLeg{
public:
    KOFixedLeg(const KOFixedLegMaker* maker,
               const InstrumentSettlementSP settle);

    // constructor without maker.
    KOFixedLeg(const FixedLegSP fix,
               const KOStubRuleSP koRule,
               const InstrumentSettlementSP settle);

    // constructor without maker w/o settle
    KOFixedLeg(const FixedLegSP fix,
               const KOStubRuleSP koRule);

    // return the critical dates.
    DateTimeArray getCritDates() const;

    // return AccrueDates.
    DateTimeArray getAccrueDates() const;

    CashFlowArraySP  getKnownCashFlows(const DateTime hitDate);

    // return simple cash flow array
	CashFlowArrayConstSP	 getCashFlows();

    // return cash flow when KO occur on hitDate
	CashFlow	  getKOCashFlow(const DateTime hitDate);

    // return KO value as of hitDate
    double        getKOValue(const DateTime hitDate, const YieldCurve* yc);

    // return non-KO (plain) value as of hitDate
	double		  getPV(const DateTime hitDate,const YieldCurve* yc);

	double		  getKOPV(const DateTime&   hitDate,
                          const YieldCurve* discount);


private:
    KOFixedLegMakerSP   maker;
    FixedLegSP          fixedLeg;
    KOSettleSP          koSettle;
    bool                hasKOSettle;

    CashFlowArraySP     knownCFLs;         // knownCFL

public: //I'd like to hide this...
    CashFlowArraySP		getKOCashFlowArray(const DateTime& hitDate,
                                           const YieldCurve* discount,
                                           const string koStubRule,
                                           const bool isReverse);

private:
	bool				getCashFlowOnKO(const DateTime hitDate,
                                        const string koStubRule,
                                        CashFlow*  cfl);

};

typedef refCountPtr<KOFixedLeg> KOFixedLegSP;

DRLIB_END_NAMESPACE

#endif
