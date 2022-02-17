//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWAverage.hpp
//
//   Description   averaging classes for equity swap.
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ESW_AVERAGE_HPP
#define EDG_ESW_AVERAGE_HPP

#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ESWAvIn : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultESWAvIn(){
        return new ESWAvIn();
    }

    ESWAvIn();

    typedef enum {NONE=0,USE_FWDS=1,USE_SPOT=2} ConsolidationType;
    ConsolidationType consolidationType;  // $unregistered

    DateTimeArray   sampleDates;
	DoubleArray		sampleLevels;
	DoubleArray		weights;
	DateTimeArray	liborFixingDates;
	DoubleArray		liborFixing;
    CStringArray    rateTypes;
//	DoubleArray		spreads;
	DoubleArray		fxFixing;
	DoubleArray		fwdLevels;
	bool			isConsolidated;
	DateTime		consoliDate;

    // unregistered
    string      ccyTreatment;


	// numSharesPerDollar: returns the number of shares, given notional = 1.  
    // should be used for notionally specified swaps.
	// i.e., E[1/sum(S_i)] where 
	//     i runs from 1 to numSamples, if so far = false.
	//     i runs from 1 to currentSample, if so far = true.
    // soFar = false should be used for computing number of shares for equity leg.
    // soFar = true should be used for computing number of shares for dividend and libor legs
    double numSharesPerDollar(
        const DateTime& valueDate, 
        CAssetSP        asset, 
        const DateTime& refDate = DateTime(0,0), 
        bool            soFar = false
        ) const;

	double numDollarsPerShare(
        const DateTime& valueDate, 
        CAssetSP        asset, 
        const DateTime& refDate = DateTime(0,0), 
        bool            soFar = false
        ) const;

    int getLength() const;
	
    double getLevel(const DateTime& valueDate, int iAveSample, CAssetSP asset) const;

	void validate( const DateTime& beginDate, 
                   const DateTime& endDate,
                   const DateTime& lastEquityAccrueStartDate,
                   const DateTime& valDate,
				   bool			   isFloating) const;


    void setFixing(const DateTime& valDate,
                   const DateTime& rollDate,
                     CAssetWrapper assetWrapper,
	                 const  Theta* shift,
                     string        ccyTreatment);

	double fractionOfShares(const DateTime& valueDate,
                            CAssetSP        asset,
                            const DateTime& refDate) const;
	void init();
	void validateHistoricalSamples(const DateTime& valDate, bool isFloating) const;

};
typedef smartPtr<ESWAvIn> ESWAvInSP;

class MARKET_DLL ESWAvOut : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultESWAvOut(){
    return new ESWAvOut();
    }

    typedef enum {NONE=0,LIBOR_WITHOUT_SPREAD=1,LIBOR_WITH_SPREAD=2,FIXED=3} accrualType;

	DateTimeArray   sampleDates;
	DoubleArray		sampleLevels;
	DoubleArray		weights;
	DoubleArray		discountFactors;
	DoubleArray		fxFixing;
	DateTimeArray	payDates;
    CStringArray    rateTypes; // $unregistered
    bool            accrueIRonFullNotional;
    int             accrueAvOutType;

    // unregistered
    DoubleArray eqGrowthFactor; // $unregistered
    string      ccyTreatment;

    ESWAvOut();
	int getLength() const;
	double getLevel(const DateTime& valueDate, int iAvSample, CAssetSP asset) const;
	void validate(const DateTime& beginDate, const DateTime& endDate, const DateTime& valDate) const;
	void setFixing(const DateTime& valDate,const DateTime& rollDate,
		CAssetWrapper assetWrapper,const Theta* shift,string ccyTreatment);
	
	double numSharesPerDollar(const DateTime& valueDate,CAssetSP asset,const DateTime& refDate,bool soFar) const;
	double numDollarsPerShare(const DateTime& valueDate,CAssetSP asset,const DateTime& refDate,bool soFar) const;
	double fractionOfShares(const DateTime& valueDate,CAssetSP asset,const DateTime& refDate) const;
    void init();

	void validateHistoricalSamples(const DateTime& valDate);
};
typedef smartPtr<ESWAvOut> ESWAvOutSP;

DRLIB_END_NAMESPACE
#endif
