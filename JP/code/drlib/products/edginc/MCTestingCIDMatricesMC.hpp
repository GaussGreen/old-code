#ifndef MC_TESTING_CID_MATRICES_MC_HPP
#define MC_TESTING_CID_MATRICES_MC_HPP

#include "edginc/MCProductClient.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenPathWeight.hpp"
#include "edginc/FXAsset.hpp"
/*
  Monte Carlo implementation of MCTestingCIDMatrices
*/

DRLIB_BEGIN_NAMESPACE

class MCTestingCIDMatrices;

/* MC view of the MCTestingCIDMatrices */
class PRODUCTS_DLL MCTestingCIDMatricesMC :   public MCProductClient
/*, public IMCStatelessProductClient*/
{
private:
	// SVs:
	//    SVDateOfDefaultSP dd; // date of default
	vector<SVExpSurvDiscFactorSP> esdf; // expected survival prob
	// FIXME: not implemented    SVRecoveryRateSP rr; // recovery rate
	SVDiscFactorSP domDf; // domestic DFs
	SVDiscFactorSP forDf; // foreign DFs
	
	// SV Gens:
	//    SVGenDateOfDefaultSP ddGen;
	vector<SVGenExpectedSurvivalDiscFactorSP> esdfGen;
	// FIXME    SVGenRecoveryRateSP rrGen;
	SVGenDiscFactorSP domDfGen;
	SVGenDiscFactorSP forDfGen;
	
	MCPath::IStateVarSP fx;
	SVGenSpotSP fxGen;
  
    SVPathWeightSP  pathW;
    SVGenPathWeightSP  pathWGen;

private:
    vector<YieldCurveConstSP> ycAsset;
    vector<ICDSParSpreadsConstSP> cdsAsset;
    vector<FXAssetConstSP> fxAsset;
    map<string,int> fxAssetIDMap;
    DateTime today;
    DateTimeArray maturityDates;
    DateTime measureDate;
    bool doLog;

	int nbAssets;
    int nbMaturityDates;
    bool needFx;
	
protected:
	/** Override default method on IMCProduct. This method is called every time
		the path generator is changed (which is, at the moment, when the
		past path generator is created, and then when the future path
		generator is created  */
	virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);
	
public:
	/** Appends 'true' (ie non derived) state variable generators
		required to the supplied collector.*/
	virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;
	
	/** equivalent to InstIntoMCProduct. Need to call parent's constructor */
//   MCTestingCIDMatricesMC( const MCTestingCIDMatrices* inst,
// 						  const SimSeriesSP& simSeries,
// 						  InstrumentSettlementSP instSettle);
	MCTestingCIDMatricesMC(IMultiMarketFactorsSP assets,
						   const vector<YieldCurveConstSP>& ycAsset,
						   const vector<ICDSParSpreadsConstSP>& cdsAsset,
						   const vector<FXAssetConstSP>& fxAsset,
						   const map<string,int>& fxAssetIDMap,
						   const YieldCurve* discount,
						   const DateTime& today,
						   const DateTimeArray& maturityDates,
						   const DateTime& measureDate,
						   bool doLog);
						 
    /** Called within the simulation loop */
	virtual void payoff(
		const IPathGenerator* pathGen,
		IMCPrices& prices);

	IMCPrices* createOrigPrices(int  nbIter,
								int  nbSubSamples,
								int  mode);

	void recordExtraOutput(CControl*     control,
						   Results*      results,
						   const IMCPrices& prices) const;

};

DRLIB_END_NAMESPACE
#endif
