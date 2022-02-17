//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimpathIInstrument.hpp
//
//   Description : SimpathIInstrument
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SIMPATHIINSTRUMENT_MC_HPP
#define EDR_SIMPATHIINSTRUMENT_MC_HPP

#include "edginc/MCProductClient.hpp"
#include "edginc/IReportResults.hpp"
#include "edginc/DateTime.hpp" 
// #include "edginc/StateVariableGen.hpp"
#include "edginc/MCWriter.hpp"
#include "edginc/ISimPathICollection.hpp"
#include "edginc/IRngSlotAssigner.hpp"

DRLIB_BEGIN_NAMESPACE

class IMultiMarketFactors;
class YieldCurve;
class MCPathGenerator;
class IMCPrices;
class IStateVariableGen;
class IStateVariableGen::IStateGen;

FORWARD_DECLARE( IStateVariableCollector );
FORWARD_DECLARE( IStateVariable );

class SimpathIInstrumentMC : public MCProductClient,
							 virtual public IReportResults,
							 virtual public IRngSlotAssigner
{
public:
    SimpathIInstrumentMC(
        const IMultiMarketFactors* mFactors,
        const DateTime& today,
        const YieldCurve* discount,
        const SimPathICollectionMap& collectionMap,
        const DateTimeArray& _timeLine,
        IMCWriterSP _mcWriter,
        const map<string, int>& asset2rng,
		bool  doPrices);

    //    void validatePop2Object(){}

    virtual void pathGenUpdated( IStateVariableGen::IStateGen* newPathGen );

    virtual void collectStateVars( IStateVariableCollectorSP svCollector ) const;

    virtual void payoff(
        const MCPathGenerator* pathGen,
        IMCPrices& prices );

    RM_Assets::IMCWriterSP getWriter() const;
    virtual IObjectSP getResults();
    virtual IMCPrices* createOrigPrices( int nbIter,
                                         int nbSubSamples,
                                         int mode );

    virtual void finalize();
	virtual int getRngIdx(int globalIdx, const string& assetName) const; // assigns RNG slots based on information in the collectionMap.
private:
    // Store pointers to information owned by the SimpathIInstrument
    vector<IStateVariableGen*> svGenSet;
    mutable SimPathICollectionMap collectionMap;

    //    const SimpathICurveInfoMap* swapCurveInfoMap;
    //    const SimpathICurveInfoMap* cdsCurveInfoMap;
    //    const SimpathISpotInfoMap* fxInfoMap;
    //    const SimpathISpotInfoMap* eqInfoMap;

    const DateTimeArray timeLine;

    // Store SP for heap allocated information owned by this class
    vector<IStateVariableSP> svSet;

    // lnDiscountFactor[tpIdx] -> list of ir assets diffused on this date
    // lnDiscountFactor[tpIdx][assetIdx] -> list of ln discount factors
    // lnDiscountFactor[tpIdx][assetIdx][0] -> discount factor from today to timeLine[tpIdx]
    // lnDiscountFactor[tpIdx][assetIdx][1] -> expected discount factor from timeLine[tpIdx] to first CDS date
    //   ...
    // lnDiscountFactor[tpIdx][assetIdx][n] -> expected discount factor from timeLine[tpIdx] to nth CDS date
    vector<vector<vector<double> > > lnDiscountFactor;

    // lnSurvivalDiscountFactor[tpIdx][assetIdx][pillarIdx] with similar meaning
    // as lnDiscountFactor (just substitute "survival discount" for "discount" above.
    vector<vector<vector<double> > > lnSurvivalDiscountFactor;

    // fxRate[tpIdx][assetIdx] (DO WE HAVE PILLAR DATES OR EXPECTED FXRATES?)
    vector<vector<double> > fxRate;


    /// for writing information to files/streams/etc
    RM_Assets::IMCWriterSP mcWriter; // represents one or a chain of writers that react on events like notifyStartPath

	// The following governs if we want to return diffusion results back to Qlib.
	// For SAMPRAS/SimDesk we definitely do not want to do this as it consumes a lot of memory and time.
	bool doPrices; /// when set, product will return results of diffusion as a 3D array [path][timepoint][asset]->diffusion_numbers; it is disabled when at least one writer is requested (BINARY|DEBUG)
    
    map<string, int> asset2rng; // map that specifies for each asset its RNG slot, so that adding or removing assets does not change random numbers of other ones.
};

DRLIB_END_NAMESPACE

#endif
