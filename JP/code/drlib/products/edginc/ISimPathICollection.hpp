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

#ifndef EDR_ISIMPATHI_COLLECTION_HPP
#define EDR_ISIMPATHI_COLLECTION_HPP

#include "edginc/Class.hpp"
#include "edginc/StateVariableGen.hpp"
#include "edginc/MCWriter.hpp"
//#include "edginc/SVGenExpectedDiscFactor.hpp"
//#include "edginc/SVGenDiscFactor.hpp"
//#include "edginc/SVGenSurvivalDiscFactor.hpp"
//#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
// #include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"

// #include "edginc/FXAsset.hpp"
#include "edginc/SVGenSpot.hpp"

//#include "edginc/SVGenExpectedEnergyFuture.hpp"
// //#include "edginc/SVGenExpectedBasisForward.hpp"
// #include "edginc/RegressionTest.hpp" // declares IReportResults
// #include "edginc/IReportResults.hpp"

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
using namespace RM_Assets;

FORWARD_DECLARE( SimpleEquity );
FORWARD_DECLARE( YieldCurve );
FORWARD_DECLARE( IStateVariableCollector );

FORWARD_DECLARE( SVGenExpectedDiscFactor );
// FORWARD_DECLARE( SVGenDiscFactor );
class SVGenDiscFactor;
typedef smartPtr<SVGenDiscFactor> SVGenDiscFactorSP;
FORWARD_DECLARE( SVDiscFactor );
FORWARD_DECLARE( SVExpectedDiscFactor );

// FORWARD_DECLARE( SVGenSurvivalDiscFactor );
class SVGenSurvivalDiscFactor;
typedef smartPtr<SVGenSurvivalDiscFactor> SVGenSurvivalDiscFactorSP;

FORWARD_DECLARE( SVGenExpectedSurvivalDiscFactor );
FORWARD_DECLARE( SVGenAggregatedSurvivalDiscFactor );

FORWARD_DECLARE( SVSurvivalDiscFactor );
FORWARD_DECLARE( SVExpSurvDiscFactor );
FORWARD_DECLARE( SVAggregatedSurvDiscFactor );

FORWARD_DECLARE( SVGenExpectedEnergyFuture );
FORWARD_DECLARE( SVExpEnergyFuture );

FORWARD_DECLARE( SVGenExpectedBasisFwdSpread );
FORWARD_DECLARE( SVExpectedBasisFwdSpread );

class IRNGSlot {
public:
    virtual int getRngSlot() const = 0; // return assets's RNG slot in SAMPRAS convention
};

class ISimPathICollection : public virtual IObject,
                            public virtual IRNGSlot
{
public:

    virtual void initAssetInfo( RM_Assets::MCWriterDataHolder& writerData ) = 0;
    virtual void initSV( IStateVariableGen::IStateGen* pathGen ) = 0;
    virtual void updateDiffusedAssetInfo( int iStep, const DateTime& date ) = 0;
    virtual void appendSVGen( IStateVariableCollectorSP svCollector ) = 0;

    //TODO: later should be done properly -- all FX accessed independently,
    //but for now... (note that both Fx and Eq are here, disguised as FX)
    virtual void setFxEqGen( SVGenSpotSP fx )
    { }
    virtual void setFxEqSV( SVGenSpot::IStateVarSP _fx )
    { }
    virtual SVGenSpotSP getFxEqGen()
    {
        return SVGenSpotSP();
    }
    virtual IAssetWritableSP getAssetInfo() const = 0;
    // virtual void getRngSlot() const will be overriden via dominance in SimpathICurveInfo
};

//typedef refCountPtr<ISimPathICollection> ISimPathICollectionSP;
DECLARE( ISimPathICollection );
// the holder for the whole collection
typedef map<string, ISimPathICollectionSP > SimPathICollectionMap;
typedef map<string, ISimPathICollectionSP >::iterator SimPathICollectionMapIterator;

struct SimpathICurveInfo : public CObject,
                            public virtual IRNGSlot
{
    static CClassConstSP const TYPE;

    string name;
    DateTime maxDiffDate;
    DateTime maxCurveMaturity;
    
    virtual int getRngSlot() const {return rngSlot;}
    
private:
    int rngSlot;  // location in random number generator

    SimpathICurveInfo();
    static IObject* defaultSimpathICurveInfo();
    static void load( CClassSP& clazz );
};

DECLARE( SimpathICurveInfo );
typedef map<string, SimpathICurveInfoSP> SimpathICurveInfoMap; // FIXME: map should be to SP ?

struct SimpathICurrencyCollection : public SimpathICurveInfo, public ISimPathICollection
{
    // constructor
    SimpathICurrencyCollection(const SimpathICurveInfo& base, const vector<long>& offsets, const YieldCurve* y, const FXAsset* f)
        : SimpathICurveInfo(base), vOffsets(offsets), fxPosition(0), yc(y), fx(f) {}

    // SV generators
    SVGenSpotSP fxGen;
    SVGenDiscFactorSP dfGen;
    vector<SVGenExpectedDiscFactorSP> vExpDFGen;

    // SV themselves
    SVGenSpot ::IStateVarSP fxSV;
    SVDiscFactorSP dfSV;
    vector<SVExpectedDiscFactorSP> vExpDFSV;

    // these might be name-specific, say G7 are detailed, the rest - sparse
    vector<long> vOffsets;

    RM_Assets::AssetInfoWritable_CurrencySP ptrAssetInfo;
    virtual IAssetWritableSP getAssetInfo() const {
        return ptrAssetInfo;
    }
    
    
    void initAssetInfo( RM_Assets::MCWriterDataHolder& writerData );
    void initSV( IStateVariableGen::IStateGen* pathGen );

    void updateDiffusedAssetInfo( int iStep, const DateTime& date );
    void appendSVGen( IStateVariableCollectorSP svCollector );


    //TODO: later should be done properly -- all FX accessed independently,
    //but for now...
    // temporary solution to the lack of individual fxSV
    long fxPosition;

    void setFxEqGen( SVGenSpotSP _fx )
    {
        fxGen = _fx;
    }
    void setFxEqSV( SVGenSpot::IStateVarSP _fx )
    {
        fxSV = _fx;
    }
    SVGenSpotSP getFxEqGen()
    {
        return fxGen;
    }

    const YieldCurve*   yc;
    const FXAsset* fx;
};

//typedef refCountPtr<SimpathICurrencyCollection> SimpathICurrencyCollectionSP;
DECLARE( SimpathICurrencyCollection );
typedef map<string, SimpathICurrencyCollectionSP> SimPathICurrencyMap;

struct SimpathICreditCollection : public SimpathICurveInfo, public ISimPathICollection
{
    // constructor
    SimpathICreditCollection(const SimpathICurveInfo& base, const vector<long>& offsets, const ICDSParSpreads* a)
        : SimpathICurveInfo(base), vOffsets(offsets), cds(a) {}

    // SV generators
    SVGenSurvivalDiscFactorSP sdfGen;
    vector<SVGenExpectedSurvivalDiscFactorSP> vExpSDFGen;
    SVGenAggregatedSurvivalDiscFactorSP aggSDFGen;

    // SV themselves
    SVSurvivalDiscFactorSP sdfSV;
    vector<SVExpSurvDiscFactorSP> vExpSDFSV;
    SVAggregatedSurvDiscFactorSP aggSDFSV;
    vector<long> vOffsets;

    RM_Assets::AssetInfoWritable_CreditSP ptrAssetInfo;
    virtual IAssetWritableSP getAssetInfo() const {
        return ptrAssetInfo;
    }
    
    
    void initAssetInfo( RM_Assets::MCWriterDataHolder& writerData );
    void initSV( IStateVariableGen::IStateGen* pathGen );
    const ICDSParSpreads* cds;

    void updateDiffusedAssetInfo( int iStep, const DateTime& date );
    void appendSVGen( IStateVariableCollectorSP svCollector );
};

//typedef refCountPtr<SimpathICreditCollection> SimpathICreditCollectionSP;
DECLARE( SimpathICreditCollection );
typedef map<string, SimpathICreditCollectionSP> SimPathICreditMap;


// Energy stuffs
struct SimpathIEnergyCollection : public SimpathICurveInfo, public ISimPathICollection
{
    // constructor
    SimpathIEnergyCollection(const SimpathICurveInfo& base, const vector<TDate>& dates, const EnergyFuturesCurve* en)
        : SimpathICurveInfo(base), vDates(dates), enrg(en) {}

    // SV generators
    vector<SVGenExpectedEnergyFutureSP> vExpFpGen;

    // SV themselves
    vector<SVExpEnergyFutureSP> vExpFpSV;

    vector<TDate> vDates;

	RM_Assets::AssetInfoWritable_EnergySP ptrAssetInfo;
    virtual IAssetWritableSP getAssetInfo() const {
        return ptrAssetInfo;
    }
    
    
    void initAssetInfo( RM_Assets::MCWriterDataHolder& writerData );
    void initSV( IStateVariableGen::IStateGen* pathGen );
    const EnergyFuturesCurve* enrg;

    void updateDiffusedAssetInfo( int iStep, const DateTime& date );
    void appendSVGen( IStateVariableCollectorSP svCollector );
};

/*typedef refCountPtr<SimpathIEnergyCollection> SimpathIEnergyCollectionSP;*/
DECLARE( SimpathIEnergyCollection );
typedef map<string, SimpathIEnergyCollectionSP> SimPathIEnergyMap;


// Basis collection
struct SimpathIBasisCollection : public SimpathICurveInfo, public ISimPathICollection
{
    // constructor
    SimpathIBasisCollection(const SimpathICurveInfo& base, const vector<long>& offsets, const IBasisIndexCurve* b)
        : SimpathICurveInfo(base), vOffsets(offsets) /*, fxPosition(0)*/ , basis(b) {}

    // SV generators
    vector<SVGenExpectedBasisFwdSpreadSP> vExpBasisFwdSprdGen;

    // SV themselves
    vector<SVExpectedBasisFwdSpreadSP> vExpBasisFwdSprdSV;

    // these might be name-specific, say G7 are detailed, the rest - sparse
    vector<long> vOffsets;

    virtual IAssetWritableSP getAssetInfo() const {
        return ptrAssetInfo;
    }
    
    void initAssetInfo( RM_Assets::MCWriterDataHolder& writerData );
    void initSV( IStateVariableGen::IStateGen* pathGen );

    RM_Assets::AssetInfoWritable_BasisSP  ptrAssetInfo; // NOTE: reuse Rates MCWriter for basis
    const IBasisIndexCurve* basis;

    void updateDiffusedAssetInfo( int iStep, const DateTime& date );
    void appendSVGen( IStateVariableCollectorSP svCollector );
};

DECLARE( SimpathIBasisCollection );
typedef map<string, SimpathIBasisCollectionSP> SimPathIBasisMap;


struct SimpathISpotInfo :   public CObject,
                            public virtual IRNGSlot
{
    static CClassConstSP const TYPE;

    string name;
    DateTime maxDiffDate;
    virtual int getRngSlot() const {return rngSlot;}
 
private:
    int rngSlot;  // location in random number generator

    SimpathISpotInfo();
    static IObject* defaultSimpathISpotInfo();
    static void load( CClassSP& clazz );
};
DECLARE( SimpathISpotInfo );
typedef map<string, SimpathISpotInfoSP> SimpathISpotInfoMap;

struct SimpathIEquityCollection : public SimpathISpotInfo, public ISimPathICollection
{
    // constructor
    SimpathIEquityCollection(const SimpathISpotInfo& base, /*const string& ccy*/const SimpleEquity* a) : SimpathISpotInfo(base), eqPosition(-1), /*yc(ccy)*/ eq(a) {}

    //TODO: later should be done properly -- all EQ accessed independently,
    //but for now...
    // temporary solution to the lack of individual eqSV
    long eqPosition;

    void setFxEqGen( SVGenSpotSP _eq )
    {
        eqGen = _eq;
    }
    void setFxEqSV( SVGenSpot::IStateVarSP _eq )
    {
        eqSV = _eq;
    }
    SVGenSpotSP getFxEqGen()
    {
        return eqGen;
    }

    // SV generators
    SVGenSpotSP eqGen;

    // SV themselves
    SVGenSpot::IStateVarSP eqSV;

    RM_Assets::AssetInfoWritable_EquitySP ptrAssetInfo;

    //string                            yc;   // pay ccy
    const SimpleEquity*                 eq;   // to be hold in AssetInfoWritable_Equity for time 0 data output, and maybe calibrated data output.

    void initAssetInfo( RM_Assets::MCWriterDataHolder& writerData );
    void initSV( IStateVariableGen::IStateGen* pathGen );

    void updateDiffusedAssetInfo( int iStep, const DateTime& date );
    void appendSVGen( IStateVariableCollectorSP svCollector );
    virtual IAssetWritableSP getAssetInfo() const {
        return ptrAssetInfo;
    }
    
};

DECLARE( SimpathIEquityCollection );

typedef map<string, SimpathIEquityCollectionSP> SimPathIEquityMap;


DRLIB_END_NAMESPACE

#endif
