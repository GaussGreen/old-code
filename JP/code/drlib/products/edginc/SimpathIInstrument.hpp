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

#ifndef EDR_SIMPATHIINST_HPP
#define EDR_SIMPATHIINST_HPP

#include "edginc/Class.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/MCWriter.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DECLARE.hpp"

#include "edginc/ISimPathICollection.hpp"

DRLIB_BEGIN_NAMESPACE
using namespace RM_Assets;

FORWARD_DECLARE(SimpleEquity);
FORWARD_DECLARE(YieldCurve);
FORWARD_DECLARE(FXAsset);
FORWARD_DECLARE(ICDSParSpreads);
FORWARD_DECLARE(EnergyFuturesCurve);
FORWARD_DECLARE(IBasisIndexCurve);

class SimpathIInstrument: public GenericNFactor,
                          virtual public IMCIntoProduct {
public:
    static CClassConstSP const TYPE;

    IMCProduct* createProduct(const MonteCarlo* model) const;
//    IMCWriterSP getIMCWriter(void) const {return stateIMCWriter;}

    static void classifyAssets(
        IMultiMarketFactorsSP   assets, // (I)
        vector<YieldCurveConstSP>& ycAsset, // (O)
        vector<FXAssetConstSP>& fxAsset,// (O)
        vector<ICDSParSpreadsConstSP>& cdsAsset,// (O)
		vector<SimpleEquityConstSP>& eqAsset,// (O)
		vector<EnergyFuturesCurveConstSP>& enrgAsset,// (O)
        vector<IBasisIndexCurveConstSP>& basisAsset,// (O)
        map<string,int>& fxAssetIDMap // (O)
        );

	static DateTimeArray calcForwardDates(  const DateTime& currentDate,
											const StringArray&  offsets,
											const DateTime& maxDate);
	
	static DateTimeArray calcForwardDates(  const DateTime& currentDate,
											const vector<long>&  offsets,
											const DateTime& maxDate);
	

private:
    SimpathIInstrument();
    static IObject* defaultSimpathIInstrument();
    static void load(CClassSP& clazz);
    void validatePop2Object() {}; // CAN ONE FINALIZE AN OBJECT HERE?

    void initializeMaps() const; // NEEDED THIS SINCE RESULT OF validatePop2Object is lost.

    // The following insertXYZ methods are not used.
    // TODO: We should call these in createProduct() to reduce the size of it...
    void insertDFSVGen(
        const DateTimeArray& timeLine,
        const vector<YieldCurveConstSP>& ycAsset,
        vector<IStateVariableGen*>& gens ) const;

    void insertExpectedDFSVGen(
        const DateTime& currentDate,
        vector<YieldCurveConstSP>& ycAsset,
        vector<IStateVariableGen*>& gens ) const;

    void insertSurvivalDFSVGen(
        const DateTimeArray& timeLine,
        const vector<ICDSParSpreadsConstSP>& cdsAsset,
        vector<IStateVariableGen*>& gens ) const;

    void insertExpectedSurvivalDFSVGen(
        const DateTime& currentDate,
        const vector<ICDSParSpreadsConstSP>& cdsAsset,
        vector<IStateVariableGen*>& gens ) const;

	// note: no need to insert spot energy SV generator.  Only expected SVs are needed.
	void insertExpectedEnergyFutureSVGen(
		const DateTime& currentDate,
		const vector<EnergyFuturesCurveConstSP>& enrgAsset,
		vector<IStateVariableGen*>& gens ) const;

    void insertExpectedBasisFwdSprdSVGen(
        const DateTime& currentDate,
        const vector<IBasisIndexCurveConstSP>& basisAsset,
        vector<IStateVariableGen*>& gens ) const {};

    // Create decorated writer and report if we want to report results via Qlib's MCPrices mechanism as well.
	std::pair<IMCWriterSP, bool> createWriter(const MCWriterDataHolder& writableData, const std::vector<TDate>& dates) const;

    // Variables for construction of SimpathIInstrument
    StringArraySP timeLineOffsets;    // 1D, 2D, etc.
    StringArraySP yieldCurveOffsets;
    StringArraySP cdsOffsets;
    StringArraySP basisCurveOffsets;
	// Note: energy curve has a fixed set of maturities, don't use offsets!

    SimpathICurveInfoArraySP swapCurveInfo;
    SimpathICurveInfoArraySP cdsCurveInfo;
    SimpathICurveInfoArraySP enrgCurveInfo;
    SimpathICurveInfoArraySP basisCurveInfo;
    SimpathISpotInfoArraySP fxInfo;
    SimpathISpotInfoArraySP eqInfo;

    // Hold data in more convenient form
    mutable SimpathICurveInfoMap swapCurveInfoMap; // $unregistered
    mutable SimpathICurveInfoMap cdsCurveInfoMap; // $unregistered
	mutable SimpathICurveInfoMap enrgCurveInfoMap; // $unregistered
    mutable SimpathICurveInfoMap basisCurveInfoMap; // $unregistered
    mutable SimpathISpotInfoMap fxInfoMap; // $unregistered
    mutable SimpathISpotInfoMap eqInfoMap; // $unregistered

    // Pointer to the Writer class to look into the diffused Quantities
//    mutable IMCWriterSP  stateIMCWriter; // sorry for that, canoot assign SP without that $unregistered
    string       filenamePrefix; // filename for the XML output
    string       format; // XML/DEBUG/BINARY
    mutable map<string, int> asset2rng; // map [asset_name]->[rng slot number]
};

DECLARE(SimpathIInstrument);


DRLIB_END_NAMESPACE

#endif
