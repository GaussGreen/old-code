//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRM.hpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------

// #include <fstream>

#include "edginc/config.hpp"
#include "edginc/MCPathConfigSRM.hpp"
#include "edginc/MCPathConfigSRMGen.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MarketDataFetcherSRM.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/SRMEquityUtil.hpp" // TO DO - why not dependency on SRMFXUtil as well?
#include "edginc/MDFUtil.hpp"
#include "edginc/QMCGenDiffusibleAsset.hpp"
#include ext_hash_set

#include "edginc/CorrelationCategory.hpp"


DRLIB_BEGIN_NAMESPACE

START_PUBLIC_ENUM_DEFINITION(EquityDiffusionStyle, "Values for EquitiesDiffusionStyle");
    ENUM_VALUE_AND_NAME(eqDEFAULT, "DEFAULT", 
                        "The default equities diffusion: currently Euler");
    ENUM_VALUE_AND_NAME(eqEULER, "EULER",
                        "The Euler scheme");
    ENUM_VALUE_AND_NAME(eqSECOND_ORDER, "SECOND_ORDER",
                        "Second order scheme");
    ENUM_VALUE_AND_NAME(eqSTRICT_SECOND_ORDER, "STRICT_SECOND_ORDER",
                        "Second order scheme with possibly less simulation dates");
    ENUM_VALUE_AND_NAME(eqMAPPING_METHOD, "MAPPING_METHOD",
                        "Mapping method")
END_ENUM_DEFINITION(EquityDiffusionStyle);

///// static fields  ////

const string MCPathConfigSRM::CALIB_NIL("NIL"); // "NIL" - no calibration against swaptns vol
const string MCPathConfigSRM::ZC_GROWTH("GROWTH");// "GROWTH" - use 'growth' zero curve for diffusion
const string MCPathConfigSRM::ZC_DISCOUNT("DISCOUNT"); // "DISCOUNT" - use 'discount' instead
const string MCPathConfigSRM::SRM_MODEL_FAMILY_NAME("SRM");

const string MCPathConfigSRM::CREDIT_MODEL_TYPE_HJM("HJM");
const string MCPathConfigSRM::CREDIT_MODEL_TYPE_CIR("CIR");
const string MCPathConfigSRM::CREDIT_MODEL_TYPE_LIBOR("LIBOR");
const string MCPathConfigSRM::CREDIT_MODEL_TYPE_CID("CID");
const string MCPathConfigSRM::CREDIT_MODEL_TYPE_CID_FULL("CID_FULL");
const string MCPathConfigSRM::CREDIT_MODEL_TYPE_CID_FAST("CID_FAST");

const string MCPathConfigSRM::RATES_MODEL_TYPE_HJM("HJM");
const string MCPathConfigSRM::RATES_MODEL_TYPE_LIBOR("LIBOR");
const string MCPathConfigSRM::RATES_MODEL_TYPE_DETERM("DETERM");

const string MCPathConfigSRM::COMMON_MARKET_FACTOR_FOR_CID("CommonMarketFactorCID");

//// Special MDF class that routes call for getComponentMarketData through
//// path config
class MCPathConfigSRM::MDF: public MarketDataFetcherSRM {
public:
    static CClassConstSP const  TYPE;
    const MCPathConfigSRM*      pathConfig; // ref - NB can't be a SP

    /** Fill fetch stochastic yield curves plus type of smile and model
    specified. */
    MDF(const string&      irCalibSmileType, 
        const string&      irCalibModelType,
        bool               getSwaptionVols,
        bool               useIRVolPair,
        const string&      volType,
        const StringArray& fxVolType,
        const string&      crCalibSmileType,
        const string&      cdsVolType,
        const MCPathConfigSRM*   pathConfig):
        MarketDataFetcherSRM(irCalibSmileType, irCalibModelType,
                             getSwaptionVols, useIRVolPair, volType, fxVolType, 
                             crCalibSmileType, cdsVolType),
        pathConfig(pathConfig) {}

    /** Having retrieved a market object from the cache (or if
        presupplied) it is necessary to ensure that it has all its
        market data. Reroute to pathConfig */
    virtual void getComponentMarketData(const IModel*        model,
                                        const MarketData*    market,
                                        MarketObjectSP       mo) const{
        pathConfig->getComponentMarketData(model, market, mo, 
                                           const_cast<MDF*>(this)); ////////
    }
};

/** Create a MarketDataFetcher which will be used by the [MonteCarlo] model
    for retrieving market data etc */
MarketDataFetcherSP MCPathConfigSRM::marketDataFetcher() const 
{
    //See MCPathConfigSRM::getComponentMarketData for the case
    //where there is more than one element in numIRFactors.

    int nFactors = numIRFactors.front();
    string irCalibModelType;
    if (nFactors != 0) {
        irCalibModelType = IRCalib::getModelType(nFactors)->getName();
    } else { 
        //can't calibrate 0 factors, but need to retrieve something.
        irCalibModelType = IRCalib::getModelType(1)->getName();
    }
    MarketDataFetcherSP mdf(
        new MDF(
            "IRCalib::Smile2Q",                    // irCalibSmileType
            irCalibModelType,  
            calibrationStyle.front() != CALIB_NIL, // getSwaptionVols?
            useIRVolPair,
            getVolType(), 
            getFxVolType(),                        // fxVolType
            "CRCalib::Smile",                      // crCalibSmileType
            getCDSVolType(),                       // cdsVolType
            this));
    dependenceMaker->modifyMarketDataFetcher(mdf);    
    return mdf;
}


//// helper method: If isoCode empty returns 0 otherwise returns
//// the index in isoCode. Throws an exception if not found
int MCPathConfigSRM::findIsoCodeIndex(const string& isoCodeToFind) const
{
    if (isoCode.empty()){
        return 0;
    }
    for (int i = 0; i < isoCode.size(); i++){
        if (isoCode[i] == isoCodeToFind){
            return i;
        }
    }
    throw ModelException("MCPathConfigSRM::findIsoCodeIndex", "No data "
                         "specified for currency "+isoCodeToFind);
}

int MCPathConfigSRM::findEnergyIndex(const string & energyName) const
{
	for (int i = 0; i < energyNames.size(); ++i) {
		if (energyName == energyNames[i])
			return i;
	}
	throw ModelException("MCPathConfigSRM::findEnergyIndex", "Energy asset "
		"with name: " + energyName + " not found!");
}

string MCPathConfigSRM::cdsVolTypeFromCreditModelType(
    const string& _creditModelType) const
{
    static const string method("MCPathConfigSRM::cdsVolTypeFromCreditModelType");
    if (_creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_HJM)
        return "FlatCDSSpotVolHJM";
    else if (_creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_CIR || 
             _creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_CID || 
             _creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_CID_FAST || 
             _creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_CID_FULL )
        return "FlatCDSSpotVol";
    else if (_creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_LIBOR)
        // credit Libor model switched off for now
        throw ModelException(method, "creditModelType "+_creditModelType
        +" is currently switched off");
        // return "FlatCDSSpotVol";
    else
        throw ModelException(method, "creditModelType "+_creditModelType
        +" is unknown. Please select one of {"+ getAllSupportedCreditModels() +"}");
}


// TODO: Should use a map for a more efficient look up...
//// returns the number of factor of the energy model
int MCPathConfigSRM::getEnergyNumFactor(const string & energyName) const
{
	return energyNbFactors[findEnergyIndex(energyName)];
}

//// returns the start date of the energy correlation instrument
// TODO: Should use a map for a more efficient look up...
const string & MCPathConfigSRM::getEnergyCorrInstrStart(const string & energyName) const
{
	return energyCorrInstrStarts[findEnergyIndex(energyName)];
}

//// returns the maturity of the energy correlation instrument
// TODO: Should use a map for a more efficient look up...
const string & MCPathConfigSRM::getEnergyCorrInstrMaturity(const string & energyName) const
{
	return energyCorrInstrMaturities[findEnergyIndex(energyName)];
}

void MCPathConfigSRM::validatePop2Object()
{
    static const string method("MCPathConfigSRM::validatePop2Object");

    if (CommandLineParams::hasParameter(CommandLineParams::FewIter)){
        numICERuns = numICERuns > 3? 2: numICERuns;
    }
    if (timePointsPerYear < 1 || 
        (matchLegacySRM3 && timePointsPerYear > DateTime::DAYS_PER_YEAR)){
        throw ModelException(method, "timePointsPerYear is < 1 or > 365");
    }
    if (calibrationStyle.empty()){
        throw ModelException(method, "No calibrationStyle specified");
    }
    if (calibrationMaturity.empty()){
        throw ModelException(method, "No calibrationMaturity specified");
    }
    if (numIRFactors.empty()){
        throw ModelException(method, "No IR factors specified");
    }
    if (zcToDiffuse.empty()){
        throw ModelException(method, "No yield curve ID's specified");
    }
    if (irModelParams.empty()){
        throw ModelException(method, "No irModelParams specified");
    }
    if (irSmileParams.empty()){
        throw ModelException(method, "No irSmileParams specified");
    }
    if (correlationSwapStart.empty()){
        throw ModelException(method, "No correlationSwapStart specified");
    }
    if (correlationSwapMat.empty()){
        throw ModelException(method, "No correlationSwapMat specified");
    }
    if (correlationSwapDCC.empty()){
        throw ModelException(method, "No correlationSwapDCC specified");
    }
    if (correlationSwapFreq.empty()){
        throw ModelException(method, "No correlationSwapFreq specified");
    }
    int numFactorsSpecified = numIRFactors.size();
    int numModelParams = irModelParams.size();
    int numSmileParams = irSmileParams.size();
    int numCalibrationStyle = calibrationStyle.size();
    int numCalibrationMaturity = calibrationMaturity.size();
    int numCalibrationMaturityCMS = calibrationMaturityCMS.size();
    int numCorrelationSwapStart = correlationSwapStart.size();
    int numCorrelationSwapMat = correlationSwapMat.size();
    int numCorrelationSwapDCC = correlationSwapDCC.size();
    int numCorrelationSwapFreq = correlationSwapFreq.size();
    int numSkipNegativeIRVols = skipNegativeIRVols.size();
    int numYCIDs = zcToDiffuse.size();
    int numFXVolModes = fxVolBootStrapMode.size();
    for (int i = 0; i < numFXVolModes; i++){
        const string& s = fxVolBootStrapMode[i];
        if (s != SRMFXDiffuse::CONSTANT_SPOT_VOL && 
            s != SRMFXDiffuse::NO_FAILURE_ALLOWED &&
            s != SRMFXDiffuse::USE_LAST_LEVEL){
            throw ModelException(method, "fxVolBootStrapMode '"+s+
                                 "'not recognises");
        }
    }                        
    int numFXCutoffs = fxCutOffLevel.size();
    int numISOCodes = isoCode.size();
    if ((numFactorsSpecified > 1 && numISOCodes != numFactorsSpecified) ||
        (numModelParams > 1 && numISOCodes != numModelParams) ||
        (numSmileParams > 1 && numISOCodes != numSmileParams) ||
        (numCalibrationStyle > 1 && numISOCodes != numCalibrationStyle)||
        (numCalibrationMaturity > 1 && 
         numISOCodes != numCalibrationMaturity) ||
        (numCorrelationSwapStart > 1 && 
         numISOCodes != numCorrelationSwapStart) ||
        (numCorrelationSwapMat > 1 && 
         numISOCodes != numCorrelationSwapMat) ||
        (numCorrelationSwapDCC > 1 && 
         numISOCodes != numCorrelationSwapDCC) ||
        (numCorrelationSwapFreq > 1 && 
         numISOCodes != numCorrelationSwapFreq) ||
        (numSkipNegativeIRVols > 1 &&
         numISOCodes != numSkipNegativeIRVols) ||
        (numYCIDs >1 && numISOCodes != numYCIDs) ||
        (numFXCutoffs >1 && numISOCodes != numFXCutoffs) ||
        (numFXVolModes >1 && numISOCodes != numFXVolModes)){
        throw ModelException(method, "When specifying more than one option"
                             " in an array, the number specified must "
                             "match the number of ISO codes supplied");
    }

    if (numCalibrationMaturityCMS == 0)
	{
		// always calibrate to 5Y column, somewhat arbitrary (would want flexibility for YCSO)
        calibrationMaturityCMS = StringArray(1, "60M");
	}

	if (numCalibrationMaturityCMS > 0  && numCalibrationMaturityCMS < numCalibrationMaturity)
	{
		throw ModelException(method, "If specifying more than one CMS maturity"
			" in an array, you should specify CMS maturities for all curves");
	}

	if ( fabs(corrLower) > 1. || fabs(corrUpper) > 1. )
	{
		throw ModelException(method, "Target swap rate correlations must lie in (-1,1) range");
	}

    if (creditModelType == CREDIT_MODEL_TYPE_CID && wCIDParametersWrapper.isEmpty() )
        throw ModelException(method, "The Credit Model Type is selected as CID, but the CIDParametersWrapper is empty.");

    if (cdsVolType.empty() )
        cdsVolType = cdsVolTypeFromCreditModelType(creditModelType);

}

/** Creates a PathGenerator used for historic dates */
MCPathGeneratorSP MCPathConfigSRM::pastPathGenerator(
    const IMCProduct* prod )
{
    const MCProductClient* prodClient =
        dynamic_cast<const MCProductClient*>(prod);
    // to be removed once everything is state variable based
    if (!prodClient){
        throw ModelException("MCPathConfigSRM::pastPathGenerator",
                             "Product does not support state variables");
    }
    return MCPathGeneratorSP(new PastPathGen(prodClient));
}

/** Creates a PathGenerator used for future dates which supports
    'random access' to the paths. That is if a path generator is
    created using savePaths = true, then subsequently another path
    generator can be created (with savePaths = false) which will
    support the pathIdx parameter in the generatePath() method, ie the
    paths can be generated in any order. */
MCPathGeneratorSP MCPathConfigSRM::futurePathGenerator(
    int                      cachingMode,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control, 
    Results*                 results,
    DateTimeArray&           simDates )
{
    static const string method("MCPathConfigSRM::futurePathGenerator");
    if (cachingMode & IMCProduct::CACHE_PRODUCT_BIT){
        throw ModelException(method, "Caching not supported yet");
    }
    return makePathGenerator(false, numPaths, pastPathGenerator, 
                             prod, control, results, simDates);
}

/** Essentially a pass through for futurePathGenerator except that the
    relevant caches are created/updated etc */
MCPathGeneratorSP MCPathConfigSRM::makePathGenerator(
    bool                               cachingRequested,
    int                                numPaths,
    const MCPathGeneratorSP&           pastPathGenerator,
    const IMCProduct*                   prod,
    Control*                           control, 
    Results*                           results,
    DateTimeArray&                     simDates )
{
    try {
        const MCProductClient* prodClient = 
            &dynamic_cast<const MCProductClient&>(*prod);
        return MCPathGeneratorSP(new MCPathConfigSRMGen(numPaths, this, 
                                                        pastPathGenerator,
                                                        prodClient,simDates,
                                                        control, results));
    } catch (exception& e){
        throw ModelException(e, "MCPathConfigSRM::makePathGenerator");
    }
}

/** Having retrieved a market object from the cache (or if presupplied) it
    is necessary to ensure that it has all its market data. */
void MCPathConfigSRM::getComponentMarketData(
    const IModel*         model,
    const MarketData*     market,
    MarketObjectSP        mo,
    MarketDataFetcherSRM* mdf) const 
{
    static const string method = "MCPathConfigSRM::getComponentMarketData";
    try {
        // use this to catch the names and types of all the market data
        // that we will model. Also allows us to skip XCB's etc should we
        // want to.
        // Watch out for IYieldCurve. Need to make sure we get the
        // right IRCalib::Model from the cache
        // This routine populates :
        // - 'yieldCurveNames' for subsequent retrieval in getMarket
        // - 'marketFactors' for subsequent retrieval of corrs in getMarket
        // - 
        if (IYieldCurve::TYPE->isInstance(mo)){
            IYieldCurve* stochYC = DYNAMIC_CAST(IYieldCurve, mo.get());
            const string& name = stochYC->getName();
            yieldCurveNames.insert(name); // save the name
            const string& isoCode =// can't ask YC for ISO code (not got market)
                market->getYieldCurveISOCode(name);
            // record the isocode in list used for correlations as well
            marketFactors.insert(MFId(isoCode, mo->getClass()));
            if (numIRFactors.size() > 1){
                int nFactors = getNumIRFactors(isoCode);
                string irCalibModelType;
                if (nFactors != 0) {
                    irCalibModelType = IRCalib::getModelType(nFactors)->getName();
                }
                else { 
                    //can't calibrate 0 factors, but need to retrieve something.
                    irCalibModelType = IRCalib::getModelType(1)->getName();
                }
                mdf->setIrCalibModelType(irCalibModelType);
                mdf->setSwaptionVolFlag(calibrateAgainstSwaptionVols(isoCode));
            }
        }

        // call parent's method on mdf
        mdf->MarketDataFetcherSRM::getComponentMarketData(model, market,mo);

        // keep track of all 'names' of market factors that we're 
        // interested in
        if (ICDSParSpreads::TYPE->isInstance(mo)){
            marketFactors.insert(MFId(mo->getName(), mo->getClass()));
        } 
		else if (EnergyFuturesCurve::TYPE->isInstance(mo)) {
			marketFactors.insert(MFId(mo->getName(), mo->getClass()));
		}
		else if (IBasisIndexCurve::TYPE->isInstance(mo)) {
			marketFactors.insert(MFId(mo->getName(), mo->getClass()));
		}
        else if (IMarketFactor::TYPE->isInstance(mo) &&
                 !IYieldCurve::TYPE->isInstance(mo)){ // YC's done above
            // Correlations only needed between vanilla ccy treatment assets
            // since those are the ones we'll be simulating
            IMarketFactor* factor = DYNAMIC_CAST(IMarketFactor, mo.get());
            if (CAsset::TYPE->isInstance(factor)) {
                CAssetConstSP asset(dynamic_cast<CAsset*>(factor));
                if (asset->getCcyTreatment() == CAsset::CCY_TREATMENT_NONE ||
                    asset->getCcyTreatment() == CAsset::CCY_TREATMENT_VANILLA) {
                    marketFactors.insert(MFId(factor->getName(), factor->getClass()));
                    
                } else {
                    // ?? Perhaps record something about any "composite" assets 

                    // This will be the DOMESTIC ccy, though cannot check that here
                    const string& domCcyName = asset->getYCName();
                    // Update marketFactors with FXAsset names for all pairs versus domccy
                    // If I understand correctly by having this AFTER the call to "mo->getMarket"
                    // it means the base ccy will have been added to yieldCurveNames above
                    for (set<string>::const_iterator iter = yieldCurveNames.begin();
                         iter != yieldCurveNames.end(); ++iter){
                        const string& ycName = *iter;
                        if (ycName != domCcyName) {
                            const string& fxAssetName = market->getFXName(ycName,
                                                                          domCcyName);
                            marketFactors.insert(MFId(fxAssetName, FXAsset::TYPE));
                        }
                    }
                }
            } 
        } 
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

const ICIDParameters* MCPathConfigSRM::getCIDParameters() const
{
    return wCIDParametersWrapper.get();
}


//// returns the calibration style for specific currency
const string& MCPathConfigSRM::getCalibrationStyle(
    const string& isoCode) const
{
    return (calibrationStyle.size() == 1? calibrationStyle.front(): 
            calibrationStyle[findIsoCodeIndex(isoCode)]);
}
    
//// returns the calibration maturity for specific currency
const string& MCPathConfigSRM::getCalibrationMaturity(
    const string& isoCode) const
{
    return (calibrationMaturity.size() == 1? calibrationMaturity.front(): 
            calibrationMaturity[findIsoCodeIndex(isoCode)]);
}

//// returns the calibration maturity for specific currency
const string& MCPathConfigSRM::getCalibrationMaturityCMS(
	const string& isoCode) const
{
	return (calibrationMaturityCMS.size() == 1? calibrationMaturityCMS.front(): 
	calibrationMaturityCMS[findIsoCodeIndex(isoCode)]);
}

int MCPathConfigSRM::getNumIRFactors(const string& isoCode) const
{
    return (numIRFactors.size() == 1? numIRFactors.front(): 
            numIRFactors[findIsoCodeIndex(isoCode)]);
}

bool MCPathConfigSRM::calibrateAgainstSwaptionVols(
    const string& isoCode) const
{
    return (getCalibrationStyle(isoCode) != CALIB_NIL);
}
        
//// returns the ir model params for specific currency
const string& MCPathConfigSRM::getIRModelParams(
    const string& isoCode) const
{
    return (irModelParams.size() == 1? irModelParams.front(): 
            irModelParams[findIsoCodeIndex(isoCode)]);
}
        
//// returns the ir smile params for specific currency
const string& MCPathConfigSRM::getIRSmileParams(
    const string& isoCode) const
{
    return (irSmileParams.size() == 1? irSmileParams.front(): 
            irSmileParams[findIsoCodeIndex(isoCode)]);
}

//// returns the cr smile params for specific credit name
const string& MCPathConfigSRM::getCRSmileParams(
    const string& nameToFind) const
{
    int idx = 0; 
    if (crSmileParams.empty()) {
        // do nthg => idx = 0
    } else {
        for (int i = 0; i < crSmileParams.size(); i++) {
            string name = diffCDSs[i]->getName();
            if (name == nameToFind) {
                idx = i;
                return (crSmileParams.size() == 1 ? crSmileParams.front() : crSmileParams[idx]);
            }
        }
    }
    // crSmileParams is non-empty, but a match was not found 
    throw ModelException("MCPathConfigSRM::getCRSmileParams", 
                         "No smile data specified for credit "+nameToFind);
    // now return the smile params
}

const string& MCPathConfigSRM::getFXVolBootStrapMode(
    const string& isoCode) const
{
    if (fxVolBootStrapMode.empty()){
        throw ModelException("MCPathConfigSRM::getFXVolBootStrapMode",
                             "No fx vol boot strap modes specified");
    }
    return (fxVolBootStrapMode.size() == 1? fxVolBootStrapMode.front(): 
            fxVolBootStrapMode[findIsoCodeIndex(isoCode)]);
}

double MCPathConfigSRM::getFXCutOffLevel(const string& isoCode) const
{
    if (fxCutOffLevel.empty()){
        throw ModelException("MCPathConfigSRM::getFXCutOffLevel",
                             "No fx cut off levels specified");
    }
    return (fxCutOffLevel.size() == 1? fxCutOffLevel.front(): 
            fxCutOffLevel[findIsoCodeIndex(isoCode)]);
}

const string& MCPathConfigSRM::getFXVolType(const string& isoCode) const
{
    if (fxVolType.empty()){
        throw ModelException("MCPathConfigSRM::getFXVolType",
                             "No fx vol type specified");
    }
    return (fxVolType.size() == 1? fxVolType.front(): 
            fxVolType[findIsoCodeIndex(isoCode)]);
}

const string& MCPathConfigSRM::getEqVolBootStrapMode(int index) const
{
    if (eqVolBootStrapMode.empty()){
        throw ModelException("MCPathConfigSRM::getEqVolBootStrapMode",
                             "No eq vol boot strap modes specified");
    }
    if (index<0) {
        throw ModelException("MCPathConfigSRM::getEqVolBootStrapMode",
                             "Attempt to access for asset not in MultiAsset");
    }
    return (eqVolBootStrapMode.size() == 1? eqVolBootStrapMode.front(): 
            eqVolBootStrapMode[index]);
}

double MCPathConfigSRM::getEqCutOffLevel(int index) const
{
    if (eqCutOffLevel.empty()){
        throw ModelException("MCPathConfigSRM::getEqCutOffLevel",
                             "No eq cut off levels specified");
    }
    if (index<0) {
        throw ModelException("MCPathConfigSRM::getEqCutOffLevel",
                             "Attempt to access for asset not in MultiAsset");
    }
    return (eqCutOffLevel.size() == 1? eqCutOffLevel.front(): 
            eqCutOffLevel[index]);
}

/** Returns GROWTH or DISCOUNT etc for supplied iso code */
const string& MCPathConfigSRM::getZCToDiffuse(const string& isoCode) const
{
    return (zcToDiffuse.size() == 1? zcToDiffuse.front(): 
            zcToDiffuse[findIsoCodeIndex(isoCode)]);
}

/** Returns correlation swap start for supplied iso code */
const string& MCPathConfigSRM::getCorrelationSwapStart(
    const string& isoCode) const
{
    return (correlationSwapStart.size() == 1? correlationSwapStart.front(): 
            correlationSwapStart[findIsoCodeIndex(isoCode)]);
}

/** Returns correlation swap maturity for supplied iso code */
const string& MCPathConfigSRM::getCorrelationSwapMat(
    const string& isoCode) const
{
    return (correlationSwapMat.size() == 1? correlationSwapMat.front(): 
            correlationSwapMat[findIsoCodeIndex(isoCode)]);
}

/** Returns correlation swap day count convention for supplied iso code */
const string& MCPathConfigSRM::getCorrelationSwapDCC(
    const string& isoCode) const
{
    return (correlationSwapDCC.size() == 1? correlationSwapDCC.front(): 
            correlationSwapDCC[findIsoCodeIndex(isoCode)]);
}

/** Returns correlation swap frequency for supplied iso code */
const string& MCPathConfigSRM::getCorrelationSwapFreq(
    const string& isoCode) const
{
    return (correlationSwapFreq.size() == 1? correlationSwapFreq.front(): 
            correlationSwapFreq[findIsoCodeIndex(isoCode)]);
}

/** Returns a VolProcessedBSIRSP for the supplied stochastic yield curve
    appropriate for the specified calibration style. Note returns null
    for calibration style NIL */
CVolProcessedSP MCPathConfigSRM::getProcessedIRVol(
    const string& isoCode,
    IYieldCurveConstSP   stochYC) const
{
    static const string method("MCPathConfigSRM::getProcessedIRVol");
    const string& calib = getCalibrationStyle(isoCode);
    if (calib == CALIB_NIL){
        return VolProcessedBSIRSP();
    }
    // get hold of the maturity
    const string& calibMat = getCalibrationMaturity(isoCode);
    if (calibMat.empty()){
        throw ModelException(method, "Calibration Maturity for "+isoCode+
                             " is blank");
    }
    ExpirySP calibExpiry(new MaturityPeriod(calibMat));
    CVolRequestSP request(
                          new SwapMaturityVolRequest(calibExpiry.get(), calib));
    // get hold of processed vol
    return CVolProcessedSP(stochYC->getProcessedVol(request.get()));
}

bool MCPathConfigSRM::skipNegIRVols(const string& isoCode) const
{
    if (skipNegativeIRVols.empty()){
        return false;
    }
    return skipNegativeIRVols[findIsoCodeIndex(isoCode)];
}
                            
/** Invoked after instrument has got its market data. Allows model to
    get any extra data required. */
void MCPathConfigSRM::getMarket(
    const IModel*      model,
    const MarketData*  market,
    IInstrumentCollectionSP instrument)
{
    static const string method("MCPathConfigSRM::getMarket");
    try {

        if (!wCIDParametersWrapper.isEmpty())
            wCIDParametersWrapper.getData(model, market);

        /** switch back mechanism */        
        DependenceMakerGaussSrm* dependenceMakerTmp = 
            dynamic_cast<DependenceMakerGaussSrm*>(dependenceMaker.get());
        if (dependenceMakerTmp) { // it could be gauss term ... :S
            if (dependenceMakerTmp->getCorrTermStructureMode()) {
                MarketObjectArraySP mo = 
                    market->GetAllDataWithType(CorrelationCategory::TYPE); 
                if (mo->size()==0) {                
                    // switch back since there are no correlation category object    
                    dependenceMakerTmp->setCorrTermStructureMode(false);
                    dependenceMaker = DependenceMakerSP(dependenceMakerTmp); 
                }
            }
        }

        /** for backwards compatibility, another check */
        if (DependenceMakerGaussTermSrm::TYPE->isInstance(dependenceMaker)) {
            MarketObjectArraySP mo = 
                market->GetAllDataWithType(CorrelationCategory::TYPE); 
            if (mo->size()==0) {
                DependenceMakerGaussTermSrm* dependenceMakerTmp = 
                    dynamic_cast<DependenceMakerGaussTermSrm*>(dependenceMaker.get());
                // switch back since there are no correlation category object
                dependenceMakerTmp->setCorrTermStructureMode(false);
                dependenceMaker = DependenceMakerSP(dependenceMakerTmp);
            }
        }

        // need to retrieve any extra IR correlations etc
        // Also form a list of assets (non-YC) we shall diffuse
        int numAssets = marketFactors.size(); // all IR+FX+EQ+CR+Enrg+Basis
        corrObjArray.clear();
        eqeqCorrObjArray.clear();
        eqeqCorrTermObjArray.clear();
        // TO DO:  This is way overkill for Sampras (ie. when strictCorr is 
        // false) since most correlations will be zero.  Is there a better
        // guess in this case?  Commented out for now.
        //corrs.reserve(numAssets*(numAssets-1)/2);
        diffAssets.clear();
        //diffAssets.reserve(numAssets);
        diffCDSs.clear();
        //diffCDSs.reserve(numAssets); // usually too big, since numAssets incl all assets TODO
		//diffEnrgs.reserve(numAssets); // TODO: overkill to reserve numAssets so many enrg curves
        diffBasisArray.clear(); 
        //diffBasisArray.reserve(numAssets);

        MarketDataFetcherSP fetcher = model->getMDF();
        for (MFHashSet::iterator iter1 = marketFactors.begin(); iter1 != marketFactors.end(); ++iter1)
        {


            // special handling - for CID. The correlation of CID's common market factor is 
            // optional, yet might be taken into account

            if (isCreditCID())
            {
                if ( market->hasCorrelationData(iter1->name, COMMON_MARKET_FACTOR_FOR_CID) )
                {
                    string corrName = market->getCorrelationName(iter1->name, COMMON_MARKET_FACTOR_FOR_CID);

                    CorrelationCommonSP corr = CorrelationCommonSP::dynamicCast(
                        model->getCorrelation(
                            corrName,
                            iter1->clazz, 
                            0,
                            CorrelationCommon::TYPE,
                            market));
                    corrObjArray.push_back(corr); 
                }
            }

            /*
             *  Note that if post-calibration (aka low-level) correlation may be input,
             *  it is possible that we have non-trivial correlation for, e. g., 
             *  (USD, USD) as a matrix
             */
            for (MFHashSet::iterator iter2 = iter1; iter2 != marketFactors.end(); ++iter2)
            {
                /** retrieve corr name from market data cache
                    AB: Why not get this from MultiAsset which also
                    collects the correlation information? */
                bool found = market->hasCorrelationData(iter1->name, iter2->name);
                if ( !found && strictCorr && iter1!=iter2 ) {
                    throw ModelException(method, 
                         "No correlation available with respect to " 
                         + iter1->name + " and " + iter2->name);
                }
                if ( found ) {
                    string corrName = market->getCorrelationName(iter1->name, iter2->name);

                    CorrelationCommonSP corr = CorrelationCommonSP::dynamicCast(
                        model->getCorrelation(
                            corrName,
                            iter1->clazz, 
                            iter2->clazz,
                            CorrelationCommon::TYPE,
                            market));
                    corrObjArray.push_back(corr); 

                    /** this is needed for corr mapping, if applicable */ 
                    if (CAsset::TYPE->isAssignableFrom(iter1->clazz) &&
                        CAsset::TYPE->isAssignableFrom(iter2->clazz)) {
                        eqeqCorrObjArray.push_back(corr);

                        /** even more work todo in case corr term structure is to be used 
                            difference to MultiAsset: assign name even to zero objects, since we need 
                            to sort object array later on ... */
                        if (MDFUtil::useCorrTerm(*fetcher)) {
                            if (market->hasData(iter1->name, CorrelationCategory::TYPE)) {
                                MarketObjectSP mo1 = 
                                    model->GetMarket(market, iter1->name, CorrelationCategory::TYPE);
                                CorrelationCategorySP category1 = CorrelationCategorySP::dynamicCast(mo1);
                                if (market->hasData(iter2->name, CorrelationCategory::TYPE)) {
                                    MarketObjectSP mo2 = 
                                        model->GetMarket(market, iter2->name, CorrelationCategory::TYPE);
                                    CorrelationCategorySP category2 = CorrelationCategorySP::dynamicCast(mo2);
                                    const string& corrTermName = 
                                        market->getCorrelationTermName(category1->getCategoryName(), category2->getCategoryName());
                                    MarketObjectSP moTerm = model->GetMarket(market, corrTermName, CorrelationTerm::TYPE);
                                    eqeqCorrTermObjArray.push_back(CorrelationTermSP::dynamicCast(moTerm));
                                    eqeqCorrTermObjArray.back()->getMarket(model, market);
                                } else {
                                    eqeqCorrTermObjArray.push_back(CorrelationTermSP(new 
                                                                                     CorrelationTerm(true)));
                                }
                            } else {
                                eqeqCorrTermObjArray.push_back(CorrelationTermSP(new CorrelationTerm(true)));
                            }
                        }
                    }
                }
                /** Ignore missing correlations if reach this point.
                    Reason : strictCorr must be false */
            }
            // YieldCurves dealt with separately below
            if (!IYieldCurve::TYPE->isAssignableFrom(iter1->clazz)) {
                if(ICDSParSpreads::TYPE->isAssignableFrom(iter1->clazz)) {
                    MarketObjectSP credit(market->GetData(model, iter1->name, iter1->clazz));
                    diffCDSs.push_back(ICDSParSpreadsSP::dynamicCast(credit));
				} 
				else if (EnergyFuturesCurve::TYPE->isAssignableFrom(iter1->clazz)) {
					MarketObjectSP enrg(market->GetData(model, iter1->name, iter1->clazz));
					diffEnrgs.push_back(EnergyFuturesCurveSP::dynamicCast(enrg));
				}
                // What to do with a basis curve
                else if (IBasisIndexCurve::TYPE->isAssignableFrom(iter1->clazz)) {
                    MarketObjectSP basis(market->GetData(model, iter1->name, iter1->clazz));
                    diffBasisArray.push_back(IBasisIndexCurveSP::dynamicCast(basis));
                }
				else {
                    MarketObjectSP asset(market->GetData(model, iter1->name, 
                                                         iter1->clazz));
                    diffAssets.push_back(CAssetSP::dynamicCast(asset));
                }
            }

        }

        /* then get any 'GROWTH' yield curves (discount ones will be in the
           instrument/market data already) */
        diffYCs.reserve(yieldCurveNames.size());
        for (set<string>::const_iterator iter = yieldCurveNames.begin();
             iter != yieldCurveNames.end(); ++iter){
            const string& ycName = *iter;
            const string& ycID = getZCToDiffuse(*iter);
            if (ycID == ZC_GROWTH){
                try{
                    MarketObjectSP yc(market->GetData(model, ycName, 
                                                      YieldCurve::TYPE));
                    diffYCs.push_back(YieldCurveSP::dynamicCast(yc));
                    diffYCs.back()->setProjectionCurve(); // switch zc
                } catch (exception& e){
                    throw ModelException(e, method, "Could not retrieve "
                                         "yield curve for "+ycName);
                }
            } else if (ycID != ZC_DISCOUNT){
                throw ModelException(method, "Unrecognized zero curve to "
                                     "diffuse (must be GROWTH or DISCOUNT)");
            }
        }    
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Given the discount YC, what YC do we diffuse */
IYieldCurveConstSP MCPathConfigSRM::getDiffusedYC(IYieldCurveConstSP discount)
{
    const string& ccy = discount->getCcy();
    for (int i = 0; i < diffYCs.size(); i++){
        if (diffYCs[i]->getCcy() == ccy){
            return diffYCs[i];
        }
    }
    // if not listed then diffusing discount
    return discount;
}

DateTimeArray MCPathConfigSRM::getCriticalDates(
    IQMCGenDiffusibleAsset* gen,
    const DateTime& start,       // likely to be "today"
    const DateTime& finish) // the latest requested date
{
    ISupportPathConfigSRM* genSRM = dynamic_cast<ISupportPathConfigSRM*>(gen);
    QLIB_VERIFY( genSRM != 0, "Failed to create SRM Util object");
    return genSRM->getSRMCriticalDates(this,start,finish);
}

void MCPathConfigSRM::createUtil( 
    IQMCGenDiffusibleAsset* gen, 
    const DateTime&        today,
    DateTimeArrayConstSP   simDates)
{
    ISupportPathConfigSRM* genSRM = dynamic_cast<ISupportPathConfigSRM*>(gen);
    QLIB_VERIFY( genSRM != 0, "Failed to create SRM Util object");
    genSRM->createSRMUtil(this,today,simDates);
}

void MCPathConfigSRM::setDiffusibleAsset(
    const DateTime&     today,              // base date of the run
    IQMCGenDiffusibleAsset* gen,
    const IPastValues*  pastValues,         // historic values
    DependenceSP        dependence)    // factorized correlations
{
    ISupportPathConfigSRM* genSRM = dynamic_cast<ISupportPathConfigSRM*>(gen);
    QLIB_VERIFY( genSRM != 0, "Failed to create SRM Util object");
    genSRM->setSRMDiffusibleAsset(today, this, pastValues, dependence);
}

MCPathConfigSRM::MCPathConfigSRM(): 
    MCPathConfigParametric(
        TYPE, 
        DependenceMakerGaussSrmSP(new DependenceMakerGaussSrm())), 
    fxVolType(1, "SRMFX::Vol"), 
    numICERuns(2),
    isCarefulRandoms(false), 
    skipNegativeIRVols(false),
    flatVolIr(1.0),
    choiceCutoff("default"),
    cutoffValue(0.5),
    numSDMaxEffRateIR(2.0), 
    numSDMinEffRateIR(2.0),
    eqProtAdjMethod("DEFAULT"),
    crSmileParams(1, "CRCalib::Smile"),
    numSDMaxEffRateCR(2.0),
    numSDMinEffRateCR(2.0),
    numSDMaxEffRateSP(2.0),
    numSDMinEffRateSP(2.0),
    deICE(true), 
    matchLegacySRM3(false), 
    offCycleICE(false),
    dependenceType("not used"),
    strictCorr(true),
    equityDiffusionStyle(eqDEFAULT),
    creditModelType(CREDIT_MODEL_TYPE_CIR),
    ratesModelType(RATES_MODEL_TYPE_HJM),
    betaCorrMax(0.8),
	corrLower(1.),
	corrUpper(1.),
	useIRVolPair(false)
{
}

IObject* MCPathConfigSRM::defaultConstructor()
{
    return new MCPathConfigSRM();
}

/** Invoked when Class is 'loaded' */
void MCPathConfigSRM::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MCPathConfigSRM, clazz);
    SUPERCLASS(MCPathConfigParametric);
    IMPLEMENTS(IMCPathConfig);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(timePointsPerYear, "time points per year");
	FIELD(isCarefulRandoms, "isCarefulRandoms"); // ?
	FIELD_MAKE_OPTIONAL(isCarefulRandoms);
    FIELD(useIRVolPair, "true: Use caplets for short end in addition to swaptions (needed for consistency with rates work).  Defaulted to false for regression tests.");
    FIELD_MAKE_OPTIONAL(useIRVolPair);
    FIELD(matchLegacySRM3, "true: match FI SRM3 numbers, "
                 "default false");
    FIELD_MAKE_OPTIONAL(matchLegacySRM3);

	// "ICE"
    FIELD(offCycleICE, "To help match certain legacy SRM3 numbers");
    FIELD_MAKE_OPTIONAL(offCycleICE);
	FIELD(numICERuns, "How many times to run the Iteratice "
		"Calibration Engine");
	FIELD_MAKE_OPTIONAL(numICERuns);
	FIELD(deICE, "True: deactive Iterative Calibration Engine");
	FIELD_MAKE_OPTIONAL(deICE); // default false

	// rates 
	FIELD(ratesModelType, "Rates model to use");
	FIELD_MAKE_OPTIONAL(ratesModelType);
	FIELD(calibrationStyle, "How to interpolat IR Vols eg CMS");
	FIELD(calibrationMaturity, "Maturity for IR vol interpolation");
	FIELD(calibrationMaturityCMS, "Optional second maturity for IR vol interpolation");
	FIELD_MAKE_OPTIONAL(calibrationMaturityCMS);
	FIELD(isoCode,"Keys different calibration data with currencies");
	FIELD_MAKE_OPTIONAL(isoCode);
	FIELD(numIRFactors, "Number of factors to use for IR");
	FIELD(irModelParams, "Type of model params to use for "
		"IR factors");
	FIELD(irSmileParams, "Type of smile params to use for "
		"IR factors");
	FIELD(correlationSwapStart, "How to correlate with other assets"
		", e.g. 10Y");
	FIELD(correlationSwapMat, "How to correlate with other assets"
		", e.g. 10Y");
	FIELD(correlationSwapDCC, "How to correlate with other assets"
		", e.g. 30/360");
	FIELD(correlationSwapFreq, "How to correlate with other assets"
		", e.g. A");
	FIELD(skipNegativeIRVols, "True: skip negative IR Vols when "
		"bootstrapping");
	FIELD_MAKE_OPTIONAL(skipNegativeIRVols);
	FIELD(flatVolIr, "flat IR vol");
	FIELD_MAKE_OPTIONAL(flatVolIr);
	FIELD(choiceCutoff, "cutoff choice");
	FIELD_MAKE_OPTIONAL(choiceCutoff);
	FIELD(cutoffValue, "cutoff values");
	FIELD_MAKE_OPTIONAL(cutoffValue);
	FIELD(numSDMaxEffRateIR, "max number of standard deviations for "
		"the Interest Rate volatility cutoff");
	FIELD_MAKE_OPTIONAL(numSDMaxEffRateIR);
	FIELD(numSDMinEffRateIR, "min number of standard deviations for "
		"the Interest Rate volatility cutoff");
	FIELD_MAKE_OPTIONAL(numSDMinEffRateIR);
	FIELD(zcToDiffuse, "Which type of zero curve to diffuse");

	// correlation parameters for rates model
	FIELD(corrLower, "");
	FIELD_MAKE_OPTIONAL(corrLower);
	FIELD(corrUpper, "");
	FIELD_MAKE_OPTIONAL(corrUpper);

    // CID parameters
    FIELD(wCIDParametersWrapper, "");
	FIELD_MAKE_OPTIONAL(wCIDParametersWrapper);

    // stratified sampling
    FIELD(stratification, "a bunch of stratas with their relative weights");
    FIELD_MAKE_OPTIONAL(stratification);

	// equity 

	FIELD(volType, "Type of vol to use");
	FIELD_MAKE_OPTIONAL(volType);  // for now 
	FIELD(eqProtAdjMethod, "How to do the quanto adjustment");
	FIELD_MAKE_OPTIONAL(eqProtAdjMethod);
	FIELD(eqVolBootStrapMode, "How to bootstrap eq vol (per equity)");
	FIELD_MAKE_OPTIONAL(eqVolBootStrapMode);
	FIELD(eqCutOffLevel, "Cut off level for EQ Vols");
	FIELD_MAKE_OPTIONAL(eqCutOffLevel);


	// fx 
	FIELD(fxVolType, "Type of fx vol to use");
	FIELD_MAKE_OPTIONAL(fxVolType);
	FIELD(fxVolBootStrapMode, "How to bootstrap fx vol (per foreign "
		"currency)");
	FIELD_MAKE_OPTIONAL(fxVolBootStrapMode);
	FIELD(fxCutOffLevel, "Cut off level for FX Vols");
	FIELD_MAKE_OPTIONAL(fxCutOffLevel);

	// credit
	FIELD(creditModelType, "Credit model to use");
    FIELD_MAKE_OPTIONAL(creditModelType);
	FIELD(crSmileParams, "Type of smile params to use for Credit");
	FIELD_MAKE_OPTIONAL(crSmileParams); // in order to make it work without credit
    FIELD(cdsVolType, "Type of cds vol to use");
    FIELD_MAKE_OPTIONAL(cdsVolType);
	FIELD(numSDMaxEffRateCR, "max number of standard deviations for "
		"the Credit volatility cutoff");
	FIELD_MAKE_OPTIONAL(numSDMaxEffRateCR);
	FIELD(numSDMinEffRateCR, "min number of standard deviations for "
		"the Credit volatility cutoff");
	FIELD_MAKE_OPTIONAL(numSDMinEffRateCR);

	// energy 
	FIELD(energyNames, "Names should match energy curve names");
	FIELD_MAKE_OPTIONAL(energyNames);
	FIELD(energyNbFactors, "Number of factors");
	FIELD_MAKE_OPTIONAL(energyNbFactors);
	FIELD(energyCorrInstrStarts, "Start label for the energy corr instrument");
	FIELD_MAKE_OPTIONAL(energyCorrInstrStarts);
	FIELD(energyCorrInstrMaturities, "Maturity label for the energy corr instrument");
	FIELD_MAKE_OPTIONAL(energyCorrInstrMaturities);

    FIELD(rateDriverTargets, "target rates for rate driver rates");
    FIELD_MAKE_OPTIONAL(rateDriverTargets);
    FIELD(rateDrivers, "rate driver rates");
    FIELD_MAKE_OPTIONAL(rateDrivers);

    FIELD(creditDriverTargets, "target credits for credit driver rates");
    FIELD_MAKE_OPTIONAL(creditDriverTargets);
    FIELD(creditDrivers, "credit driver rates");
    FIELD_MAKE_OPTIONAL(creditDrivers);

    FIELD(energyDriverTargets, "target energies for energy driver rates");
    FIELD_MAKE_OPTIONAL(energyDriverTargets);
    FIELD(energyDrivers, "energy driver rates");
    FIELD_MAKE_OPTIONAL(energyDrivers);

	// correlations
	FIELD(dependenceType, "not used");
	FIELD_MAKE_OPTIONAL(dependenceType);
	FIELD(strictCorr, "FALSE means to ignore missing correlations.  "
                 "TRUE would give error.  Default value is true");
    FIELD_MAKE_OPTIONAL(strictCorr);
    FIELD_NO_DESC(corrObjArray);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrObjArray);
	FIELD(betaCorrArray, "Array of Beta Correlations");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(betaCorrArray); // TODO: is it good enough or we need "Transient" here ?
    FIELD_NO_DESC(betaCorrMax);
    FIELD_MAKE_OPTIONAL(betaCorrMax); 
	FIELD_NO_DESC(eqeqCorrObjArray);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(eqeqCorrObjArray);
    FIELD_NO_DESC(eqeqCorrTermObjArray);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(eqeqCorrTermObjArray);

    FIELD_NO_DESC(diffYCs);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(diffYCs);
    FIELD_NO_DESC(diffAssets);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(diffAssets);
    FIELD_NO_DESC(diffCDSs);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(diffCDSs);
	FIELD_NO_DESC(diffEnrgs);
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(diffEnrgs);
	FIELD_NO_DESC(diffBasisArray);
	FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(diffBasisArray);

    //mapping method, Euler, etc
    FIELD(equityDiffusionStyle, "equity diffusion discretisation type");
    FIELD_MAKE_OPTIONAL(equityDiffusionStyle);


}

CClassConstSP const MCPathConfigSRM::TYPE = 
    CClass::registerClassLoadMethod(
        "MCPathConfigSRM", typeid(MCPathConfigSRM), load);

/**
 * Specialised MonteCarlo Model implementation for SRM --- required because of
 * limitations to polymorphism allowed in IMS.
 */

class MonteCarloSRMDefault: public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloSRMDefault():MonteCarlo(TYPE) {
        useStateVars = true;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(MonteCarloSRMDefault, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloSRM);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultMonteCarloSRM(){
        return new MonteCarloSRMDefault();
    }
};

CClassConstSP const MonteCarloSRMDefault::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloSRMDefault", typeid(MonteCarloSRMDefault), load);

// symbol (referenced by MonteCarloLib.cpp) to ensure file gets linked in
bool MCPathConfigSRMLoad() {
    return MCPathConfigSRM::TYPE != NULL &&
           MonteCarloSRMDefault::TYPE != NULL;
}

DRLIB_END_NAMESPACE
