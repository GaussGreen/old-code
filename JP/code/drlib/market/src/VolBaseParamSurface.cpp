//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBaseParamSurface.hpp
//
//   Description : Implements the tweak dispatch for Parametrized Vol Surfaces
//                 sensitivities
//
//   Author      : Jean-Noël Juston
//
//   Date        : 01 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/PDFParamLNStrike.hpp"
#include "edginc/VolProcessedBSParam.hpp"
#include "edginc/VolProcessedDVFParam.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/VolProcessedDispatch.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/VegaPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

CVolBaseParamSurface::~CVolBaseParamSurface(){}

/** Returns name of vol */
string CVolBaseParamSurface::getName() const{
    return name;
}

CClassConstSP const CVolBaseParamSurface::TYPE =
CClass::registerClassLoadMethod(
    "VolBaseParamSurface", typeid(CVolBaseParamSurface), load);

void CVolBaseParamSurface::load(CClassSP& clazz){
    REGISTER(CVolBaseParamSurface, clazz);
    SUPERCLASS(CVolBase);
    IMPLEMENTS(IPDFCalculator);
    IMPLEMENTS(IRestorableWithRespectTo<VolParallel>);
    IMPLEMENTS(VegaMatrix::IRestorableShift);
    IMPLEMENTS(IRestorableWithRespectTo<VolPointwise>);
    IMPLEMENTS(RootTimeVega::IRestorableShift);
    IMPLEMENTS(VegaSkewParallel::IRestorableShift);
    IMPLEMENTS(VegaSkewPointwise::IRestorableShift);
    IMPLEMENTS(VolLevel::Shift);
    IMPLEMENTS(VolParallelShift::Shift);
    IMPLEMENTS(VolBenchmarkShift::Shift);
    IMPLEMENTS(PowerVega::Shift);
    IMPLEMENTS(Theta::IShift);
    IMPLEMENTS(VolRelativeShift::IShift);
    IMPLEMENTS(VolAbsoluteShift::IShift);
    FIELD(name, "Vol identifier");
    FIELD(volSurfaceForBackbone, "vol backbone");
    FIELD_MAKE_TRANSIENT(volSurfaceForBackbone);

    // register the flavours of vol request that we support
    VolProcessedDispatch::registerVolMethod(
        &CVolBaseParamSurface::getProcessedVolTime);
    VolProcessedDispatch::registerVolMethod(
        &CVolBaseParamSurface::getProcessedVolLN);
    VolProcessedDispatch::registerVolMethod(
        &CVolBaseParamSurface::getProcessedVolDVF);
}

CVolBaseParamSurface::CVolBaseParamSurface(const CClassConstSP& clazz):
    CVolBase(clazz) {}

    CVolBaseParamSurface::CVolBaseParamSurface(const CClassConstSP& clazz, const string& name):
    name(name), CVolBase(clazz) {}

//// constructor (eg see VolSpline on how to calculate volParam)
CVolBaseParamSurface::CVolBaseParamSurface(const CClassConstSP&    clazz,
                                           const VolSurface&       volSurface,
                                           const CVolParamSP&      volParam):
    CVolBase(clazz), name(volSurface.getName()),
    volSurfaceForBackbone(copy(&volSurface)),
    myVolParam(volParam) {}

//// copies over the fields as appropriate
IObject* CVolBaseParamSurface::clone() const{
    // start by calling parent's clone method
    CVolBaseParamSurface& paramSurf = 
        dynamic_cast<CVolBaseParamSurface&>(*CVolBase::clone());
    // we can just copy the SP's over since the CVolParam doesn't contain 
    // anything that changes
    paramSurf.myVolParam = myVolParam;
    paramSurf.previousVols = previousVols;
    return &paramSurf;
}

void CVolBaseParamSurface::getMarket(const IModel*     model, 
                                     const MarketData* market,
                                     const string&     name){
    try{
        if (!volSurfaceForBackbone){
            volSurfaceForBackbone = VolSurfaceSP(VolSurfaceSP::dynamicCast(
                market->GetData(name,VolSurface::TYPE)));
        }
        /* we now have to manually call getMarket on our realVolBase here.
           Normally this would be done for us, but we are bypassing the
           normal route here (normally routed through Model but that
           would result in us going in circles) */
        // Dummy try/catch to avoid vc6.opt crash
		try {
            volSurfaceForBackbone->getMarket(model, market);
        }
        catch (...) { throw; }

        // give derived class chance to do any caching (or alternatively
        // could call getVolParam() here)
        // Dummy try/catch to avoid vc6.opt crash
		try {
            buildCache();
        } catch (...) { throw; }

        myVolParam = CVolParamSP(createVolParam());
    } catch (exception& e){
        throw ModelException(e, "CVolBaseParamSurface::getMarket");
    }
}

void CVolBaseParamSurface::getMarket(const IModel*     model, 
                                     const MarketData* market){
    getMarket(model, market, getName());
}

/** Returns the parameterised vol - needed to calls can be made directly
    to ComputeImpVol (which we want to keep out of this interface) by,
    for example, the pdfCalculator. 
    Note that either getMarket() must have been called or this object needs
    to have been constructed with a supplied vol surface for this method
    to work */
CVolParamConstSP CVolBaseParamSurface::getVolParam() {
    if (!myVolParam){
        if (!volSurfaceForBackbone){
            throw ModelException("CVolBaseParamSurface::getVolParam", "No "
                                 "market data available");
        }
        // give derived class chance to do any caching
        // Dummy try/catch to avoid vc6.opt crash
        try { buildCache(); } catch (...) { throw; }
        myVolParam = CVolParamSP(createVolParam());
    }
    return myVolParam;
}


/** returns a constant reference to surface to be used for the backbone */
const VolSurface* CVolBaseParamSurface::getBackboneSurface() const{
    return volSurfaceForBackbone.get();
}

void CVolBaseParamSurface::buildCache() {} // default: do nothing

void CVolBaseParamSurface::buildCache(VolSurfaceSP volSurfaceForBackbone) {} // default: do nothing

static bool isLNRequest(const CVolRequest* volRequest){
    bool lnRequest;
    CClassConstSP clazz = volRequest->getClass();
    // for speed - specific test for common ones
    if (LinearStrikeTSVolRequest::TYPE == clazz ||
        LinearStrikeVolRequest::TYPE == clazz){
        lnRequest = true;
    } else if (CVolRequestDVF::TYPE == clazz){
        lnRequest = false;
    } else if (CVolRequestDVF::TYPE->isInstance(volRequest)){
        lnRequest = false;
    } else if (CVolRequestLN::TYPE->isInstance(volRequest)){
        lnRequest = true;
    } else {
        throw ModelException("CVolBaseParamSurface::getProcessedVol",
                             "VolRequest of type "+ clazz->getName()+
                             " not supported");
    }
    return lnRequest;
}
//// handles VolRequestTime
IVolProcessed* CVolBaseParamSurface::getProcessedVolTime(
    const VolRequestTime* volRequest,
    const CAsset*         asset) const{
    return VolRequestTime::createVolProcessed(
        getName(), volSurfaceForBackbone->getTimeMetric());
}

//// handles LNRequest
IVolProcessed* CVolBaseParamSurface::getProcessedVolLN(
    const CVolRequestLN* volRequest,
    const CAsset*        asset) const{
    return new CVolProcessedBSParam(this, myVolParam, volRequest, asset);
}

//// handles RequestDVF
CVolProcessed* CVolBaseParamSurface::getProcessedVolDVF(
    const CVolRequestDVF* volRequest,
    const CAsset*         asset) const{
    return new CVolProcessedDVFParam(this, myVolParam, volRequest, 
                                     asset,volSurfaceForBackbone.get());
}

/** Creates processed vol - CVolProcessedBS or DVF */
IVolProcessed* CVolBaseParamSurface::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      asset) const{
    // we delegate the work of doing a 'double dispatch' to VolProcessedDispatch
    // Here 'double dispatch' means choosing the method based upon the type
    // of the vol and the request
    return VolProcessedDispatch::dispatch(this, volRequest, asset);
}

/** Creates Struck processed vol - CVolProcessedBS or DVF */
CVolProcessed* CVolBaseParamSurface::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    static const string routine("CVolBaseParamSurface::getProcessedVol"
                                " (Struck)");
    bool lnRequest = isLNRequest(volRequest);
    if (lnRequest){
        // cf getProcessedVol() above - struck case not optimised (yet)
        const CVolRequestLN& request =
            static_cast<const CVolRequestLN&>(*volRequest);
        return new CVolProcessedBSParam(this, myVolParam, &request, eqAsset, 
                                        fxAsset, eqFXCorr);
    } else {
        const CVolRequestDVF& request = 
            static_cast<const CVolRequestDVF&>(*volRequest);
        return new CVolProcessedDVFParam(this, myVolParam, &request, 
                                         eqAsset,volSurfaceForBackbone.get(), 
                                         fxAsset, eqFXCorr);
    }
}

bool CVolBaseParamSurface::genericNonRestorableShift(ITweakID* shift){
    const DateTime& valueDate = volSurfaceForBackbone->getBaseDate();
    TimeMetricConstSP timeMetric = volSurfaceForBackbone->getTimeMetric();
    myVolParam = CVolParamSP(new CVolParamTweak(myVolParam, valueDate, 
                                                timeMetric, shift));
    return true;
}    

bool CVolBaseParamSurface::genericRestorableShift(ITweakID* shift){
    previousVols.push_back(myVolParam);
    return genericNonRestorableShift(shift);
}

void CVolBaseParamSurface::genericRestore(){
    if (previousVols.empty()){
        throw ModelException("CVolBaseParamSurface::genericRestore",
                             "Incorrect call to restore");
    }
    myVolParam = previousVols.back();
    previousVols.pop_back();
}

/** Shifts the backbone by supplied shift */
bool CVolBaseParamSurface::shiftBackbone(SensControl* shift){
    IObjectSP surfAsObject(volSurfaceForBackbone);
    shift->shift(surfAsObject);
    // then update the derived class
    // Dummy try/catch to avoid vc6.opt crash
    try { buildCache(); } catch (...) { throw; }
    return true;
}

/** VegaParallel */
string CVolBaseParamSurface::sensName(const VolParallel*) const{
    return getName();
}
TweakOutcome CVolBaseParamSurface::sensShift(const PropertyTweak<VolParallel>& tweak){
    AbstractPropertyTweakHypothesisSP vp(
        new PropertyTweakHypothesis<VolParallel>(tweak));
    return TweakOutcome(
        tweak.coefficient,
        genericRestorableShift(auto_ptr<ITweakOptID>(vp->asTweakOptID()).get()));
}
void CVolBaseParamSurface::sensRestore(const PropertyTweak<VolParallel>& tweak){
    genericRestore();
}


/** VegaMatrix Interface */
string CVolBaseParamSurface::sensName(VegaMatrix* shift) const{
    return getName();
}

/** Returns the array of expiries that need to be tweaked for this vol */
ExpiryArrayConstSP CVolBaseParamSurface::sensExpiries(VegaMatrix* shift) const{
    return volSurfaceForBackbone->sensExpiries(shift);
}

/** VegaMatrix Interface */
bool CVolBaseParamSurface::sensShift(VegaMatrix* shift){
    return genericRestorableShift(shift);
}

void CVolBaseParamSurface::sensRestore(VegaMatrix* shift){
    genericRestore();
}

/** VegaPointwise Interface */
string CVolBaseParamSurface::sensName(const VolPointwise*) const{
    return getName();
}
ExpiryWindowArrayConstSP CVolBaseParamSurface::sensQualifiers(const VolPointwise*) const {
    return volSurfaceForBackbone->sensQualifiers((VolPointwise*)0);
}
TweakOutcome CVolBaseParamSurface::sensShift(const PropertyTweak<VolPointwise>& shift){
    AbstractPropertyTweakHypothesisSP vp(
        new PropertyTweakHypothesis<VolPointwise>(shift));
    return TweakOutcome(
        shift.coefficient,
        genericRestorableShift(auto_ptr<ITweakOptID>(vp->asTweakOptID()).get()));
}
void CVolBaseParamSurface::sensRestore(const PropertyTweak<VolPointwise>&){
    genericRestore();
}

/** RootTimeVega Interface */
string CVolBaseParamSurface::sensName(RootTimeVega* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(RootTimeVega* shift){
    return genericRestorableShift(shift);
}

void CVolBaseParamSurface::sensRestore(RootTimeVega* shift){
    genericRestore();
}

/** VegaSkewParallel Interface */
string CVolBaseParamSurface::sensName(VegaSkewParallel* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(VegaSkewParallel* shift){
    return genericRestorableShift(shift);
}

void CVolBaseParamSurface::sensRestore(VegaSkewParallel* shift){
    genericRestore();
}

/** VegaSkewPointwise Interface */
string CVolBaseParamSurface::sensName(VegaSkewPointwise* shift) const{
    return getName();
}
ExpiryArrayConstSP CVolBaseParamSurface::sensExpiries(
    VegaSkewPointwise* shift) const{
    return volSurfaceForBackbone->sensExpiries(shift);
}
bool CVolBaseParamSurface::sensShift(VegaSkewPointwise* shift){
    return genericRestorableShift(shift);
}
void CVolBaseParamSurface::sensRestore(VegaSkewPointwise* shift){
    genericRestore();
}

/** VolLevel */
string CVolBaseParamSurface::sensName(VolLevel* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(VolLevel* shift){
    return genericNonRestorableShift(shift);
}

/** VolParallelShift */
string CVolBaseParamSurface::sensName(VolParallelShift* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(VolParallelShift* shift) {
    return genericNonRestorableShift(shift);
}

/** VolBenchmarkShift Interface */
string CVolBaseParamSurface::sensName(VolBenchmarkShift* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(VolBenchmarkShift* shift){
    // need to store all the backbone expiries on the shift
    shift->cacheExpiries(ExpiryWindow::expiries(
        volSurfaceForBackbone->sensQualifiers((VolPointwise*)0)));

    return genericNonRestorableShift(shift); 
}

/** PowerVega Shift Interface */
string CVolBaseParamSurface::sensName(PowerVega* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(PowerVega* shift){
    return genericNonRestorableShift(shift);
}

/** theta */
bool CVolBaseParamSurface::sensShift(Theta* shift){
    return shiftBackbone(shift);
}

/** VolRelativeShift Interface */
string CVolBaseParamSurface::sensName(VolRelativeShift* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(VolRelativeShift* shift){
    return genericNonRestorableShift(shift); 
}

/** VolAbsoluteShift Interface */
string CVolBaseParamSurface::sensName(VolAbsoluteShift* shift) const{
    return getName();
}
bool CVolBaseParamSurface::sensShift(VolAbsoluteShift* shift){
    return genericNonRestorableShift(shift); 
}

PDFCalculator* CVolBaseParamSurface::getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const {
    static const string method("CVolBaseParamSurface::getPDFCalculator");
    try {
        const DateTime& valDate = volSurfaceForBackbone->getBaseDate();
        // temporay hack until fwd starting for vol params is sorted out
        const PDFRequestLNStrike* lnRequest =
            dynamic_cast<const PDFRequestLNStrike*>(request);
        if (!lnRequest){
            throw ModelException(method, "No support for requests of type "+
                                 request->getClass()->getName());
        }
        return new PDFParamLNStrike(
            valDate, CAssetConstSP::attachToRef(asset),
            volSurfaceForBackbone->getTimeMetric(),
            CVolBaseConstSP::attachToRef(this),
            myVolParam, 
            PDFRequestLNStrikeConstSP(copyIfRef(lnRequest)));
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// returns the current parameterised vol
CVolParamSP CVolBaseParamSurface::getMyVolParam(){
    return myVolParam;
}

/** addin to calculate implied vols */
class ImpliedVolAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    CVolBaseWrapper    vol;
    DoubleArray        strikes;
    DateTimeArray      dates;
    string             volType;
    CMarketDataSP      market;
    double             shiftSize;     // for vega
    ExpirySP           expiryToShift; // null => parallel shift
    bool               fwdStarting;
    DateTime           startDate;     // if fwd starting
    double             spot;          // if fwd starting
    double             fwdAtStart;    // if fwd starting
    double             strike;        // for vega matrix
    DoubleArraySP      sensStrikes;   // for vega matrix
    

    // Return a double matrix at requested strikes and maturities
    static IObjectSP run(ImpliedVolAddin* params){
        return params->calcImpliedVols();
    }

    IObjectSP calcImpliedVols(){
        static const string method = "ImpliedVolAddin::calcImpliedVols";
        if (!vol){
            // Ideally like to instantiate a Model that will select what
            // type of vol to get - however not in right directory.
            // Need to manually pull the vol out of the cache since we would
            // need the model to specify the type of vol to get and there
            // may be many BS vols in the cache
            CClassConstSP volClass = CClass::forName(volType);
            vol.setObject(market->GetData(vol.getName(), volClass));
        }
        // create a non pricing model in order to get the market data
        NonPricingModel dummyModel;
        vol->getMarket(&dummyModel, market.get());
        // convert expiries into dates
        // check dates are increasing and not empty
        DateTime::ensureIncreasing(dates, "Benchmark dates", true);
        if (strikes.empty()){
            throw ModelException(method, "No strikes specified");
        }

        CVolBaseSP theVol(copy(vol.get())); // avoid altering input data
        if (!CVolBaseParamSurface::TYPE->isInstance(theVol.get())){
            throw ModelException(method, "Addin only works for parameterised "
                                 "vols derived from CVolBaseParamSurface");
        }

        if (!expiryToShift){
            // do vega
            PropertyTweakHypothesis<VolParallel>(shiftSize).applyTo(theVol);
        } else if (Maths::isZero(strike)){
            // get hold of expiries from CVolBaseParamSurface
            CVolBaseParamSurface& paramVolBase = 
                dynamic_cast<CVolBaseParamSurface&>(*theVol);

            ExpiryWindowSP window(
                ExpiryWindow::find(paramVolBase.sensQualifiers((VolPointwise*)0),
                                   expiryToShift));

            PropertyTweakHypothesis<VolPointwise>(
                shiftSize, OutputNameConstSP(), window).applyTo(theVol);
        }
        else {
            CVolBaseParamSurface& paramVolBase = 
                dynamic_cast<CVolBaseParamSurface&>(*theVol);
            VegaMatrixSP sensTmp = VegaMatrixSP(new VegaMatrix(shiftSize));
            ExpiryArrayConstSP expiries = 
                paramVolBase.sensExpiries(sensTmp.get());

            DoubleArraySP allStrikes;

            if (sensStrikes.get()) {
                allStrikes = sensStrikes;
            }
            else {
                allStrikes = DoubleArraySP(new DoubleArray(1));
                (*allStrikes)[0] = strike;
            }
            //DoubleArraySP allStrikes(new DoubleArray(1));
            //(*allStrikes)[0] = strike;

            bool gotK = false;
            int strikeIdx;

            for (strikeIdx = 0; 
                 strikeIdx < allStrikes->size() && !gotK; strikeIdx++) {
                gotK = Maths::equals(strike, (*allStrikes)[strikeIdx]);
            }
            strikeIdx--;
            if (!gotK) {
                throw ModelException(method,
                                     "couldn't find strike " + 
                                     Format::toString(strike) + 
                                     " in sensitive strikes");
            }

            int expiryIdx = expiryToShift->search(expiries.get());

            SensControlPerNameSP(new VegaMatrix(shiftSize, 
                                                expiryIdx,
                                                expiries,
                                                strikeIdx,
                                                allStrikes))->
                    findAndShift(theVol, OutputNameConstSP());
        }

        CVolBaseParamSurface& paramVolBase = 
            dynamic_cast<CVolBaseParamSurface&>(*theVol);
        // get our parameterised vol
        CVolParamSP paramVol(paramVolBase.getMyVolParam());
        // then set up parameters for the computeImpVol function call
        vector<int>  latticeShape(dates.size(), strikes.size());
        CLatticeDouble strikesAsLattice(latticeShape);
        for (int i = 0; i < dates.size(); i++){
            for (int j = 0; j < strikes.size(); j++){
                strikesAsLattice[i][j] = strikes[j];
            }
        }
        CLatticeDouble impliedVols(vector<int>(dates.size(), strikes.size()));
        // do calc
        if (!fwdStarting){
            paramVol->ComputeImpVol(&paramVolBase, strikesAsLattice,
                                    dates, impliedVols);
        } else {
            // use surface to get value date and time metric
            const VolSurface* surf = paramVolBase.getBackboneSurface();
            CVolParam::FwdStart fwdStart(surf->getBaseDate(),
                                         startDate, 
                                         surf->getTimeMetric(),
                                         spot, fwdAtStart);
            paramVol->computeFwdStartImpVol(&paramVolBase, fwdStart,
                                            strikesAsLattice, 
                                            false, // absolute strikes
                                            dates,
                                            impliedVols);
        }
        // then map output to double matrix:
        // create space for returned object
        CDoubleMatrixSP output(new DoubleMatrix(strikes.size(), dates.size()));
        DoubleMatrix&  matrix = *output;
        for (int j = 0; j < strikes.size(); j++){
            for (int k = 0; k < dates.size(); k++){
                matrix[j][k] = impliedVols[k][j];
            }
        }
        return output;
    }
    

    ImpliedVolAddin(): CObject(TYPE), volType(IVolatilityDVF::TYPE->getName()),
                       market(new MarketData()), 
                       shiftSize(0), fwdStarting(false), strike(0.0) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ImpliedVolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObj);
        FIELD(vol, "Volatility to generate surface from");
        FIELD(strikes, "Strikes to calculate vols on");
        FIELD(dates, "Dates to calculate vols on");
        FIELD(volType, "Type of the vol in the market data cache");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(market, "Market data cache");
        FIELD_MAKE_OPTIONAL(market);
        FIELD(shiftSize, "Size of vega shift to apply");
        FIELD_MAKE_OPTIONAL(shiftSize);
        FIELD(expiryToShift, "Specific benchmark to shift else parallel");
        FIELD_MAKE_OPTIONAL(expiryToShift);
        FIELD(fwdStarting, "Should fwd start adjustment be made");
        FIELD_MAKE_OPTIONAL(fwdStarting);
        FIELD(startDate, "Inst start date if fwd starting");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(spot, "Asset spot price if fwd starting");
        FIELD_MAKE_OPTIONAL(spot);
        FIELD(fwdAtStart, "Asset spot at start date if fwd starting");
        FIELD_MAKE_OPTIONAL(fwdAtStart);
        FIELD(strike, "strike for VEGA_MATRIX");
        FIELD_MAKE_OPTIONAL(strike);
        FIELD(sensStrikes, "sens strikes for VEGA_MATRIX");
        FIELD_MAKE_OPTIONAL(sensStrikes);
        
        Addin::registerClassObjectMethod("CALC_IMPLIED_VOLS",
                                         Addin::MARKET,
                                         "Calculates implied vols from "
                                         "supplied parameterised vol",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)run);
    }
    
    static IObject* defaultObj(){
        return new ImpliedVolAddin();
    }
 
};

CClassConstSP const ImpliedVolAddin::TYPE = CClass::registerClassLoadMethod(
    "ImpliedVolAddin", typeid(ImpliedVolAddin), load);

DRLIB_END_NAMESPACE
