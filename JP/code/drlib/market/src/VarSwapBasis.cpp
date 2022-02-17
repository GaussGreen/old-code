//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VarSwapBasis.cpp
//
//   Description : 
//
//   Author      : Zhijiang Huang
//
//   Date        : August 15, 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VarSwapBasis.hpp"
#include "edginc/ProtAsset.hpp"
#include "edginc/Model.hpp"
#include "edginc/DeltaToStrike.hpp"
#include "edginc/CashSettlePeriod.hpp"

DRLIB_BEGIN_NAMESPACE

VarSwapBasisSP VarSwapBasis::fetchBasis(const IModel*     model,
                                        const MarketData* market,
                                        const string&     name) {
    static const string routine = "VarSwapBasis::fetchBasis";

    VarSwapBasisSP basis;
    
    try {
        MarketDataFetcherSP fetcher = model->getMDF();
        
        // Try fetching the basis
        MarketObjectSP mo = fetcher->fetch(market, name, TYPE, model);
        if(!mo) {
            // Fetcher disallows basis i.e. non-VarSwap model
            basis = VarSwapBasisSP(   );
            return basis;
        } else {
            basis = VarSwapBasisSP::dynamicCast(mo);
            return basis;
        }
    } catch (exception&) {
        
        // Return NULL for safety
        basis = VarSwapBasisSP(   );
        return basis;
    }
}


const string VarSwapBasis::NONE            = "none";
const string VarSwapBasis::DEFAULT         = "default";
const string VarSwapBasis::SKEW_CUTOFF     = "skew_cutoff";
const string VarSwapBasis::PRICE_CUTOFF    = "price_cutoff";
const string VarSwapBasis::BASIS_ONLY      = "basis_only";

const int VarSwapBasis::NONE_ID         = 0;
const int VarSwapBasis::SKEW_CUTOFF_ID  = 1;
const int VarSwapBasis::PRICE_CUTOFF_ID = 2;
const int VarSwapBasis::BASIS_ONLY_ID   = 3;


////////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef BLACK_DELTA
/** THIS OBJECT MUST NOT PERSIST AS IT USES SHALLOW COPIES OF MARKET DATA.
    Implementation of Vanilla::DeltaImpliedStrike hidden here.
    Local market data are shallow copies. */
class VarSwapBasis::DeltaImpliedStrike {
public: 
    /** Full constructor */
    DeltaImpliedStrike(const DateTime&             valueDate,
                       const DateTime&             maturityDate,
                       bool                        isCall,
                       const CAsset*               asset,
                       const YieldCurve*           discount,
                       double                      tgtDelta):
    asset(asset),
    isCall(isCall),
    valueDate(valueDate),
    maturityDate(maturityDate),
    tgtDelta(tgtDelta) { 

        static const string routine = "VarSwapBasis::DeltaImpliedStrike::DeltaImpliedStrike"; 
                    
        try {
            if(!asset) {
                throw ModelException("Internal error: asset is NULL");
            }

            if(!discount) {
                throw ModelException("Internal error: yield curve is NULL");
            }

            if(valueDate > maturityDate) {
                throw ModelException("Value date " + 
                                     valueDate.toString() +
                                     " is after maturity date "+
                                     maturityDate.toString());
            }

            // Spot
            spot = asset->getSpot();
            
            // Compute forward
            if(Asset::IQuanto::TYPE->isInstance(asset)) {
                fwd = dynamic_cast<const Asset::IQuanto*>(asset)->unadjustedFwdValue(maturityDate);
            } else {
                fwd = asset->fwdValue(maturityDate);
            }
            
            // Compute PV factor
            // Ideally we would like to have the settlement here
            pv = discount->pv(valueDate, maturityDate);

            // Setup a VolRequest once            
            volReq = LinearStrikeTSVolRequestSP(new LinearStrikeTSVolRequest(0.0, valueDate, maturityDate, false));


        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    // calc strike implied by delta
    double calcStrike(double lowerStrike,
                      double upperStrike,
                      double strikeAbsAcc) const {
        static const string routine = "VarSwapBasis::DeltaImpliedStrike::calcStrike";
        try {
            if (Maths::isZero(tgtDelta)) {
                if(isCall) {
                    throw ModelException(routine, "Zero delta for call option is not allowed.");
                } else {
                    return 0.0;
                }
            }
            // validate 0 <= lowerStrike < upperStrike
            if (Maths::isNegative(lowerStrike)
                || !Maths::isPositive(upperStrike - lowerStrike)){
                throw ModelException(routine, "we should have 0.0 <= lowerStrike < upperStrike; got "
                                              + Format::toString(lowerStrike)
                                              + " and "
                                              + Format::toString(upperStrike)
                                              + ", respectively");
            }
            // bracket the root
            double lowerLogStrike = log(lowerStrike);
            double upperLogStrike = log(upperStrike);
            DeltaFunc deltaFunc(this, &DeltaImpliedStrike::calcDelta);
            DeltaDiffFunc deltaDiffFunc(deltaFunc, tgtDelta);
            try{
                ZBrac_bracket(deltaDiffFunc,
                              lowerLogStrike,
                              upperLogStrike);
            } catch (exception& e) {
                throw ModelException::addTextToException(e, 
                                                         "Failed to bracket root a maturity "
                                                         + maturityDate.toString()
                                                         + " and delta level "
                                                         + Format::toString(tgtDelta));
            }            
            lowerStrike = exp(lowerLogStrike);
            upperStrike = exp(upperLogStrike);
            // find the root using ZBrent
            double avgStrike = 0.5 * (lowerStrike + upperStrike);
            double logStrikeAbsAcc = strikeAbsAcc / avgStrike;
            double logStrike;
            try{
                logStrike = ZBrent_solve(deltaDiffFunc,
                                         lowerLogStrike,
                                         upperLogStrike,
                                         logStrikeAbsAcc);
            } catch (exception& e) {
                throw ModelException::addTextToException(e, 
                                                         "Failed to solve for root a maturity "
                                                         + maturityDate.toString()
                                                         + " and delta level "
                                                         + Format::toString(tgtDelta));
            }            
            return exp(logStrike);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    /** Disabled default constructor */
    DeltaImpliedStrike();
    
    // calc delta
    double calcDelta(double logStrike) const{
        static const string routine = "VarSwapBasis::DeltaImpliedStrike::calcDelta";
        try {
            double          strike = exp(logStrike);
            
            // Compute var
            volReq->setStrike(strike);
            CVolProcessedBSSP vol(asset->getProcessedVol(volReq.get()));
            double            var = vol->CalcVar(valueDate, maturityDate);

            //The Black formula is only for spot starting option, need to modify to accommodate forward starting
            return Black::delta(isCall, 
                                spot,
                                fwd,
                                strike,
                                pv,
                                var);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // Typedefs for Functors
    typedef const DeltaImpliedStrike* Ptr;
    typedef double (DeltaImpliedStrike::* Func1DConstPtr)(double) const;
    typedef MemFuncWrapper<Ptr, Func1DConstPtr> DeltaFunc;
    typedef FuncShiftWrapper<DeltaFunc> DeltaDiffFunc;

    
    // MANDATORY FIELDS
    const CAsset*               asset;              //!< Shallow copy of asset
    bool                        isCall;
    DateTime                    valueDate;
    DateTime                    maturityDate;
    double                      tgtDelta;           //!< Positive for call, negative for put

    // TRANSIENT FIELDS
    double                      spot;
    double                      fwd;
    double                      pv;
    LinearStrikeTSVolRequestSP  volReq;
};

#else

class VarSwapBasis::DeltaImpliedStrike{};

#endif


////////////////////////////////////////////////////////////////////////////////////////////////////


VarSwapBasis::VarSwapBasisProc::VarSwapBasisProc(const VarSwapBasis* basis,
                                                 const CAsset* asset,
                                                 const VolSurface* backbone,
                                                 const VarSwapBasis::VarSwapBasisRequest* req,
                                                 const string& tenorTime):
CObject(TYPE), tweakDeltaCutoff(basis->tweakDeltaCutoff), methodId(basis->methodId) {
    
    static const string routine = "VarSwapBasis::VarSwapBasisProc::VarSwapBasisProc";

    try {
        // Copy TimeMetric
        metric = TimeMetricSP(copy(backbone->getTimeMetric().get()));
        
        if(methodId == VarSwapBasis::NONE_ID) {
            // Do something quickly and return
            return;
        }
        
        if(!asset) {
            throw ModelException("Internal error: asset is NULL");
        }

        if(!backbone) {
            throw ModelException("Internal error: Vol backbone is NULL");
        }
    
        int tenorTimeInt = DateTime::timeConvert(tenorTime);

        // Get absolute strke corresponding to deltaCutoff at different maturities
        // and setup strikeCutoffInterp
        const DateTimeArray&   dates = backbone->getDates();
        int          numStrikeBMs = dates.size();
        DateTime     baseDate = backbone->getBaseDate();
        DoubleArray  strikeTradYears(numStrikeBMs);
        DoubleArray  strikes(numStrikeBMs, 0.0);

        double spot         = asset->getSpot();
        double lowerStrike  = spot * 0.1;
        double upperStrike  = spot * 10.;
        double strikeAbsAcc = spot * 0.001;

        int iBM;
        for (iBM = 0; iBM<numStrikeBMs; iBM++) {
            // Old Delta
            strikeTradYears[iBM] = metric->yearFrac(baseDate,dates[iBM]);

#ifdef BLACK_DELTA            
            DeltaImpliedStrike tmp(baseDate,
                                   dates[iBM],
                                   false,
                                   asset,
                                   req->getAssetYC().get(),
                                   basis->deltaCutoff);

            double blackDelta = tmp.calcStrike(lowerStrike,
                                               upperStrike,
                                               strikeAbsAcc);

            strikes[iBM] = blackDelta;
#else
        
            // New Delta
            string className("CVanilla::DeltaImpliedStrikeMakerLN");
            IDeltaToStrikeMakerSP maker(dynamic_cast<IDeltaToStrikeMaker*>(CClass::forName(className)->newInstance()));
            if(!maker) {
                throw ModelException(routine, "Failed to create CVanilla::DeltaImpliedStrikeMakerLN");
            }

            InstrumentSettlementSP instSetllement(new CashSettlePeriod(0));

            IDeltaToStrikeMaker::IDeltaToStrikeConstSP calc(maker->make(
                baseDate,
                dates[iBM],
                false,
                asset,
                req->getAssetYC().get(),
                instSetllement.get(),
                Delta::DEFAULT_SHIFT,
                -basis->deltaCutoff,
                "VolPreferred",
                true));
            
            double modelDelta = calc->calcStrike(lowerStrike, upperStrike, strikeAbsAcc);
            strikes[iBM] = modelDelta;
#endif
        }

        strikeCutoffInterp = LinearInterpolantNonVirtualSP(new 
            LinearInterpolantNonVirtual(strikeTradYears, strikes));

        // setup volBasisInterp
        const StringArray&  basisBMs = basis->volBasisBMs;
        const DoubleArray&  basisVol = basis->volBasis;
        int                 numBasisBMs = basisBMs.size();
        DateTimeArray       basisDates(numBasisBMs);
        DoubleArray         basisTradYears(numBasisBMs);

        for (iBM = 0; iBM<numBasisBMs; iBM++) {
            MaturityPeriod  period(basisBMs[iBM]);
            DateTime date(period.toDate(baseDate).getDate(), tenorTimeInt);
            basisDates[iBM] = date;
            basisTradYears[iBM] = metric->yearFrac(baseDate,date);
        }

        DateTime::ensureIncreasing(basisDates, "Basis dates", false);
        if(basisDates[0] == baseDate) {
            throw ModelException("First basis date is identical to value date. Has to be greater");
        }
        
        volBasisInterp = LinearInterpolantNonVirtualSP(new 
            LinearInterpolantNonVirtual(basisTradYears, basisVol));
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


double VarSwapBasis::VarSwapBasisProc::interpVolBasis(double time) const {
    if(methodId != VarSwapBasis::NONE_ID) {
        // Interpolate
        return volBasisInterp->value(time);
    } else {
        return 0.0;
    }
}


double VarSwapBasis::VarSwapBasisProc::interpPriceCutoff(double time) const {
    if(methodId == VarSwapBasis::PRICE_CUTOFF_ID) {
        return strikeCutoffInterp->value(time);
    } else {
        return 0.0;
    }
}


double VarSwapBasis::VarSwapBasisProc::interpSkewCutoff(double time) const {
    if(methodId == VarSwapBasis::SKEW_CUTOFF_ID) {
        return strikeCutoffInterp->value(time);
    } else {
        return 0.0;
    }
}


void VarSwapBasis::VarSwapBasisProc::recordRequests(double tradYear, Control* control, CResults* results) {
    static const string routine("VarSwapBasis::VarSwapBasisProc::recordRequests");
    try {
        if(control && results) {
            OutputRequest* request = control->requestsOutput(OutputRequest::VAR_SWAP_VOL_BASIS);
            if (request) {
                double volBasis = interpVolBasis(tradYear);
                results->storeRequestResult(request, volBasis);
            }

            request = control->requestsOutput(OutputRequest::VAR_SWAP_CUTOFF);
            if (request) {
                double cutoff = 0.0;
                if(methodId != NONE_ID) {
                    cutoff = strikeCutoffInterp->value(tradYear);
                }
                results->storeRequestResult(request, cutoff);
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


string VarSwapBasis::VarSwapBasisProc::getName() const {
    static const string routine = "VarSwapBasis::VarSwapBasisProc::getName";
    throw ModelException(routine, "Method not supported");
}
        

double VarSwapBasis::VarSwapBasisProc::calcTradingTime(const DateTime &date1,
                                                       const DateTime &date2) const {
    return metric->yearFrac(date1, date2);
}


TimeMetricConstSP VarSwapBasis::VarSwapBasisProc::GetTimeMetric() const {
    return metric;
}


void VarSwapBasis::VarSwapBasisProc::load(CClassSP& clazz) {
    REGISTER(VarSwapBasis::VarSwapBasisProc, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
    EMPTY_SHELL_METHOD(defaultVarSwapBasisProc);
    FIELD(tweakDeltaCutoff, "Tweak Delta Cutoff");
    FIELD(methodId, "Cutoff Methodology");
    FIELD(strikeCutoffInterp, "Strike Interpolant");
    FIELD(volBasisInterp, "Vol Basis Interpolant");
}


IObject* VarSwapBasis::VarSwapBasisProc::defaultVarSwapBasisProc() {
    return new VarSwapBasisProc();
}


VarSwapBasis::VarSwapBasisProc::VarSwapBasisProc(): CObject(TYPE) {}


CClassConstSP const VarSwapBasis::VarSwapBasisProc::TYPE =
CClass::registerClassLoadMethod("VarSwapBasis::VarSwapBasisProc", typeid(VarSwapBasis::VarSwapBasisProc), load);


/////////////////////////////////////////////////////////////////////////////////////////

string VarSwapBasis::VarSwapBasisProcError::getName() const {
    static const string routine = "VarSwapBasis::VarSwapBasisProcError::getName";
    throw ModelException(routine, "Method not supported");
}
        

double VarSwapBasis::VarSwapBasisProcError::calcTradingTime(const DateTime &date1,
                                                            const DateTime &date2) const {
    static const string routine = "VarSwapBasis::VarSwapBasisProcError::calcTradingTime";
    throw ModelException(routine, "Method not supported");

}


TimeMetricConstSP VarSwapBasis::VarSwapBasisProcError::GetTimeMetric() const {
    static const string routine = "VarSwapBasis::VarSwapBasisProcError::GetTimeMetric";
    throw ModelException(routine, "Method not supported");

}


/** Invoked when Class is 'loaded' */
void VarSwapBasis::VarSwapBasisProcError::load(CClassSP& clazz) {
    REGISTER(VarSwapBasis::VarSwapBasisProcError, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IVolProcessed);
}


VarSwapBasis::VarSwapBasisProcError::VarSwapBasisProcError(const ModelException& e): CObject(TYPE), e(e) {}


const ModelException& VarSwapBasis::VarSwapBasisProcError::getException() const {
    return e;
}


CClassConstSP const VarSwapBasis::VarSwapBasisProcError::TYPE =
CClass::registerClassLoadMethod("VarSwapBasis::VarSwapBasisProcError", typeid(VarSwapBasis::VarSwapBasisProcError), load);


/////////////////////////////////////////////////////////////////////////////////////////


VarSwapBasis::VarSwapBasisRequest::VarSwapBasisRequest(const DateTimeArray& dates, YieldCurveConstSP yc):
CVolRequest(TYPE), dates(dates), yc(copy(yc.get())) {}
        

const DateTimeArray& VarSwapBasis::VarSwapBasisRequest::getDates() const {
    return dates;
}


YieldCurveConstSP VarSwapBasis::VarSwapBasisRequest::getAssetYC() const {
    return yc;
}


void VarSwapBasis::VarSwapBasisRequest::load(CClassSP& clazz) {
    REGISTER(VarSwapBasis::VarSwapBasisRequest, clazz);
    SUPERCLASS(CVolRequest);
    FIELD(dates, "Relevant maturities");
    FIELD(yc, "Asset's domestic ccy YC");
}


CClassConstSP const VarSwapBasis::VarSwapBasisRequest::TYPE =
CClass::registerClassLoadMethod("VarSwapBasis::VarSwapBasisRequest", typeid(VarSwapBasis::VarSwapBasisRequest), load);


/////////////////////////////////////////////////////////////////////////////////////////


VarSwapBasis::VarSwapBasis(): 
MarketObject(TYPE), 
tweakDeltaCutoff(true),
deltaCutoff(0.0), 
method(NONE), 
volBasisBMs(StringArray(1,"1Y")), 
volBasis(DoubleArray(1,0.0)){}


string VarSwapBasis::getName() const {
    return name;
}


double VarSwapBasis::interpVolBasis(double time) const {
    if(basisProc.get()) {
        return basisProc->interpVolBasis(time);
    } else {
        return 0.0;
    }
}


double VarSwapBasis::interpPriceCutoff(double time) const {
    if(basisProc.get()) {
        return basisProc->interpPriceCutoff(time);
    } else {
        return 0.0;
    }
}


double VarSwapBasis::interpSkewCutoff(double time) const {
    if(basisProc.get()) {
        return basisProc->interpSkewCutoff(time);
    } else {
        return 0.0;
    }
}


VarSwapBasis::VarSwapBasisProc* VarSwapBasis::setup(const CAsset* asset, 
                                                    const VolSurface* backbone, 
                                                    const VarSwapBasis::VarSwapBasisRequest* req, 
                                                    const string& tenorTime) {
    // Remove implementation and create a new one
    basisProc = VarSwapBasisProcSP(   );    
    basisProc = VarSwapBasisProcSP(new VarSwapBasisProc(this, asset, backbone, req, tenorTime));
    return basisProc.get();
}


void VarSwapBasis::validatePop2Object() {
    static const string routine("VarSwapBasis::validatePop2Object");
    try{
        // Same size for arrays
        if (volBasisBMs.size() != volBasis.size()){
            throw ModelException(routine,
                                 "the number of basis benchmarks should be equal to that of basis");
        }

        // OK if the sizes are zero let's convert them to size 1 just to be on the safe side
        if (volBasisBMs.size() == 0){
            volBasisBMs = StringArray(1, "1Y"); 
            volBasis    = DoubleArray(1, 0.0);
        }

        // Map DEFAULT to concrete methodology
        if (CString::equalsIgnoreCase(method, DEFAULT)){
            method = SKEW_CUTOFF;
        }

        // Map string to id and do some validation
        if(CString::equalsIgnoreCase(method, NONE)) {
            methodId = NONE_ID;
        } else if(CString::equalsIgnoreCase(method, SKEW_CUTOFF)) {
            methodId = SKEW_CUTOFF_ID;
        } else if(CString::equalsIgnoreCase(method, PRICE_CUTOFF)){ 
            methodId = PRICE_CUTOFF_ID;
        } else if(CString::equalsIgnoreCase(method, BASIS_ONLY)){
            methodId = BASIS_ONLY_ID;
            if(!Maths::isZero(deltaCutoff)) {
                throw ModelException(routine, 
                                     "Delta cutoff must be zero when using BASIS_ONLY");
            }
        } else {
            throw ModelException(routine,
                                "method must be DEFAULT, NONE, SKEW_CUTOFF, PRICE_CUTOFF or BASIS_ONLY");
        }

        // Verify deltaCutoff and change its sign internally
        if (Maths::isNegative(deltaCutoff) || !Maths::isPositive(1.0 - deltaCutoff)){
            throw ModelException(routine,
                                 "deltaCutoff must be between 0.0 (included) and 1.0 (excluded)");
        }
    } catch(exception& e){
        throw ModelException(e, routine);
    }
}


/** Invoked when Class is 'loaded' */
void VarSwapBasis::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VarSwapBasis, clazz);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultVarSwapBasis);
    FIELD(name, "VarSwapBasis identifier");
    FIELD(tweakDeltaCutoff, "Recompute delta cutoff in tweaks");
    FIELD_MAKE_OPTIONAL(tweakDeltaCutoff);
    FIELD(deltaCutoff, "Delta cutoff level");
    FIELD_MAKE_OPTIONAL(deltaCutoff);
    FIELD(method, "Cutoff methodology");
    FIELD_MAKE_OPTIONAL(method);
    FIELD(volBasisBMs, "Vol Basis BMs");
    FIELD_MAKE_OPTIONAL(volBasisBMs);
    FIELD(volBasis, "VolBasis");
    FIELD_MAKE_OPTIONAL(volBasis);
    FIELD(basisProc, "Basis processed");
    FIELD_MAKE_TRANSIENT(basisProc);
    FIELD(methodId, "");
    FIELD_MAKE_TRANSIENT(methodId);
    // by default don't get var swap basis
    MarketDataFetcher::setDefaultRetrievalMode(TYPE, false, NULL);
}


IObject* VarSwapBasis::defaultVarSwapBasis(){
    return new VarSwapBasis();
}


VarSwapBasis::~VarSwapBasis() {}


IObject* VarSwapBasis::clone() const {
    return MarketObject::clone();
}


CClassConstSP const VarSwapBasis::TYPE =
CClass::registerClassLoadMethod("VarSwapBasis", typeid(VarSwapBasis), load);


/* external symbol to allow class to be forced to be linked in */
bool VarSwapBasisLinkIn(){
    return true;
}
    
DRLIB_END_NAMESPACE
