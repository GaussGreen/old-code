//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : ProxyVol.cpp
//
//   Description : Volatility generated from a set of proxies (a composite vol)
//                 This used to be inside the Fund class but was factored out
//                 for wider use (i.e. with assets like ContangoCommodity)
//
//                  Vol surface is built from proxies as if we're doing an XCB
//                  with an optional spread applied followed by an optional skew
//                  Finally a minVol is used to floor the resultant vols
//
//   Author      : Ian S Stares
//
//   Date        : 27 September 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ProxyVol.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/StruckAsset.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/VolRequestDVF.hpp"
#include "edginc/CompositeVol.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/PDFParamLNStrike.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/ContangoRhoParallel.hpp"

DRLIB_BEGIN_NAMESPACE

// so Fund can build one inside itself (for now)
ProxyVol::ProxyVol(string&             name,
                   YieldCurveWrapper&  yc,
                   CAssetWrapperArray& proxies,
                   CStringArray&       ccyTreatments,
                   DoubleArray&        weights,
                   HolidayWrapper&     marketHols,
                   DateTime&           baseDate,
                   DoubleMatrix&       correlations,
                   TimeMetricSP&       timeMetric,
                   ExpiryArraySP&      spreadDates,
                   DoubleArray&        volSpread,
                   DoubleArray&        volSkew,
                   double              minVol,
                   bool                useCorrelationTerm) :
        CVolBase(TYPE), name(name), yc(yc), proxies(proxies), 
        ccyTreatments(ccyTreatments), weights(weights), marketHols(marketHols),
        baseDate(baseDate), correlations(correlations), timeMetric(timeMetric), 
        spreadDates(spreadDates), volSpread(volSpread), volSkew(volSkew), 
        minVol(minVol), useCorrelationTerm(useCorrelationTerm), 
        useVegaMatrix(false), useSkewParallel(false), useSkewPointwise(false), pdfBoundaryProb(CVolProcessedBS::LEGACY_PDF_BOUNDARY_PROB) {
    validatePop2Object();
}

/** Pull out the component assets & correlations from the market data */
void ProxyVol::getMarket(const IModel* model, const MarketData* market){
    static const string method("ProxyVol::getMarket");
    try {
        int numAssets = proxies.size();
        /* correlations: if we are explicitly supplied with just the numbers,
            then corr obj were built in validatePop2Object */

        if (!correlations.empty()){
            int pos = 0;
            for (int i = 0; i < proxies.size(); i++) {
                for (int j = i + 1; j < proxies.size(); j++, pos++) {
#if 0
                    corrObjects[pos]->configureForSensitivities(
                        proxies[i]->getClass(), proxies[j]->getClass());
#else
                    // force it to be EQ-EQ
                    corrObjects[pos]->configureForSensitivities(
                        CAsset::TYPE, CAsset::TYPE);
#endif
                    corrObjects[pos]->getMarket(model, market);
                }
            }
        } else {
            int pos = 0;
            corrObjects.resize((numAssets * numAssets - numAssets)/2);
            for (int i = 0; i < numAssets; i++) {                
                for (int j = i + 1; j < numAssets; j++, pos++) {
                    string corrName = 
                        market->getCorrelationName(proxies[i].getName(),
                                                       proxies[j].getName());
                    CorrelationBaseSP corr(model->getCorrelation(corrName,
#if 0
                        proxies[i]->getClass(),proxies[j]->getClass(),
#else
                        CAsset::TYPE, CAsset::TYPE, // force it to be EQ-EQ
#endif
                        Correlation::TYPE,market));
                     corrObjects[pos] = CorrelationSP::dynamicCast(corr);                    
                }
            }
        }

        if(useCorrelationTerm) {
        /* now: all corr category and corr term objects should be available.
            if not, then we wont fail, but assign a dummy CorrelationTerm */
            
            // allocate space
            corrTermArray.resize((numAssets * numAssets - numAssets)/2);

            int pos=0, i,j;
            for (i = 0; i < numAssets; i++) {
                if (market->hasData(proxies[i].getName(), CorrelationCategory::TYPE)) {
                    MarketObjectSP mo1 = 
                        market->GetData(proxies[i].getName(), CorrelationCategory::TYPE);
                    CorrelationCategorySP category1 = CorrelationCategorySP::dynamicCast(mo1);
                    for (j = i+1; j < numAssets; j++, pos++) {
                        if (market->hasData(proxies[j].getName(), CorrelationCategory::TYPE)) {
                            MarketObjectSP mo2 = 
                                market->GetData(proxies[j].getName(), CorrelationCategory::TYPE);
                            CorrelationCategorySP category2 = CorrelationCategorySP::dynamicCast(mo2);

                            const string& corrTermName = 
                                market->getCorrelationTermName(category1->getCategoryName(), category2->getCategoryName());
                            MarketObjectSP moTerm = market->GetData(corrTermName, CorrelationTerm::TYPE);
                            corrTermArray[pos] = CorrelationTermSP::dynamicCast(moTerm);
                            corrTermArray[pos]->getMarket(model, market);
                        } else {
                            // fill with dummy data
                            corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                        }
                    }
                } else {
                    for (j = i+1; j < proxies.size(); j++, pos++) {
                        // fill with dummy data
                        corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                    }
                }
            }
        }   

        // get the rest of the data we need from the cache
        market->GetReferenceDate(baseDate);
        yc.getData(model, market);
        for (int i = 0; i < proxies.size(); i++){
            CAsset::getAssetMarketData(model, market, ccyTreatments[i],
                                       yc.getName(), proxies[i]);
        }
        marketHols.getData(model, market);
        timeMetric->getMarket(model, market);
    } 
    catch (exception& e){
        throw ModelException(e, method, "Failed for fund " + getName());
    }
}

/** Validation */
void ProxyVol::validatePop2Object(){
    static const string method("ProxyVol::validatePop2Object");
    try {
        int numAssets = proxies.size();

        if (numAssets < 1) {
            throw ModelException(method,
                                 "Proxy Vols must have at least one component");
        }
        if (Maths::isNegative(minVol)) {
            throw ModelException(method,
                                 "minimum vol (" + Format::toString(minVol) + 
                                 ") is negative");
        }

        if (spreadDates->size() != volSpread.size() ||
            spreadDates->size() != volSkew.size()) {
            throw ModelException(method,
                                 "spread dates, vol spreads & skew spreads "
                                 "have inconsistent lengths");
        }

        // correlations: if we are explicitly supplied with just the numbers,
        // then generate correlation objects 
        if (!correlations.empty()){
            // validation of correlations
            if (correlations.numRows() != numAssets ||
                correlations.numCols() != numAssets){
                throw ModelException(method,
                    "Correlations must be an n x n matrix "
                    "with n = number of proxies");
            }
            correlations.checkSymmetric();
            // First allocate space
            corrObjects.resize((numAssets * numAssets - numAssets)/2);
            // then generate correlation objects
            int pos = 0, i,j;
            for (i = 0; i < numAssets; i++) {
                for (j = i + 1; j < numAssets; j++, pos++) {
                    // the Correlation constructor reorders the names into
                    // alphabetical order if no name supplied. For backwards
                    // compatibility we keep the original order
                    string corrName(proxies[i].getName() + "_" + 
                                    proxies[j].getName());
                    corrObjects[pos] = 
                        CorrelationSP(new Correlation(corrName,
                                                      proxies[i].getName(),
                                                      proxies[j].getName(),
                                                      correlations[i][j]));                    
                }
            }
        }

        // do our weights sum to 1 ?
        //AssetUtil::checkWeights(weights, proxies.size());
    }
    catch (exception& e){
        throw ModelException(e, method, "Failed for proxy vol " + getName());
    }
}

/* Returns a reference to the internal weights to use for combining assets */
const CDoubleArray& ProxyVol::weightsRef(const CAsset* asset) const {

    bool emptyWeights = adjWeights.empty();
    try {
        if (emptyWeights){
            adjWeights = DoubleArray(proxies.size());
        }
        // always recalc - spots may have changed
        double spot = asset->getSpot();
        for (int i = 0; i < proxies.size(); i++){
            adjWeights[i] = spot * weights[i]/proxies[i]->getSpot();
        }
        return adjWeights;
    }
    catch (exception& e){
        throw ModelException(e, "ProxyVol::weightsRef", "Failed to "
                             "calculate weights for ProxyVol "+getName());
    } 
}

/** Returns name of vol */
string ProxyVol::getName() const{
    return name;
}

//** calculate the vol interp scale factor for smile */
double ProxyVol::volInterpSmileScale(
    double assetSpot,   /* (I) the current asset spot price */
    int    i) const    /* (I) index of asset in compAsset */
{
    double spot = proxies[i]->getSpot();
    double scale = spot/assetSpot;
    return scale;
}

CVolRequestLNArraySP ProxyVol::getComponentVolRequests(
    const CVolRequest* volRequest,
    const CAsset* asset) const
{
    static const string  routine("ProxyVol::getComponentVolRequests");

    if (!CVolRequestLN::TYPE->isInstance(volRequest)){
        throw ModelException(routine,
                             "Only LN vol requests currently supported");
    }
    const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
        *(IObject*)(volRequest));
    
    // ensure no negative forward variance
    request.checkNegativeFwdVarDisallowed("ProxyVol::getComponentVolRequests");

    CVolRequestLNConstSP origInterp(
        CVolRequestLNConstSP::attachToRef(&request));

    // create an array of vol requests - one for each asset
    CVolRequestLNArraySP interps(new CVolRequestLNArray(proxies.size()));
    // cache the fund spot price
    double assetSpot = asset->getSpot();

    for (int i = 0; i < proxies.size(); i++) {
        /** copy each of the interps so that we can scale as needed */
        (*interps)[i] = CVolRequestLNSP(origInterp.clone());
        double scale = volInterpSmileScale(assetSpot, i);
        (*interps)[i]->scale(scale);
    }
    return interps;
}

// do we need to add any skew shift ?
bool ProxyVol::needToSkew() const {
    bool shift = false;
    for (int i = 0; !shift && i < spreadDates->size(); i++) {
        shift = !Maths::equals(1.0, volSkew[i]);
    }
    return shift;
}

// build a surface (pre-skew shift) with ATM & 90% strikes
VolSurfaceSP ProxyVol::volSkewATM(const CAsset* asset) const {
    static const string method = "ProxyVol::volSkewATM";
    try {
        double          spot = asset->getSpot();
        DoubleArray     strikes(2);

        strikes[0] = spot * VegaSkewParallel::IShift::SKEW_NORMALISATION;
        strikes[1] = spot;

        CDoubleMatrixSP matrix(interpVolMatrix(spreadDates, strikes, asset));

        VolSurfaceSP vol = VolSurfaceSP(new VolSurface(asset->getName(), 
                                                       timeMetric.get(), 
                                                       strikes, 
                                                       *matrix,
                                                       spreadDates.get(), 
                                                       baseDate));
        return vol;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// what's the 90-100 skew at a given expiry ?
double ProxyVol::measureSkew(Expiry* expiry, VolSurface* vol, const CAsset* asset) const {
    static const string method = "ProxyVol::measureSkew";
    try {
        double spot = asset->getSpot();
        double skew;
        DateTime skewDate = expiry->toDate(baseDate);

        CVolProcessedBSSP volATM(vol->getProcessedVol(spot));
        CVolProcessedBSSP vol90(vol->getProcessedVol(spot * 
                                                     VegaSkewParallel::IShift::SKEW_NORMALISATION));


        double vatm = volATM->CalcVol(baseDate, skewDate);
        double v90  = vol90->CalcVol(baseDate, skewDate);

        //skew = vatm - v90;
        skew = v90 - vatm;
        return skew;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** where are we trying to interpolate */
DoubleArraySP ProxyVol::interpLevels(const CVolRequest* volRequest, 
                                     const CAsset* asset) const {
    static const string method = "ProxyVol::interpLevels";
    try {
        if (!CVolRequestLN::TYPE->isInstance(volRequest)){
            throw ModelException(method,
                                 "Only LN vol requests currently supported");
        }
        const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
            *(IObject*)(volRequest));
        DoubleArraySP levels(new DoubleArray(0));
        DoubleArraySP unique(new DoubleArray(0));
        double        spot = asset->getSpot();
        request.getSensitiveStrike(spot, levels);
        sort(levels->begin(), levels->end());

        // strip duplicates
        if (!levels->empty()) {
            unique->push_back((*levels)[0]);
        }

        for (int i = 1; i < levels->size(); i++) {
            if ((*levels)[i] > (*levels)[i-1] + DBL_EPSILON) {
                unique->push_back((*levels)[i]);
            }
        }

        return unique;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** build vol surface */
VolSurfaceSP ProxyVol::volSurface(const CVolRequest* volRequest,
                                  const CAsset* asset) const {
    static const string method = "ProxyVol::volSurface";
    try {
        DoubleArray strikes = *(interpLevels(volRequest, asset));

        return volSurface(strikes, asset);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

VolSurfaceSP ProxyVol::volSurface(const DoubleArray& strikes,
                                  const CAsset* asset) const {
    static const string method = "ProxyVol::volSurface";
    try {
        CDoubleMatrixSP matrix(interpVolMatrix(spreadDates, strikes, asset));
        int i;
        int j;

        // is there a skew shift to add ?
        if (needToSkew()) {
            double spot = asset->getSpot();

            // build a surface with 90/100 strikes
            VolSurfaceSP skewVol = volSkewATM(asset);

            for (i = 0; i < spreadDates->size(); i++) {
                double skew = measureSkew((*spreadDates)[i].get(), 
                                          skewVol.get(), asset);

                for (j = 0; j < strikes.size(); j++) {
                    (*matrix)[j][i] += (volSkew[i] - 1.0) * 
                        VegaSkewParallel::skewShift(skew, strikes[j], spot);

                     // enforce minimum vol
                    (*matrix)[j][i] = Maths::max(minVol, (*matrix)[j][i]); 
                }
            }
        }

        // enforce minimum vol
        for (i = 0; i < spreadDates->size(); i++) {
            for (j = 0; j < strikes.size(); j++) {               
                (*matrix)[j][i] = Maths::max(minVol, (*matrix)[j][i]); 
            }
        }

        // if we're tweaking, add on any tweak spreads
        if (!sensSpread.empty()) {
            for (i = 0; i < sensSpread.size(); i++) {
                for (j = 0; j < strikes.size(); j++) {
                    (*matrix)[j][i] += sensSpread[i];
                }
            }
        }

        VolSurfaceSP vol = VolSurfaceSP(new VolSurface(asset->getName(), 
                                                       timeMetric.get(), 
                                                       strikes, 
                                                       *matrix,
                                                       spreadDates.get(), 
                                                       baseDate));
        return vol;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DoubleArray ProxyVol::interpVolCurve(
    const CVolRequest* volRequest,
    CompositeVol*      compVol,
    const CAsset*      asset) const // can be null for single asset
{
    static const string routine("ProxyVol::interpVolCurve");
    try {
        const DoubleArray& volWeights = weightsRef(asset);
        // get the vol requests for the proxies
        CVolRequestLNArraySP interps = getComponentVolRequests(volRequest, asset);
        
        CVolProcessed* rawVol = 0;

        // scale the composite vol by the total weight - lazy way to give
        // a vol scale (as well as shift/skew) without changing the interface
        double totalWeight = 0.0;
        for (int j = 0; j < weights.size(); j++) {
            totalWeight += weights[j];
        }

        // if only one asset, use its vol directly, 
        // otherwise build a composite vol
        if (proxies.size() == 1) {
            rawVol = proxies[0]->getProcessedVol((*interps)[0].get());
        } else {
            CorrTermDataSP data = CorrelationTerm::getCorrelationTermSqueezesAndExpiries(
                baseDate, proxies.size(),
                corrObjects, corrTermArray);

            int numComponents = proxies.size();
            DoubleMatrix corrMatrix(numComponents, numComponents);
            int pos = 0;
            for (int i = 0; i < proxies.size(); i++) {
                corrMatrix[i][i] = 1.0;
                for (int j = i + 1; j < proxies.size(); j++, pos++) {
                    if (Correlation::TYPE->isInstance(corrObjects[pos].get()))
                    {
                        CorrelationSP p = CorrelationSP::dynamicCast(corrObjects[pos]);
                        corrMatrix[i][j] = p->getCorrelation();
                        corrMatrix[j][i] = corrMatrix[i][j];
                    }
                    else
                    {
                        // TO DO: handle the case the correlation is post calibrated
                        throw ModelException(routine, "correlation is calibrated");
                    }
                }
            }

            data->correlationMatrix = CDoubleMatrixSP(new DoubleMatrix(corrMatrix));

            rawVol = compVol->getProcessedVol(volRequest,
                                                   asset,
                                                   false,  // no fwd start
                                                   baseDate,
                                                   asset->getName(),
                                                   timeMetric.get(),
                                                   (*interps),
                                                   volWeights,
                                                   *data->correlationMatrix, 
                                                   *data->benchmarkTermExpiries,
                                                   0.0,
                                                   0.0,
                                                   *data->shortTermSqueeze, 
                                                   *data->shortTermExpiries,
                                                   *data->longTermSqueeze,
                                                   *data->longTermExpiries);
        }
        CVolProcessedSP rawVolSP(rawVol);

        // cast to the type of vol we're expecting
        CVolProcessedBSSP rawVolBS = CVolProcessedBSSP::dynamicCast(rawVolSP);
        if (!rawVolBS){
            throw ModelException(routine, "Not a Black Scholes vol");
        }
         
        // convert the expiries to dates for the vol calculations
        DateTimeArray bmDates(spreadDates->size());
        int i;
        for (i = 0; i < bmDates.size(); i++) {
            bmDates[i] = (*spreadDates)[i]->toDate(baseDate);
        }
        // interpolate at each of the spreadDates for the vol spread
        DoubleArray interpVols(bmDates.size());
        for (i = 0; i < interpVols.size(); i++) {
            interpVols[i] = rawVolBS->CalcVol(baseDate, bmDates[i]);
            // scale by total weight
            interpVols[i] *= totalWeight;
            // add that spread
            interpVols[i] += volSpread[i];
       }
        return interpVols;
    } 
    catch (exception& e){
        throw ModelException(e, routine,
                             "failed to build composite vol " +
                             getName());
    }
}

CDoubleMatrixSP ProxyVol::interpVolMatrix(
    const ExpiryArraySP benchmarks,
    const DoubleArray&  strikes,
    const CAsset* asset) const 
{
    static const string method = "ProxyVol::interpVolMatrix";
    try {
        int numBM = benchmarks->size();
        int numK  = strikes.size();

        CDoubleMatrixSP matrix(new DoubleMatrix(numK, numBM));

        // convert the expiries to dates for the vol calculations
        DateTimeArray bmDates(numBM);
        int i;
        int j;
        for (i = 0; i < numBM; i++)
        {
            bmDates[i] = (*benchmarks)[i]->toDate(baseDate);
        }
        // cache asset info across multiple vol calculations
        auto_ptr<CompositeVol> compositeVol(proxies.size() > 1?
                                            (new CompositeVol(proxies)): 0);

        for (i = 0; i < numK; i++) {
            // interpolate at each strike
            CVolRequestSP volRequest(new LinearStrikeVolRequest(strikes[i],
                                                                baseDate,
                                                                baseDate,
                                                                false));

            // interpolate the composite vol
            DoubleArray vols(interpVolCurve(volRequest.get(), 
                                            compositeVol.get(), asset));
                          
            for (j = 0; j < numBM; j++) {
                (*matrix)[i][j] = vols[j];
            }
        }
        return matrix;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Generate a Proxy vol surface using 'default' strikes */
VolSurfaceSP ProxyVol::defaultVolSurface(const CAsset* asset) const{
    if (pdfStrikes.empty()){
        // to avoid instability in the greeks reuse the same strikes
        // from the original pricing run. We assume that this method
        // will be called from the pricing run if it is called during
        // the tweaks
        pdfStrikes = CVolProcessedBS::defaultStrikes(baseDate, asset, this);
    }
    // then create a vol surface with these strike
    VolSurfaceSP myVolSurface(volSurface(pdfStrikes, asset));
    // are we doing vega skew ?
    applyVegaSkew(myVolSurface.get());
    return myVolSurface;
}

/** if we're building an actual VolSurface then need to handle case
    where we're operating under a vega skew sensitivity and need to apply 
    this skew shift to our generated surface*/
void ProxyVol::applyVegaSkew(VolSurface* vol) const {
    // are we doing vega skew ?
    if (useSkewParallel) {
        vol->sensShift(skewParallel.get());
    }
    if (useSkewPointwise) {
        vol->sensShift(skewPointwise.get());
    }
}    

/** Combines market and instrument data together to give a Processed Vol */
CVolProcessed * ProxyVol::getProcessedVol(const CVolRequest* volRequest, 
                                          const CAsset* asset) const
{
    static const string routine("ProxyVol::getProcessedVol");
    try {
        if (CVolRequestDVF::TYPE->isInstance(volRequest)){
            if (useVegaMatrix) {
                // this shouldn't happen. Certainly not clear what we should do
                throw ModelException(routine, "Internal error - doing "
                                     "strikewise vega with local vol");
            }
            // to support local vol we choose a set of strikes and build a
            // vol surface for those strikes. We then spline that surface 
            // and use that spline to support the local vol.
            // Create a vol surface with these strike
            VolSurfaceSP myVolSurface(defaultVolSurface(asset));
            // Pass volRequest to surface (it will spline surface)
            return myVolSurface->getProcessedVol(volRequest, asset);
        }

        // are we doing vega matrix ?
        if (useVegaMatrix) {
            VolSurfaceSP derivedVol = buildVegaMatrixVol(asset);
            return derivedVol->getProcessedVol(volRequest, asset);
        }

        VolSurfaceSP vol(volSurface(volRequest, asset));

        // are we doing vega skew ?
        applyVegaSkew(vol.get());

        return VolSurface::getProcessedVol(vol, volRequest, asset);
    } 
    catch (exception& e){
        throw ModelException(e, routine,
                             "failed to build composite vol for proxyVol " +
                             name);
    }
}

/** Combines market and instrument data together to give a
    Processed Vol. Here the processed volatility is a processed
    struck volatility ie it reflects the combination of this
    asset together with the supplied FX asset and the
    correlation between this CVolBase and the vol of the
    FX. Note that the struckAsset is indeed the struckAsset cf
    'this' which is the non struck asset */
CVolProcessed* ProxyVol::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      struckAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    static const string routine("ProxyVol::getProcessedVol");
    try {
        // Get hold of asset inside struck asset
        CAssetConstSP plainAsset;
        if (StruckAsset::TYPE->isInstance(struckAsset)) {
            StruckAssetConstSP struck(dynamic_cast<const StruckAsset*>(struckAsset));
            plainAsset = struck->getPlainAsset();
        } else if (StruckEquity::TYPE->isInstance(struckAsset)) {
            StruckEquityConstSP struck(dynamic_cast<const StruckEquity*>(struckAsset));
            plainAsset = struck->getPlainAsset();
        }

        // are we doing vega matrix ? 
        if (useVegaMatrix) {
            VolSurfaceSP derivedVol = buildVegaMatrixVol(plainAsset.get());
            return derivedVol->getProcessedVol(volRequest, struckAsset, 
                                            fxAsset, eqFXCorr);
        }

        // Need to scale the volRequest into our units **
        if (!CVolRequestLN::TYPE->isInstance(volRequest)){
            throw ModelException(routine, "Only LN vol supported");
        }
        const CVolRequestLN& lnRequest = 
            dynamic_cast<const CVolRequestLN&>(*volRequest);
        CVolRequestLNSP volRequestCopy(copy(&lnRequest)); 
        double fxSpot = fxAsset->getSpot();
        volRequestCopy->scale(1.0/fxSpot);

        // then can build our surface
        VolSurfaceSP vol(volSurface(volRequestCopy.get(), plainAsset.get()));

        // are we doing vega skew ?
        applyVegaSkew(vol.get());

        return vol->getProcessedVol(volRequest, struckAsset, 
                                    fxAsset, eqFXCorr);
    } 
    catch (exception& e){
        throw ModelException(e, routine,
                             "failed to build composite vol for " + name);
    }
}
   
/** handles the delta shift size adjustment */
void ProxyVol::adjustDeltaShiftSize(ShiftSizeCollector* collector,
                                    const string assetName,
                                    double spot) const {
    try {
        for (int i = 0; i < proxies.size(); i++) {
            proxies[i]->accept(collector);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "ProxyVol::adjustDeltaShiftSize",
                             "failed for vol " + getName());
    } 
}

/** Shifts the object using given shift (see Theta::Shift)*/
bool ProxyVol::sensShift(Theta* shift) 
{
    try {
        pdfStrikes.clear(); // use fresh strikes when doing theta
        baseDate = shift->rollDate(baseDate); 
    }
    catch (exception& e) {
        throw ModelException(e, "ProxyVol::sensShift (theta)"); 
    }
    return true; // tweak stuff inside
}


/*returns sensitive strikes for a given vol request and asset*/
void ProxyVol::getSensitiveStrikes(const CAsset* asset,
                               const CVolRequest* volRequest,
                               OutputNameConstSP outputName,
                               const SensitiveStrikeDescriptor& sensStrikeDesc,
                               DoubleArraySP sensitiveStrikes) const
{
    static const string method = "ProxyVol::getSensitiveStrikes";
    try {
        // get the vol requests for the proxies
        CVolRequestLNArraySP interps = getComponentVolRequests(volRequest, asset);
        
        for (int i = 0; i < proxies.size(); i++) {
            // let the component take care of its sensitive strikes
            proxies[i]->getSensitiveStrikes((*interps)[i].get(),
                                                outputName, 
                                                sensStrikeDesc,
                                                sensitiveStrikes);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method); 
    }        
}

/** record forwards at maturity */
void ProxyVol::recordFwdAtMat(OutputRequest*  request,
                          CResults*       results,
                          const DateTime& maturityDate) const
{
    // record forward at maturity for each proxy
    for (int i = 0; i < proxies.size() ; ++i)
    {
        proxies[i]->recordFwdAtMat(request, results, maturityDate);
    }
}

void ProxyVol::acceptCollector(const ProxyVol* proxyVol, StartDateCollector* collector)
{
    // check start dates for each component
    for (int i = 0; i < proxyVol->proxies.size() ; ++i)
    {
        proxyVol->proxies[i]->accept(collector);
    }
}

void ProxyVol::acceptNameCollector(const ProxyVol* proxyVol, AssetNameCollector* collector)
{
//    collector->assetNameValidate(asset->getName());
    try {
        // and check myself internally
        AssetNameCollectorSP visitor(new AssetNameCollector());
        
        for (int i = 0; i < proxyVol->proxies.size() ; ++i)
        {
            proxyVol->proxies[i]->accept(visitor.get());
        }
    }
    catch (exception& e) {
        throw ModelException(e, "ProxyVol::acceptNameCollector",
                             "failed for vol " + proxyVol->getName());           
    }
}

void ProxyVol::acceptFutureCollector(const ProxyVol* proxyVol,
                                 FutureExpiryCollector* collector)
{
    // check future expiry date for each component
    for (int i = 0; i < proxyVol->proxies.size() ; ++i)
    {
        proxyVol->proxies[i]->accept(collector);
    }
}

void ProxyVol::acceptValueDateCollector(const ProxyVol* proxyVol, 
                                    CValueDateCollector* collector)
{
    collector->valueDateValidate(proxyVol->baseDate, proxyVol->getName());
}

void ProxyVol::acceptImntCcy(const ProxyVol*        proxyVol,
                         AssetCcyCollector* collector)
{
    collector->currencyValidate(proxyVol->yc->getCcy(),
                                proxyVol->getName());
    
    try {
        for (int i = 0; i < proxyVol->proxies.size(); ++i) {
            IObjectConstSP component(
                IObjectConstSP::attachToRef(proxyVol->proxies[i].get()));
            
            AssetCcyCollector::
                validateAllCurrencies(component, proxyVol->yc->getCcy());
        }
    }
    catch (exception& e) {
        throw ModelException(e, "ProxyVol::acceptImntCcy",
                             "failed for vol " + proxyVol->getName());
    }
}

void ProxyVol::acceptHoliday(const ProxyVol*       proxyVol,
                         HolidayCollector* collector)
{
    collector->setHoliday(proxyVol->marketHols.getSP());
}

/** Returns name identifying vol for vega parallel */
string ProxyVol::sensName(const VolParallel*) const{
    return getName();
}

/** Shifts the object using given shift */
TweakOutcome ProxyVol::sensShift(const PropertyTweak<VolParallel>& tweak){
    static const string method("ProxyVol::sensShift");
    if (!Maths::isZero(tweak.coefficient)){
        // simply store shift so can apply when create proxy vol
        sensSpread = DoubleArray(spreadDates->size());
        for (int i = 0; i < spreadDates->size(); i++) {
            sensSpread[i] = tweak.coefficient;
        }
    }
    return TweakOutcome(tweak.coefficient, false);
}

/** Restores the object to its original form */
void ProxyVol::sensRestore(const PropertyTweak<VolParallel>& tweak) {
    // simply zap cache of shifts
    sensSpread.clear();
}

/** Returns name identifying vol for vega pointwise */
string ProxyVol::sensName(const VolPointwise*) const{
    return getName();
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryWindowArrayConstSP ProxyVol::sensQualifiers(const VolPointwise*) const{
    return ExpiryWindow::series(spreadDates);
}

/** Shifts the object using given shift */
TweakOutcome ProxyVol::sensShift(const PropertyTweak<VolPointwise>& shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        if (!Maths::isZero(shift.coefficient)) {
            // simply store shift so can apply when create proxy vol
            sensSpread = DoubleArray(spreadDates->size());

            for (int i = 0; i < spreadDates->size(); i++) {
                sensSpread[i] = 0.0;
            }

            int expiryIdx = shift.qualifier->expiry->search(spreadDates.get());

            sensSpread[expiryIdx] = shift.coefficient;
        }

        return TweakOutcome(shift.coefficient, false);
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Restores the object to its original form */
void ProxyVol::sensRestore(const PropertyTweak<VolPointwise>&) {
    // simply zap cache of shifts
    sensSpread.clear();
}

/** Returns name identifying vol for vega matrix */
string ProxyVol::sensName(VegaMatrix* shift) const{
    return getName();
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP ProxyVol::sensExpiries(VegaMatrix* shift) const{
    return spreadDates;
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VegaMatrix* shift) {
    static const string method = "ProxyVol::sensShift";
    try {
        // retrieve shift information from MatrixShift 
        const ExpiryArray& shiftExpiries = *(shift->getExpiries());
        int                expiryIdx     = shift->getExpiryIdx();

        // check that the expiry to bump matches ours (it always should since
        // we supplied the expiries 
        if (expiryIdx >= spreadDates->size() ||
            !shiftExpiries[expiryIdx]->equals((*spreadDates)[expiryIdx].get())){
            throw ModelException(method, "Maturity "+
                                 shiftExpiries[expiryIdx]->toString()+
                                 " unexpected");
        }

        // copy the shift, we'll build it later once we've got the asset
        vegaMatrix = VegaMatrixSP(copy(shift));
        useVegaMatrix = true;
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

VolSurfaceSP ProxyVol::buildVegaMatrixVol(const CAsset* asset) const {
    static const string method = "ProxyVol::buildVegaMatrixVol";
    try {
        // retrieve shift information from MatrixShift 
        const DoubleArray& shiftAxis = *(vegaMatrix->getXAxisValues());

        VolSurfaceSP derivedVol = volSurface(shiftAxis, asset);

        derivedVol->sensShift(vegaMatrix.get());

        return derivedVol;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns name identifying vol for root time vega */
string ProxyVol::sensName(RootTimeVega* shift) const{
    return getName();
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(RootTimeVega* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            // simply store shift so can apply when create fund vol
            sensSpread = DoubleArray(spreadDates->size());

            for (int i = 0; i < spreadDates->size(); i++) {
                DateTime bm = (*spreadDates)[i]->toDate(baseDate);
                double   years = timeMetric->yearFrac(baseDate, bm);

                sensSpread[i] = shift->rtVegaShift(years);
            }
        }
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

/** Restores the object to its original form */
void ProxyVol::sensRestore(RootTimeVega* shift) {
    // simply zap cache of shifts
    sensSpread.clear();
}

string ProxyVol::sensName(VegaSkewParallel* shift) const{
    return getName();
}

bool ProxyVol::sensShift(VegaSkewParallel* shift){
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        // only bother if non zero 
        useSkewParallel = true;
        skewParallel = VegaSkewParallelSP(copy(shift));
    }
    return false; // no more tweaking required here
}

string ProxyVol::sensName(VegaSkewPointwise* shift) const{
    return getName();
}

ExpiryArrayConstSP ProxyVol::sensExpiries(VegaSkewPointwise* shift) const{
    return spreadDates;
}


/** Apply a skew shift to a specific benchmark. Vols for all strikes are
 * shifted by -shift * log(strike/spot)/log(1.1) */
bool ProxyVol::sensShift(VegaSkewPointwise* shift){
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        // only bother if non zero 
        useSkewPointwise = true;
        skewPointwise = VegaSkewPointwiseSP(copy(shift));
    }
    return false; // no more tweaking required here
}

/** Returns name identifying vol */
string ProxyVol::sensName(VolParallelShift* shift) const{
    return getName();
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VolParallelShift* shift){
    static const string method("ProxyVol::sensShift");
    double shiftSize = shift->getShiftSize();
    if (!Maths::isZero(shiftSize)){
        // bump the vol shift
        for (int i = 0; i < spreadDates->size(); i++) {
            volSpread[i] += shiftSize;
        }
    }
    return false;
}

/** Returns name identifying vol */
string ProxyVol::sensName(VolBenchmarkShift* shift) const{
    return getName();
}


/** Shifts the object using given shift */
bool ProxyVol::sensShift(VolBenchmarkShift* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            bool found = false;
            // bump the vol shift
            for (int i = 0; i < spreadDates->size(); i++) {
                if (!found && shift->expiryEquals((*spreadDates)[i].get())) {
                    volSpread[i] += shiftSize;
                    found = true;
                }
            }

            if (!found) {
                throw ModelException(routine, "benchmark not found "
                                     " on proxy vol surface "+ getName());
            }
        }
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

/** Returns name identifying vol for power vega */
string ProxyVol::sensName(PowerVega* shift) const{
    return getName();
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(PowerVega* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {

            for (int i = 0; i < spreadDates->size(); i++) {
                DateTime bm = (*spreadDates)[i]->toDate(baseDate);
                double   years = timeMetric->yearFrac(baseDate, bm);

                volSpread[i] += shift->powerShift(years);
            }
        }
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

// Returns name identifying yield curve for rho parallel **
// only want rho to 'fund' ccy 
string ProxyVol::sensName(const RateParallel* shift) const {
    return yc->getCcy();
}

// Shifts the object using given shift **
TweakOutcome ProxyVol::sensShift(const PropertyTweak<RateParallel>& shift){
    PropertyTweakHypothesis<RateParallel> hyp(shift);

    for (int i = 0; i < proxies.size(); i++) {
        hyp.applyTo(proxies[i].getSP());
    }

    return TweakOutcome(shift.coefficient, true); // fund has rho
}

/** Returns the name of the yield curve - used to determine whether 
    to tweak the object 
   only want rho to 'fund' ccy */ 
string ProxyVol::sensName(const RatePointwise* shift) const{
    return yc->getCcy();
}
   
/* Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryWindowArrayConstSP ProxyVol::sensQualifiers(const RatePointwise* shift) const {
    IObjectConstSP world = IObjectConstSP(yc.get());
    return RiskProperty<RatePointwise>().subjectQualifiers(
        world, OutputName::SP(yc->getCcy()));
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the proxies within the object
    which implements this interface */
TweakOutcome ProxyVol::sensShift(const PropertyTweak<RatePointwise>& shift) {
    PropertyTweakHypothesis<RatePointwise> hyp(shift);
    for (int i = 0; i < proxies.size(); i++) {
        hyp.applyTo(proxies[i].getSP());
    }
    return TweakOutcome(shift.coefficient, true); // fund has rho
}

/** Returns name identifying vol for VolRelativeShift */
string ProxyVol::sensName(VolRelativeShift* shift) const{
    return getName();
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VolRelativeShift* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        // shift is ABSOLUTE based on a RELATIVE shift to the vol 
        // interpolated at "getLevel" (e.g. ATM vol)

        for (int i = 0; i < spreadDates->size(); i++) {
            double shiftSize = shift->shiftSize(baseDate,(*spreadDates)[i]->toDate(baseDate));
            volSpread[i] += shiftSize;
        }
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

/** Returns name identifying vol for VolAbsoluteShift */
string ProxyVol::sensName(VolAbsoluteShift* shift) const{
    return getName();
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VolAbsoluteShift* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        for (int i = 0; i < spreadDates->size(); i++) {
            double shiftSize = shift->shiftSize(baseDate,
                                                (*spreadDates)[i]->toDate(baseDate));
            volSpread[i] += shiftSize;
        }
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

// Appends the market data names inside supplied proxy to the passed names array
void ProxyVol::proxySensNames(SensControlPerName* sens, const Asset* proxy, OutputNameArray& sensNames) const {
    static const string method = "ProxyVol::proxySensNames";
    try {
        // For struck/protected proxies will get name(s) inside plain asset aswell as name inside
        // FXAsset
        // For vol proxies will get names inside these vol proxies unless sensitivity 
        // is suppressed for proxies (most cases).
        // For xcb proxies will get names inside xcb proxies.
        OutputNameArrayConstSP names(sens->names(proxy));
        
        for (int i=0; i<names->size(); i++) {
            sensNames.push_back((*names)[i]);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Returns the expiries of the target market data (in sens) found in the fund proxies.
ExpiryArrayConstSP ProxyVol::proxySensExpiries(VectorShift *sens) const {
    static const string method = "ProxyVol::proxySensExpiries";
    try {
        OutputNameArray sensNames(0);
        for (int i = 0; i < proxies.size(); i++) {
            // Collect names so that can get expiries given a match
            proxySensNames(sens, proxies[i].get(), sensNames);
            for (int j=0; j < sensNames.size(); j++) {
                if (sens->getMarketDataName()->equals(sensNames[j].get())) {
                    return sens->getExpiries(proxies[i].get());
                }
            }
            sensNames.clear();
        }  
        return ExpiryArrayConstSP(new ExpiryArray(0));
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns true if this object matches the supplied name with 
    respect to the DeltaProxy sensitivity */
bool ProxyVol::sensNameMatches(DeltaProxy*       shift,
                           const OutputName& name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameArray deltaNames(0);
        Delta delta(Delta::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&delta, proxies[i].get(), deltaNames);
            for (int j=0; j < deltaNames.size(); j++) {
                if (name.equals(deltaNames[j].get())) {
                    return true;
                }
            }
            deltaNames.clear();
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/*Appends the name(s) of this object with respect to
    the DeltaProxy sensitivity to the supplied list. */
void ProxyVol::sensAppendName(DeltaProxy*      shift,
                          OutputNameArray& namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        Delta delta(Delta::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&delta, proxies[i].get(), namesList);
        }  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }        
}

/* Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the proxies within the object
    which implements this interface */
bool ProxyVol::sensShift(DeltaProxy* shift) {
    OutputNameArray deltaNames(0);
    // Use a new delta in proxySensNames() because toTweak proxies passed from DeltaProxy
    // to delta in makeDeltaShift may interfere with name collection
    DeltaSP dummyDelta(new Delta(Delta::DEFAULT_SHIFT));

    for (int i = 0; i < proxies.size(); i++) {
        // Collect names so that can set initialValue given a match
        proxySensNames(dummyDelta.get(), proxies[i].get(), deltaNames);

        for (int j=0; j < deltaNames.size(); j++) {
            if (shift->marketDataNameMatches(deltaNames[j])) {
                TweakOutcome outcome = PropertyTweakHypothesis<Spot>(
                    shift->getShiftSize(), shift->getMarketDataName()).
                        applyTo_TweakOutcome(proxies[i].getSP());
                shift->setInitialValue(outcome.oldValue());
            }
        }
        deltaNames.clear();
    }
    return false;
}
   
/** Returns true if this object matches the supplied name with 
    respect to the VegaProxyParallel sensitivity */
bool ProxyVol::sensNameMatches(VegaProxyParallel* shift,
                           const OutputName&  name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        for (int i = 0; i < proxies.size(); i++) {
            OutputNameArrayConstSP vegaNames(
                RiskProperty<VolParallel>().subjectNames(proxies[i].getSP()));

            for (int j=0; j < vegaNames->size(); j++) {
                if (name.equals((*vegaNames)[j].get())) {
                    return true;
                }
            }
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Appends the name(s) of this object with respect to
    the VegaProxyParallel sensitivity to the supplied list. */
void ProxyVol::sensAppendName(VegaProxyParallel* shift,
                          OutputNameArray&   namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        for (int i = 0; i < proxies.size(); i++) {
            OutputName::append(namesList,
                               RiskProperty<VolParallel>().subjectNames(
                                   proxies[i].getSP()));
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/* Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the proxies within the object
    which implements this interface */
bool ProxyVol::sensShift(VegaProxyParallel* shift) {
    PropertyTweakHypothesis<VolParallel> v(
        shift->getShiftSize(), shift->getMarketDataName());

    for (int i = 0; i < proxies.size(); i++) {
        v.applyTo(proxies[i].getSP());
    }
    return false;
}

/** Returns true if this object matches the supplied name with 
    respect to the VegaProxyPointwise sensitivity */
bool ProxyVol::sensNameMatches(VegaProxyPointwise* shift,
                           const OutputName&   name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        for (int i = 0; i < proxies.size(); i++) {
            OutputNameArrayConstSP vegaNames(
                RiskProperty<VolParallel>().subjectNames(
                    proxies[i].getSP()));

            for (int j=0; j < vegaNames->size(); j++) {
                if (name.equals((*vegaNames)[j].get())) {
                    return true;
                }
            }
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Appends the name(s) of this object with respect to
    the VegaProxyPointwise sensitivity to the supplied list. */
void ProxyVol::sensAppendName(VegaProxyPointwise* shift,
                          OutputNameArray&   namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        for (int i = 0; i < proxies.size(); i++) {
            OutputName::append(namesList,
                               RiskProperty<VolPointwise>().subjectNames(
                                   proxies[i].getSP()));
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
} 

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP ProxyVol::sensExpiries(VegaProxyPointwise* shift) const{
    static const string method = "ProxyVol::sensExpiries";
    try {
        for (int i = 0; i < proxies.size(); ++i) {
            OutputNameArrayConstSP names(
                RiskProperty<VolPointwise>().subjectNames(proxies[i].getSP()));
            for (int n = 0; n < names->size(); ++n) {
                if ((*names)[n]->equals(shift->getMarketDataName().get())) {
                    return ExpiryWindow::expiries(
                        RiskProperty<VolPointwise>().subjectQualifiers(
                            proxies[i].getSP(), (*names)[n]));
                }
            }
        }

        return ExpiryArrayConstSP(new ExpiryArray());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VegaProxyPointwise* proxy){

    static const string routine = "ProxyVol::sensShift";
    try {
        PropertyTweakHypothesis<VolPointwise> tweak(proxy->getShiftSize(),
                                                    proxy->getMarketDataName(),
                                                    proxy->getExpiryWindow());

        for (int i = 0; i < proxies.size(); i++) {
            tweak.applyTo(proxies[i].getSP());
        }
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

/** Returns true if this object matches the supplied name with 
    respect to the sensitivity */
bool ProxyVol::sensNameMatches(VegaProxyMatrix*  shift,
                           const OutputName& name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameArray vegaNames(0);
        VegaMatrix vega(VegaMatrix::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), vegaNames);
            for (int j=0; j < vegaNames.size(); j++) {
                if (name.equals(vegaNames[j].get())) {
                    return true;
                }
            }
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Appends the name(s) of this object with respect to
    the sensitivity to the supplied list. */
void ProxyVol::sensAppendName(VegaProxyMatrix* shift,
                          OutputNameArray& namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        VegaMatrix vega(VegaMatrix::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), namesList);
        }  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
} 

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP ProxyVol::sensExpiries(VegaProxyMatrix* shift) const{
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameConstSP  name = shift->getMarketDataName();
        VegaMatrix         vega(shift->getShiftSize());
        vega.setMarketDataName(name);
        OutputNameArray vegaNames(0);
        for (int i = 0; i < proxies.size(); i++) {
            // Collect names so that can get expiries given a match
            proxySensNames(&vega, proxies[i].get(), vegaNames);
            for (int j=0; j < vegaNames.size(); j++) {
                if (vega.getMarketDataName()->equals(vegaNames[j].get())) {
                    return vega.getExpiries(proxies[i].get());
                }
            }
            vegaNames.clear();
        }  
        return ExpiryArrayConstSP(new ExpiryArray(0));
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VegaProxyMatrix* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        VegaMatrixSP vega(VegaMatrix::fromProxy(shift));
        OutputNameConstSP name = vega->getMarketDataName();
        for (int i = 0; i < proxies.size(); i++) {
            vega->findAndShift(proxies[i].getSP(), name);
        }
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

/** Returns true if this object matches the supplied name with 
    respect to the sensitivity */
bool ProxyVol::sensNameMatches(VegaSkewProxyParallel* shift,
                           const OutputName&      name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameArray vegaNames(0);
        VegaSkewParallel vega(VegaSkewParallel::DEFAULT_SHIFT);

        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), vegaNames);

            for (int j=0; j < vegaNames.size(); j++) {
                if (name.equals(vegaNames[j].get())) {
                    return true;
                }
            }
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Appends the name(s) of this object with respect to
    the sensitivity to the supplied list. */
void ProxyVol::sensAppendName(VegaSkewProxyParallel* shift,
                          OutputNameArray&       namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        VegaSkewParallel vega(VegaSkewParallel::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), namesList);
        }  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the proxies within the object
    which implements this interface */
bool ProxyVol::sensShift(VegaSkewProxyParallel* shift) {
    VegaSkewParallel vega(shift->getShiftSize());
    OutputNameConstSP name = shift->getMarketDataName();

    for (int i = 0; i < proxies.size(); i++) {
        vega.findAndShift(proxies[i].getSP(), name);
    }
    return false;
}

/** Returns true if this object matches the supplied name with 
    respect to the sensitivity */
bool ProxyVol::sensNameMatches(VegaSkewProxyPointwise* shift,
                           const OutputName&       name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameArray vegaNames(0);
        VegaSkewPointwise vega(VegaSkewPointwise::DEFAULT_SHIFT);

        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), vegaNames);

            for (int j=0; j < vegaNames.size(); j++) {
                if (name.equals(vegaNames[j].get())) {
                    return true;
                }
            }
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/*Appends the name(s) of this object with respect to
    the sensitivity to the supplied list. */
void ProxyVol::sensAppendName(VegaSkewProxyPointwise* shift,
                          OutputNameArray&        namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        VegaSkewPointwise vega(VegaSkewPointwise::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), namesList);
        }  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
} 

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP ProxyVol::sensExpiries(VegaSkewProxyPointwise* shift) const{
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameConstSP  name = shift->getMarketDataName();
        VegaSkewPointwise  vega(shift->getShiftSize());
        vega.setMarketDataName(name);
        return proxySensExpiries(&vega);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
}

/** Shifts the object using given shift */
bool ProxyVol::sensShift(VegaSkewProxyPointwise* shift){
    static const string routine = "ProxyVol::sensShift";
    try {
        VegaSkewPointwiseSP vega(VegaSkewPointwise::fromProxy(shift));
        OutputNameConstSP name = vega->getMarketDataName();
        for (int i = 0; i < proxies.size(); i++) {
            vega->findAndShift(proxies[i].getSP(), name);
        }
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
    return false;
}

/** Returns true if this object matches the supplied name with 
    respect to the sensitivity */
bool ProxyVol::sensNameMatches(RootTimeVegaProxy* shift,
                           const OutputName&  name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameArray vegaNames(0);
        RootTimeVega vega(RootTimeVega::DEFAULT_SHIFT);

        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&vega, proxies[i].get(), vegaNames);

            for (int j=0; j < vegaNames.size(); j++) {
                if (name.equals(vegaNames[j].get())) {
                    return true;
                }
            }
        }

        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Appends the name(s) of this object with respect to
    the sensitivity to the supplied list. */
void ProxyVol::sensAppendName(RootTimeVegaProxy* shift,
                          OutputNameArray&   namesList) const {
    static const string method = "ProxyVol::sensAppendName";
    try {
        RootTimeVega vega(RootTimeVega::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
           proxySensNames(&vega, proxies[i].get(), namesList);
        }  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the proxies within the object
    which implements this interface */
bool ProxyVol::sensShift(RootTimeVegaProxy* shift) {
    RootTimeVega vega(shift->getShiftSize());
    OutputNameConstSP name = shift->getMarketDataName();

    for (int i = 0; i < proxies.size(); i++) {
        vega.findAndShift(proxies[i].getSP(), name);
    }
    return false;
}

/** implementation of DeltaSurfaceProxy::IShift interface */
bool ProxyVol::sensNameMatches(DeltaSurfaceProxy* shift,
                           const OutputName&  name) const {
    static const string method = "ProxyVol::sensNameMatches";
    try {
        OutputNameArray deltaNames(0);
        Delta delta(Delta::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&delta, proxies[i].get(), deltaNames);
            for (int j=0; j < deltaNames.size(); j++) {
                if (name.equals(deltaNames[j].get())) {
                    return true;
                }
            }
            deltaNames.clear();
        }  
        return false;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void ProxyVol::sensAppendName(DeltaSurfaceProxy* shift,
                          OutputNameArray&   namesList) const{
    static const string method = "ProxyVol::sensAppendName";
    try {
        Delta delta(Delta::DEFAULT_SHIFT);
        for (int i = 0; i < proxies.size(); i++) {
            proxySensNames(&delta, proxies[i].get(), namesList);
        }  
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

bool ProxyVol::sensShift(DeltaSurfaceProxy* shift) {
    DeltaSurfaceSP delta(shift->makeDeltaSurfaceShift());
    OutputNameConstSP name = delta->getMarketDataName();
    OutputNameArray deltaNames(0);

    // Use a new delta in proxySensNames() because toTweak proxies passed from DeltaProxy
    // to delta in makeDeltaShift may interfere with name collection
    DeltaSP dummyDelta(new Delta(Delta::DEFAULT_SHIFT));

    for (int i = 0; i < proxies.size(); i++) {
        AssetSP proxy = proxies[i].getSP();
        proxySensNames(dummyDelta.get(), proxies[i].get(), deltaNames);
        for (int j=0; j < deltaNames.size(); j++) {
            if (shift->marketDataNameMatches(deltaNames[j])) {
                delta->findAndShift(proxy, name);
                shift->setInitialValue(delta->getInitialValue());
            }
        }
        deltaNames.clear();
    }
    return false;
}

static CFieldConstSP proxiesField;

/** used to turn off proxy greeks */
bool ProxyVol::recurse(const CFieldConstSP& field,
                   const CClassConstSP& targetClass) const {

    // this gets called as part of the tweaking and allows us to specify 
    // whether the fields within this class should be tweaked or not. The
    // target class indicates what is being shifted
    if (field == proxiesField) {
        if (targetClass == MuParallel::Shift::TYPE                     ||
            targetClass == MuParallel::RestorableShift::TYPE           ||
            targetClass == MuPointwise::IShift::TYPE                   ||
            targetClass == MuPointwise::IRestorableShift::TYPE         ||
            targetClass == MuSpecial::IShift::TYPE                     ||
            targetClass == MuSpecial::IRestorableShift::TYPE           ||
            targetClass == BorrowParallelShift::Shift::TYPE            ||
            targetClass == RhoBorrowParallel::IShift::TYPE             ||
            targetClass == RhoBorrowParallel::IRestorableShift::TYPE   ||
            targetClass == RhoBorrowPointwise::IShift::TYPE            ||
            targetClass == RhoBorrowPointwise::IRestorableShift::TYPE  ||
            targetClass == ContangoRhoParallel::Shift::TYPE ||
            targetClass == ITweakableWithRespectTo<Spot>::TYPE         ||
            targetClass == DeltaSurface::Shift::TYPE                   ||
            targetClass == ITweakableWithRespectTo<VolParallel>::TYPE                   ||
            targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE   ||
            targetClass == VegaMatrix::IShift::TYPE                    ||
            targetClass == VegaSkewParallel::IShift::TYPE              ||
            targetClass == VegaSkewPointwise::IShift::TYPE             ||

            targetClass == VolParallelShift::IShift::TYPE              ||
#if 0
            // unclear if we want to change these too
            targetClass == VolBenchmarkShift::IShift::TYPE             ||
            targetClass == PowerVega::IShift::TYPE                     ||
            targetClass == VolRelativeShift::IShift::TYPE              ||
            targetClass == VolAbsoluteShift::IShift::TYPE              ||
#endif
            targetClass == RootTimeVega::IShift::TYPE                  ||
            targetClass == ITweakableWithRespectTo<RateParallel>::TYPE ||
            targetClass == ITweakableWithRespectTo<RatePointwise>::TYPE ||
            targetClass == IRestorableWithRespectTo<RatePointwise>::TYPE) {
            return false;
        }
    }

    return true;
}

PDFCalculator* ProxyVol::getPDFCalculator(const PDFRequest* request,
                                          const CAsset*     asset) const {
    static const string method("ProxyVol::getPDFCalculator");
    if (!PDFRequestLNStrike::TYPE->isInstance(request)){
        throw ModelException(method, "Only LN pdf supported");
    }
    if (useVegaMatrix) {
        // this shouldn't happen. Certainly not clear what we should do
        throw ModelException(method, "Internal error - doing "
                             "strikewise vega with local vol");
    }
    const PDFRequestLNStrike& lnRequest = 
        dynamic_cast<const PDFRequestLNStrike&>(*request);
    PDFRequestLNStrikeConstSP lnRequestSP(copyIfRef(&lnRequest));
    double spot = asset->getSpot();
    // for performance we choose some strikes, calculate a vol surface for
    // these. Then spline the vol surface and use the PDFParamLNStrike class
    // to do the work. This avoids umpteen calculations of fwds and vols.
    
    // Step 1 - choose strikes. (done in defaultVolSurface())
    // Step 2 - create a vol surface with these strike
    VolSurfaceSP myVolSurface(defaultVolSurface(asset));
    // Step 3 - create spline
    smartPtr<VolSpline> volSpline(new VolSpline(*myVolSurface, spot));
    // Step 4 - use PDFParamLNStrike class
    // need a CVolParam though
    CVolParamConstSP volParam(volSpline->getVolParam());
    smartConstPtr<VolSpline> volSplineConst(volSpline);
    return (new PDFParamLNStrike(baseDate, CAssetConstSP::attachToRef(asset), 
                                 timeMetric, volSplineConst, 
                                 volParam, lnRequestSP));
}


double ProxyVol::getPDFBoundaryProb() const
{
    return pdfBoundaryProb;
}


    
/** for reflection */
ProxyVol::ProxyVol(): CVolBase(TYPE), minVol(0.0), useVegaMatrix(false), 
    useSkewParallel(false), useSkewPointwise(false), useCorrelationTerm(false), pdfBoundaryProb(CVolProcessedBS::LEGACY_PDF_BOUNDARY_PROB) {}
    
class ProxyVolHelper{
public:
    // Invoked when Class is 'loaded' **
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ProxyVol, clazz);
        SUPERCLASS(CVolBase);
        IMPLEMENTS(IVolTypeSensitiveStrikes);
        IMPLEMENTS(ObjectIteration::IOverride);
        IMPLEMENTS(IPDFCalculator);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IRestorableWithRespectTo<VolParallel>);
        IMPLEMENTS(IRestorableWithRespectTo<VolPointwise>);
        IMPLEMENTS(VegaMatrix::IShift);
        IMPLEMENTS(RootTimeVega::IRestorableShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VolParallelShift::Shift);
        IMPLEMENTS(VolBenchmarkShift::Shift);
        IMPLEMENTS(PowerVega::Shift);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(ITweakableWithRespectTo<RateParallel>);
        IMPLEMENTS(ITweakableWithRespectTo<RatePointwise>);
        IMPLEMENTS(DeltaProxy::Shift);
        IMPLEMENTS(VegaProxyParallel::Shift);
        IMPLEMENTS(VegaProxyPointwise::IShift);
        IMPLEMENTS(VegaProxyMatrix::IShift);
        IMPLEMENTS(VegaSkewProxyParallel::IShift);
        IMPLEMENTS(VegaSkewProxyPointwise::IShift);
        IMPLEMENTS(RootTimeVegaProxy::IShift);
        IMPLEMENTS(DeltaSurfaceProxy::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(VolAbsoluteShift::IShift);
        IMPLEMENTS(IPDFBoundaryProb);

        EMPTY_SHELL_METHOD(defaultProxyVol);
        FIELD(name, "Vol identifier");
        FIELD(yc, "Yield curve for asset's base currency");
        FIELD(proxies, "Assets whose vols are used to build the ProxyVol");
        FIELD(ccyTreatments, "Currency treatments");
        FIELD(weights, "Percentage of each asset");
        FIELD(marketHols, "Market hols for ProxyVol");
        FIELD(timeMetric, "Trading time for composite vol");
        FIELD(spreadDates, "Benchmarks for vol spread");
        FIELD(volSpread, "Vol spread for composite vol");
        FIELD(volSkew, "Vol skew for composite vol");
        FIELD(minVol, "minimum fund vol");
        FIELD(useCorrelationTerm, "useCorrelationTerm");
        FIELD(baseDate, "today");
        FIELD_MAKE_OPTIONAL(baseDate);
        FIELD(correlations, "Asset-Asset Correlations");
        FIELD_MAKE_OPTIONAL(correlations);

        FIELD(adjWeights, "cached weights");
        FIELD_MAKE_TRANSIENT(adjWeights);

        // for tweaking
        FIELD(sensSpread, "for fund vol greeks");
        FIELD_MAKE_TRANSIENT(sensSpread);
        FIELD(useVegaMatrix, "useVegaMatrix");
        FIELD_MAKE_TRANSIENT(useVegaMatrix);
        FIELD(vegaMatrix, "for use with VegaMatrix");
        FIELD_MAKE_TRANSIENT(vegaMatrix);
        FIELD(useSkewParallel, "useSkewParallel");
        FIELD_MAKE_TRANSIENT(useSkewParallel);
        FIELD(skewParallel, "skewParallel");
        FIELD_MAKE_TRANSIENT(skewParallel);
        FIELD(useSkewPointwise, "useSkewPointwise");
        FIELD_MAKE_TRANSIENT(useSkewPointwise);
        FIELD(skewPointwise, "skewPointwise");
        FIELD_MAKE_TRANSIENT(skewPointwise);     

        FIELD(pdfStrikes, "pdfStrikes");
        FIELD_MAKE_TRANSIENT(pdfStrikes);
        FIELD(corrTermArray, "corrTermArray");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrTermArray);
        FIELD(corrObjects, "corrObjects");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrObjects);

        FIELD(pdfBoundaryProb, "limits implied pdf");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(pdfBoundaryProb);

        ClassSetAcceptMethod(ProxyVol::acceptCollector);
        ClassSetAcceptMethod(ProxyVol::acceptNameCollector);
        ClassSetAcceptMethod(ProxyVol::acceptFutureCollector);
        ClassSetAcceptMethod(ProxyVol::acceptValueDateCollector);
        ClassSetAcceptMethod(ProxyVol::acceptImntCcy);
        ClassSetAcceptMethod(ProxyVol::acceptHoliday);

        // look up field for use on recurse
        proxiesField = clazz->getDeclaredField("proxies");
    }

    static IObject* defaultProxyVol(){
        return new ProxyVol();
    }
};

CClassConstSP const ProxyVol::TYPE = CClass::registerClassLoadMethod(
    "ProxyVol", typeid(ProxyVol), ProxyVolHelper::load);

bool ProxyVolLoad() {
    return (ProxyVol::TYPE != NULL);
}

DRLIB_END_NAMESPACE

