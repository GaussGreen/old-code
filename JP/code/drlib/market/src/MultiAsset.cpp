//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiAsset.cpp
//
//   Description :
//
//   Date        : July 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MultiAsset.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/CorrelationCategory.hpp"
#include "edginc/MDFUtil.hpp"
#include "edginc/CDSParSpreads.hpp"


DRLIB_BEGIN_NAMESPACE

/** This is not an asset, but captures multiple factors data and
    implements a choice of views of them. */

void MultiAsset::checkFactorRange(const char* routine, int iFactor) const {
    if (iFactor < 0 || iFactor >= assets.size()) {
        throw ModelException(routine,
                             "Factor index (" +Format::toString(iFactor) +
                             ") out of range [0," +
                             Format::toString(assets.size()) + ")");
    }
}

/** Returns the number of factors in this collection (this does not include
    MarketFactors which are contained within top level MarketFactors) */
int MultiAsset::numFactors() const {
    return assets.size();
}

//// for backwards compatibility
int MultiAsset::NbAssets() const {
    return assets.size();
}

//// Returns the correlations between the immediate factors.
//// This is numFactors() x numFactors()
CDoubleMatrixConstSP MultiAsset::factorsCorrelationMatrix() const {
    static const char* method = "MultiAsset::factorsCorrelationMatrix";
    // create it on the fly
    CDoubleMatrixSP theMatrix(new DoubleMatrix(assets.size(),
                              assets.size()));
    DoubleMatrix& matrix = *theMatrix;
    int pos = 0;
    for (int i = 0; i < assets.size(); i++) {
        matrix[i][i] = 1.0; // set diagonal element
        for (int j = i + 1; j < assets.size(); j++, pos++) {
            //// TO DO:  determine what to do in Sampras case
            //// when strictCorr is false and returning full
            //// correlation matrix with zeros for missing entries
            //// would eat up to much memory.
            //// JD: temporarily return zeros, to avoid an exception.
            if (!strictCorr && !corrObjects[pos]) { 
                matrix[i][j] = 0.0;
                matrix[j][i] = matrix[i][j];
            }
            else if (Correlation::TYPE->isInstance(corrObjects[pos].get())) {
                CorrelationSP p = CorrelationSP::dynamicCast(corrObjects[pos]);
                matrix[i][j] = p->getCorrelation();
                matrix[j][i] = matrix[i][j];
            } else {
                // TO DO: handle the case the correlation is post calibrated
                throw ModelException(method, "correlation is calibrated");
            }
        }
    }
    return theMatrix;
}

/** Returns the name of the specified 'market factor' */
string MultiAsset::getName(int index) const {
    checkFactorRange("MultiAsset::getName", index);
    return assets[index]->getName();
}

/** Returns the MarketFactor with the specifed index */
IMarketFactorConstSP MultiAsset::getFactor(int index) const {
    checkFactorRange("MultiAsset::getFactor", index);
    return assets[index].getSP();
}

/** Pull out the component assets & correlations from the market data */
void MultiAsset::getMarket(const IModel* model, const MarketData* market) {
    static const string method("MultiAsset::getMarket");
    try {
        int numAssets = assets.size();
        // do the easy bit first
        yc.getData(model, market);
        for (int i = 0; i < numAssets; i++) {
            IMarketFactor::getMarketData(model, market, ccyTreatments[i],
                                         yc.getName(), assets[i]);
        }
        market->GetReferenceDate(valueDate);
        /** correlations: if we are explicitly supplied with just the numbers,
            then corr obj were built in validatePop2Object */
        if (!correlations.empty()) {
            // No need to worry about strictCorr flag since memory
            // issues with 0 entries would have surfaced earlier
            // when the correlations matrix was being put together.
            int pos = 0;
            for (int i = 0; i < assets.size(); i++) {
                for (int j = i + 1; j < assets.size(); j++, pos++) {
#if 0
                    corrObjects[pos]->configureForSensitivities(
                        assets[i]->getClass(), assets[j]->getClass());
#else
                    // force it to be EQ-EQ
                    corrObjects[pos]->configureForSensitivities(
                        CAsset::TYPE, CAsset::TYPE);
#endif

                    corrObjects[pos]->getMarket(model, market);
                }
            }
        } else {
            // TO DO:  comment out next line and use push_back
            int pos = 0;
            corrObjects.resize((numAssets * numAssets - numAssets)/2);
            for (int i = 0; i < assets.size(); i++) {
                for (int j = i + 1; j < assets.size(); j++, pos++) {
                    string corrName;
                    bool found =
                        market->hasCorrelationData(
                            assets[i].getName(),
                            assets[j].getName() );
                    if ( !found && strictCorr ) {
                        throw ModelException(method,
                                             "No correlation available with respect to "+
                                             assets[i].getName()+" and "+assets[j].getName());
                    }
                    if ( found ) {
                        corrName = market->getCorrelationName(
                                       assets[i].getName(),
                                       assets[j].getName() );
                        // then get hold of the actual object
                        CorrelationBaseSP corrBase = model->getCorrelation(
                                                         corrName,
#if 0
                                                         assets[i]->getClass(),
                                                         assets[j]->getClass(),
#else
                                                         // force it to be EQ-EQ
                                                         CAsset::TYPE, CAsset::TYPE,
#endif
                                                         CorrelationCommon::TYPE,
                                                         market);
                        corrObjects[pos] = CorrelationCommonSP::dynamicCast(corrBase);
                    } // The following would take up too much memory for sampras (couple of gigs).
                    // Instead, when strictCorr=false, we ignore
                    // missing correlations
                    /*else {
                          corrName = assets[i].getName()+"_"+assets[j].getName();
                          CorrelationSP corr(
                              new Correlation( 
                                  corrName,
                                  assets[i].getName(),
                                  assets[j].getName(),
                                  0 ) );
                          corrObjects[pos] = corr;
                      }
                      */
                }
            }
        }
        // ask the model whether we want to have correlation term structure
        MarketDataFetcherSP fetcher = model->getMDF();
        if (MDFUtil::useCorrTerm(*fetcher)) {
            /** all corr category and corr term objects should be available.
                if not, then we wont fail, but assign zero correlationTerm objects */

            // allocate space
            corrTermArray.resize((numAssets * numAssets - numAssets)/2);
            skipFwdCorrelation.resize((numAssets * numAssets - numAssets)/2);

            int pos=0, i,j;
            for (i = 0; i < assets.size(); i++) {
                if (market->hasData(assets[i].getName(), CorrelationCategory::TYPE)) {
                    MarketObjectSP mo1 =
                        model->GetMarket(market, assets[i].getName(), CorrelationCategory::TYPE);
                    CorrelationCategorySP category1 = CorrelationCategorySP::dynamicCast(mo1);
                    if (!category1) {
                        throw ModelException(method, "Required correlation category object not retrievable");
                    }
                    if (!CString::equalsIgnoreCase(category1->getCategoryName(),
                                                   CorrelationCategory::SKIP_FWD_CORRELATION)) {
                        for (j = i+1; j < assets.size(); j++, pos++) {
                            if (market->hasData(assets[j].getName(), CorrelationCategory::TYPE)) {
                                MarketObjectSP mo2 =
                                    model->GetMarket(market, assets[j].getName(), CorrelationCategory::TYPE);
                                CorrelationCategorySP category2 = CorrelationCategorySP::dynamicCast(mo2);
                                if (!category2) {
                                    throw ModelException(method, "Required correlation category object not retrievable");
                                }

                                if (!CString::equalsIgnoreCase(category2->getCategoryName(),
                                                               CorrelationCategory::SKIP_FWD_CORRELATION)) {
                                    /** both assets have correlation category, which is NOT "NO_FWD_COOR"
                                        => assign proper CorrelationTerm object, which has to be in the market
                                            set skipFwdCorrelation to false */
                                    const string& corrTermName =
                                        market->getCorrelationTermName(category1->getCategoryName(), category2->getCategoryName());
                                    MarketObjectSP moTerm = model->GetMarket(market, corrTermName, CorrelationTerm::TYPE);
                                    corrTermArray[pos] = CorrelationTermSP::dynamicCast(moTerm);
                                    if(!corrTermArray[pos]) {
                                        throw ModelException(method, "Required corr term object not retrieveable");
                                    }
                                    corrTermArray[pos]->getMarket(model, market);
                                    skipFwdCorrelation[pos] = false;
                                } else {
                                    /** asset j has correlation category "NO_FWD_CORR"
                                        => assign zero object & set skipFwdCorrelation to true */
                                    corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                                    skipFwdCorrelation[pos] = true;
                                }
                            } else {
                                /** asset j does not have a correlation cateogory
                                    => assign zero object & set skipFwdCorrelation to false */
                                corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                                skipFwdCorrelation[pos] = false;
                            }
                        }
                    } else {
                        /** asset i has correlation category "NO_FWD_CORR"
                            => assign zero objects & set skipFwdCorrelation to true */
                        for (j = i+1; j < assets.size(); j++, pos++) {
                            corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                            skipFwdCorrelation[pos] = true;
                        }
                    }
                } else {
                    /** asset i does not have a correlation category
                        => assign zero objects & set skipFwdCorrelation to false */
                    for (j = i+1; j < assets.size(); j++, pos++) {
                        // assign zero object
                        corrTermArray[pos] = CorrelationTermSP(new CorrelationTerm(true));
                        skipFwdCorrelation[pos] = false;
                    }
                }
            }
        }
        if (MDFUtil::useCorrSkew(*fetcher)) {
                /** all corr category and corr skew objects should be available.
                    if not, then we wont fail, but assign zero correlationTerm objects */
                
                // allocate space
                localCorrSqueezeIndexArray.resize(numAssets);                
                StringArray regionCollection(1);
                regionCollection[0] = "";
                
                localCorrSqueezeArray.resize(1); // first element is always a dummy
                localCorrSqueezeArray[0] = LocalCorrSqueezeSP(new LocalCorrSqueeze(true));
                
                for (int iAsset = 0; iAsset < assets.size(); iAsset++) {                    
                    if (market->hasData(assets[iAsset].getName(), CorrelationCategory::TYPE)) {
                        MarketObjectSP mo = 
                            model->GetMarket(market, assets[iAsset].getName(), CorrelationCategory::TYPE);
                        CorrelationCategorySP category = CorrelationCategorySP::dynamicCast(mo);
                        string categoryName = category->getCategoryName(); // this is the name of the squeeze object
                        if (categoryName.empty()) {
                            throw ModelException(method, "Entry categoryName for CorrelationCategory "
                                + category->getName() + " missing, but expected.");
                        }
                        int iRegion;
                        for (iRegion=0; iRegion<regionCollection.size(); iRegion++) {
                            if (CString::equalsIgnoreCase(categoryName, regionCollection[iRegion])) {
                                localCorrSqueezeIndexArray[iAsset]=iRegion;
                                break;
                            }
                        }
                        if (iRegion == regionCollection.size()) {
                            localCorrSqueezeIndexArray[iAsset] = iRegion;
                            
                            regionCollection.resize(iRegion+1);
                            regionCollection[iRegion] = categoryName;
                            
                            localCorrSqueezeArray.resize(iRegion+1);

                            MarketObjectSP moSqueeze = 
                                model->GetMarket(market, categoryName, LocalCorrSqueeze::TYPE);
                            localCorrSqueezeArray[iRegion] = LocalCorrSqueezeSP::dynamicCast(moSqueeze); 
                            if (!localCorrSqueezeArray[iRegion]) {
                                throw ModelException(method, "Required local corr squeeze obj not retrievable");
                            }
                            localCorrSqueezeArray[iRegion]->getMarket(model, market); 
                        }
                    } else {
                        localCorrSqueezeIndexArray[iAsset] = 0;
                    }
                }
            }                
    } catch (exception& e) {
        throw ModelException(e, "MultiAsset::getMarket");
    }
}

/** Validation */
void MultiAsset::validatePop2Object() {
    static const string method("MultiAsset::validatePop2Object");
    try {
        int numAssets = assets.size();

        /** correlations: if we are explicitly supplied with just the numbers,
            then generate correlation objects */
        if (!correlations.empty()) {
            // validation of correlations
            if (correlations.numRows() != numAssets ||
                    correlations.numCols() != numAssets) {
                throw ModelException(method,
                                     "Correlations must be an n x n matrix "
                                     "with n = number of assets in MultiAsset");
            }
            correlations.checkSymmetric();
            // First allocate space
            // No need to worry about strictCorr flag since memory
            // issues with 0 entries would have surfaced earlier
            // when the correlations matrix was being put together.
            corrObjects.resize((numAssets * numAssets - numAssets)/2);
            // then generate correlation objects
            int pos = 0, i,j;
            for (i = 0; i < assets.size(); i++) {
                for (j = i + 1; j < assets.size(); j++, pos++) {
                    // the Correlation constructor reorders the names into
                    // alphabetical order if no name supplied. For backwards
                    // compatibility we keep the original order
                    string corrName(assets[i].getName() + "_" +
                                    assets[j].getName());
                    corrObjects[pos] =
                        CorrelationSP(new Correlation(corrName,
                                                      assets[i].getName(),
                                                      assets[j].getName(),
                                                      correlations[i][j]));
                }
            }
        }

        if (ccyTreatments.size() != assets.size()) {
            throw ModelException(method,
                                 "Not enough ccyTreatments supplied");
        }

        // validate everthing else .... TBD
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Maps onto assetCrossValidate for each constituent
void MultiAsset::crossValidate(const DateTime&       startDate,
                               const DateTime&       valueDate,
                               const YieldCurve*     discCcy,
                               const CInstrument*    instrument) const {
    static const string method("MultiAsset::crossValidate");

    // Check also consistent ccy between instrument and MultiAsset itself
    if (discCcy->getCcy() != yc->getCcy()) {
        throw ModelException(method, "Mismatch between ccy of instrument ("+
                             discCcy->getCcy() +
                             ") and MultiAsset (" + yc->getCcy() + ")");
    }

    for( int i = 0; i < assets.size(); i++) {
        // need to work out correct behaviour for non CAssets
        const IMarketFactor* mFactor = assets[i].get();
        if (CAsset::TYPE->isInstance(mFactor)) {
            CAsset* asset = STATIC_CAST(CAsset, mFactor);
            AssetUtil::assetCrossValidate(asset,
                                          startDate>valueDate,/* fwdStarting
                                                                                                               flag */
                                          startDate,
                                          valueDate,
                                          discCcy,
                                          instrument);
        }
    }
}


/** Returns a list of factor indexes indicating which asset is
    sensitive to combination of supplied SensControl and
    OutputName. If crossAssetSensitivities is true then this
    only covers sensitivities such as phi which can live at
    this level otherwise these are excluded. name must not be
    NULL */
IntArray MultiAsset::getSensitiveFactors(
    const IPerNameSensitivity* sens,
    bool                       crossAssetSensitivities,
    const OutputName*          name) const {
    IntArray assetIdxs;
    if (!crossAssetSensitivities) {
        // loop over assets
        for (int i = 0; i < assets.size(); i++) {
            // get list of what this is sensitive to
            OutputNameArrayConstSP names(sens->allNames(assets[i].get()));
            for (int j = 0; j < names->size(); j++) {
                if ((*names)[j]->equals(name)) {
                    assetIdxs.push_back(i);
                    break;
                }
            }
        }
    } else if (sens->shiftInterface()->
               isAssignableFrom(CorrelationTerm::TYPE) && (corrTermArray.size() > 0) ) {
        // Check the CorrelationTerm array
        int pos = 0;
        for (int i = 0; i < assets.size(); i++) {
            for (int j = i + 1; j < assets.size(); j++, pos++) {
                string name1 = corrTermArray[pos]->getName();
                if (name->equals(corrTermArray[pos]->getName())) {
                    assetIdxs.push_back(i);
                    assetIdxs.push_back(j);
                }
            }
        }
    } else if (sens->shiftInterface()->
               isAssignableFrom(LocalCorrSqueeze::TYPE) && (localCorrSqueezeArray.size() > 1) ) {
        int pos = 0;
        for (int i = 0; i < assets.size(); i++) {
            for (int j = i + 1; j < assets.size(); j++, pos++) {
                assetIdxs.push_back(i);
                assetIdxs.push_back(j);
            }
        }
    } else if (sens->shiftInterface()->isAssignableFrom(CorrSwapBasisAdj::TYPE) ||
               sens->shiftInterface()->isAssignableFrom(CorrSwapSamplingAdj::TYPE) ) {
        int pos = 0;
        for (int i = 0; i < assets.size(); i++) {
            for (int j = i + 1; j < assets.size(); j++, pos++) {
                assetIdxs.push_back(i);
                assetIdxs.push_back(j);
            }
        }
    } else if (sens->shiftInterface()->isAssignableFrom(Correlation::TYPE)) {
        // check the correlations (for Phi, FXPhi etc)
        int pos = 0;
        for (int i = 0; i < assets.size(); i++) {
            for (int j = i + 1; j < assets.size(); j++, pos++) {
                // TO DO:  Not sure how strictCorr would apply here
                // The code could crash when strictCorr is false
                // and a correlation was not provided since then
                // the associated corrObjects[pos] would be null .
                if (name->equals(corrObjects[pos]->getName()) &&
                        corrObjects[pos]->isSensitiveTo(sens)) {
                    assetIdxs.push_back(i);
                    assetIdxs.push_back(j);
                }
            }
        }
    }
    // In case we have the same object used twice ...
    // sort then
    sort(assetIdxs.begin(), assetIdxs.end());
    // remove duplicates
    vector<int>::iterator newEnd = unique(assetIdxs.begin(),
                                          assetIdxs.end());
    assetIdxs.resize(newEnd - assetIdxs.begin());
    return assetIdxs;
}

/** Returns a list of asset indexes indicating which asset is
    sensitive to supplied SensControl. If crossAssetSensitivities
    is true then this only covers sensitivities such as phi which
    can live at this level otherwise these are excluded. The
    relevant name of what is being tweaked must be stored in the
    SensControl */
IntArray MultiAsset::getSensitiveFactors(
    const Sensitivity* sens,
    bool               crossAssetSensitivities) const {
    return getSensitiveAssets(sens, crossAssetSensitivities);
}

/** Returns a list of asset indexes indicating which asset is
    sensitive to supplied SensControl. If crossAssetSensitivities
    is true then this only covers sensitivities such as phi which
    can live at this level otherwise these are excluded. The
    relevant name of what is being tweaked must be stored in the
    SensControl */
IntArray MultiAsset::getSensitiveAssets(
    const Sensitivity* sens,
    bool               crossAssetSensitivities) const {
    // are we some sort of combination of tweaks
    const ITwoSidedDeriv* twoSidedDeriv =
        dynamic_cast<const ITwoSidedDeriv*>(sens);
    if (!twoSidedDeriv) {
        // is it a sens control per name?
        const IPerNameSensitivity* sensCtrl =
            dynamic_cast<const IPerNameSensitivity*>(sens);
        // and look at the name
        if (!sensCtrl || !sensCtrl->getMarketDataName()) {
            /* if it's not a SensControl we're a bit stuck
               (but this should be typically covered by
               ITwoSidedDeriv since most 'composite' greeks are
               done using another control) */
            // if no name - probably a theta type of
            // shift. Return everything
            IntArray sensAssets(numFactors());
            for (int i = 0; i < sensAssets.size(); i++) {
                sensAssets[i] = i;
            }
            return sensAssets;
        }
        // use low level function
        return getSensitiveFactors(sensCtrl, crossAssetSensitivities,
                                   sensCtrl->getMarketDataName().get());
    } else {
        const ScalarShiftArray& shifts =
            twoSidedDeriv->getComponentShifts();
        // loop through shifts merging array together
        IntArray sensAssets;
        for (int i = 0; i < shifts.size(); i++) {
            OutputNameConstSP cmptName(shifts[i]->getMarketDataName());
            // get for this component
            IntArray moreAssets(
                getSensitiveFactors(shifts[i].get(),
                                    crossAssetSensitivities,
                                    cmptName.get()));
            // merge them in
            sensAssets.insert(sensAssets.end(),
                              moreAssets.begin(), moreAssets.end());
            // sort then
            sort(sensAssets.begin(), sensAssets.end());
            // remove duplicates
            vector<int>::iterator newEnd =
                unique(sensAssets.begin(), sensAssets.end());
            sensAssets.resize(newEnd - sensAssets.begin());
        }
        return sensAssets;
    }
}

/** Shifts the object using given shift (see Theta::Shift)*/
bool MultiAsset::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    } catch (exception& e) {
        throw ModelException(e, "Fund::sensShift (theta)");
    }
    return true; // tweak stuff inside
}


/** record forwards at maturity for each factor */
void MultiAsset::recordFwdAtMat(OutputRequest*  request,
                                Results*        results,
                                const DateTime& maturityDate) const {
    for (int i = 0; i < assets.size() ; ++i) {
        if (IGeneralAsset::TYPE->isInstance(assets[i].get())) {
            // unclear exactly what methods we want on MarketFactor
            // at the moment
            IGeneralAsset* genAsset = STATIC_CAST(IGeneralAsset,
                                                  assets[i].get());
            genAsset->recordFwdAtMat(request, results, maturityDate);
        }
    }
}

// returns ccy treatment
string MultiAsset::getCcyTreatment(int iFactor) const {
    return ccyTreatments[iFactor];
}

MultiAsset::MultiAsset(YieldCurveWrapper yc,
                       IMarketFactorWrapperArray assets,
                       CStringArray ccyTreatments)
        : CObject(TYPE), yc(yc), assets(assets), ccyTreatments(ccyTreatments),
        useCorrelationTermNotUsed(false), strictCorr(true) {}


/* for reflection */
MultiAsset::MultiAsset() : CObject(TYPE), assets(0), ccyTreatments(0),
        useCorrelationTermNotUsed(false), strictCorr(true) {}

/** Invoked when Class is 'loaded' */
void MultiAsset::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MultiAsset, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMultiMarketFactors);
    EMPTY_SHELL_METHOD(defaultMultiAsset);
    IMPLEMENTS(Theta::Shift);
    FIELD(yc, "MultiAsset's currency");
    FIELD(assets, "Assets in MultiAsset");
    FIELD(ccyTreatments, "Currency treatments");
    FIELD(correlations, "Factor-Factor Correlations");
    FIELD_MAKE_OPTIONAL(correlations);
    FIELD(corrTermArray, "corrTermArray");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrTermArray);
    FIELD(localCorrSqueezeArray, "localCorrSqueezeArray");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(localCorrSqueezeArray);        
    FIELD(localCorrSqueezeIndexArray, "localCorrSqueezeIndexArray");
    FIELD_MAKE_TRANSIENT(localCorrSqueezeIndexArray);
    FIELD(corrObjects, "Correlations");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrObjects);
    FIELD(valueDate, "value date");
    FIELD_MAKE_TRANSIENT(valueDate);
    FIELD(skipFwdCorrelation, "for some underlyings, we only use spot correlation");
    FIELD_MAKE_TRANSIENT(skipFwdCorrelation);
    FIELD_USING_ALIAS(useCorrelationTermNotUsed, useCorrelationTerm, "DO NOT USE");
    FIELD_MAKE_OPTIONAL(useCorrelationTerm);
    FIELD_USING_ALIAS(shortTermSpreadNotUsed, shortTermSpread, "DO NOT USE");
    FIELD_MAKE_OPTIONAL(shortTermSpread);
    FIELD_USING_ALIAS(longTermSpreadNotUsed, longTermSpread, "DO NOT USE");
    FIELD_MAKE_OPTIONAL(longTermSpread);
    FIELD_USING_ALIAS(corrTermArrayNotUsed, corrTermArray, "DO NOT USE");
    FIELD_MAKE_OPTIONAL(corrTermArray);
    FIELD(strictCorr, "FALSE means to ignore missing correlations.  "
          "TRUE would give error.  Default value is true");
    FIELD_MAKE_OPTIONAL(strictCorr);
}

IObject* MultiAsset::defaultMultiAsset() {
    return new MultiAsset();
}









/** View of MultiAsset as IMultiFactors */
MultiAsset::AsMultiFactors::AsMultiFactors(const CAssetArray& assets,
        MultiAssetConstSP  multiAsset):
CObject(TYPE), assets(assets), multiAsset(multiAsset) {}

// MultiFactors interface. Here this is a degenerate case with
// NbFactors==NbAssets
int MultiAsset::AsMultiFactors::NbAssets() const {
    return assets.size();
}

int MultiAsset::AsMultiFactors::numFactors() const {
    return assets.size();
}

/** Returns the name of the specified 'market factor' */
string MultiAsset::AsMultiFactors::getName(int index) const {
    if (!assets[index])
        return multiAsset->getName(index); // just do this?
    else
        return assets[index]->getName();
}

/** Returns the MarketFactor with the specifed index */
IMarketFactorConstSP MultiAsset::AsMultiFactors::getFactor(int index) const {
    if (!assets[index])
        return multiAsset->getFactor(index); // just do this?
    else
        return assets[index];
}

/** Pull out the component assets & correlations from the market data */
void MultiAsset::AsMultiFactors::getMarket(const IModel* model, const MarketData* market) {
    // assume all done already
}


/** Validate that the market data and instrument are consistent with
    each other etc */
void MultiAsset::AsMultiFactors::crossValidate(const DateTime&       startDate,
        const DateTime&       valueDate,
        const YieldCurve*     discCcy,
        const CInstrument*    instrument) const {
    multiAsset->crossValidate(startDate, valueDate, discCcy, instrument);
}

void MultiAsset::AsMultiFactors::checkAssetRange(const char* routine, int iAsset) const {
    if (iAsset<0 || iAsset>=NbAssets()) {
        throw ModelException(routine,
                             "Asset index (" +Format::toString(iAsset) +
                             ") out of range [0," +
                             Format::toString(NbAssets()) + ")");
    }
}

string MultiAsset::AsMultiFactors::assetGetName(int iAsset) const {
    checkAssetRange("MultiAsset::assetGetName", iAsset);

    if (!assets[iAsset])
        return multiAsset->getName(iAsset); // just do this?
    else
        return assets[iAsset]->getName();
}
string MultiAsset::AsMultiFactors::assetGetTrueName(int iAsset) const {
    checkAssetRange("MultiAsset::assetGetTrueName", iAsset);

    if (!assets[iAsset])
        return multiAsset->getName(iAsset); // just do this? , no getTrueName()
    else
        return assets[iAsset]->getTrueName();
}
string MultiAsset::AsMultiFactors::assetGetCcyTreatment(int iAsset) const {
    checkAssetRange("MultiAsset::assetGetCcyTreatment", iAsset);

    if (!assets[iAsset])
        return multiAsset->getCcyTreatment(iAsset); // just do this?
    else
        return assets[iAsset]->getCcyTreatment();
}
double MultiAsset::AsMultiFactors::assetGetSpot(int iAsset) const {
    checkAssetRange("MultiAsset::assetGetSpot", iAsset);

    if (!assets[iAsset]) {
        throw ModelException("assetGetSpot() is invliad call on " + multiAsset->getName(iAsset));
        return 0.0;
    } else
        return assets[iAsset]->getSpot();
}
double MultiAsset::AsMultiFactors::assetFwdValue(int             iAsset,
        const DateTime& date) const {
    checkAssetRange("MultiAsset::assetFwdValue", iAsset);

    if (!assets[iAsset]) {
        throw ModelException("assetFwdValue() is invliad call on " + multiAsset->getName(iAsset));
        return 0.0;
    } else
        return assets[iAsset]->fwdValue(date);
}
void MultiAsset::AsMultiFactors::assetFwdValue(int                  iAsset,
        const DateTimeArray& dates,
        CDoubleArray&        result) const {
    checkAssetRange("MultiAsset::assetFwdValue", iAsset);

    if (!assets[iAsset])
        throw ModelException("assetFwdValue() is invliad call on " + multiAsset->getName(iAsset));
    else
        assets[iAsset]->fwdValue(dates, result);
}

int MultiAsset::AsMultiFactors::NbFactors() const {
    return assets.size(); // for now at least
}
void MultiAsset::AsMultiFactors::checkFactorRange(const char* routine, int iFactor) const {
    if (iFactor<0 || iFactor>=NbFactors()) {
        throw ModelException(routine,
                             "Factor index ("+Format::toString(iFactor)+
                             ") out of range [0," +
                             Format::toString(NbFactors()) + ")");
    }
}
string MultiAsset::AsMultiFactors::factorGetName(int iFactor) const {
    checkFactorRange("MultiAsset::factorGetName", iFactor);

    if (!assets[iFactor])
        return multiAsset->getName(iFactor); // just do this?
    else
        return assets[iFactor]->getName();
}
CVolProcessed * MultiAsset::AsMultiFactors::factorGetProcessedVol(
    int                iFactor,
    const CVolRequest* volRequest) const {
    checkFactorRange("MultiAsset::factorGetProcessedVol", iFactor);

    if (!assets[iFactor]) {
        throw ModelException("factorGetProcessedVol() is invliad call on " + multiAsset->getName(iFactor));
        return 0;
        //return multiAsset->factorGetProcessedVol(iFactor, volRequest); // just do this?
    } else
        return assets[iFactor]->getProcessedVol(volRequest);
}

double MultiAsset::AsMultiFactors::factorGetSpot(int iFactor) const {
    checkFactorRange("MultiAsset::factorGetSpot", iFactor);

    if (!assets[iFactor]) {
        throw ModelException("factorGetSpot() is invliad call on " + multiAsset->getName(iFactor));
        return 0.0;
    } else
        return assets[iFactor]->getSpot();
}
double MultiAsset::AsMultiFactors::factorFwdValue(int             iFactor,
        const DateTime& date) const {
    checkFactorRange("MultiAsset::factorFwdValue", iFactor);

    if (!assets[iFactor]) {
        throw ModelException("factorFwdValue() is invliad call on " + multiAsset->getName(iFactor));
        return 0.0;
    } else
        return assets[iFactor]->fwdValue(date);
}

void MultiAsset::AsMultiFactors::factorFwdValues(int                  iFactor,
        const DateTimeArray& dates,
        CDoubleArray&        result) const {
    checkFactorRange("MultiAsset::factorFwdValue", iFactor);

    if (!assets[iFactor])
        throw ModelException("factorFwdValue() is invliad call on " + multiAsset->getName(iFactor));
    else
        assets[iFactor]->fwdValue(dates, result);
}
CDoubleMatrixConstSP MultiAsset::AsMultiFactors::factorsCorrelationMatrix() const {
    return multiAsset->factorsCorrelationMatrix();
}

BoolArray MultiAsset::AsMultiFactors::getSkipFwdCorrArray() const {
    return multiAsset->skipFwdCorrelation;
}
// CorrelationArray getCorrObjArray() const {
CorrelationCommonArray MultiAsset::AsMultiFactors::getCorrObjArray() const {
    return multiAsset->corrObjects;
}
CorrelationTermArray MultiAsset::AsMultiFactors::getCorrTermObjArray() const {
    return multiAsset->corrTermArray;
}
TimeMetricArray MultiAsset::AsMultiFactors::getTimeMetricArray() const {
    TimeMetricArray timeMetricArray(assets.size());
    for(int iAsset=0; iAsset<assets.size(); iAsset++ ) {
        ATMVolRequestSP volRequest(new ATMVolRequest());
        CVolProcessedBSSP volBS(assets[iAsset]->getProcessedVol(volRequest.get()));
        timeMetricArray[iAsset] = TimeMetricSP::constCast(volBS->GetTimeMetric());
    }
    return timeMetricArray;
}

/** record forwards at maturity for each factor */
void MultiAsset::AsMultiFactors::recordFwdAtMat(OutputRequest*  request,
        CResults*       results,
        const DateTime& maturityDate) const {
    multiAsset->recordFwdAtMat(request, results, maturityDate);
}

void MultiAsset::AsMultiFactors::CollapseFactors(
    int             len,                 // length of each array
    const double**  factorLevels,        // [NbFactors] arrays
    double**        assetLevels) const {
    for(int iAsset=0; iAsset< assets.size(); iAsset++) {
        for(int iStep=0; iStep<len; iStep++) {
            assetLevels[iAsset][iStep] = factorLevels[iAsset][iStep];
        }
    }
}

/** Using supplied outputName to select the relevent asset, adds
    sensitive strikes to the sensitiveStrikes parameter */
void MultiAsset::AsMultiFactors::assetGetSensitiveStrikes(
    int                               iAsset,
    const CVolRequest*                volRequest,
    const OutputNameConstSP&          outputName,
    const SensitiveStrikeDescriptor&  sensStrikeDesc,
    const DoubleArraySP&              sensitiveStrikes) const {
    checkAssetRange("MultiAsset::assetGetSensitiveStrikes", iAsset);
    assets[iAsset]->getSensitiveStrikes(volRequest, outputName,
                                        sensStrikeDesc,
                                        sensitiveStrikes);
}

/** Returns a reference to the i-th Asset */
const CAsset& MultiAsset::AsMultiFactors::getAsset(int iAsset) const {
    checkAssetRange("MultiAsset::getAsset", iAsset);
    return *assets[iAsset].get();
}

IntArray MultiAsset::AsMultiFactors::getSensitiveFactors(
    const Sensitivity* sens,
    bool               crossAssetSensitivities) const {
    return multiAsset->getSensitiveFactors(sens, crossAssetSensitivities);
}

/** Returns a list of asset indexes indicating which asset is
    sensitive to supplied SensControl. If crossAssetSensitivities
    is true then this only covers sensitivities such as phi which
    can live at this level otherwise these are excluded. The
    relevant name of what is being tweaked must be stored in the
    SensControl */
IntArray MultiAsset::AsMultiFactors::getSensitiveAssets(
    const Sensitivity* sens,
    bool               crossAssetSensitivities) const {
    return multiAsset->getSensitiveFactors(sens, crossAssetSensitivities);
}

/** Returns PDF calculator for iAsset */
PDFCalculator* MultiAsset::AsMultiFactors::assetPdfCalculator(const PDFRequest* request,
        int               iAsset) const {
    checkAssetRange("MultiAsset::assetPdfCalculator", iAsset);
    return assets[iAsset]->pdfCalculator(request);
}

/** Invoked when Class is 'loaded' */
void MultiAsset::AsMultiFactors::load(CClassSP& clazz) {
    // not visible to EAS/spreadsheet
    REGISTER(AsMultiFactors, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IMultiFactors);
    FIELD(assets, "Factors as Assets");
    FIELD(multiAsset, "Underlying MultiAsset");
}



IMultiFactors* MultiAsset::asMultiFactors() const {
    CAssetArray theAssets(numFactors());
    for (int i = 0; i < theAssets.size(); i++) {
        IMarketFactorConstSP mFactor = assets[i].getSP();
        if (YieldCurve::TYPE->isInstance(mFactor) ||
            CDSParSpreads::TYPE->isInstance(mFactor)) {
            theAssets[i] = CAssetSP(   ); 
            // empty pointer for YieldCurve and CDSParSpreads types
        } else if (!CAsset::TYPE->isInstance(mFactor)) {
            throw ModelException("MultiAsset::asMultiFactors", "Cannot fulfill"
                                 " request for IMultiFactors view of object "
                                 "since asset "+mFactor->getName()+" is of "
                                 "type "+mFactor->getClass()->getName());
        } else
            theAssets[i] = CAssetSP::dynamicCast(
                               IMarketFactorSP::constCast(mFactor));
    }
    return new AsMultiFactors(theAssets, MultiAssetConstSP(this));
}


CClassConstSP const MultiAsset::AsMultiFactors::TYPE =
    CClass::registerClassLoadMethod(
        "MultiAsset::AsMultiFactors", typeid(MultiAsset::AsMultiFactors), load);

CClassConstSP const MultiAsset::TYPE = CClass::registerClassLoadMethod(
                                           "MultiAsset", typeid(MultiAsset), load);


bool MultiAssetLinkIn() {
    return (MultiAsset::TYPE != NULL);
}


DRLIB_END_NAMESPACE
