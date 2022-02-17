//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : MarketDataFetcherCIS.cpp
//
//   Description : Helper class for models that use CDS Basket.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MarketDataFetcherCIS.hpp"
#include "edginc/CreditIndexBasis.hpp"
#include "edginc/CreditIndex.hpp"
#include "edginc/AdjustedCDSParSpreads.hpp"
#include "edginc/AdjustedCDSPSwithTweaking.hpp"

DRLIB_BEGIN_NAMESPACE

MarketDataFetcherCIS::~MarketDataFetcherCIS()
{
    initialise();
}


/** default constructor */
MarketDataFetcherCIS::MarketDataFetcherCIS(): 
    MarketDataFetcherCDS(), 
    useIndexBasis(true),
    calculateIndexBasis(true){
    initialise();
}


/** Will fetch stochastic yield curves plus type of smile and model
    specified. Will probably need some sort of vol choice once we know
    what we're doing */
MarketDataFetcherCIS::MarketDataFetcherCIS(bool getSwaptionVols): 
    MarketDataFetcherCDS(getSwaptionVols), 
    useIndexBasis(true),
    calculateIndexBasis(true){
    initialise();
}


/** Will fetch stochastic yield curves plus type of smile and model
    specified. Also accept an extra parameter indicating whether to 
    fetch the index basis or not. */
MarketDataFetcherCIS::MarketDataFetcherCIS(bool getSwaptionVols,
                                           bool calcIndexBasis,
                                           bool useIndexBasis):
    MarketDataFetcherCDS(getSwaptionVols),
    useIndexBasis(useIndexBasis){
    initialise();
    if (useIndexBasis) {
        calculateIndexBasis = calcIndexBasis;
    }
    else {
        // Since we do not want to use the index basis, set calculateIndexBasis 
        // so that we do not use the index basis provided in the market data
        // (we will set a meaningfull "dummy" instead)
        calculateIndexBasis = true;
    }
}


/** Will fetch stochastic yield curves plus type of smile and model
    specified. Will probably need some sort of vol choice once we know
    what we're doing */
MarketDataFetcherCIS::MarketDataFetcherCIS(
          bool                   getSwaptionVols,
          ICreditIndexMapConstSP indexMap,
          const string&          indexMapTreatment,
          bool                   calculateIndexBasis) :
    MarketDataFetcherCDS(getSwaptionVols),
    useIndexBasis(true),
    calculateIndexBasis(calculateIndexBasis),
    indexMap(indexMap),
    indexMapTreatment(indexMapTreatment){
    initialise();
}

/** configure MarketDataFetcher properly */
void MarketDataFetcherCIS::initialise(){
    // There are a number of combinations of par spread curves and
    // their sub-classes that require care to be taken when determining
    // exactly what we can apply index basis to.
    // This logic will need to be kept up-to-date as new forms of par
    // spread curves are implemented.
    // With the hiearchy at the time of writing it is clear that
    // -- we do not want to adjust a par spread curve that forms
    //    part of a CreditIndex description
    //    insideCreditIndex will be > 0 in this case
    // -- we do not want to adjust a par spread curve contained within
    //    an adjusted curve,but we do want to adjust the adjusted
    //    curve itself
    //    insideAdjustedCDS will be > 0 in this case
    // -- conversely, we do not want to adjust a Quanto, but we do want
    //    to adjust the par spread within
    //    insideQuantoCDS will be > 0 in this case
    //    and since a quanto is also a cds, insideCDSCurve is also
    //    incremented, in which case, if the counters are the same
    //    we know that we are dealing with the par spread curve just
    //    inside the Quanto, and no deeper

    /* Put in a rule to say get AdjustedCDSPSwithTweaking::TYPE if
       asked for a ICDSParSpreads. */
    setRetrievalMode(ICDSParSpreads::TYPE, true, 
                     AdjustedCDSPSwithTweaking::TYPE);
    /* Then put in rules to say get type
       as asked for when inside a CreditIndex, an AdjustedCDSParSpreads,
       an AdjustedCDSPSwithTweaking. For quanto, we never fetch a
       quanto curve so we should never be here with a quanto so we
       can never adjust it. */
    setRetrievalMode(ICDSParSpreads::TYPE, CreditIndex::TYPE, true, NULL);
    setRetrievalMode(ICDSParSpreads::TYPE, AdjustedCDSParSpreads::TYPE,
                     true, NULL);
    setRetrievalMode(ICDSParSpreads::TYPE, AdjustedCDSPSwithTweaking::TYPE,
                     true, NULL);
}

/** When we get an ICDSParSpreads object, see if we should specialise the
    type - normally we would do this within the MarketDataFetcher
    framework but here we have hijacked it so we know when to get the
    adjusted curve */
CClassConstSP MarketDataFetcherCIS::getCDSPSType(
    CClassConstSP requestedType) const{
    /* One potential issue with the approach in initialise() is that
       when we actually ask for the [Par Spread] curve we lose any
       ability to ask for anything other than a ICDSParSpreads.  Hence
       this method was put in as a marker as where it could be done if
       needed. */ 
    return requestedType;
}

/** Pulls out index basis objects from the cache if required */
MarketObjectSP MarketDataFetcherCIS::fetch(const MarketData*    market,
                                           const string&        name,
                                           const CClassConstSP& type,
                                           const IModel*        model) const 
{
    static const string method = "MarketDataFetcherCIS::fetch";
   
    try {
        if (CreditIndexBasis::TYPE->isAssignableFrom(type)) {
            if (!useIndexBasis) {
                 // return a "dummy" indexBasis
                return MarketObjectSP(new CreditIndexBasis(false));
            }
            else if (calculateIndexBasis) {
                // This is a CreditIndexBasis and the calculateIndexBasis flag
                // indicates that we should not fetch it -> return null
                return MarketObjectSP();
            }
        }
        // else
        // find out what type we should actually get
        CClassConstSP typeWanted = 
            stochasticYCFix(market, name, getTypeToRetrieve(type));
        if (typeWanted == AdjustedCDSPSwithTweaking::TYPE){
            // intercept this and do the actual adjustment. Start by retrieving
            // the ICDSParSpreads curve - see if we should alter type
            typeWanted = getCDSPSType(type);
            MarketObjectSP mo(market->GetObjectData(name, typeWanted, model));
            // then create the adjusted curve
            return applyIndexBasis(mo);
        } 
        return !typeWanted? 
            MarketObjectSP():
            market->GetObjectData(name, typeWanted, model);
    } catch (exception& e) {
        throw ModelException(e, 
                             "MarketDataFetcherCIS::fetch", 
                             "Failed to get market data " + name + 
                             " of type " + type->getName());
    }
}


/** ensures correct data is retrieved for quanto curves. Models should redirect
    calls to getComponentMarketData to this method */
void MarketDataFetcherCIS::getComponentMarketData(const IModel*     model,
                                                  const MarketData* market,
                                                  MarketObjectSP    mo) const 
{
    if (CreditIndexBasis::TYPE->isInstance(mo) && calculateIndexBasis){
        const CreditIndexBasis* cib = STATIC_CAST(CreditIndexBasis, mo.get());
        if (cib->used()){
            // This is a CreditIndexBasis and the
            // calculateIndexBasis flag indicates that we should
            // not fetch it, therefore we should not have tried to
            // get its components - unless the
            // haveSetDummyIndexBasis flag is set, meaning that
            // what we have fetched is just a "dummy" because the
            // index basis is not used at all.  Although being in
            // this method it is not obvious, a detailed
            // description of the actual issue is included in the
            // exception raised
            throw ModelException(
                "MarketDataFetcherCIS::getComponentMarketData",
                "Conflicting information: The model is configured "
                "either not to use the index basis at all, or to "
                "compute it internally, but a pre-computed basis "
                "is actually embedded in the indexBasisWrapper.");
        }
    }
    // Just pass to parent
    MarketDataFetcherCDS::getComponentMarketData(model, market, mo);
}


// Locally apply index basis if applicable
MarketObjectSP MarketDataFetcherCIS::applyIndexBasis(
    const MarketObjectSP mo) const
{
    static const string method = "MarketDataFetcherCIS::applyIndexBasis";

    try {
        //in the presence of a CreditIndexMap
        //replace any names in the portfolio with index basis adjusted names
        MarketObjectSP moAdj;

        if (!indexMap) {            
            moAdj = mo; //just return the object as is
        }
        else {
            //actual behaviour is dependent upon the indexMapTreatment input
            if (indexMapTreatment == CreditIndexBasis::APPLY_TO_NONE) {
                moAdj = mo; //do not modify
            }
            else if (indexMapTreatment == 
                     CreditIndexBasis::APPLY_TO_MAPPED_NAMES_ONLY)
            {
                //get the credit index associated to this name from the map
                CreditIndexBaseConstSP index(indexMap->getIndex(mo->getName()));

                if (!index) {
                    moAdj = mo; //do not modify
                }
                else {
                    moAdj = adjustCurve(mo, index); //reset the object to return
                }
            }
            else if (indexMapTreatment == CreditIndexBasis::APPLY_TO_ALL_NAMES){
                //get the credit index irrespective of this name
                CreditIndexBaseConstSP index = indexMap->getIndex();

                if (!index) {
                    //map supports for than one index at once, so cant apply
                    throw ModelException(
                        method, 
                        "Supplied index map is not appropriate "
                        "for adjusting all names. The map must "
                        "define a unique index.");
                }
                else {
                    //reset the object to return
                    moAdj = adjustCurve(mo, index);
                }
            }
            else {
                throw ModelException(method, 
                                     "Unknown form of indexMapTreatment");
            }
        }

        return moAdj;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


// Factor out the actual application of index basis
MarketObjectSP MarketDataFetcherCIS::adjustCurve(
    const MarketObjectSP         mo,
    const CreditIndexBaseConstSP index) const
{
    static const string method = "MarketDataFetcherCIS::adjustCurve";

    //specialise the market object
    BootstrappableCDSParSpreadsSP toAdj = 
        BootstrappableCDSParSpreadsSP(
            dynamic_cast<BootstrappableCDSParSpreads*>(mo.get()));

    if (!toAdj) {
        throw ModelException(method,
                             "Failed to create bootstrappable curve.");
    }

    //create the adjusted name
    AdjustedCDSParSpreadsSP spdAdj = index->adjustCurveNtwk(toAdj);

    return spdAdj;
}


DRLIB_END_NAMESPACE
