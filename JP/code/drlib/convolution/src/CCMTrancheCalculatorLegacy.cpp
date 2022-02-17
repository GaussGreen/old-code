//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CCMTrancheCalculatorLegacy.hpp"
#include "edginc/CCMLossCalculatorBase.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CCMLossUnit.hpp"
#include "edginc/CCMTrancheUtils.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/CcmOnlyParameters.hpp"
#include "edginc/RflOnlyParameters.hpp"
#include "edginc/PortfolioName.hpp" // for the betaTweak function
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
// avoid source file for single method
ITrancheLossCalculatorLegacy::~ITrancheLossCalculatorLegacy(){}
ITrancheLossCalculatorLegacy::ITrancheLossCalculatorLegacy(){}


/** Hash function (needed by hash_map template) */
struct HashUtil{
    size_t operator()(CCMTrancheCalculatorLegacySP c) const {
        return c->hashCode();;
    }
};

/** Equal function (needed by hash_map template) */
struct EqualUtil {
    bool operator()(CCMTrancheCalculatorLegacySP c1, CCMTrancheCalculatorLegacySP c2) const {
        return c1->equal(c2.get());
    }
};

/** The hash table used for caching (CDOKey, pair<DoubleArraySP, DoubleArraySP>) pairs */
typedef hash_map<CCMTrancheCalculatorLegacySP, pair<DoubleArraySP, DoubleArraySP>, HashUtil, EqualUtil> CDOCacheHashtable;

/** */
class CDOCache {
public:
    /** */
    static const pair<DoubleArraySP, DoubleArraySP> get(
                const CCMTrancheCalculatorLegacySP key,
                const Control*               control,
                const CreditMetricsModel*    model,
                const int                    timepoint,
                const int                    cacheSet) {
        
        // Calls "find(key)" only once to improve performance
        CDOCacheHashtable::iterator iter = hashtable.find(key);
        if (iter == hashtable.end()) {
            // Key not found in the hash table : need to compute the
            // densities and insert the (CDOKey, pair<DoubleArraySP, DoubleArraySP>)
            // pair in the hash table

            // First test if the cache is not full
            if (hashtable.size() >= MAX_SIZE) {
                // Cache is full, so we simply empty it !
                // (might do something more clever here, but not sure...)
                hashtable.clear();
            }
            
            // Here is the (potentially) time consuming operation !
            key->convolutionNoDensityCache(control,model,timepoint,cacheSet);
            pair<DoubleArraySP, DoubleArraySP> densities = key->getDensities();

            // Insert the new (key, data) pair into hashtable
            hashtable[key] = densities;
/* this functionality is bypassed for the time being
            if (model->allowQuickGreeks())
            {
                if (control->isPricing())
                {
                    hashtable[key] = densities;
                }
            }
            else
            {
                hashtable[key] = densities;
            }
*/
            return densities;
        } else {
            return iter->second;
        }
    }
private:

    // The cache itself
    static CDOCacheHashtable hashtable;
    
    // Maximum number of entries in the cache
    static const unsigned int MAX_SIZE;
};

/** The cache itself */
CDOCacheHashtable CDOCache::hashtable(MAX_SIZE);

/** Maximum number of entries in the cache */
const unsigned int CDOCache::MAX_SIZE = 100;

/** Builds calculator skipping all variant data ie data that changes
    across time */
CCMTrancheCalculatorLegacy::CCMTrancheCalculatorLegacy(
        bool                           doCreditMetrics, /* (I) */
        bool                           doRFL,           /* (I) */
        const bool                     recoverNotional, /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        double                         lossUnit,
        CounterPartyCreditConstSP      counterParty): 
    CObject(TYPE), 
    lossUnit(lossUnit)
{
    try {
        pastLoss = tranche->portfolioLoss();

        if (doRFL) {
            CCMTrancheCalculatorLegacy::createBasketInfoRFL(tranche, basketInfo);
            cpty = CCMTrancheCalculatorLegacy::createCptyInfoRFL(counterParty);
        }
        else {
            CCMTrancheCalculatorLegacy::createBasketInfo(tranche, basketInfo);
            cpty = CCMTrancheCalculatorLegacy::createCptyInfo(counterParty);
        }

        CCMTrancheCalculatorLegacy::populateBasketInfoCM(tranche, 
                                                         recoverNotional, 
                                                         basketInfo);
        CCMTrancheCalculatorLegacy::populateCptyInfoCM(counterParty, cpty);


        if (!doCreditMetrics) {
            // do additional stuff for CCM
            populateBasketInfoCCM(tranche, basketInfo);
            populateCptyInfoCCM(counterParty, *cpty);
        }
        if (doRFL) {
            // do additional stuff for RFL
            populateBasketInfoRFL(tranche, basketInfo);
            populateCptyInfoRFL(counterParty, *cpty);
        }

        int n = basketInfo.size(); // for ease
        
        /*  renormalize notionals and allocate density */
        ntl.resize(n);
        for (int i = 0; i < n; ++i) {
            basketInfo[i]->ntl /= lossUnit;
            ntl[i] = basketInfo[i]->ntl;
        }
        
        CCMLossUnit::calcLossDensityBounds(ntl, maxl, maxs);
        
        density.resize(maxs+maxl+1);
        densityCond.resize(maxs+maxl+1);
    } 
    catch (exception& e){
        throw ModelException(e, "CCMTrancheCalculatorLegacy::"
                             "CCMTrancheCalculatorLegacy");
    }
}



/** Need to move this class to its own file. Probably split into a
    CreditMetricsTrancheCalculator which CCMTrancheCalculatorLegacy would
    derive from. This setup is the same as the setup in CDO.cpp but
    takes a ConvolutionProduct (the goal is to scrap the one in
    CDO.cpp and break the dependency between it and CDO) */
void CCMTrancheCalculatorLegacy::setupCM(
    const DoubleArray&        survivalProb, /* name survival proba */
    double                    counterPartyProb, /* name survival proba */
    const DoubleArray&        betaOverride) // optional per name override
{
    bool useOverride = !betaOverride.empty();
    for (int i = 0; i < basketInfo.size(); i++){
        basketInfo[i]->survival = survivalProb[i];
        if (useOverride){
            basketInfo[i]->beta = betaOverride[i];
        }
    }
    if (cpty.get()){
        cpty->survival = counterPartyProb;
    }
}

/** Same as setupCM but for CCM rather than Credit Metrics */
void CCMTrancheCalculatorLegacy::setupCCM(
    const DoubleArray&        survivalProb, /* name survival proba */
    const DoubleArray&        floorSurvivalProb, /* from 'senior' curve */
    double                    counterPartyProb, /* name survival proba */
    const DoubleArray&        betaOverride)  // optional per name override
{
    bool useOverride = !betaOverride.empty();
    // just copy over probabilities
    for (int i = 0; i < basketInfo.size(); i++){
        basketInfo[i]->survival = survivalProb[i];
        basketInfo[i]->pdep = floorSurvivalProb[i];
        if (useOverride){
            basketInfo[i]->beta = betaOverride[i];
        }
    }
    if (cpty.get()){
        cpty->survival      = counterPartyProb;
    }
}

/** creates basketInfo field for RFL. */
void CCMTrancheCalculatorLegacy::createBasketInfoRFL(
    CreditTrancheLossConfigConstSP  tranche,         /* (I) */
    CCMConvolution::NameParamArray& basketInfo)      /* (O) */
{
    int numNames = tranche->numInnerLossConfigs();
    basketInfo.clear();
    basketInfo.reserve(numNames);
    // the whole approach here is pretty weak since all but one of the 
    // parameters is time independent yet we copy everything over each time.
    // Should review all of this - not sure how necessary it really is
    for (int i = 0; i < numNames; ++i) {
        CCMConvolution::RFLNameParamSP nm(new CCMConvolution::RFLNameParam());
        //add to the basket
        basketInfo.push_back(nm);
    }
}

/** creates basketInfo field for non RFL. */
void CCMTrancheCalculatorLegacy::createBasketInfo(
    CreditTrancheLossConfigConstSP  tranche,         /* (I) */
    CCMConvolution::NameParamArray& basketInfo)      /* (O) */
{
    int numNames = tranche->numInnerLossConfigs();
    basketInfo.clear();
    basketInfo.reserve(numNames);
    // the whole approach here is pretty weak since all but one of the 
    // parameters is time independent yet we copy everything over each time.
    // Should review all of this - not sure how necessary it really is
    for (int i = 0; i < numNames; ++i) {
        CCMConvolution::NameParamSP nm(new CCMConvolution::NameParam());
        //add to the basket
        basketInfo.push_back(nm);
    }
}

/** creates basketInfo field for non RFL. */
CCMConvolution::NameParamSP CCMTrancheCalculatorLegacy::createCptyInfo(
    CounterPartyCreditConstSP counterParty)      /* (I) */
{
    CCMConvolution::NameParamSP cpty;
    if (!!counterParty) {
        cpty.reset(new CCMConvolution::NameParam());
    }
    return cpty;
}

/** creates basketInfo field for non RFL. */
CCMConvolution::NameParamSP CCMTrancheCalculatorLegacy::createCptyInfoRFL(
    CounterPartyCreditConstSP counterParty)       /* (I) */
{
    CCMConvolution::NameParamSP cpty;
    if (!!counterParty) {
        cpty.reset(new CCMConvolution::RFLNameParam());
    }
    return cpty;
}


/** populates basketInfo field. Same as CDO::getBasketInfo but is on this
    class and operates through ConvolutionProduct and is for CreditMetrics -
    see comments above about splitting up class into 2.*/
void CCMTrancheCalculatorLegacy::populateBasketInfoCM(
    CreditTrancheLossConfigConstSP  tranche,         /* (I) */
    const bool                      recoverNotional, /* (I) */
    CCMConvolution::NameParamArray& basketInfo)      /* (M) */
{
    int numNames = tranche->numInnerLossConfigs();
    // the whole approach here is pretty weak since all but one of the 
    // parameters is time independent yet we copy everything over each time.
    // Should review all of this - not sure how necessary it really is
    for (int i = 0; i < numNames; ++i) {
        CCMConvolution::NameParamSP nm = basketInfo[i]; // for ease
        // NB Non credit metrics params are set to zero in constructor
        nm->beta          = tranche->nameBeta(i);
        tranche->nameLGDRanges(i, nm->lgdNotional, nm->lgdFloor, nm->lgdCap);
        nm->pdep          = 1.0;
        SingleCreditAssetConstSP asset(tranche->nameAsset(i));
        nm->nameId        = asset->getName();
        nm->ntl           = tranche->nameDefaulted(i)? 
            0.0: tranche->nameNotional(i);
        if (recoverNotional) { //we invert the recovery rate
            nm->R = 1 - tranche->nameRecovery(i);
        } 
        else {
            nm->R = tranche->nameRecovery(i);
        }
    }
}

/** populates cpty field. Same as CDO::getCptyInfo but is on this
    class and operates through ConvolutionProduct and is for CreditMetrics -
    see comments above about splitting up class into 2.*/
void CCMTrancheCalculatorLegacy::populateCptyInfoCM(
    CounterPartyCreditConstSP   counterParty, /* (I) */
    CCMConvolution::NameParamSP cpty)         // (M)  
{
    if (!!counterParty) {
        cpty->beta          = counterParty->getBeta();
        cpty->pdep          = 1.0;
        cpty->nameId        = counterParty->getName();
        cpty->R             = counterParty->getRecovery();
    }
}

/** populates existing basketInfo parameter with additional inputs needed by 
    CCM */
void CCMTrancheCalculatorLegacy::populateBasketInfoCCM(
    CreditTrancheLossConfigConstSP  tranche,     /* (I) */
    CCMConvolution::NameParamArray& basketInfo)  // (M)  
{
    int numNames = basketInfo.size();
    for (int i = 0; i < numNames; ++i) {
        CCMConvolution::NameParamSP nm = basketInfo[i]; // for ease
        SingleCreditAssetConstSP asset(tranche->nameAsset(i));
        CcmOnlyParametersConstSP ccmp(
            CCMLossCalculatorBase::getCCMEngineParameters(asset)); 
        nm->indep         = ccmp->getIndependenceFactor();
        nm->betaRec       = ccmp->getBetaR();
        nm->cataRecFactor = ccmp->getCatastrophicRecoveryFactor();        
        nm->lambda        = ccmp->getRecDispersion();
        nm->qM            = ccmp->getQM();

        // tweak the beta historical if a beta Tweak is specified
        if (ccmp->getBetaTweak() != 0.0) {
            // use the beta tweak from BC just for not duplicating the code
            nm->beta = PortfolioName::betaTweak(nm->beta, ccmp->getBetaTweak());
        }
    }
}

/** populates cpty field with additional info needed for CCM - needs 
    refactoring */
void CCMTrancheCalculatorLegacy::populateCptyInfoCCM(
    CounterPartyCreditConstSP  counterParty, /* (I) */
    CCMConvolution::NameParam& cpty)         // (M)  
{
    if (!!counterParty) {
        SingleCreditAssetConstSP asset(counterParty->getAsset());
        CcmOnlyParametersConstSP ccmp(
            CCMLossCalculatorBase::getCCMEngineParameters(asset)); 
        cpty.indep         = ccmp->getIndependenceFactor();
        cpty.cataRecFactor = ccmp->getCatastrophicRecoveryFactor();        
        cpty.qM            = ccmp->getQM();
    
        // tweak the beta historical if a beta Tweak is specified
        if (ccmp->getBetaTweak() != 0.0) {
            // use the beta tweak from BC just for not duplicating the code
            cpty.beta = PortfolioName::betaTweak(cpty.beta, ccmp->getBetaTweak());
        }
    }
}


/** populates existing basketInfo parameter with additional inputs needed by 
    CCM */
void CCMTrancheCalculatorLegacy::populateBasketInfoRFL(
    CreditTrancheLossConfigConstSP  tranche,      /* (I) */
    CCMConvolution::NameParamArray& basketInfo)   // (M)  
{
    int numNames = basketInfo.size();
    for (int i = 0; i < numNames; ++i) {
        CCMConvolution::RFLNameParamSP nm(
            dynamic_cast<CCMConvolution::RFLNameParam*>(basketInfo[i].get()));

        SingleCreditAssetConstSP asset(tranche->nameAsset(i));
        RflOnlyParametersConstSP rflp(
            CreditMetricsLossCalculatorBase::getRFLEngineParameters(asset)); 

        // populate the RFL parameters into RFLNameParam as a CONST pointer
        nm->rflParam = rflp;
    }
}


/** populates cpty field with additional info needed for RFL - needs 
    refactoring */
void CCMTrancheCalculatorLegacy::populateCptyInfoRFL(
    CounterPartyCreditConstSP  counterParty, /* (I) */
    CCMConvolution::NameParam& cpty)         // (M)  
{
    if (!!counterParty) {
        SingleCreditAssetConstSP asset(counterParty->getAsset());
        RflOnlyParametersConstSP rflp(
            CreditMetricsLossCalculatorBase::getRFLEngineParameters(asset)); 

        CCMConvolution::RFLNameParamSP rflCpty(
            dynamic_cast<CCMConvolution::RFLNameParam*>(&cpty));

        // populate the RFL parameters into RFLNameParam as a CONST pointer
        rflCpty->rflParam = rflp;    
    }
}


CCMTrancheCalculatorLegacy::CCMTrancheCalculatorLegacy():CObject(TYPE) {}

CCMTrancheCalculatorLegacy::CCMTrancheCalculatorLegacy(CCMTrancheCalculatorLegacy* copy):CObject(TYPE)
{
    basketInfo = *CCMConvolution::NameParamArraySP((dynamic_cast<CCMConvolution::NameParamArray*>(copy->basketInfo.clone())));

    if (copy->cpty.get() != 0)
    {
        cpty    = CCMConvolution::NameParamSP( copy->cpty.clone());
    }

    maxl = copy->maxl;
    maxs = copy->maxs;
    pastLoss = copy->pastLoss;
    lossUnit = copy->lossUnit;
    ntl = copy->ntl;
    density = copy->density;
    densityCond = copy->densityCond;
}

CCMTrancheCalculatorLegacy::~CCMTrancheCalculatorLegacy() {}

pair<DoubleArraySP, DoubleArraySP> CCMTrancheCalculatorLegacy::getDensities() const {
    DoubleArraySP first(new DoubleArray(density));
    DoubleArraySP second(new DoubleArray(densityCond));
    return pair<DoubleArraySP, DoubleArraySP>(first, second);
}
    
int CCMTrancheCalculatorLegacy::hashCode() const {
    int hCode = 0;
    for (int i = 0; i < basketInfo.size(); ++i) {
        hCode ^= basketInfo[i]->hashCode();
    }
    if (cpty.get() != 0) {

        hCode ^= cpty->hashCode();
    }
    return hCode;
}

bool CCMTrancheCalculatorLegacy::equal(CCMTrancheCalculatorLegacy* c ) {
    bool isEqual = true;
    if (basketInfo.size() != c->basketInfo.size()) {
        return false;
    }
    for (int i = 0; i < basketInfo.size(); ++i) {
        isEqual = isEqual && c->basketInfo[i]->equals(basketInfo[i]);
    }
    if (cpty.get() != 0) {
        isEqual = isEqual && c->cpty->equals(cpty);
    }
    isEqual = isEqual && maxl == c->maxl;
    isEqual = isEqual && maxs == c->maxs;
    isEqual = isEqual && pastLoss == c->pastLoss;
    isEqual = isEqual && lossUnit == c->lossUnit;
    return isEqual;
}

IObject* CCMTrancheCalculatorLegacy::defaultCCMTrancheCalculatorLegacy()
{
    return new CCMTrancheCalculatorLegacy();
}

void CCMTrancheCalculatorLegacy::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CCMTrancheCalculatorLegacy, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCCMTrancheCalculatorLegacy);

    //transient fields
    FIELD( pastLoss,    "");
    FIELD( basketInfo,  "");
    FIELD       ( cpty,        "");
    FIELD( ntl,         "");
    FIELD( density,     "");
    FIELD( densityCond, "");
    FIELD( maxs,        "");
    FIELD( maxl,        "");
    FIELD( lossUnit,    "");
}

CCMConvolutionSP CCMTrancheCalculatorLegacy::convolutionNoCache()
{
    static const string method = "CCMTrancheCalculatorLegacy::convolutionNoDensityCache";
    try
    {
        CCMConvolutionSP conv = CCMConvolutionSP(
            new CCMConvolution(
                false, //use fast
                basketInfo,
                maxl, 
                maxs));

        conv->calcLossDensityRecM4();

        return conv;

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** 
 * function calculating expected loss given strikes 
 * can be invoked directly or as a callback by london floor adjustment function
 * @return 0 for success, -1 if the first failed, -2 if the conditional loss failed
 */
void CCMTrancheCalculatorLegacy::loss(
    double        K1,    /* (I) lower strike      */
    double        K2,    /* (I) upper strike      */
    double       &L,     /* (O) tranche loss amt  */
    double       &Lcond) /* (O) tranche loss amt cond on cpty surviving */
    const
{
    L = CCMTrancheUtils::calcExpectedTrancheLoss(K1,
                                                 K2,
                                                 pastLoss,
                                                 lossUnit,
                                                 density,
                                                 maxs,
                                                 maxl);

    Lcond = CCMTrancheUtils::calcExpectedTrancheLoss(K1,
                                                     K2,
                                                     pastLoss,
                                                     lossUnit,
                                                     densityCond,
                                                     maxs,
                                                     maxl);
}

void CCMTrancheCalculatorLegacy::convolution(const Control*            control,
                                             const CreditMetricsModel* model,
                                             const int                 timepoint,
                                             const int                 cacheSet)
{
    CCMTrancheCalculatorLegacySP key(new CCMTrancheCalculatorLegacy(this));
    key->density.resize(maxs+maxl+1);
    key->densityCond.resize(maxs+maxl+1);
    
    pair<DoubleArraySP, DoubleArraySP> densities = CDOCache::get(key,
                                                                 control,
                                                                 model,
                                                                 timepoint,
                                                                 cacheSet);
    int n = densities.first->size();
    for (int i = 0; i < n; ++i) {
        density[i] = (*densities.first)[i];
        densityCond[i] = (*densities.second)[i];
    }
}

void CCMTrancheCalculatorLegacy::convolutionNoDensityCache(
    const Control*            control,
    const CreditMetricsModel* model,
    const int                 timepoint,
    const int                 cacheSet)
{
    static const string method = "CCMTrancheCalculatorLegacy::convolutionNoDensityCache";
    try
    {
        density.resize(maxs+maxl+1);
        densityCond.resize(maxs+maxl+1);

        CCMConvolutionSP convo;
        convo = convolutionNoCache();
        convo->getDensities(cpty, density, densityCond);

/* this functionality is bypassed for the time being
   ConvolutionCacheSP cache = model->getConvolutionCache();

   //can we make use of previous convolution results ?
   if (control->isPricing())
   {
   // no - we need to convolute
   convo = convolutionNoCache();

   //and store the result
   if (model->allowQuickGreeks())
   {
   CCMTrancheCalculatorLegacySP thisAsSP = CCMTrancheCalculatorLegacySP(this);
   cache->save(cacheSet, timepoint, thisAsSP, convo);
   }
   }
   else
   {
   //maybe
        
   //does the model permit it
   if (model->allowQuickGreeks())
   {
   CCMTrancheCalculatorLegacySP untweakedState;
   CCMConvolutionSP       untweakedResult;

   cache->get(cacheSet,timepoint,untweakedState,untweakedResult);

   //ensure we got something back
   if (!untweakedState || !untweakedResult)
   {
   //something went wrong, so convolute from scratch
   //but dont store
   convo = convolutionNoCache();
   }
   else
   {
   //are we performing a candidate tweak ?
   //i.e. one which results in a change to candidate parameters
   //in a single name

   int nmIdx = NameParam::singleTweakDiff(basketInfo,
   untweakedState->getBasketInfo());
   if (nmIdx != -1)
   {
   //yes - reuse the baseline result
   if (nmIdx == 0)
   {
   //actually no difference, so just return the baseline result
   convo = untweakedResult;
   }
   else
   {
   //deep copy the baseline result
   convo = CCMConvolutionSP(untweakedResult.clone());

   //reconvolute this name
   convo->reconvolute((untweakedState->getBasketInfo())[nmIdx-1],
   basketInfo[nmIdx-1]);
   }
   }
   else
   {
   //no - we need to convolute, but NOT store
   // Here is the time consuming operation !
   convo = convolutionNoCache();
   }
   }
   }
   else
   {
   //no - we need to convolute, but NOT store
   // Here is the time consuming operation !
   convo = convolutionNoCache();
   }
   }
   //now get the densities
   convo->getDensities(cpty, density, densityCond);
*/
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const CCMTrancheCalculatorLegacy::TYPE = CClass::registerClassLoadMethod(
    "CCMTrancheCalculatorLegacy", typeid(CCMTrancheCalculatorLegacy), load);

DEFINE_TEMPLATE_TYPE(CCMTrancheCalculatorLegacyArray);


DRLIB_END_NAMESPACE
