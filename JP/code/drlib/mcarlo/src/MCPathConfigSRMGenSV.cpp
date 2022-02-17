//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRMGenSV.cpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MCPathConfigSRMGenSV.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/SRMGenDiffusibleAsset.hpp"
#include "edginc/SVGenExpectedSpot.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/IRngSlotAssigner.hpp"

#include <functional>
#include <algorithm>
#include <typeinfo>

USING_DRLIB_NAMESPACE

// as this convenient SGI STL extension is not in STLport
struct my_select2nd {
    template <typename K, typename V>
    const V& operator() (const std::pair<K,V> &p) { return p.second; }
};

/************************************************************************/
/*  Non-asset-type specific set of functions                            */
/************************************************************************/

// constructor
SV::SV(MCPathConfigSRM*       pathConfig,
       const MCProductClient* prodClient): // (M)
    multiFactors(prodClient->getMultiMarketFactors()),
    pastValues(prodClient->getMCPastValues()),
//    storedNumIRFactors(-1),
    storedNumAllFactors(-1)
{
    static const string method("SV::SV");

    // get hold of payoff yield curve (create mapIRGen map)
    domYC.reset(prodClient->getDiscount());
    saveDiscountYC(pathConfig, domYC);
    domISOCode = domYC->getCcy();

    // Check that all FX assets have the domestic ccy as their base ccy
    validateBaseFXCurrencies(pathConfig->diffAssets, domYC);
     
    // Collect state variables from product and categorize them
    StateVariableCollectorSP svCollector(new StateVariableCollector());
    prodClient->collectStateVars(svCollector); // collect them
    // and process them. 

    IElemStateVariableGenArray stateVarGenArray(svCollector->getElemStateVarGens());
    spots = filterStateVars<const SVGenSpot>(stateVarGenArray); // modifies stateVarGenArray

    // here: need a very special handling for the SVGenPathWeight variables

    bool isCID = pathConfig->isCreditCID();

    // Build map of diffusing assets and connect to Multifactors where appropriate
    extractFromMultiFactors(pathConfig, domISOCode, pathConfig->diffAssets, spots); 
    // Build map for credit ... 
    if (isCID)
        buildPerCreditMapCID(pathConfig, pathConfig->diffCDSs);
    else
        buildPerCreditMap(pathConfig, pathConfig->diffCDSs);
	// Build map for energy ...
	buildPerEnergyMap(pathConfig, pathConfig->diffEnrgs);
    // Build map for basis currencies
    buildPerBasisMap(pathConfig, pathConfig->diffBasisArray);

    //setup currency treatment for all gens
    setupAssetCurrencyTreatment();

    size_t numSVHandled = processAllSVGens(svCollector);


    // then get hold of the YC's to diffuse/order the ccys 
    orderedCcys.resize(mapIRGen.size());
    int pos = 1; // pos = 0 reserved for discount ccy

	if (mapIRGen.find(domISOCode) == mapIRGen.end()) // the cycle below implies _there is_ domISOCode or we get memory corruption
	    throw ModelException(method, (string)"Domestic Currency "+ domISOCode +" is not in the diffused list!");

    for (map<string, QMCGenDiffusibleIRSP>::iterator iter = mapIRGen.begin();
         iter != mapIRGen.end(); ++iter){
        // retrieve numFactors
        int nFactors = pathConfig->getNumIRFactors(iter->first);
        SRMGenDiffusibleIRDeterm* pSRMGenDiffusibleIRDet = dynamic_cast<SRMGenDiffusibleIRDeterm*>(iter->second.get());
        QLIB_VERIFY((nFactors == 0 && pSRMGenDiffusibleIRDet != NULL) 
            || (nFactors > 0 && pSRMGenDiffusibleIRDet == NULL), 
            "Number of interest rate factors for ccy " + iter->first +
            " is inconsistent with the rates model type " + pathConfig->ratesModelType);       
        iter->second->nFactors = nFactors;
 
        // retrieve yc
        iter->second->ycs.second = 
            pathConfig->getDiffusedYC(iter->second->ycs.first); // ycs.first = ptr to dicounting yc; ycs.second = ptr to diffusion yc
        // then populate array - put discount ccy first (match SRM3)
        if (iter->first == domISOCode){
            orderedCcys[0] = iter->second;
        } else {
            orderedCcys[pos] = iter->second;
            pos++;
        }
    }


    transform(mapFXGen.begin(),     mapFXGen.end(),     back_inserter(orderedFXs),      my_select2nd());    // For FX
    transform(mapEquityGen.begin(), mapEquityGen.end(), back_inserter(orderedEquities), my_select2nd());    // For EQ
    // TODO: make sure we do not do this for CID credits, as the order is very important and assigned somewhere else
    if (!isCID)
       transform(mapCreditGen.begin(), mapCreditGen.end(), back_inserter(orderedCredits),  my_select2nd());    // For CR
    transform(mapEnergyGen.begin(), mapEnergyGen.end(), back_inserter(orderedEnergys),  my_select2nd());    // For Energy
    transform(mapBasisGen.begin(),  mapBasisGen.end(),  back_inserter(orderedBasis),    my_select2nd());    // For Basis
    // sort the vector of QMCGenDiffusibleEnergySPs in increasing tier level order:
    DiffEnrgGenSPComparator diffEnrgGenSPCompare;
    sort(orderedEnergys.begin(), orderedEnergys.end(), diffEnrgGenSPCompare);

    for(size_t i=0; i<orderedFXs.size(); ++i)  
        orderedFXs[i]->multiSpots = &spots;

    setByAssetPosition(); // initialize byAssetPosition
    // create vector sv->byAssetPosition in the order of:
	// IR, FX, EQ, CR and Enrg
    assignRandomIdx(prodClient); // create vector sv->byAssetPosition in the order of: IR, FX, EQ, CR, Enrg, Basis
}

/** We check that all FX asset have the dom ccy as its baseCcy. 
  * If the risk ccy was the dom ccy, then the dom ccy would get a nonzero pQMCGenDiffusibleFX, 
  * which would make some domestic assets think they were Quantoed! 
  */
void SV::validateBaseFXCurrencies(const CAssetArray& diffAssets, YieldCurveConstSP domYC)
{
    string domCcy = domYC->getCcy();
    for (int i = 0;i<diffAssets.size();i++) {
        if (FXAsset::TYPE->isInstance(diffAssets[i])){
            FXAssetSP fxAsset = FXAssetSP::dynamicCast(diffAssets[i]);
            YieldCurveConstSP ycRisk = fxAsset->getRiskCcy();
            YieldCurveConstSP ycBase = fxAsset->getBaseCcy();
            string riskCcy = ycRisk->getCcy();
            string baseCcy = ycBase->getCcy();
            string name = fxAsset->getName();
            QLIB_VERIFY(baseCcy == domCcy, 
                "The FX asset " + name + " should have " + domCcy + 
                " as its base currency. Its base currency is currently " + 
                baseCcy);
        }
    }
}

void SV::setByAssetPosition()
{
    byAssetPosition.clear();

    for (size_t i = 0; i<orderedCcys.size(); ++i)
       byAssetPosition.push_back(orderedCcys[i]);

    for (size_t i = 0; i<orderedPureJumps.size(); ++i)  // so far only CID -- but other jump-diffusions might use it later
        byAssetPosition.push_back(orderedPureJumps[i]);

    for (size_t i = 0; i<orderedFXs.size(); ++i)
        byAssetPosition.push_back(orderedFXs[i]);
    
    for (size_t i = 0; i<orderedEquities.size(); ++i)
        byAssetPosition.push_back(orderedEquities[i]);
    
    for (size_t i = 0; i<orderedCredits.size(); ++i)
        byAssetPosition.push_back(orderedCredits[i]);

    for (size_t i = 0; i<orderedBasis.size(); ++i)
        byAssetPosition.push_back(orderedBasis[i]);
    
    for (size_t i = 0; i<orderedEnergys.size(); ++i)
        byAssetPosition.push_back(orderedEnergys[i]);
}

// helper template to transform map<string, smartPtr> into map<ptr, string> 
template <typename X, typename Y, typename Z>
void flip(const map<X, Y>& in, map<Z, X>& out)
{
	for(typename  map<X, Y>::const_iterator it = in.begin(); it != in.end(); ++it)
		out[it->second.get()]=it->first;
}

/** Now arrange all the assets in some order and assign all the randomIdx. 
    For regression purposes we continue using ordering like in original SRM3:
    IR, FX, EQ, CR and Enrg */

void SV::assignRandomIdx(const MCProductClient* prodClient)
{
	IRngSlotAssignerSP defAssigner;
	const IRngSlotAssigner * slotAssigner = dynamic_cast<const IRngSlotAssigner*>(prodClient); // check if product knows how to assign name -> rngSlot

	if (slotAssigner == NULL)
	{
		defAssigner = IRngSlotAssignerSP(new RngSlotAssignerLinear);
		slotAssigner = defAssigner.get();
	}
	
    int globalIdx = 0;

	map<IQMCGenDiffusibleAsset*, string> asset2name;
//     for(map<string, QMCGenDiffusibleIRSP>::iterator it=mapIRGen.begin();
//         it != mapIRGen.end();
//         ++it)
//         asset2name[it->second.get()]=it->first;

	flip(mapIRGen, asset2name);
	flip(mapFXGen, asset2name);
	flip(mapEquityGen, asset2name);
	flip(mapCreditGen, asset2name);
	flip(mapEnergyGen, asset2name);
	flip(mapBasisGen, asset2name);
	
    
    for (size_t i = 0; i != byAssetPosition.size(); ++i)
    {
        /** 
         * There are two strategies to assign rnadom slots to assets.
         * The default one is to assign slots as they needed without any gaps.
         * This works in many cases, but sometimes we want to guarantee that same assets get same slots regardless of the new assets (ex: qSRM).
         * Then, we expect that slot numbers are passed from the outside, and we merely use them.
         * 
         * returns:
         *    -1 -- asset should be excluded as we exhausted available slots for this kind of assets
         *  >= 0 -- the assigned rngSlot for this asset.
         * We pass in: 
         * - how many slots were used by assets of this kind
         * - asset's name
         * - number of factors (to see if we fit the 
         * */

		if (byAssetPosition[i]->numFactors() > 0) {
			string assetName = asset2name[byAssetPosition[i].get()]; // assets without name will be resolved to ""
			int rngIdx = slotAssigner->getRngIdx(globalIdx, assetName);
			QLIB_VERIFY(rngIdx >= 0, "Random slots expected to be positive");
			byAssetPosition[i]->setRandomIdx(rngIdx); // rngIdx is += nfactors
			globalIdx = std::max(globalIdx, rngIdx); // keep track of the #of factors used
		}
#if 0
        cerr << "Asset [" << i << "] (" << asset2name[byAssetPosition[i].get()] << ") has " << byAssetPosition[i]->numFactors() << " factors and assigned slot " << byAssetPosition[i]->getRandomIdx() << endl;
#endif
    }

    /** the total number of factors -- i.e. the size or Random noise matrix to produce */
    storedNumAllFactors = globalIdx;
    verifyRandomIdx();

}
    
/// returns iff random number slots do not intersect otherwise throws exception
void SV::verifyRandomIdx(void) {
    map<int, size_t> slot2pos; // map to map rngSlot to byAssetPosition 
    for(size_t i = 0; i != byAssetPosition.size(); ++i)
        if (byAssetPosition[i]->numFactors() > 0) {
            int slot = byAssetPosition[i]->getRandomIdx();
            QLIB_VERIFY(slot2pos.find(slot) == slot2pos.end(),
                "ERROR: Slot " + Format::toString(slot) + " is used by both assets " +  Format::toString(i) + " and " + Format::toString(slot2pos[slot]));
            
            slot2pos[slot]=i;
        }
    // check that slots do not intersect        
    int nextAvailable = 0;
    map<int, size_t>::iterator prevIt;
    for(map<int, size_t>::iterator it = slot2pos.begin();
        it != slot2pos.end();
        prevIt=it, ++it)
    {
        QLIB_VERIFY(it->first >= 0, "Negative slots are not allowed");
        QLIB_VERIFY(it->first >= nextAvailable,"ERROR: Asset " + Format::toString( it->second) + " assigned slot " + Format::toString( it->first) + " that intersects with the slot assigned to asset " + Format::toString( prevIt->second) +" that was assigned slot " + Format::toString( prevIt->first));
        nextAvailable = it->first+byAssetPosition[it->second]->numFactors();
    }
}
/** Creates a QMCGenDiffusibleIRSP according to the Rates model type */
QMCGenDiffusibleIRSP SV::diffusibleIRGenSPFromModelType(string ratesModelType){
    static const string method("SV::diffusibleIRGenSPFromModelType");
    if (ratesModelType == MCPathConfigSRM::RATES_MODEL_TYPE_HJM)
        return QMCGenDiffusibleIRSP(new SRMGenDiffusibleIRHJM());
    else if (ratesModelType == MCPathConfigSRM::RATES_MODEL_TYPE_LIBOR)
        return QMCGenDiffusibleIRSP(new SRMGenDiffusibleIRLibor());
    else if (ratesModelType == MCPathConfigSRM::RATES_MODEL_TYPE_DETERM)
        return QMCGenDiffusibleIRSP(new SRMGenDiffusibleIRDeterm());
    else
        throw ModelException(method, "ratesModelType "+ratesModelType
        +" is unknown. Please select one of {"+
            MCPathConfigSRM::RATES_MODEL_TYPE_HJM+", "+
            MCPathConfigSRM::RATES_MODEL_TYPE_DETERM+", "+
            MCPathConfigSRM::RATES_MODEL_TYPE_LIBOR+"}");
}

/** Saves this discount yield curve in a map against the isocode of the ccy */
void SV::saveDiscountYC(MCPathConfigSRM* pathConfig,
                        YieldCurveConstSP discYC){
    static const string method("SV::saveDiscountYC");
    QMCGenDiffusibleIRSP& ccy = mapIRGen[discYC->getCcy()];
    if (!ccy.get()) {
        ccy = diffusibleIRGenSPFromModelType(pathConfig->getRatesModelType());
    }
    ccy->name = discYC->getCcy();

    ccy->ycs.first = IYieldCurveConstSP::dynamicCast(discYC); // save discount
}


/** Creates a QMCGenDiffusibleCreditSP according to the Credit model type */
QMCGenDiffusibleCreditSP SV::diffusibleCreditGenSPFromModelType(string creditModelType, MCPathConfigSRM* pathConfig){
    static const string method("SV::diffusibleCreditGenSPFromModelType");
    if (creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_HJM)
        return QMCGenDiffusibleCreditSP(new SRMGenDiffusibleCreditHJM());
    else if (creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_CIR)
        return QMCGenDiffusibleCreditSP(new SRMGenDiffusibleCreditCIR());
    else if (creditModelType == MCPathConfigSRM::CREDIT_MODEL_TYPE_LIBOR)
        // credit Libor model switched off for now
        throw ModelException(method, "creditModelType "+creditModelType
        +" is currently switched off");
        // return QMCGenDiffusibleCreditSP(new SRMGenDiffusibleCreditLibor());
    else
        throw ModelException(method, "creditModelType "+creditModelType
        +" is unknown. Please select one of {"+ pathConfig->getAllSupportedCreditModels() +"}");
}

/** Saves CDS curve in a map against the name of the credit */
void SV::buildPerCreditMap(MCPathConfigSRM*     pathConfig,
                           ICDSParSpreadsArray& diffCDSs) {
    static const string method("SV::buildPerCreditMap");
    /** extend the mapCreditGen map (ie add this credit name) and fill in basic information 
        similar to what is done in extractFromMultiFactors */ 
    numCr = 0;
    for (int i = 0; i<pathConfig->diffCDSs.size(); i++) {
        ICDSParSpreadsConstSP cdsCurve = pathConfig->diffCDSs[i];
        string creditName = cdsCurve->getName();
        QMCGenDiffusibleCreditSP& perCreditData = mapCreditGen[creditName]; // extend the map (by checking whether this credit name already exists)
        if (!perCreditData.get()) // build one if necessary  (alas, should always be)
            perCreditData = diffusibleCreditGenSPFromModelType(pathConfig->getCreditModelType(), pathConfig);
        perCreditData->name = creditName; // useful once these are no longer in a map
        perCreditData->ccy = cdsCurve->getCcy(); 
        perCreditData->cds = diffCDSs[i]; // used only to create util!
        numCr++;
    }
    // update any xccy info needed for credit 
    for (map<string, QMCGenDiffusibleCreditSP>::iterator iter = mapCreditGen.begin();
         iter != mapCreditGen.end(); ++iter)
    {
        QMCGenDiffusibleIRSP pccy = mapIRGen[iter->second->ccy];
        iter->second->pQMCGenDiffusibleIR = pccy.get();

        if (iter->second->ccy != domISOCode) {
             // should exist already
            if (!pccy->pQMCGenDiffusibleFX || !pccy->pQMCGenDiffusibleFX->fxAsset.get()) {
                throw ModelException(method,
                                     "Incomplete FX data!");
            }
        }
    }

}

// this is a specialized version that builds CID-specific set of objects
// because it differs from "a bunch of single-name diffusions" approach
/** Saves CDS curve in a map against the name of the credit for CID model */
void SV::buildPerCreditMapCID(MCPathConfigSRM*     pathConfig,
                              ICDSParSpreadsArray& diffCDSs) {
    static const string method("SV::buildPerCreditMapCID");
    /** extend the mapCreditGen map (ie add this credit name) and fill in basic information 
        similar to what is done in extractFromMultiFactors */ 

    //  if we are here, then the requested model is CID and we will do only CID, i.e. 
    //  no single-name non-CID diffusions are allowed for the time being.
    QMCGenDiffusibleIRSP pDomCcy = mapIRGen[domISOCode];

    // (1) create generator of common markets, not in the name-based map
    QMCGenDiffusibleCreditCIDSP commonMarket(new QMCGenDiffusibleCreditCID(true));
    commonMarket->name = MCPathConfigSRM::COMMON_MARKET_FACTOR_FOR_CID; // useful once these are no longer in a map
    commonMarket->ccy = domISOCode;  // hard-coded for the time being
    commonMarket->pQMCGenDiffusibleIR = pDomCcy.get();
    commonMarket->isFullMC = true;   // no matter what model - we use full diffusion for common markets
    orderedCredits.push_back(commonMarket);

    // (2) create generator of pure jumps, not in the name-based map
    QMCGenPureJumpSP pureJumps(new QMCGenPureJump);
    pureJumps->name = "Pure Poisson Generator for CID";
    orderedPureJumps.push_back(pureJumps);
    // TODO: need to store it as a class member


    const ICIDParameters* pCIDParameters = pathConfig->getCIDParameters();
    QLIB_VERIFY(pCIDParameters != NULL, 
        "Cannot initialize CID jumps for credit -- CIDParameters are not defined in MCPathConfig.");

    for (int i = 0; i<pathConfig->diffCDSs.size(); i++) {
        ICDSParSpreadsConstSP cdsCurve = pathConfig->diffCDSs[i];
        string creditName = cdsCurve->getName();


        // check that the name exists in the CID Parameters
        if (! pCIDParameters->getSingleNameRecord(creditName).get())
            continue; // for the time being we just continue

        // weight of common market
        double A = pCIDParameters->getSingleNameRecord(creditName)->A;

        // getting recovery rates
        double catastrophicRecovery = pCIDParameters->getSingleNameRecord(creditName)->cataR;
        double recoveryRate = cdsCurve->getRecovery();

        // (3) create jump contributions per name, not in the name-based map
        QMCGenDiffusibleCreditCIDJumpsSP cidJumps(new QMCGenDiffusibleCreditCIDJumps);
        cidJumps->name = creditName;
        cidJumps->ccy = domISOCode; // hard-coded for now
        cidJumps->pQMCGenDiffusibleIR = pDomCcy.get();
        cidJumps->pJumpSource = pureJumps;
        cidJumps->isFullMC = pathConfig->isFullCreditMC();
        orderedCredits.push_back(cidJumps);

        // (4) create idiosyncratic diffusion gen per name, not in the name-based map
        QMCGenDiffusibleCreditCIDSP cidDiffuse(new QMCGenDiffusibleCreditCID(false));
        cidDiffuse->name = creditName;
        cidDiffuse->ccy = domISOCode; // hard-coded for now
        cidDiffuse->pQMCGenDiffusibleIR = pDomCcy.get();
        cidDiffuse->isFullMC = pathConfig->isFullCreditMC();
        orderedCredits.push_back(cidDiffuse);

        // (5) create composites, finally: the only thing in the name-based map (for attaching SV)
        QMCGenDiffusibleCreditCompositeSP pComposite(new QMCGenDiffusibleCreditComposite);
        pComposite->name = creditName; // useful once these are no longer in a map
        pComposite->ccy = domISOCode;
        pComposite->pQMCGenDiffusibleIR = pDomCcy.get();
        pComposite->regularRecoveryRate = recoveryRate;    // setting up the recovery rate
        pComposite->catastrophicRecoveryRate = catastrophicRecovery;
        pComposite->addComponent(cidDiffuse, 1.0);  // idiosyncratic part
        pComposite->addComponent(cidJumps, 1.0);    // jumps
        pComposite->addComponent(commonMarket, A);  // common market part
        orderedCredits.push_back(pComposite);

        QMCGenDiffusibleCreditSP& perCreditData = mapCreditGen[creditName]; // extend the map (by checking whether this credit name already exists)
        perCreditData = pComposite;

    }

    numCr = orderedCredits.size();
}

void SV::buildPerEnergyMap(
	MCPathConfigSRM * pathConfig,
	EnergyFuturesCurveArray & diffEnrgs
	)
{
	const string method = "SV::buildPerEnergyMap()";

	// loop thru each energy curve, add them to the map and do some book keeping
    numEnrg = 0;
	for (int i = 0; i < diffEnrgs.size(); ++i) 
	{
		EnergyFuturesCurveConstSP energyCurve = diffEnrgs[i];
		const string & energyName = energyCurve->getName();
        EnergyUnderlyerConstSP energyUnderlyer = energyCurve->getEnergyUnderlyer();

		// add an entry to the energy map 
		QMCGenDiffusibleEnergySP & perEnergyData = mapEnergyGen[energyName];
		if (perEnergyData.get()) { // check if we have duplicate curve
			throw ModelException(method, 
				"Duplicate energy curve with name: " + energyName);
		} else {
			perEnergyData = QMCGenDiffusibleEnergySP(new SRMGenDiffusibleEnergy());
		}

        // 1. record energy curve and its name
        // note: all the sv generators reference to the same future curve 
        perEnergyData->futureCurve = energyCurve;
        perEnergyData->name = energyName;
        perEnergyData->modelType = energyUnderlyer->getUnderlyerType();
        perEnergyData->isTier1 = energyUnderlyer->getParentName().empty() ? true : false;

        // 2. record model type and number of factor
        if (perEnergyData->isTier1) 
            perEnergyData->nFactors = pathConfig->getEnergyNumFactor(energyName);	// get the number of factor
        else 
            perEnergyData->nFactors = 1; // tier 2 spread hard-coded to have 1 factor (for now)

        // 3. record pricing ccy data
        EnergyUnderlyerConstSP enrgUnderlyer = perEnergyData->futureCurve->getEnergyUnderlyer();
        perEnergyData->ccy = enrgUnderlyer->getPricingCurrency();
        //QMCGenDiffusibleIRSP pccy = mapIRGen[perEnergyData->ccy];
        map<string, QMCGenDiffusibleIRSP>::iterator pccyIt = mapIRGen.find(perEnergyData->ccy);
        if (pccyIt == mapIRGen.end())
            throw ModelException(method, "Ccy data not found for energy name: " + energyName);

        QMCGenDiffusibleIRSP pccy = pccyIt->second;
        perEnergyData->pQMCGenDiffusibleIR = pccy.get();

        // 4. get fx data from pricing ccy for doing quanto adjustment
        if (perEnergyData->ccy != domISOCode) {
            // should exist aleady
            if (!pccy->pQMCGenDiffusibleFX || !pccy->pQMCGenDiffusibleFX->fxAsset.get()) {
                throw ModelException(method, 
                    "FX data not found for energy name: " + energyName);
            }
        }
        ++numEnrg;
	}

    // now go over each energies again to set up tree relationships
    for (map<string, QMCGenDiffusibleEnergySP>::iterator it = mapEnergyGen.begin(); it != mapEnergyGen.end(); ++it) 
    {
        QMCGenDiffusibleEnergySP & perEnergyData = it->second;
        if (!perEnergyData->isTier1) { 
            // is tier 2
            string parentName = perEnergyData->futureCurve->getEnergyUnderlyer()->getParentName();
            map<string, QMCGenDiffusibleEnergySP>::iterator itParent = mapEnergyGen.find(parentName);
            
            // check if parent exist
            if (itParent == mapEnergyGen.end())
                throw ModelException(method, "Parent '" + parentName + 
                    "' not found for tier2 energy: " + perEnergyData->name);
            
            // TODO: check energy tree structural integrity (e.g. circular references)
            // hook up parent energy
            perEnergyData->pParentQMCGenDiffusibleEnergy = itParent->second.get();
        }
    }

    // now go over the energys one last time to compute tier level
    for (map<string, QMCGenDiffusibleEnergySP>::iterator it = mapEnergyGen.begin(); it != mapEnergyGen.end(); ++it) 
    {
        QMCGenDiffusibleEnergySP & perEnergyData = it->second;
        if (!perEnergyData->isTier1) // is higher tier:
            perEnergyData->setTierLevel();
    }
}

void SV::buildPerBasisMap(
	MCPathConfigSRM * pathConfig,
	IBasisIndexCurveArray & diffBasis
	)
{
	const string method = "SV::buildPerBasisMap()";

	// loop through each basis curve, add them to the map and do some book keeping
	for (int i = 0; i < diffBasis.size(); ++i) 
	{
		IBasisIndexCurveConstSP basisCurve = diffBasis[i];
		const string & name = basisCurve->getName();

		// add an entry to the energy map 
		QMCGenDiffusibleBasisSpreadSP & perBasisData = mapBasisGen[name];
		if (perBasisData.get())  // check if we have duplicate curve
			throw ModelException(method, "Duplicate basis curve with name: " + name);
		else
			perBasisData = QMCGenDiffusibleBasisSpreadSP(new SRMGenDiffusibleBasisSpread());
		
		// 1. record basis curve and its name
		// note: all the sv generators reference to the same future curve 
		perBasisData->basisCurve = basisCurve;
		perBasisData->name = name;

		// 2. record pricing ccy and fx data for quanto adjustment (if necessary)
		perBasisData->ccy = basisCurve->getRefCurve()->getName();
		QMCGenDiffusibleIRSP pccy = mapIRGen[perBasisData->ccy];
		perBasisData->pQMCGenDiffusibleIR = pccy.get();

		// verify that fx data from pricing ccy is available
		if (perBasisData->ccy != domISOCode) 
        {
			// should exist already
			if (!pccy->pQMCGenDiffusibleFX || !pccy->pQMCGenDiffusibleFX->fxAsset.get()) 
				throw ModelException(method, "FX data not found for basis curve name: " + name);
			
		}
	}
}

/** FYI (as of 08-Nov-2006): 
  * CAsset::makeAssetUsingCache() is where the Prot or Struck assets are created. 
  * Struck/Prot/Vanilla is also checked in MCPathConfigSRM::getComponentMarketData, 
  * and this is another possible place to put the logic currently in setupCurrencyTreatment
  */
void SV::setCurrencyTreatment(IQMCGenDiffusibleAssetSP perAssetData, IMarketFactorConstSP factor)
{
    CurrencyTreatment ccyT = ccyVanilla;
    if (CAsset::IStruck::TYPE->isInstance(*factor)) {
        ccyT = ccyStruck;
    }
    else if (CAsset::IQuanto::TYPE->isInstance(*factor)) {
        ccyT = ccyProtected;
    }
    if (ccyT != ccyVanilla) {
        perAssetData->setCurrencyTreatment(ccyT);
    }
}

void SV::setupAssetCurrencyTreatment()
{
    for (int i = 0; i < multiFactors->numFactors(); i++){
        // get hold of the actual asset
        IMarketFactorConstSP factor = multiFactors->getFactor(i);
        const string& truename = factor->getName();

        //Currently only Equities can be 'struck', and CUPS adjustment is handled
        //slightly differently for other assets. Uncomment if and when other assets
        //support different currency treatments via multifactors (probably never!!)

        if (CAsset::TYPE->isInstance(factor)) {
            CAssetConstSP asset(CAssetConstSP::dynamicCast(factor));
            const string& truename = asset->getTrueName();
            if (mapEquityGen.find(truename) != mapEquityGen.end()) {
                setCurrencyTreatment(mapEquityGen[truename], factor);
            } /* else if (mapFXGen.find(truename) != mapFXGen.end()) {
                setCurrencyTreatment(mapFXGen[truename], factor);
            } else if (mapCreditGen.find(truename) != mapCreditGen.end()) {
                setCurrencyTreatment(mapCreditGen[truename], factor);
            } else if (mapEnergyGen.find(truename) != mapCreditGen.end()) {
                setCurrencyTreatment(mapEnergyGen[truename], factor);
            } else if (mapBasisGen.find(truename) != mapBasisGen.end()) {
                setCurrencyTreatment(mapBasisGen[truename], factor);
            } */
        }
    }
}

// Build 'perAsset' map and extend 'mapIRGen' map.
// Use compositeAssets and diffAssets to connect to multiFactors/SVs
void SV::extractFromMultiFactors(MCPathConfigSRM* pathConfig,
                                 string         domISOCode,
                                 CAssetArray&   diffAssets,
                                 vector<const SVGenSpot*> &spots){
    static const string method("SV::extractFromMultiFactors");
    // Counts .... for SVs. Though of course we're wrong since properly
    // only need SV (for asset) for an entry in MultiFactors.
    numFX = 0;
    numEq = 0;

    // Build the 'perAsset' map. 
    for(int j=0; j<diffAssets.size(); j++) {
        string assetName = diffAssets[j]->getName();

        if (FXAsset::TYPE->isInstance(diffAssets[j])){ 
            QMCGenDiffusibleFXSP& perAssetData = mapFXGen[assetName]; // extend the map
            if (!perAssetData.get()) // build one if necessary  (alas, should always be)
                perAssetData = QMCGenDiffusibleFXSP(new SRMGenDiffusibleFX());

            perAssetData->name = assetName; // useful once these are no longer in a map
            FXAssetConstSP fx(FXAssetConstSP::dynamicCast(diffAssets[j]));
            perAssetData->fxAsset = fx;

            // FX is special 'cos it connects with mapIRGen
            //  get hold of the foreign yield curve
            YieldCurveConstSP forYC(fx->getRiskCcy());
            perAssetData->ccy = forYC->getCcy();
            QMCGenDiffusibleIRSP& perCcyData = mapIRGen[perAssetData->ccy]; // this extends mapIRGen too 
            if (!perCcyData.get()) { // build one if necessary  (alas, should always be)
                perCcyData = diffusibleIRGenSPFromModelType(pathConfig->getRatesModelType());
            }

            perAssetData->multiSpots = &spots;
            // introducing bidirectional link of ccy(i.e. ir) and fx 
            perCcyData->pQMCGenDiffusibleFX    = perAssetData.get();
            perAssetData->pForQMCGenDiffusibleIR = perCcyData.get();
            saveDiscountYC(pathConfig, forYC); // need to record foreign YC
            numFX++;
        } else {
            QMCGenDiffusibleEquitySP& perAssetData = mapEquityGen[assetName]; // extend the map
            if (!perAssetData.get()) // build one if necessary  (alas, should always be)
                perAssetData = QMCGenDiffusibleEquitySP(new SRMGenDiffusibleEquity(pathConfig->getEquityDiffusionStyle()));

            perAssetData->multiSpots = &spots;
            
            perAssetData->name = assetName; // useful once these are no longer in a map
            perAssetData->ccy = AssetUtil::assetCcy(diffAssets[j].get());
            perAssetData->eqAsset = diffAssets[j];
            numEq++;
        }
    }
    // update any xccy info needed for equity
    for (map<string, QMCGenDiffusibleEquitySP>::iterator aiter = mapEquityGen.begin();
         aiter != mapEquityGen.end(); ++aiter)
    {

        QMCGenDiffusibleIRSP pccy = mapIRGen[aiter->second->ccy];
        aiter->second->pQMCGenDiffusibleIR = pccy.get();

        if (aiter->second->ccy != domISOCode) {
             // should exist already
            if (!pccy->pQMCGenDiffusibleFX || !pccy->pQMCGenDiffusibleFX->fxAsset.get()) {
                throw ModelException(method,"Incomplete FX data!");
            }
        }
    }


    // Identify in multiFactors - this can tell us whether the instrument has explicit dates
    // at which a diffused quantity is to be accessed. Otherwise we can simply diffuse it
    // and not publish the path

    //Also: need to identify currency protection type for equity assets by examining the multifactors.
    //E.g. 'struck' or 'protected':
    for (int i = 0; i < multiFactors->numFactors(); i++){
        // get hold of the actual asset
        IMarketFactorConstSP factor = multiFactors->getFactor(i);
        const string& fname = factor->getName();
        if (mapFXGen.find(fname) != mapFXGen.end()) {
            // exists - update it
            mapFXGen[fname]->assetIdx = i;
        } else if (mapEquityGen.find(fname) != mapEquityGen.end()) {
            // exists - update it
            mapEquityGen[fname]->assetIdx = i;
        } else if (IYieldCurve::TYPE->isInstance(factor)){
            IYieldCurveConstSP yc(IYieldCurveConstSP::dynamicCast(factor));
            const string& ccy = yc->getCcy();
            if (ccy != domISOCode){
                throw ModelException(method, "Yield curves not in "
                                     "domestic currency not supported");
            }
        } else if (CAsset::TYPE->isInstance(factor)) {
            // The assetIdx allows fetching SV-required sim dates, so connect for composite assets too
            CAssetConstSP asset(CAssetConstSP::dynamicCast(factor));
            const string& truename = asset->getTrueName();
            if (mapFXGen.find(truename) != mapFXGen.end()) {
                mapFXGen[truename]->assetIdx = i;
            } else if (mapEquityGen.find(truename) != mapEquityGen.end()) {
                mapEquityGen[truename]->assetIdx = i;
            }
		} /*else if (EnergyFuturesCurve::TYPE->isInstance(factor)) {
            // For some reason we get CDSParSpreads among multiFactors, which are not handled above.
            // FIXME disable this check for now and clarify what we should do
            //throw ModelException(method, "Unsupported market factor of "
            //                     "type "+factor->getClass()->getName());
			EnergyFuturesCurveConstSP engAsset(EnergyFuturesCurveConstSP::dynamicCast(factor));
			const string& name = engAsset->getName();
			if (mapFXGen.find(name) != mapFXGen.end()) {
				mapFXGen[name].assetIdx = i;
			} else if (mapEquityGen.find(name) != mapEquityGen.end()) {
				mapEquityGen[name].assetIdx = i;
			}
		}*/ 
		/*else {
            throw ModelException(method, "Unsupported market factor of "
                                 "type "+factor->getClass()->getName());
        }*/
    }

 
    // two auxiliary objects to allow lookups by AssetIdx for EQ and FX
    for(map<string, QMCGenDiffusibleFXSP>::iterator aiter = mapFXGen.begin();
        aiter != mapFXGen.end(); ++aiter)
        {
            
            //The situation is that there are ccy protected equities,
            //and so the FX needs to be diffused (to get the quanto
            //adjustment) but there is no explicit mention of them in
            //the multiAsset. This latter means they have no index,
            //but there should be no need of such an index, since
            //there should be no SV in the instrument explicitly
            //referring to the FX. -- Stewart N.

            if (aiter->second->assetIdx != -1) 
                mFXbyId[aiter->second->assetIdx] = aiter->second;
        }

    for(map<string, QMCGenDiffusibleEquitySP>::iterator biter = mapEquityGen.begin();
             biter != mapEquityGen.end(); ++biter)
        {

            if (biter->second->assetIdx != -1)
                mEQbyId[biter->second->assetIdx] = biter->second;
            else
                throw ModelException(method, "EQ asset has assetIdx == -1 (i.e. was not found in MultiFactor)");
        }

}

// int SV::processAllSVGens(StateVariableCollectorSP    svCollector)
// {
//     IElemStateVariableGenArray stateVarGenArray(
//                                                 svCollector->getElemStateVarGens());
//     for(size_t i=0; i<stateVarGenArray.size(); ++i)
//         stateVarGenArray[i]->attachSVGen(this);
//     return stateVarGenArray.size();
// 
// }

void SV::processSVGen(const SVGenDiscFactor* df)
{
    string name = df->getYieldCurve()->getCcy();
    QLIB_VERIFY(mapIRGen.count(name), 
        "Attempted to create a state variable while the underlying asset ("+name+") is not initialized in the QMC engine");
    mapIRGen[name]->dfs.first.push_back(df);
}
void SV::processSVGen(const SVGenExpectedDiscFactor* edf)
{
    string name = edf->getYieldCurve()->getCcy();
    QLIB_VERIFY(mapIRGen.count(name), 
        "Attempted to create a state variable while the underlying asset ("+name+") is not initialized in the QMC engine");
    mapIRGen[name]->dfs.second.push_back(edf);
}
void SV::processSVGen(const SVGenSurvivalDiscFactor* sdf)
{
    string name = sdf->getCDSParSpreadCurve()->getName();
    QLIB_VERIFY(mapCreditGen.count(name), 
        "Attempted to create a state variable while the underlying asset ("+name+") is not initialized in the QMC engine");
    mapCreditGen[name]->sdf.first.push_back(sdf);
}
void SV::processSVGen(const SVGenExpectedSurvivalDiscFactor* esdf)
{
    string name = esdf->getName();
    QLIB_VERIFY(mapCreditGen.count(name), 
        "Attempted to create a state variable while the underlying asset ("+name+") is not initialized in the QMC engine");
    mapCreditGen[name]->sdf.second.push_back(esdf);
}

void SV::processSVGen(const SVGenAggregatedSurvivalDiscFactor* asdf)
{
    string name = asdf->getName();
    QLIB_VERIFY(mapCreditGen.count(name), 
        "Attempted to create a state variable while the underlying asset ("+name+") is not initialized in the QMC engine");
    mapCreditGen[name]->sdf.third.push_back(asdf);
}

void SV::processSVGen(const SVGenExpectedSpot* es)
{
    // complication: need to find out whether this is equity of fx:
    // the "assetIdx" hack existed for some time, I have merely moved it around.
    // the old comment said that lookup by name sometimes failed - TODO: see what's wrong

    int assetIdx = es->getAssetIdx();
    if (mFXbyId.count(assetIdx))
        mFXbyId[assetIdx]->fxExpSpots.push_back(es);
    else if (mEQbyId.count(assetIdx))
        mEQbyId[assetIdx]->expSpots.push_back(es);
    else
        throw ModelException("SV::processSVGen(const SVGenExpectedSpot*)", 
            "Missing data - no asset with id "+Format::toString(assetIdx)+" defined!");
}
void SV::processSVGen(const SVGenExpectedEnergyFuture* eEnrgFutureGen)
{
	string energyName = eEnrgFutureGen->getName();
    QLIB_VERIFY(mapEnergyGen.count(energyName), 
        "Attempted to create a state variable while the underlying asset ("+energyName+") is not initialized in the QMC engine");
	mapEnergyGen[energyName]->expFpGens.push_back(eEnrgFutureGen);
}
void SV::processSVGen(const SVGenExpectedBasisFwdSpread* s)
{
    string basisName = s->getName();
    QLIB_VERIFY(mapBasisGen.count(basisName), 
        "Attempted to create a state variable while the underlying asset ("+basisName+") is not initialized in the QMC engine");
    mapBasisGen[basisName]->expSPGens.push_back(s);
}
void SV::processSVGen(const SVGenDateOfDefault* dod)
{
    string name = dod->getCDSCurveName();
    QLIB_VERIFY(mapCreditGen.count(name), 
        "Attempted to create a state variable while the underlying asset ("+name+") is not initialized in the QMC engine");
    mapCreditGen[name]->sdf.dod.push_back(dod);
}
void SV::processSVGen(const SVGenSpot* s)
{
    s->getSimSeries(); // TODO: see what this function actually does
}

void SV::processSVGen(const SVGenPathWeight* s)
{
    pathWeightSVGen.push_back(s);
}



void SV::processSVGen(const IElemStateVariableGen* s)
{
    // this is a fallback -- if SVCollector collected something we do not know how to process

    throw ModelException("processSVGen(const IElemStateVariableGen* )",
       string( "unable to process request to add unidentified SV Generators to the simulation. Generator's type is " )+ typeid(s).name());
}

int SV::processAllSVGens(StateVariableCollectorSP    svCollector)
{
    IElemStateVariableGenArray stateVarGenArray(
            svCollector->getElemStateVarGens());
    for(size_t i=0; i<stateVarGenArray.size(); ++i)
        stateVarGenArray[i]->attachSVGen(this);
    return stateVarGenArray.size();
}
