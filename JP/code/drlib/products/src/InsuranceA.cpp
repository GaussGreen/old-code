
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceA.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 common features for GMDB/GMAB/GMWB
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/Random.hpp"

#include "edginc/InsuranceA.hpp"
#include "edginc/VolRequestTime.hpp"


DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////
// returns the array of all sampling dates the instrument will ever need
// excluding (possibly) the ref level dates (handled in GenericNFBase)
const DateTimeArray InsuranceA::samplingDates() const {
    return monitorDates;
}

void InsuranceA::Validate() {
    static const string routine("InsuranceA::Validate");


    // checks that there is at least one age group
    if (ageAtStart.size() <= 0)
        throw ModelException(routine, "ageAtStart array size must be > 0.");

    // checks that the number of ages and the number of age percentages are consistent
    if (ageAtStart.size() != ageGpSizeAtStart.size())
        throw ModelException(routine, "ageAtStart and ageGpSizeAtStart arrays must be the same size.");

    // computes the sum of age percentages
    double sum = 0.0;
    int iGroup;
    for (iGroup = 0 ; iGroup < ageGpSizeAtStart.size() ; iGroup++)
    {
        sum += ageGpSizeAtStart[iGroup];
    }

    // checks that the sum of age percentages is 100%
    if (!Maths::equals(1.0, sum))
        throw ModelException(routine, "The sum of the age percentages has to be 100%.");

    //check the sum of featureGpSizeAtStart is 100%
    sum = 0.0;
    for (iGroup = 0; iGroup < featureGpSizeAtStart.size(); iGroup++){
        sum += featureGpSizeAtStart[iGroup];
    }

    if (!Maths::equals(1.0, sum)){
        throw ModelException(routine, "The sum of the different feature group percentages has to be 100%.");
    }

    if ((lapseRateMultiplierType != 1) && (lapseRateMultiplierType != 2)){
        throw ModelException(routine, "Wrong lapseRateMultiplierType! should be 1 or 2 ");
    }

    //should provide at least 1 monitor date
    if (monitorDates.size()==0){
        throw ModelException(routine, "At lease 1 monitor Date should be provided!");
    }

    // these tests have to be here due to monitor dates for this structure is monthly, to review!!!

    // computes the maximum age among the groups
    int maxAge = 0;
    int maxAgeIndex = 0;
    for (int iGroup = 0 ; iGroup < ageGpSizeAtStart.size() ; iGroup++){
        if (ageAtStart[iGroup] > maxAge){
            maxAge = max(maxAge, ageAtStart[iGroup]);
            maxAgeIndex = iGroup;
        }
    }

    //monitor date can be monthly or yearly
    //GMWB1 is for yearly, and GMWB2 is monthly, but, GMWB2 should also work for yearly
    //to review
    int dealTermInY = monitorDates[0].yearFrac(monitorDates[monitorDates.size()-1]) + 1; 

    // checks that the mortality rates are provided up to the last simulation date
    //if (dieOption && (mortalityRate.size() < maxAge + monitorDates.size()))
    if (dieOption && (mortalityRate.size() < maxAge + dealTermInY))
        throw ModelException(routine, "the mortality rates for the last simulation dates are missing.");

    // checks that the lapse rates are provided up to the last simulation date
    //if (lapseOption && (lapseRate.size() < monitorDates.size()))            
    if (lapseOption && (lapseRate.size() < dealTermInY))            
        throw ModelException(routine, "the lapse rates for the last simulation dates are missing.");

    // check that we've got as many weights as there are assets 
    // and if % weights that they sum to 100%
    AssetUtil::checkWeights(weights, assets->NbAssets());

    // validate dates are not empty - order is handled by SimSeries
    if (monitorDates.empty()) {
        throw ModelException(routine, "No monitoring Dates supplied!");
    }		
    
    //to review, do we really need this check?
    // checks the first sample date is on or before value date

    //if (spotSamples->getFirstDate() > valueDate)
//    if (monitorDates[0] > valueDate)            
//        throw ModelException(routine, "there must be one sample date in the past (or = valueDate).");

}

/** Invoked when Class is 'loaded' */
void InsuranceA::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet

    REGISTER(InsuranceA, clazz);
    SUPERCLASS(GenericNFBase);
    IMPLEMENTS(IMCIntoProduct);
    IMPLEMENTS(LastSensDate);   
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultInsuranceA);
    FIELD(weights,          "Weights");
    FIELD(monitorDates, "sample dates ");
    FIELD(feeRateJPM, "fee rate as % of the account value, charged by JPM");
    FIELD(feeRateIns, "fee rate as % of the account value, charged by Insurance company");
    FIELD_MAKE_OPTIONAL(feeRateIns);

    FIELD(ageAtStart, "age of the people in the age group at start date");
    FIELD(ageGpSizeAtStart, "the percentage of people in the corresponding age group");

    //optional inputs    
    FIELD(featureGpSizeAtStart, "the % of people in the group with different feature.");
    FIELD_MAKE_OPTIONAL(featureGpSizeAtStart);

    //mortality
    FIELD(dieOption, "if yes, the mortality rates table is taken into account");
    FIELD(mortalityRate, "mortality rate as function of the client age");
    //lapse
    FIELD(lapseOption, "if yes, the lapse rate table is taken into account");
    FIELD(lapseRate, "contract lapse rate as function of the time from start date");
    
    //optional inputs
    FIELD(isDynamicLapseRate, "if yes, the lapse rate is dynamic (i.e.determined with respect to ITM)");
    FIELD_MAKE_OPTIONAL(isDynamicLapseRate);
    FIELD(lapseRateMultiplierType, "1 power ft; 2: AV/Benefits, both are active only if Benefit/AV >=1");
    FIELD_MAKE_OPTIONAL(lapseRateMultiplierType);
    FIELD(withdrawMultiplier, "additional multiplier for lapse rate: depends on the current status of withdrawal.");
    FIELD_MAKE_OPTIONAL(withdrawMultiplier);
    FIELD(noWithdrawMultiplier, "additional multiplier for lapse rate: depends on the current status of withdrawal.");
    FIELD_MAKE_OPTIONAL(noWithdrawMultiplier);
}

//////////////////////////////////////////////////////////////                                                                                                                                                                                                                                                                                                                                                                                                                        
//////// product class //////////

/* MC product class for insurance annuity -- state variable version*/

//////////////////////////////////////////////////////////////////////////////////////////


/** Override default method on IMCProduct. This method is called every time
    the path generator is changed (which is, at the moment, when the
    past path generator is created, and then when the future path
    generator is created  */
void InsuranceAProdSV::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
    static const string routine = "InsuranceAProdSV::pathGenUpdated";

    try {
        spotSV = spotGen->getSpotSV(newPathGen);
        refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
        matDfSV = matDfGen->getSVDiscFactor(newPathGen);
        dfSV = dfGen->getSVDiscFactor(newPathGen);
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** equivalent to InstIntoMCProduct. Need to call parent's constructor */
InsuranceAProdSV::InsuranceAProdSV(const InsuranceA* inst,
                                 const SimSeriesSP& simSeries) :
        MCProductClient(inst->assets.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        inst->refLevel,
                        simSeries,
                        inst->pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())),
        matDfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
        inst->instSettle, simSeries->getLastDate())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
        inst->instSettle, simSeries->getAllDates() /*inst->monitorDates, */)){  // to check????
        
    // create and initialize uniform random generator
    uniRand = RandUniformDefaultSP(new RandUniformDefault());

    // compute year fraction from Idx-1 to Idx
    //used to convert the annual rates such as mortality rate, fees...
	yearFrac.resize(inst->monitorDates.size(), 0.0);

    VolRequestTimeSP volTimeReq(new VolRequestTime()); // vol request for time
	CVolProcessedSP volTimeProcessed(getMultiFactors()->factorGetProcessedVol(0, volTimeReq.get())); 
	TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric(); // time metric
    
    int iDate;    
    for (iDate = 1 ; iDate < inst->monitorDates.size() ; iDate++){
        yearFrac[iDate] = timeMetric->yearFrac(inst->monitorDates[iDate-1], 
                            inst->monitorDates[iDate]);

//      if (inst->valueDate <= inst->monitorDates[iDate])
//            yearFrac[iDate] = timeMetric->yearFrac(inst->monitorDates[iDate-1], 
//                                inst->monitorDates[iDate]);
//        else{
//            yearFrac[iDate] = 999999.9; // this number should never be used
//        }

//        yearFrac[iDate] = 1.0;
	}

    // initialize path value
    pathValue = 0.0;
    pathValueSoFar = 0.0;

    // initialize present value of the fees
    feeValue = 0.0;
    feeValueSoFar = 0.0;

    int nberAgeGp = inst->ageGpSizeAtStart.size();        
    int nberFeatureGp = inst->featureGpSizeAtStart.size();

    int iGroup;
    int iWD;

    // resize and initialize account value
    accountValue.resize(nberFeatureGp);
    accountValueSoFar.resize(nberFeatureGp);
    for (iWD = 0 ; iWD < nberFeatureGp ; iWD++){
        accountValueSoFar[iWD].resize(nberAgeGp, 1.0);
        accountValue[iWD].resize(nberAgeGp,1.0);
    }
    
    groupSizeSoFar.resize(nberFeatureGp);
    groupSize.resize(nberFeatureGp);

    for (iWD = 0 ; iWD < nberFeatureGp ; iWD++){
        groupSize[iWD].resize(nberAgeGp);
        groupSizeSoFar[iWD].resize(nberAgeGp);
    }
    
    for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
        for (iWD = 0 ; iWD < nberFeatureGp ; iWD++){
            groupSize[iWD][iGroup] = inst->ageGpSizeAtStart[iGroup] 
                                    * inst->featureGpSizeAtStart[iWD]; //featureGpSizeAtStart[0] = 1 except deterministicWDType =2 or 3

            groupSizeSoFar[iWD][iGroup] = groupSize[iWD][iGroup];
        }
    }            

    benefit.resize(nberFeatureGp);
    benefitSoFar.resize(nberFeatureGp);

    for (iWD = 0 ; iWD < nberFeatureGp ; iWD++){
        benefit[iWD].resize(nberAgeGp, 1.0);
        benefitSoFar[iWD].resize(nberAgeGp, 1.0);
    }
}

//double InsuranceAProdSV::computeBasket(const IPathGenerator*  pathGen, const int iStep) {
double InsuranceAProdSV::computeBasket( const int iStep) {

     //const SVGenSpot::Path& pathAsset = assetSV->path(0);

    double baskValue = 0.0;

    for (int iAsset=0; iAsset< inst->assets->NbAssets(); iAsset++) {
        const SVPath& path = spotSV->path(iAsset);

        //baskValue += inst->weights[iAsset] * pathGen->Path(iAsset, 0/*iPath*/)[iStep];
        baskValue += inst->weights[iAsset] * path[iStep];
    }
    return baskValue;
}

/** compute fees (supposed to be paid at the begining of the period) */
/** return this period fees */
double InsuranceAProdSV::computeFees(double* fees, const int& timeFromStart, const int& iWDGroup){
    return 0.0;
}

/** update benefits */
void InsuranceAProdSV::updatebenefit(DoubleArraySP withdrawals, const int& iWDGroup){

}

// decide if the client steps up or not and update benefit
void InsuranceAProdSV::updateStepUp(const int& timeFromStart, const int& iWDGroup){
}


// update number of clients by taking into account the lapse rate
void InsuranceAProdSV::updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup) {
}

// update number of clients by taking into account the mortality rate
void InsuranceAProdSV::updateDeadClients(const int& timeFromStart, const int& iWDGroup) {
}


// compute coupons taken by clients
void InsuranceAProdSV::computeWithdrawals(const int& timeFromStart,                                                 
                                                DoubleArraySP withdrawals,
                                                const int& iWDGroup
                                                 ) {
        
}


// update account value
void InsuranceAProdSV::updateAccountValue(DoubleArraySP withdrawals, 
                                                 bool alreadySettled,
                                                 const int& iWDGroup) {
        
}

void InsuranceAProdSV::payoff(const IPathGenerator*  pathGen, IMCPrices& prices) {
}

CVolRequestLNArray InsuranceAProdSV::getVolInterp(const IMCPathGenerator* pathGen,
                                int                     iAsset) const {

    static const string method = "InsuranceAProd::getVolInterp";

    try
    {
        // shouldn't use pathGen, since it's set to 0 in the state variable version
        //double interpLevel = pathGen->refLevel(iAsset, 0);

        // one interp level/path per asset here        
        CVolRequestLNArray reqarr(1);

        const IRefLevel* refLevel = getRefLevel();
        const DateTime&  startDate = refLevel->getAllDates().front();     
        const DateTime&  today = getToday();        
        bool             fwdStarting = startDate.isGreater(today);
        const DateTime&  lastSimDate = getSimSeries()->getLastDate();

        double interpLevel = refLevelSV->refLevel(iAsset);

        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,   // to review, at refLevel
                                                                    startDate,
                                                                    lastSimDate,
                                                                    fwdStarting));
        return reqarr;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }	
}


// control variate is done here
void InsuranceAProdSV::recordExtraOutput(Control* control, Results* results, const IMCPrices& prices) const {

    const InsuranceAPrices& myprices = static_cast<const InsuranceAPrices&>(prices);

//mv pv to payoff
    // pv from last date to value date
//    double discountFactor = inst->instSettle->pv(inst->valueDate,  
//                 inst->monitorDates[inst->monitorDates.size()-1],
//                 inst->discount.get(),
//                 0/*inst->asset.get()*/);

    if (control && control->isPricing()) {
        OutputRequest* request =
        control->requestsOutput(OutputRequest::OPTION_PRICE);
        if (request) {
            double optionPrice = myprices.getOptionPrice();
            // need to present value here
            //optionPrice *= discountFactor;
            results->storeRequestResult(request, optionPrice);
        }
        request = control->requestsOutput(OutputRequest::FEE_PRICE);
        if (request) {
            double feePrice = myprices.getFeePrice();
            // need to present value here
            //feePrice *= discountFactor;
            results->storeRequestResult(request, feePrice);
        }
    }
}


//////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface - the
    implementation of this is below */

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* InsuranceA::createProduct(const MonteCarlo* model) const {
    static const string routine = "InsuranceA::createProduct";
    try {
        if(model->stateVarUsed()) {
            // we need to create a SimSeries object which says which assets need
            // which dates to be simulated
            SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                            one */
            simSeries->addDates(monitorDates);

            // State variables
            return new InsuranceAProdSV(this, simSeries);
        }else{
            throw ModelException(routine, "only useStateVars = True is supported for GMWB products.");
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

CClassConstSP const InsuranceA::TYPE = CClass::registerClassLoadMethod(
    "InsuranceA", typeid(InsuranceA), InsuranceA::load);

// force linker to include this file (avoid having header file) */
extern bool InsuranceALoad()
{
    return true && InsuranceA::TYPE;
}

DRLIB_END_NAMESPACE
