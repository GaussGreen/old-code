
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAGMWB2.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 Guaranteed Minimum Withdrawal Benefit (GMWB)
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/InsuranceAGMWB2.hpp"
#include "edginc/VolRequestTime.hpp"

DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////

void InsuranceAGMWB2::Validate() {
    static const string routine("InsuranceAGMWB2::Validate");

    InsuranceAGMWB::Validate();


    //lapseRateMultiplierType not used here
    
    // checks that the inputs are consistent with step up feature
    if (isStepUp){
        for (int i = 0; i < stepUpITMLevels.size(); i++){
            if (stepUpITMLevels[i] < 1.0)
                throw ModelException(routine, "ITM level for step up has to be >= 1.0.");        
        }
    }

    int nberFeatureGp = featureGpSizeAtStart.size();

    //stepUpITMLevels; WDITMLevels;  WDTime; WDPercent;  have to be the same size as featureGpSizeAtStart
    if ( ( nberFeatureGp!=stepUpITMLevels.size()) ||
            ( nberFeatureGp!=wDITMLevels.size()) ||
            ( nberFeatureGp!=wDTime.size()) ||
            ( nberFeatureGp!=wDPercent.size()) ){
            throw ModelException(routine,
                                    "stepUpITMLevels'size ( = " + Format::toString(stepUpITMLevels.size()) +
                                    "), WDITMLevels' size ( = "  + Format::toString(wDITMLevels.size())     +
                                    "), WDTime' size( = "        + Format::toString(wDTime.size())          +
                                    "), WDPercent's size ( = "   + Format::toString(wDPercent.size())       +
                                    "), featureGpSizeAtStart's size (= " + Format::toString(nberFeatureGp)+
                                    " have to be the same size. \n" );
    }
  
    //WDTime has to be increasing
//    for (int j = 0; j < wDTime.size() - 1; j++){
//        if ( wDTime[j+1] < wDTime[j]  ){
//            throw ModelException(routine, "withdrawal time shouldn't be an decreasing function!");
//        }                
//    }

    //WD % >= 0
    for (int j = 0; j < wDPercent.size() ; j++){
        if (Maths::isNegative(wDPercent[j])){
            throw ModelException(routine, "withdrawal % shouldn't be negative!");
        }
        if (wDPercent[j] > maxGMWBPercentage){
            throw ModelException(routine, "Amount withdrawn % ( = " + Format::toString(wDPercent[j]) +
                "), for each group can't be bigger than the maxGMWBPercentage ( = " + Format::toString(maxGMWBPercentage) +
                ")!");
        }
    }
}
   
/** Invoked when Class is 'loaded' */
void InsuranceAGMWB2::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(InsuranceAGMWB2, clazz);
    SUPERCLASS(InsuranceAGMWB);
    IMPLEMENTS(IMCIntoProduct);
    IMPLEMENTS(LastSensDate);   
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultInsuranceAGMWB2);

    // Guaranteed Minimum Withdrawal Benefit feature

    //step up
    FIELD(stepUpITMLevels, "Step up ITM levels for different feature group");
    FIELD_MAKE_OPTIONAL(stepUpITMLevels);

    //withdraw
    // Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature 

    FIELD(wDTime, "Withdrawal time for different feature group");
    FIELD_MAKE_OPTIONAL(wDTime);

    FIELD(wDPercent, "Withdrawal % for different feature group");
    FIELD_MAKE_OPTIONAL(wDPercent);

    FIELD(maturity, "maturity of the traded product");
    FIELD_MAKE_OPTIONAL(maturity);

    FIELD(longVol, "use flat longvol beyond maturity");
    FIELD_MAKE_OPTIONAL(longVol);
}


///////////// product class ////////////////////////

////////////////////////////////////////////////////////////////////////////////
//state variable
//////////////////////////////////////////////////////////////////////////
///////////// product class ////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/* MC product class for insurance annuity */

/** equivalent to InstIntoMCProduct. */
InsuranceAGMWB2ProdSV::InsuranceAGMWB2ProdSV(const InsuranceAGMWB2* instGMWB2,
                        const SimSeriesSP& simSeries) :
                InsuranceAGMWBProdSV(instGMWB2, simSeries),  instGMWB2(instGMWB2){  //??????

//to review!!!!!!!!
    //overide the yearFrac defined in the parent class by 1/12
    //since GMWB2 is always monthly 
//    VolRequestTimeSP volTimeReq(new VolRequestTime()); // vol request for time
//  CVolProcessedSP volTimeProcessed(getMultiFactors()->factorGetProcessedVol(0, volTimeReq.get())); 
//  TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric(); // time metric    
    for (int iDate = 1 ; iDate < instGMWB2->monitorDates.size() ; iDate++){
        yearFrac[iDate] = 1.0 / 12.0;
    }

    int nberAgeGp = instGMWB2->ageGpSizeAtStart.size();        
    int nberFeatureGp = instGMWB2->featureGpSizeAtStart.size();

    //overide the initial value for lastMaxGMWBWithdrawal...
    for (int iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
        for (int iWD = 0 ; iWD < nberFeatureGp ; iWD++){
            lastMaxGMWBWithdrawal[iWD][iGroup]=  instGMWB2->maxGMWBPercentage*yearFrac[1];  //account value = 1.0 initially, to review
            lastMaxGMWBWithdrawalSoFar[iWD][iGroup]= instGMWB2->maxGMWBPercentage*yearFrac[1];           
        }
    }

    //may start withdraw later
    isWithdraw.resize(nberFeatureGp);
    isWithdrawSoFar.resize(nberFeatureGp);

    for (int iWD = 0 ; iWD < nberFeatureGp ; iWD++){
        isWithdraw[iWD].resize(nberAgeGp,false);
        isWithdrawSoFar[iWD].resize(nberAgeGp,false);
    }
}

// decide if the client steps up or not and update benefit
void InsuranceAGMWB2ProdSV::updateStepUp(const int& timeFromStart, const int& iWDGroup){
    int iGroup;
    int nberAgeGp = instGMWB2->ageGpSizeAtStart.size();

    // isWaitPeriod flag was initialized to be true
    if (isWaitPeriod[iWDGroup] && (timeFromStart >= instGMWB2->initialWaitPeriod)) {
        isWaitPeriod[iWDGroup] = false;

        for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
            // authorize step-up from the end of the initial waiting period
            lastStepUp[iWDGroup][iGroup] = instGMWB2->stepUpPeriod;
        }
    }

    // step-ups only supported by GMWB benefit
    if (instGMWB2->isStepUp) {
        for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
            if ((!isWaitPeriod[iWDGroup] 
                && (lastStepUp[iWDGroup][iGroup] >= instGMWB2->stepUpPeriod))
                && (!Maths::equals(0.0, benefit[iWDGroup][iGroup]))
                && (accountValue[iWDGroup][iGroup] > instGMWB2->stepUpITMLevels[iWDGroup] * benefit[iWDGroup][iGroup]) 
                ){
                benefit[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup];
                //yearFrac from timeFromStart -1 to timeFromStart
                // convert annual wd rate  to this period wd rate
                double maxWDPerc = instGMWB2->maxGMWBPercentage 
                                                * yearFrac[timeFromStart];

                lastMaxGMWBWithdrawal[iWDGroup][iGroup] =
                    max(lastMaxGMWBWithdrawal[iWDGroup][iGroup], 
                        maxWDPerc * benefit[iWDGroup][iGroup]);

                //temporary, to remove, just for check purpose, should never enter
                if (benefit[iWDGroup][iGroup] != accountValue[iWDGroup][iGroup]){
                    throw ModelException( "InsuranceAGMWB2Prod::updateStepUp, account value diff from benefit");        
                }
                // set lastStepUp[iGroup] to -1 since lastStepUp[iGroup]
                // is going to be incremented at the end of the loop
                // possible step up >= stepUpPeriod
                //lastStepUp[iWDGroup][iGroup] = -1;
                lastStepUp[iWDGroup][iGroup] = 0;
            }


//            if (!isWaitPeriod[iWDGroup] && (lastStepUp[iWDGroup][iGroup] >= instGMWB2->stepUpPeriod)){/// to review
//                if (!Maths::equals(0.0, benefit[iWDGroup][iGroup])){                                                        
//                    // if step ups are static, compare to stepUpITMLevel
//                    if (accountValue[iWDGroup][iGroup] > (instGMWB2->stepUpITMLevels[iWDGroup] 
//                                                        * benefit[iWDGroup][iGroup]) ){
//                        benefit[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup];
//                        //yearFrac from timeFromStart -1 to timeFromStart
//                        // convert annual wd rate  to this period wd rate
//                        double maxWDPerc = instGMWB2->maxGMWBPercentage 
//                                            * yearFrac[timeFromStart];
//
//                        lastMaxGMWBWithdrawal[iWDGroup][iGroup] =
//                            max(lastMaxGMWBWithdrawal[iWDGroup][iGroup], 
//                                maxWDPerc * accountValue[iWDGroup][iGroup]);
//                        // set lastStepUp[iGroup] to -1 since lastStepUp[iGroup]
//                        // is going to be incremented at the end of the loop
//                        //lastStepUp[iWDGroup][iGroup] = 0;  //to review
//                        lastStepUp[iWDGroup][iGroup] = -1;
//                    }
//                }
//            }
        }
    }

    for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
        ++lastStepUp[iWDGroup][iGroup];
    }
}

// update number of clients by taking into account the mortality rate
void InsuranceAGMWB2ProdSV::updateDeadClients(const int& timeFromStart, const int& iWDGroup) {

    if (instGMWB2->dieOption && (timeFromStart != 0)){
        for (int iGroup = 0 ; iGroup < instGMWB2->ageGpSizeAtStart.size() ; iGroup++){   
            // current age of the group
            //timeFromStart is counted as monthly, to review ????
            int currentAge = (timeFromStart-1) / 12.0;

            //to review????, below only 
            // below only use 1 mortality rate for 0 to 11 months
            //int currentAge = instGMWB2->monitorDates[0].yearFrac(instGMWB2->monitorDates[timeFromStart]); 

            currentAge += instGMWB2->ageAtStart[iGroup] ;
            
            // mortality rate depending on the current age of the group
            //yearFrac from timeFromStart-1 to timeFromStart
            double mRForThisPeriod = yearFrac[timeFromStart]
                                    * instGMWB2->mortalityRate[currentAge];
//            double mRForThisPeriod = yearFrac[timeFromStart]
//                                    * instGMWB->mortalityRate[currentAge-1];

            groupSize[iWDGroup][iGroup] *= (1.0 - mRForThisPeriod);
        }
    }
}

// update number of clients by taking into account the lapse rate
void InsuranceAGMWB2ProdSV::updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup) {

    if (instGMWB2->lapseOption && (timeFromStart != 0)) {// timeFromStart = 0 at start date

        for (int iGroup = 0 ; iGroup < instGMWB2->ageGpSizeAtStart.size(); iGroup++){
            double ITMAv_ben = 0.0;
            if (Maths::isPositive(benefit[iWDGroup][iGroup])){
                ITMAv_ben = accountValue[iWDGroup][iGroup] 
                        / benefit[iWDGroup][iGroup] /instGMWB2->adjustmentAV;
            }
            double lapseRateMultiplier = Maths::min(1.0, ITMAv_ben);
            // base lapse rate (the same for all clients)
            // base lapse rate depending on the time from start date only

            //timeFromStart is counted as monthly, to review ????
            int timeFromStartInY = (timeFromStart==0)? 0 : (timeFromStart-1) / 12.0; 

            //to review?????, should use 1 mortality rate for 0 to 12 month, 
            // below only use 1 lapse rate for 0 to 11 months
            //int timeFromStartInY = instGMWB2->monitorDates[0].yearFrac(instGMWB2->monitorDates[timeFromStart]); 

            //double lapseRate = instGMWB2->lapseRate[timeFromStart-1]; 
            double lapseRate = instGMWB2->lapseRate[timeFromStartInY]; 

            // yearFrac from timeFromStart-1 to timeFromStart
            // convert annual lapse Rate to this period
            lapseRate *= yearFrac[timeFromStart];

            double floorLapseRate = 0.0;
            lapseRate = max(lapseRateMultiplier * lapseRate, floorLapseRate);
            groupSize[iWDGroup][iGroup] *= (1.0 - lapseRate);
        }
    }
}

// compute coupons taken by clients
void InsuranceAGMWB2ProdSV::computeWithdrawals(const int& timeFromStart,                                                 
                                                DoubleArraySP withdrawals,
                                                const int& iWDGroup
                                                 ) {        

    for (int iGroup = 0 ; iGroup < instGMWB2->ageGpSizeAtStart.size();  iGroup++){       

        //isWithdraw is flag indicating if the withdrawal takes place or not
        //isWithdraw  is initialized to false
        //timeToWithdrawal depends on ITMLevel, if > stepUpITMLevels, it may start WD earlier
        if (!isWithdraw[iWDGroup][iGroup]){//haven't start to withdraw
            decideWithdrawal(timeFromStart, withdrawals, iGroup, iWDGroup);                
        }                                            

        // for GMWB, withdrawals expressed as a percentage of the initial amount
        double isWD = (isWithdraw[iWDGroup][iGroup] ? 1.0 : 0.0);
        if (timeFromStart == 0) {
            (*withdrawals)[iGroup] = 0.0;
        }else{
            double wDAmt = ((instGMWB2->maxGMWBPercentage ==0.0) ? 0.0 : 
                Maths::min(1.0, instGMWB2->wDPercent[iWDGroup]/instGMWB2->maxGMWBPercentage) );

            wDAmt = wDAmt*lastMaxGMWBWithdrawal[iWDGroup][iGroup];
            (*withdrawals)[iGroup] = isWD * max(min(wDAmt, benefit[iWDGroup][iGroup]), 0.0) ;
        }        
    }
}

void InsuranceAGMWB2ProdSV::decideWithdrawal(const int& timeFromStart, 
                                                DoubleArraySP withdrawals,
                                                const int& iGroup,
                                                const int& iWDGroup 
                                               ){

    double ITMLevel;

    //make "jump" decision based on ITM =AV/Ben
    //if Ben <=0, the trade stop
    if (Maths::isPositive(benefit[iWDGroup][iGroup])){
        ITMLevel = accountValue[iWDGroup][iGroup] / benefit[iWDGroup][iGroup];
        if (ITMLevel < instGMWB2->wDITMLevels[iWDGroup]){ // start WD
            isWithdraw[iWDGroup][iGroup] = true;                          
        }else{
            //check WDTime
            if (timeFromStart >= instGMWB2->wDTime[iWDGroup]){
                //start WD
                isWithdraw[iWDGroup][iGroup] = true;                          
            }
        }
    }
}

void InsuranceAGMWB2ProdSV::calcPayoff( const int iWDGroup, const int db_pathCount){
    static const string routine("InsuranceAGMWB2Prod::calcPayoff");
    try {

        //take the first asset
        //just to get beginIdx and endIdx
        const SVPath& pathAsset1 = spotSV->path(0); //
        bool  doingPast = pathAsset1.doingPast();
        int beginIdx = pathAsset1.begin();
        int endIdx = pathAsset1.end();

        int iGroup;
        int nberAgeGp = instGMWB2->ageGpSizeAtStart.size();        

        DateTime settlementDate; // settlement date for the current start date anniversary
        DateTime valueDate = instGMWB2->valueDate; // value date

        double indexValue; // index value (equal to 1$ at start date)
        DoubleArraySP currentWithdrawals = DoubleArraySP(new DoubleArray(nberAgeGp)); 

        for (int Idx = beginIdx ; Idx < endIdx ; Idx++){
            double JPMFees_Idx = 0.0; // for debug purpose

            //for debug output
            // only output one stock path for now, to review!!
            //Debug_outputBegin(pathGen->doingPast(), instGMWB2->db_ithAsset, iWDGroup, db_pathCount, Idx);
            Debug_outputBegin(doingPast, instGMWB2->db_ithAsset, iWDGroup, db_pathCount, Idx);

            // update index value
            double baskV0 = computeBasket(0);
            double baskVIdx = computeBasket(Idx);
            double baskVIdxMinus1 = (Idx== 0) ? 1.0 : computeBasket(Idx-1);

            indexValue = 1.0 * (baskVIdx / baskV0);
            //indexValue = 1.0 * (path[Idx] / path[0]);

            // update settlement date
            settlementDate = instGMWB2->instSettle->settles(instGMWB2->monitorDates[Idx], 0);
            bool alreadySettled = valueDate.isGreater(settlementDate);

            // update account value (before coupon is paid)
            for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
                if (Idx != 0) {// start date
//                    accountValue[iWDGroup][iGroup] *= 
//                        Maths::max(path[Idx] / path[Idx-1]- instGMWB2->feeRateIns, 0.0);                    
                    //to review
                    // yearFrac from timeFromStart-1 to timeFromStart
                    // convert fees Rate to this period
                    double feesIdx = yearFrac[Idx] * instGMWB2->feeRateIns;

                    accountValue[iWDGroup][iGroup] *= 
                        (Maths::max(baskVIdx / baskVIdxMinus1 - feesIdx, 0.0));                    
                }
            }

            // compute fees (supposed to be paid at the begining of the period)
            if (!alreadySettled && (Idx !=0)){
                JPMFees_Idx = computeFees(&feeValue, Idx, iWDGroup);
            }

            //update step up for GMWB
            updateStepUp(Idx, iWDGroup);

            // compute coupons taken by clients
            // calculated as percentage of the initial amount for GMWB
            computeWithdrawals(Idx, currentWithdrawals, iWDGroup);

            // update account value
            updateAccountValue(currentWithdrawals, alreadySettled, iWDGroup);

            // update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)
            updatebenefit(currentWithdrawals, iWDGroup);

            // remove the clients which have lapsed during the year
            // additional lapse rate Multiplier may depend on the current status of withdrawal
            updateLapseClients(Idx, currentWithdrawals, iWDGroup); 

            // remove the clients which are dead during the year
            // (we assume that the clients die at the end of the year)
            updateDeadClients(Idx, iWDGroup);

            if (!alreadySettled && (Idx != 0)){
                // add value of Guaranteed Minimum Death Benefit (GMDB) option
                for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
                    if (Maths::isNegative(accountValue[iWDGroup][iGroup])){
                        double df = dfSV->path()[Idx];
                        pathValue +=  df *  // the fees are forward valued to last date
                            groupSize[iWDGroup][iGroup] * (-accountValue[iWDGroup][iGroup]);


                        //for debug output
                        Debug_outputMid(iWDGroup, iGroup, db_pathCount, Idx);
                        accountValue[iWDGroup][iGroup] = 0.0;
                    }
//                    // if the GMWB benefit is inferior to 0.0, we set the account value at 0.0
//                    if (!Maths::isPositive(benefit[iWDGroup][iGroup]))
//                        accountValue[iWDGroup][iGroup] = 0.0;
                }
            }

            //for debug output
            // only output one stock path for now, to review!!
            Debug_outputEnd(doingPast, currentWithdrawals, JPMFees_Idx, instGMWB2->db_ithAsset, iWDGroup, db_pathCount, Idx);
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}
//////////////////////////////////////////////////////////////////////////
//end state variable
////////////////////////////////////////

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* InsuranceAGMWB2::createProduct(const MonteCarlo* model) const {

    static const string routine = "InsuranceAGMWB::createProduct";
    try {
        if(model->stateVarUsed()) {
            // we need to create a SimSeries object which says which assets need
            // which dates to be simulated
            SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                        one */
            simSeries->addDates(monitorDates);

            // State variables
            return new InsuranceAGMWB2ProdSV(this, simSeries);
        }else{
            throw ModelException(routine, "only useStateVars = True is supported for GMWB products.");
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

CClassConstSP const InsuranceAGMWB2::TYPE = CClass::registerClassLoadMethod(
    "InsuranceAGMWB2", typeid(InsuranceAGMWB2), InsuranceAGMWB2::load);

// force linker to include this file (avoid having header file) */
extern bool InsuranceAGMWB2Load()
{
    return true && InsuranceAGMWB2::TYPE;
}

DRLIB_END_NAMESPACE
