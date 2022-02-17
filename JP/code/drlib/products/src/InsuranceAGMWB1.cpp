
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAGMWB1.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 Guaranteed Minimum Withdrawal Benefit (GMWB)
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/InsuranceAGMWB1.hpp"
#include "edginc/VolRequestTime.hpp"


DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////

void InsuranceAGMWB1::Validate() {
    static const string routine("InsuranceAGMWB1::Validate");

    InsuranceAGMWB::Validate();

//    // these tests have to be here due to monitor dates for this structure is yearly, to review!!!
//    // computes the maximum age among the groups
//    int maxAge = 0;
//    int maxAgeIndex = 0;
//    for (int iGroup = 0 ; iGroup < ageGpSizeAtStart.size() ; iGroup++){
//        if (ageAtStart[iGroup] > maxAge){
//            maxAge = max(maxAge, ageAtStart[iGroup]);
//            maxAgeIndex = iGroup;
//        }
//    }
//
//    // remember for this deal, the monitor date is monthly
//    // checks that the mortality rates are provided up to the last simulation date
//    if (dieOption && (mortalityRate.size() < maxAge + monitorDates.size()))
//        throw ModelException(routine, "the mortality rates for the last simulation dates are missing.");
//
//    // checks that the lapse rates are provided up to the last simulation date
//    if (lapseOption && (lapseRate.size() < monitorDates.size()))            
//        throw ModelException(routine, "the lapse rates for the last simulation dates are missing.");
//

    if ((deterministicWDType != 1) && (deterministicWDType != 2) && (deterministicWDType != 3)){
        throw ModelException(routine, "Wrong deterministicWDType! should be 1 or 2 or 3! ");
    }
    
    // checks that the inputs are consistent with step up feature
    if (isStepUp){
        if (!isDynamicStepUp && stepUpITMLevel < 0.0)
            throw ModelException(routine, "ITM level for step up has to be >= 0.0.");
    }

    int nbfeatureGpSizeAtStart = featureGpSizeAtStart.size();
    int nbWDITMLevels = wDITMLevels.size();
    int i;
                            
    if (deterministicWDType == 2){ //withdraw only depends on Policy Year
        if (nbWDITMLevels > 1){
            throw ModelException(routine, "If withdrawal depends only on policy year not on ITM level, then WDITMLevels should be have only 1 col.");
        }
    }

    if ((deterministicWDType == 2 || deterministicWDType == 3)) {
        //withdraw depends on Policy Year and/or ITM
        if (adjustmentAV==0.0){
            throw ModelException(routine, "adjustment Buffer for AV can't be 0.0 ");
        }

        for (i = 0; i < nbWDITMLevels - 1; i++){
            if (wDITMLevels[i] > wDITMLevels[i+1] ){
                throw ModelException(routine, "ITM levels for withdrawal should be increasing ");
            }                
        }
        
        if (startWDYears.get()){
            //check size
            int nbCols = startWDYears->numCols();
            if (nbCols != nbWDITMLevels){
                throw ModelException(routine,
                                        Format::toString("Nb of columns (%d) in startWDYears should be equal to nb of WDITMLevels (%d)",
                                                        nbCols,
                                                        nbWDITMLevels));
            }

            int nbRows = startWDYears->numRows();
            if (nbRows != nbfeatureGpSizeAtStart){
                throw ModelException(routine,
                                        Format::toString("Nb of rows (%d) in startWDYears should be equal to nb of featureGpSizeAtStart (%d)",
                                                        nbRows,
                                                        nbfeatureGpSizeAtStart));
            }

            //decrease on the dimension of nbWDITMLevels
            int j;
            for (i = 0; i < nbfeatureGpSizeAtStart; i++){
                for (j = 0; j < nbWDITMLevels - 1; j++){
                    if ( (*startWDYears)[j][i] < (*startWDYears)[j+1][i]  ){
                        throw ModelException(routine, "start withdraw year should be an decreasing function of ITM levels!");
                    }                
                }
            }

            //increase on the dimension of nbfeatureGpSizeAtStart
            for (j = 0; j < nbWDITMLevels ; j++){
                for (i = 0; i < nbfeatureGpSizeAtStart-1; i++){
                    if ( (*startWDYears)[j][i] > (*startWDYears)[j][i+1]  ){
                        throw ModelException(routine, "start withdraw year should be an increasing function!");
                    }                
                }
            }

            startWDYears->checkNonNegative();
        }

        // if no weights provided, 
        if (!startWDYears){
            int nbfeatureGpSizeAtStart = featureGpSizeAtStart.size();
            int nbWDITMLevels   = wDITMLevels.size();
            startWDYears = CDoubleMatrixSP(new CDoubleMatrix(nbfeatureGpSizeAtStart, nbWDITMLevels));
            int ifeatureGpSizeAtStart = 0;
            for (; ifeatureGpSizeAtStart < nbfeatureGpSizeAtStart; ++ifeatureGpSizeAtStart){
                int iWDITMLevels = 0;
                for (; iWDITMLevels < nbWDITMLevels; ++iWDITMLevels){
                    (*startWDYears)[ifeatureGpSizeAtStart][iWDITMLevels] = 0.0;
                }
            }
        }
    }else{
        //set the value to 100%
        featureGpSizeAtStart.resize(1,1.0 );

        if ((DEBUG_OUTPUT) && (db_ithWDGroup > 1)){
            throw ModelException(routine, "only 1 withdraw group supported!");
        }
    }
}
   
/** Invoked when Class is 'loaded' */
void InsuranceAGMWB1::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(InsuranceAGMWB1, clazz);
    SUPERCLASS(InsuranceAGMWB);
    IMPLEMENTS(IMCIntoProduct);
    IMPLEMENTS(LastSensDate);   
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultInsuranceAGMWB1);

    // Guaranteed Minimum Withdrawal Benefit feature

    //step up
    FIELD(isDynamicStepUp, "if yes, step ups are dynamic");
    FIELD_MAKE_OPTIONAL(isDynamicStepUp);
    FIELD(stepUpITMLevel, "ITM level for step up");
    FIELD_MAKE_OPTIONAL(stepUpITMLevel);

    //withdraw
    FIELD(isDynamicUtilization, "if yes, the withdrawal rate in GMWB benefit is dynamic "
                                "(i.e.determined with respect to ITM)");
    FIELD_MAKE_OPTIONAL(isDynamicUtilization);

    // Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature 
    FIELD(deterministicWDType, "withdraw starts at 1: 0; 2 or 3: different withdraw group may start WD at 0 or some time in the futur or depend on ITM");
    FIELD_MAKE_OPTIONAL(deterministicWDType);

    FIELD(startWDYears, "diff groups start to withdrawal at diff. policy year and ITM levels ");
    FIELD_MAKE_OPTIONAL(startWDYears);
}


///////////// product class ////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/* MC product class for insurance annuity */
//state variable version

/** equivalent to InstIntoMCProduct. */
InsuranceAGMWB1ProdSV::InsuranceAGMWB1ProdSV(const InsuranceAGMWB1* instGMWB1,
                        const SimSeriesSP& simSeries) :
                InsuranceAGMWBProdSV(instGMWB1, simSeries),  instGMWB1(instGMWB1){  //??????
          

    //overide the yearFrac defined in the parent class by 1
    //since GMWB1 is always Annual 
    VolRequestTimeSP volTimeReq(new VolRequestTime()); // vol request for time
	CVolProcessedSP volTimeProcessed(getMultiFactors()->factorGetProcessedVol(0, volTimeReq.get())); 
	TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric(); // time metric    
    for (int iDate = 1 ; iDate < instGMWB1->monitorDates.size() ; iDate++){
        yearFrac[iDate] = 1.0;
	}

    int nberGroup = instGMWB1->ageGpSizeAtStart.size();        
    //int nbfeatureGpSizeAtStart = 1;
    int nbfeatureGpSizeAtStart = instGMWB1->featureGpSizeAtStart.size();

    int iWD ;

    //overide the initial value for lastMaxGMWBWithdrawal...
    for (int iGroup = 0 ; iGroup < nberGroup ; iGroup++){
        for (iWD = 0 ; iWD < nbfeatureGpSizeAtStart ; iWD++){
            lastMaxGMWBWithdrawal[iWD][iGroup]=  instGMWB1->maxGMWBPercentage*yearFrac[1];  //account value = 1.0 initially, to review
            lastMaxGMWBWithdrawalSoFar[iWD][iGroup]= instGMWB1->maxGMWBPercentage*yearFrac[1];           
        }
    }

    if (  (instGMWB1->deterministicWDType == 2) //withdraw depends on policy year
        || (instGMWB1->deterministicWDType == 3)  ) //withdraw depends on policy year and ITM levels
    {
        //nbfeatureGpSizeAtStart = instGMWB1->featureGpSizeAtStart.size();
        //may start withdraw later

        isWithdraw.resize(nbfeatureGpSizeAtStart);
        isWithdrawSoFar.resize(nbfeatureGpSizeAtStart);
        for (iWD = 0 ; iWD < nbfeatureGpSizeAtStart ; iWD++){
            isWithdraw[iWD].resize(nberGroup,false);
            isWithdrawSoFar[iWD].resize(nberGroup,false);
        }
    }else{
        //isWithdraw.resize(nbfeatureGpSizeAtStart, true); //start withdraw from the beginning
        isWithdraw.resize(nbfeatureGpSizeAtStart);
        isWithdrawSoFar.resize(nbfeatureGpSizeAtStart);
        for (iWD = 0 ; iWD < nbfeatureGpSizeAtStart ; iWD++){
            isWithdraw[iWD].resize(nberGroup,true);
            isWithdrawSoFar[iWD].resize(nberGroup,true);
        }
    }
}

// decide if the client steps up or not and update benefit
void InsuranceAGMWB1ProdSV::updateStepUp(const int& timeFromStart, const int& iWDGroup){
    int iGroup;
    int nberGroup = instGMWB1->ageGpSizeAtStart.size();
    
    if (isWaitPeriod[iWDGroup] && (timeFromStart >= instGMWB1->initialWaitPeriod)) {
        isWaitPeriod[iWDGroup] = false;

        for (iGroup = 0 ; iGroup < nberGroup ; iGroup++){
            // authorize step-up from the end of the initial waiting period
            lastStepUp[iWDGroup][iGroup] = instGMWB1->stepUpPeriod;
        }
    }

    // step-ups only supported by GMWB benefit
    if (instGMWB1->isStepUp) {
        for (iGroup = 0 ; iGroup < nberGroup ; iGroup++){
            if (!isWaitPeriod[iWDGroup] && (lastStepUp[iWDGroup][iGroup] >= instGMWB1->stepUpPeriod)){/// to review
                if (!Maths::equals(0.0, benefit[iWDGroup][iGroup])){                                                        
                    // if step ups are static, compare to stepUpITMLevel
                    if (!instGMWB1->isDynamicStepUp){
                        if (accountValue[iWDGroup][iGroup] >= instGMWB1->stepUpITMLevel * benefit[iWDGroup][iGroup]){
                            benefit[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup];

                            //yearFrac from timeFromStart -1 to timeFromStart
                            // convert annual wd rate  to this period wd rate
                            double maxWDPerc = instGMWB1->maxGMWBPercentage 
                                                            * yearFrac[timeFromStart];

                            lastMaxGMWBWithdrawal[iWDGroup][iGroup] =
                                max(lastMaxGMWBWithdrawal[iWDGroup][iGroup], 
                                    maxWDPerc * accountValue[iWDGroup][iGroup]);
                            // set lastStepUp[iGroup] to -1 since lastStepUp[iGroup]
                            // is going to be incremented at the end of the loop
                            //lastStepUp[iWDGroup][iGroup] = 0;  //to review
                            lastStepUp[iWDGroup][iGroup] = -1;
                        }
                    }
                    // if step ups are dynamic, use formula for step up utilization
                    else{
                        // compute the step up ratio
                        double stepUpRatio = 0.0; // ratio of clients who steps up
                        if (!Maths::isZero(benefit[iWDGroup][iGroup])){
                            stepUpRatio = 2.0 * ((accountValue[iWDGroup][iGroup] / benefit[iWDGroup][iGroup]) - 1.0);
                            stepUpRatio = min(max(stepUpRatio, 0.0), 1.0);
                        }

                        // determine if the group is going to step up or not
                        double uniform; // uniform variable
                        uniRand->draw(1, &uniform);
                        bool stepUpUtilization = (uniform < stepUpRatio); // flag indicating if the group steps up or not

                        if (stepUpUtilization){
                            benefit[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup];

                            //yearFrac from timeFromStart -1 to timeFromStart
                            // convert annual wd rate  to this period wd rate
                            double maxWDPerc = instGMWB1->maxGMWBPercentage 
                                                            * yearFrac[timeFromStart];

                            lastMaxGMWBWithdrawal[iWDGroup][iGroup] =
                                max(lastMaxGMWBWithdrawal[iWDGroup][iGroup], 
                                    maxWDPerc * accountValue[iWDGroup][iGroup]);
                            // set lastStepUp[iGroup] to -1 since lastStepUp[iGroup]
                            // is going to be incremented at the end of the loop
                            lastStepUp[iWDGroup][iGroup] = -1;                                }
                    }
                }
            }
        }
    }

    for (iGroup = 0 ; iGroup < nberGroup ; iGroup++){
        ++lastStepUp[iWDGroup][iGroup];
    }
}

// update number of clients by taking into account the mortality rate
void InsuranceAGMWB1ProdSV::updateDeadClients(const int& timeFromStart, const int& iWDGroup) {

    if (instGMWB1->dieOption && (timeFromStart != 0)){
        for (int iGroup = 0 ; iGroup < instGMWB1->ageGpSizeAtStart.size() ; iGroup++){   
            // current age of the group
            int currentAge = instGMWB1->ageAtStart[iGroup] + timeFromStart;
            
            // mortality rate depending on the current age of the group
            //yearFrac from timeFromStart-1 to timeFromStart
            double mRForThisPeriod = yearFrac[timeFromStart]
                                    * instGMWB1->mortalityRate[currentAge-1];
            groupSize[iWDGroup][iGroup] *= (1.0 - mRForThisPeriod);
        }
    }
}

// update number of clients by taking into account the lapse rate
void InsuranceAGMWB1ProdSV::updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup) {

    double lapseRate = 0.0; // base lapse rate (the same for all clients)
    double dynamicLapseRate = 0.0; // same as base lapse rate if no dynamic assumptions
    double lapseRateMultiplier = 0.0; // lapse rate multiplier if dynamic lapse rate
    double floorLapseRate = 0.0; //instGMWB1->isDynamicLapseRate ? 0.02 : 0.0; // floor for lapse rate if dynamic lapse rate

    if (instGMWB1->lapseOption && (timeFromStart != 0)) {// timeFromStart = 0 at start date
        lapseRate = instGMWB1->lapseRate[timeFromStart-1]; // base lapse rate depending on the time from start date only
        // yearFrac from timeFromStart-1 to timeFromStart
        // convert annual lapse Rate to this period
        lapseRate *= yearFrac[timeFromStart];

        double ITMLevel; // in the money level

        int iGroup;
        for (iGroup = 0 ; iGroup < instGMWB1->ageGpSizeAtStart.size(); iGroup++){
            // ITM level for GMAB GMWB and GMDB benefits if dynamic assumptions for lapse rate
            if (instGMWB1->isDynamicLapseRate){
                if (instGMWB1->lapseRateMultiplierType == 1){ //old version     

                    floorLapseRate = 0.02; // floor for lapse rate if dynamic lapse rate
                    // account value = 0.0, ITM level is set to 1.0
                    if (!Maths::isPositive(accountValue[iWDGroup][iGroup])){
                        ITMLevel = 1.0;
                    }
                    // account value != 0.0 and GMWB benefit
                    else{
                        ITMLevel = (benefit[iWDGroup][iGroup] / accountValue[iWDGroup][iGroup]) - 1.0;
                    }

                    // the lapse rate at a floor set to 2% (when lapse rate is dynamic only)
                    lapseRateMultiplier = pow(max((1.0 - 2.0 * max(0.0, ITMLevel)), 0.000001), 1.5);
                }else if((instGMWB1->lapseRateMultiplierType == 2) ){ //only available for GMWB

                    //floorLapseRate = 0.0; //reset floor to zero for this contract

                    double ITMAv_ben = 0.0;

                    if (Maths::isPositive(benefit[iWDGroup][iGroup])){
                        ITMAv_ben = accountValue[iWDGroup][iGroup] 
                                / benefit[iWDGroup][iGroup] /instGMWB1->adjustmentAV;
                    }
                    lapseRateMultiplier = Maths::min(1.0, ITMAv_ben);

                    //applying additional multiplier depends on the current status of withdrawal.
                    if ((*withdrawals)[iGroup] > 0.0){ //current withdraw status
                        lapseRateMultiplier *= instGMWB1->withdrawMultiplier; //only active if ITM > 0
                    }else{
                        lapseRateMultiplier *= instGMWB1->noWithdrawMultiplier; //only active if ITM > 0
                    }
                }else{

                }
            }
            // ITM level set to 0.0 for GMAB GMWB and GMDB benefits if no dynamic assumptions for lapse rate
            else{
                lapseRateMultiplier = 1.0; //just use the basic lapse rate
            }

            dynamicLapseRate = max(lapseRateMultiplier * lapseRate, floorLapseRate);
            groupSize[iWDGroup][iGroup] *= (1.0 - dynamicLapseRate);
        }
    }
}


// compute coupons taken by clients
void InsuranceAGMWB1ProdSV::computeWithdrawals(const int& timeFromStart,                                                 
                                                DoubleArraySP withdrawals,
                                                const int& iWDGroup
                                                 ) {        
    int iGroup;
    for (iGroup = 0 ; iGroup < instGMWB1->ageGpSizeAtStart.size();  iGroup++){       
         //for GMWB, different group withdrawal.

        //isWithdraw is flag indicating if the withdrawal takes place or not
        // always true when no dynamic assumption for utilization
        if (instGMWB1->isDynamicUtilization) {//stochastic withdraw
            // if the current GMWB benefit is not zero
            // determine if withdrawal happens using a uniform variable
            if (!Maths::isZero(accountValue[iWDGroup][iGroup])){
                // compute the in the money level
                double ITMLevel = (benefit[iWDGroup][iGroup] / accountValue[iWDGroup][iGroup]) - 1.0;
                // compute the dynamic utilization factor
                double utilizationFactor = 0.4 * pow(max((1.0 + 2.0 * ITMLevel), 0.0), 1.5);
                utilizationFactor = min(max(utilizationFactor, 0.2), 1.0);

                double uniform; // uniform variable
                uniRand->draw(1, &uniform);

                isWithdraw[iWDGroup][iGroup] = (uniform < utilizationFactor);
            }
            // if the current GMWB benefit is not zero, the withdrawal always happens
            else{
                isWithdraw[iWDGroup][iGroup] = true;
            }
        }else{//deterministic withdraw, but may depends on ITM 
            double timeToWithdrawal;
                            
            switch (instGMWB1->deterministicWDType) {
                case 1:{
                    //withdraw from beginning
                }        
                break;
                case 2: { //timeToWithdrawal didn't depend on ITMLevel, iWDITMLevels ==0
                    if (!isWithdraw[iWDGroup][iGroup]){//haven't start to withdraw
                        timeToWithdrawal = (*instGMWB1->startWDYears)[0][iWDGroup];  
                        if (timeFromStart >= timeToWithdrawal) {
                            isWithdraw[iWDGroup][iGroup] = true; //start withdrawal
                        }else{
                            isWithdraw[iWDGroup][iGroup] = false; //not start withdrawal yet
                        }
                    }
                }
                break;
                case 3:{//timeToWithdrawal depends on ITMLevel, iWDITMLevels may switch to higher value, it may start WD earlier

                    //looking for index such that timeToWithdrawal >= startWDYears)[iWDITMLevels][iWDGroup]; 
                    if (!isWithdraw[iWDGroup][iGroup]){//haven't start to withdraw

                        calcDeterministicWDType3(timeFromStart, withdrawals, iGroup, iWDGroup);
                            
                    }else{//already start withdraw, so don't check cond. anymore
                        //do nothing
                    }                                            
                }
            }                             
        }

        // for GMWB, withdrawals expressed as a percentage of the initial amount
        (*withdrawals)[iGroup] = (isWithdraw[iWDGroup][iGroup] ? 1.0 : 0.0) * 
            ((timeFromStart == 0) ? 0.0 : max(min(lastMaxGMWBWithdrawal[iWDGroup][iGroup], benefit[iWDGroup][iGroup]), 0.0)) ;
        
    }
}

void InsuranceAGMWB1ProdSV::calcDeterministicWDType3(const int& timeFromStart, 
                                                DoubleArraySP withdrawals,
                                                const int& iGroup,
                                                const int& iWDGroup 
                                               ){

    double ITMLevel;

    if ( (Maths::isNegative(accountValue[iWDGroup][iGroup])) || (Maths::isZero(accountValue[iWDGroup][iGroup]))){
        ITMLevel = 100000000000000.0;
    }else{
        ITMLevel = instGMWB1->adjustmentAV * (benefit[iWDGroup][iGroup] / accountValue[iWDGroup][iGroup]);
    }      

    int iFirst;
    int nbWDITMLevel = instGMWB1->wDITMLevels.size();

    //looking for iFirst
    iFirst = nbWDITMLevel -1;
    while ((iFirst >= 0) && (timeFromStart >= (*instGMWB1->startWDYears)[iFirst][iWDGroup])){
        iFirst -= 1;
    }
    iFirst = (iFirst == nbWDITMLevel -1)? iFirst: (iFirst + 1);

    double timeStartWD = (*instGMWB1->startWDYears)[iFirst][iWDGroup];
    if(timeFromStart >= timeStartWD){ //check ITM cond only if timeFromStart >= scheduled WD year
        //check the ITM levels
        if (iFirst ==0){
            isWithdraw[iWDGroup][iGroup] = true;                          
        }else{
            isWithdraw[iWDGroup][iGroup] = (ITMLevel >= instGMWB1->wDITMLevels[iFirst-1]);
        }
    }


    //looking for iLast
//        int iLast = nbWDITMLevel -1;
//        for (int i = nbWDITMLevel -1; i>=0; i--){
//            if (( (*instGMWB1->startWDYears)[i][iWDGroup] <= timeFromStart) 
//                && ((*instGMWB1->startWDYears)[i][iWDGroup] > (*instGMWB1->startWDYears)[iLast][iWDGroup] ) ){
//                iLast = i;
//            }
//        }
//                            
//        double timeStartWD = (*instGMWB1->startWDYears)[iFirst][iWDGroup];
//        if(timeFromStart >= timeStartWD){ //check ITM cond only if timeFromStart >= scheduled WD year
//            //check the ITM levels
//            if (iFirst ==0){
//                isWithdraw[iWDGroup][iGroup] = (ITMLevel < instGMWB1->WDITMLevels[iLast]);                          
//            }else{
//                isWithdraw[iWDGroup][iGroup] = ((ITMLevel < instGMWB1->WDITMLevels[iLast]) 
//                                && ((ITMLevel >= instGMWB1->WDITMLevels[iFirst-1])));
//            }
//        }

}

void InsuranceAGMWB1ProdSV::calcPayoff(const int iWDGroup, const int db_pathCount){

    static const string routine("InsuranceAGMWB1Prod::calcPayoff");
    try {

//        const double* path = pathGen->Path(0,0); // access path
//        int beginIdx = pathGen->begin(0);
//        int endIdx = pathGen->end(0);

        //take the first asset
        const SVPath& pathAsset1 = spotSV->path(0); //
        bool  doingPast = pathAsset1.doingPast();
        int beginIdx = pathAsset1.begin();
        int endIdx = pathAsset1.end();

        int iGroup;
        int nberGroup = instGMWB1->ageGpSizeAtStart.size();        

        DateTime settlementDate; // settlement date for the current start date anniversary
        DateTime valueDate = instGMWB1->valueDate; // value date

        double indexValue; // index value (equal to 1$ at start date)
        DoubleArraySP currentWithdrawals = DoubleArraySP(new DoubleArray(nberGroup)); 

        for (int Idx = beginIdx ; Idx < endIdx ; Idx++){

            double JPMFees_Idx = 0.0; // for debug output purpose
            //for debug output
            //only output one stock path for now, to review!!
            Debug_outputBegin(doingPast, instGMWB1->db_ithAsset, iWDGroup, db_pathCount, Idx);

            // update index value
            double baskV0 = computeBasket(0);
            double baskVIdx = computeBasket(Idx);
            double baskVIdxMinus1 = (Idx== 0) ? 1.0 : computeBasket(Idx-1);

            indexValue = 1.0 * (baskVIdx / baskV0);
            //indexValue = 1.0 * (path[Idx] / path[0]);

            // update settlement date
            settlementDate = instGMWB1->instSettle->settles(instGMWB1->monitorDates[Idx], 0);
            bool alreadySettled = valueDate.isGreater(settlementDate);

            // update account value (before coupon is paid)
            for (iGroup = 0 ; iGroup < nberGroup ; iGroup++){
                if (Idx != 0) // start date
                    //accountValue[iWDGroup][iGroup] *= (path[Idx] / path[Idx-1]);
                    accountValue[iWDGroup][iGroup] *= (baskVIdx / baskVIdxMinus1);
            }

            // compute coupons taken by clients
            // calculated as percentage of the initial amount for GMWB
            computeWithdrawals(Idx, currentWithdrawals, iWDGroup);

            //update step up for GMWB
            updateStepUp(Idx, iWDGroup);

            // update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)
            updatebenefit(currentWithdrawals, iWDGroup);
            // update account value
            updateAccountValue(currentWithdrawals, alreadySettled, iWDGroup);

            // compute fees (supposed to be paid at the begining of the period)
            if (!alreadySettled && (Idx !=0)){
                JPMFees_Idx = computeFees(&feeValue, Idx, iWDGroup);
            }

            // remove the clients which have lapsed during the year
            // additional lapse rate Multiplier may depend on the current status of withdrawal
            updateLapseClients(Idx, currentWithdrawals, iWDGroup); 

            if (!alreadySettled && (Idx != 0)){
                // add value of Guaranteed Minimum Death Benefit (GMDB) option
                for (iGroup = 0 ; iGroup < nberGroup ; iGroup++){
                    if (Maths::isNegative(accountValue[iWDGroup][iGroup])){                        
                        double df = dfSV->path()[Idx];
                        pathValue +=  df *  // the fees are discount to value date
                            groupSize[iWDGroup][iGroup] * (-accountValue[iWDGroup][iGroup]);

                        // for output
                        Debug_outputMid(iWDGroup, iGroup, db_pathCount, Idx);

                        accountValue[iWDGroup][iGroup] = 0.0;
                    }
                    // if the GMWB benefit is inferior to 0.0, we set the account value at 0.0
                    if (!Maths::isPositive(benefit[iWDGroup][iGroup]))
                        accountValue[iWDGroup][iGroup] = 0.0;

                }
            }

            // remove the clients which are dead during the year
            // (we assume that the clients die at the end of the year)
            updateDeadClients(Idx, iWDGroup);

            //for debug output
            //only output one stock path for now, to review!!
            Debug_outputEnd(doingPast, currentWithdrawals, JPMFees_Idx, instGMWB1->db_ithAsset, iWDGroup, db_pathCount, Idx);
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


//////////////////////////////////////////////////////////////////////////

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* InsuranceAGMWB1::createProduct(const MonteCarlo* model) const {
    static const string routine = "InsuranceAGMWB::createProduct";
    try {
        if(model->stateVarUsed()) {
            // we need to create a SimSeries object which says which assets need
            // which dates to be simulated
            SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                        one */
            simSeries->addDates(monitorDates);

            // State variables
            return new InsuranceAGMWB1ProdSV(this, simSeries);
        }else{
            throw ModelException(routine, "only useStateVars = True is supported for GMWB products.");
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

CClassConstSP const InsuranceAGMWB1::TYPE = CClass::registerClassLoadMethod(
    "InsuranceAGMWB1", typeid(InsuranceAGMWB1), InsuranceAGMWB1::load);

// force linker to include this file (avoid having header file) */
extern bool InsuranceAGMWB1Load()
{
    return true && InsuranceAGMWB1::TYPE;
}

DRLIB_END_NAMESPACE


