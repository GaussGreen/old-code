
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAGMWB.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 Guaranteed Minimum Withdrawal Benefit (GMWB)
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/InsuranceAGMWB.hpp"

DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////

void InsuranceAGMWB::Validate() {
    static const string routine("InsuranceAGMWB::Validate");

    InsuranceA::Validate();

    if (adjustmentAV==0.0){
        throw ModelException(routine, "adjustment Buffer for AV can't be 0.0 ");
    }
    
    // checks that the inputs are consistent with step up feature
    if (isStepUp){
        if (initialWaitPeriod < 0)
            throw ModelException(routine, "if step ups allowed, the initial waiting period has to be >= 0.");

        if (stepUpPeriod <= 0)
            throw ModelException(routine, "if step ups allowed, the step up period has to be >= 1.");
    }
                            
    if (DEBUG_OUTPUT) {
        if ((db_ithWDGroup < 1) || (db_ithAgeGroup < 1)){
            throw ModelException(routine, "output index for groups >= 1!");
        }

        if (  (db_ithAsset + 1 > assets->NbAssets()) || (db_ithAsset < 0) ){
            throw ModelException(routine, "output iAsset >= 0 and < nbAssets!");
        }
    }    
}
   
/** Invoked when Class is 'loaded' */
void InsuranceAGMWB::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(InsuranceAGMWB, clazz);
    SUPERCLASS(InsuranceA);
    IMPLEMENTS(IMCIntoProduct);
    IMPLEMENTS(LastSensDate);   
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultInsuranceAGMWB);

    // Guaranteed Minimum Withdrawal Benefit feature
    //withdraw
    FIELD(maxGMWBPercentage, "maximal percentage of initial amount for annual withdrawal");

    FIELD(adjustmentAV, "buffer to adjust the account value for lapse or withdrawal");
    FIELD_MAKE_OPTIONAL(adjustmentAV);

    //step up
    FIELD(isStepUp, "is there a step up feature");
    FIELD_MAKE_OPTIONAL(isStepUp);
    FIELD(initialWaitPeriod, "initial waiting period");
    FIELD_MAKE_OPTIONAL(initialWaitPeriod);
    FIELD(stepUpPeriod, "step up period");
    FIELD_MAKE_OPTIONAL(stepUpPeriod);

    // Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature 
    FIELD(wDITMLevels, "Withdrawal ITM levels for different feature group");
    FIELD_MAKE_OPTIONAL(wDITMLevels);

    //Debug output helper
    FIELD(DEBUG_OUTPUT , "yes, output some debug info");
    FIELD_MAKE_OPTIONAL(DEBUG_OUTPUT);

    FIELD(DEBUG_OUTPUT_ESTIMATE , "yes, output path and some other info for estimation purpose");
    FIELD_MAKE_OPTIONAL(DEBUG_OUTPUT_ESTIMATE);

    FIELD(db_ithPath, "ith Path to be output");
    FIELD_MAKE_OPTIONAL(db_ithPath);

    FIELD(db_ithAgeGroup , "ith Age Group to be output");
    FIELD_MAKE_OPTIONAL(db_ithAgeGroup);

    FIELD(db_ithWDGroup, "ith withdraw group to be output");
    FIELD_MAKE_OPTIONAL(db_ithWDGroup);

    FIELD(db_ithAsset, "ith asset path to be output");
    FIELD_MAKE_OPTIONAL(db_ithAsset);
}


///////////// product class ////////////////////////

//////////////////////////////////////////////////////////////////////////
//state variable
////////////////////////////////////////////////////////////////////

/** equivalent to InstIntoMCProduct. */
InsuranceAGMWBProdSV::InsuranceAGMWBProdSV(const InsuranceAGMWB* instGMWB,
                        const SimSeriesSP& simSeries) :
                InsuranceAProdSV(instGMWB, simSeries),  instGMWB(instGMWB){  //??????
          
    int nberAgeGp = instGMWB->ageGpSizeAtStart.size();        
    int nberFeatureGp = instGMWB->featureGpSizeAtStart.size();

    int iWD ;

    // initialize isWaitPeriod flag
    isWaitPeriod .resize(nberFeatureGp, true);
    isWaitPeriodSoFar.resize(nberFeatureGp, true);

    // resize and initialize benefits and maximum withdrawal for GMWB
    // initialize and resize age group sizes
    lastMaxGMWBWithdrawal.resize(nberFeatureGp);
    lastMaxGMWBWithdrawalSoFar.resize(nberFeatureGp);

    // resize and initialize last step-up flag
    lastStepUp.resize(nberFeatureGp);
    lastStepUpSoFar.resize(nberFeatureGp);

    for (iWD = 0 ; iWD < nberFeatureGp ; iWD++){
        lastMaxGMWBWithdrawal[iWD].resize(nberAgeGp, instGMWB->maxGMWBPercentage*yearFrac[1]);  //account value = 1.0 initially, to review
        lastMaxGMWBWithdrawalSoFar[iWD].resize(nberAgeGp, instGMWB->maxGMWBPercentage*yearFrac[1]);           

//        lastStepUpSoFar[iWD].resize(nberAgeGp,-1);
//        lastStepUp[iWD].resize(nberAgeGp,-1);
        lastStepUpSoFar[iWD].resize(nberAgeGp,0);
        lastStepUp[iWD].resize(nberAgeGp,0);

    }
    
    //for Debug output
    if (instGMWB->DEBUG_OUTPUT){            
        db_pathCount = -1; //0;
        int size = instGMWB->monitorDates.size();

        db_index.resize(size);

        db_accountValueBeg.resize(size); 
        db_benefitBeg.resize(size);
        db_groupSizeBeg.resize(size);

        db_accountValue.resize(size); 
        db_benefit.resize(size);
        db_currentWithdrawals.resize(size);
        db_lastStepUp.resize(size);
        db_groupSize.resize(size);
        db_payoff.resize(size);
        db_JPMFees.resize(size);

        db_negativeAV.resize(size);

        if (instGMWB->DEBUG_OUTPUT_ESTIMATE){            

            db_Mindex.resize(instGMWB->monitorDates.size());            
            for(int i = 0; i <instGMWB->monitorDates.size(); i++ ){
                db_Mindex[i].resize(instGMWB->db_maxIte, 0.0);
            }

            db_benefitAtT.resize(instGMWB->db_maxIte, 0.0);
            db_accountValueAtT.resize(instGMWB->db_maxIte, 0.0);
            db_nonLapseAtT.resize(instGMWB->db_maxIte, 0.0);
        }
    }
}

/** compute fees (supposed to be paid at the begining of the period) */
/** return this period fees */
double InsuranceAGMWBProdSV::computeFees(double* fees, const int& timeFromStart, const int& iWDGroup){

    double feeForThisPeriod = 0.0;
    for (int iGroup = 0 ; iGroup < instGMWB->ageGpSizeAtStart.size(); iGroup++){
        //JPM only collects Fees if the benefit is not zero
        if (Maths::isPositive(benefit[iWDGroup][iGroup])){
            // yearFrac from timeFromStart-1 to timeFromStart
            // convert fees Rate to this period
            feeForThisPeriod = 
                        instGMWB->feeRateJPM * yearFrac[timeFromStart] ;

            feeForThisPeriod *= groupSize[iWDGroup][iGroup] 
                                * accountValue[iWDGroup][iGroup] ;

//            *fees +=  feeForThisPeriod 
//                        * growthFactor[timeFromStart]; // the fees are forward valued to last date
            double df = dfSV->path()[timeFromStart];
            *fees +=  feeForThisPeriod * df; // the fees are discount to value date
        }
    }
    return feeForThisPeriod;  //only for output purpose
}

// update number of clients by taking into account the mortality rate
//void InsuranceAGMWBProdSV::updateDeadClients(const int& timeFromStart, const int& iWDGroup) {
//
//    if (instGMWB->dieOption && (timeFromStart != 0)){
//        for (int iGroup = 0 ; iGroup < instGMWB->ageGpSizeAtStart.size() ; iGroup++){   
//            // current age of the group
//            int currentAge = instGMWB->ageAtStart[iGroup] + timeFromStart;
//            // mortality rate depending on the current age of the group
//            // yearFrac from timeFromStart-1 to timeFromStart
//            double mRForThisPeriod = yearFrac[timeFromStart]
//                                    * instGMWB->mortalityRate[currentAge-1];
//            groupSize[iWDGroup][iGroup] *= (1.0 - mRForThisPeriod);
//        }
//    }
//}

// update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)
void InsuranceAGMWBProdSV::updatebenefit(DoubleArraySP withdrawals, const int& iWDGroup) {

    for (int iGroup = 0 ; iGroup < instGMWB->ageGpSizeAtStart.size();  iGroup++){
        benefit[iWDGroup][iGroup] -= (*withdrawals)[iGroup];

        if (Maths::isNegative(benefit[iWDGroup][iGroup])){
            benefit[iWDGroup][iGroup] = 0.0;
            //normally, we should never come here.
        }
    }
}

// update account value
void InsuranceAGMWBProdSV::updateAccountValue(DoubleArraySP withdrawals, 
                                                 bool alreadySettled,
                                                 const int& iWDGroup) {        

    for (int iGroup = 0 ; iGroup < instGMWB->ageGpSizeAtStart.size();   iGroup++){        
        // for GMWB, withdrawals expressed as a percentage of the initial amount
        accountValue[iWDGroup][iGroup] -= (*withdrawals)[iGroup];
        if (alreadySettled){
            // if the account value is inferior to 0 in the past, we set it at 0
            if (Maths::isNegative(accountValue[iWDGroup][iGroup]))
                accountValue[iWDGroup][iGroup] = 0.0;

            // if the GMWB benefit is inferior or equal to 0.0, we set the account value at 0.0
            if (!Maths::isPositive(benefit[iWDGroup][iGroup]))
                accountValue[iWDGroup][iGroup] = 0.0;
        }        
    }
}

void InsuranceAGMWBProdSV::payoff(const IPathGenerator*  pathGen, IMCPrices& prices) {
    static const string routine("InsuranceAGMWBProdSV::payoff");
    try {
        InsuranceAPrices& myprices = static_cast<InsuranceAPrices&>(prices);

        int iGroup;
        int iWDGroup ;

        // set a flag isDoingPast
        //bool doingPast = pathGen->doingPast();

        const SVPath& pathAsset1 = spotSV->path(0); //
        bool  doingPast = pathAsset1.doingPast();

        //int nberAgeGp = groupSize.size();
        int nberAgeGp = instGMWB->ageGpSizeAtStart.size();        
        int nberFeatureGp = instGMWB->featureGpSizeAtStart.size(); //1 if is not GMWB and deterministicWDType 2 or 3

        // preservation of the past
        if (!doingPast){
            for (iWDGroup = 0; iWDGroup < nberFeatureGp; iWDGroup++){

                isWaitPeriod[iWDGroup] = isWaitPeriodSoFar[iWDGroup]; // isWaitingPeriod flag

                for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
                    pathValue = 0.0; //pathValueSoFar; // path value
                    feeValue = 0.0; //feeValueSoFar; // present value of the fees
                    accountValue[iWDGroup][iGroup] = accountValueSoFar[iWDGroup][iGroup]; // current account value
                    lastStepUp[iWDGroup][iGroup] = lastStepUpSoFar[iWDGroup][iGroup]; // nber of years since the last step-up for the different groups
                    benefit[iWDGroup][iGroup] = benefitSoFar[iWDGroup][iGroup]; // garanteed minimum withdrawal benefit for the different age groups                                
                    lastMaxGMWBWithdrawal[iWDGroup][iGroup] = lastMaxGMWBWithdrawalSoFar[iWDGroup][iGroup]; // last maximum GMWB withdrawal
                    groupSize[iWDGroup][iGroup] = groupSizeSoFar[iWDGroup][iGroup];
                    isWithdraw[iWDGroup][iGroup]=isWithdrawSoFar[iWDGroup][iGroup];
                }
            }
        }

        //GMWB, loop on different withdrawal groups

        //debug output, count the nb of simulation        
        if (instGMWB->DEBUG_OUTPUT){            
            db_pathCount += 1;
        }

        for (iWDGroup = 0; iWDGroup < nberFeatureGp; iWDGroup++){
            calcPayoff(iWDGroup, db_pathCount);
        }

        // preservation of the past 
        if (doingPast){
            for (iWDGroup = 0; iWDGroup < nberFeatureGp; iWDGroup++){

                isWaitPeriodSoFar[iWDGroup] = isWaitPeriod[iWDGroup]; // isWaitingPeriod flag

                for (iGroup = 0 ; iGroup < nberAgeGp ; iGroup++){
                    pathValueSoFar = pathValue; // path value
                    feeValueSoFar = feeValue; // present value of the fees
                    accountValueSoFar[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup]; // current account value;
                    benefitSoFar[iWDGroup][iGroup] = benefit[iWDGroup][iGroup]; // garanteed minimum withdrawal benefit for the different age groups
                    lastMaxGMWBWithdrawalSoFar[iWDGroup][iGroup] = lastMaxGMWBWithdrawal[0][iGroup]; // last maximum GMWB withdrawal
                    groupSizeSoFar[iWDGroup][iGroup] = groupSize[iWDGroup][iGroup];
                    isWithdrawSoFar[iWDGroup][iGroup]= isWithdraw[iWDGroup][iGroup];
                    lastStepUpSoFar[iWDGroup][iGroup] =lastStepUp[iWDGroup][iGroup]; // nber of years since the last step-up for the different groups
                }
            }
        }

        if (!doingPast){
            // Discount
            double dF = 1; //matDfSV->firstDF();
            myprices.add((instGMWB->notional * pathValue *dF), InsuranceAPrices::OPTION_PRICE); // option price
            myprices.add((instGMWB->notional * feeValue * dF), InsuranceAPrices::FEE_PRICE); // fee price
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void InsuranceAGMWBProdSV::Debug_outputBegin(const bool doingPast, 
                                          const int iAsset,   //0 if mono
                                          const int iWDGroup, 
                                          const int db_pathCount,
                                          const int Idx){   
    static const string routine("InsuranceAGMWBProdSV::Debug_outputBegin");
    try {

        //for debug output
        if (instGMWB->DEBUG_OUTPUT ){ // not doing past

            int db_iWDG = (instGMWB->DEBUG_OUTPUT )? instGMWB->db_ithWDGroup - 1 : 0;
            const SVPath& path = spotSV->path(iAsset);

            //const double* path = pathGen->Path(iAsset,0); // access path
            
            if (doingPast || ( (db_pathCount == instGMWB->db_ithPath) && (db_iWDG ==iWDGroup) ) ){
                int db_iAgeG = instGMWB->db_ithAgeGroup - 1;
                db_accountValueBeg[Idx] = accountValue[db_iWDG][db_iAgeG];
                db_benefitBeg[Idx] = benefit[db_iWDG][db_iAgeG];
                db_groupSizeBeg[Idx] = groupSize[db_iWDG][db_iAgeG];
                db_index[Idx] = path[Idx]; //output path for tests purpose
            }

            if (instGMWB->DEBUG_OUTPUT_ESTIMATE ){ // not doing past
                if(doingPast || db_pathCount < instGMWB->db_maxIte){
                    db_Mindex[Idx][db_pathCount] = path[Idx]; //output path for tests purpose
                }else{
                    throw ModelException(routine, "nb of iteration is bigger than maximum iteration allowed for DEBUG_OUTPUT_ESTIMATE.");
                }
            }
        }
    }catch (exception& e){
        throw ModelException(e, routine);
    }
}

void InsuranceAGMWBProdSV::Debug_outputMid(const int iWDGroup,
                                        const int iGroup,
                                        const int db_pathCount,
                                        const int Idx){

    static const string routine("InsuranceAGMWBProdSV::Debug_outputMid");
    try {
        // for output
        if (instGMWB->DEBUG_OUTPUT){
            int db_iWDG = (instGMWB->DEBUG_OUTPUT )? instGMWB->db_ithWDGroup - 1 : 0;

            if (( db_pathCount == instGMWB->db_ithPath) && (db_iWDG ==iWDGroup) ){
                db_negativeAV[Idx] = accountValue[iWDGroup][iGroup];
                db_payoff[Idx] = groupSize[iWDGroup][iGroup] * (-accountValue[iWDGroup][iGroup]);
            }
        }
    }catch (exception& e){
        throw ModelException(e, routine);
    }
}

void InsuranceAGMWBProdSV::Debug_outputEnd(const bool doingPast,  
                                        const DoubleArraySP withdrawals,
                                        const double JPMFees_Idx,
                                        const int iAsset,   //0 if mono
                                        const int iWDGroup,
                                        const int db_pathCount,
                                        const int Idx){

    static const string routine("InsuranceAGMWBProdSV::Debug_outputMid");
    try {
        //for debug output
        if (instGMWB->DEBUG_OUTPUT ){ // not doing past
            int db_iWDG = (instGMWB->DEBUG_OUTPUT )? instGMWB->db_ithWDGroup - 1 : 0;
            if ( doingPast || (( db_pathCount == instGMWB->db_ithPath) && (db_iWDG ==iWDGroup)) ){
                int db_iAgeG = instGMWB->db_ithAgeGroup - 1;
                db_accountValue[Idx] = accountValue[db_iWDG][db_iAgeG];
                db_benefit[Idx] = benefit[db_iWDG][db_iAgeG];
                //db_isWithdraw[Idx] = isWithdraw[db_ithWDGroup][db_ithAgeGroup]
                db_currentWithdrawals[Idx] = (*withdrawals)[db_iAgeG]; 
                db_lastStepUp[Idx]  = lastStepUp[db_iWDG][db_iAgeG];
                db_groupSize[Idx] = groupSize[db_iWDG][db_iAgeG];
                db_JPMFees[Idx] = JPMFees_Idx;
            }

            if(instGMWB->DEBUG_OUTPUT_ESTIMATE){
                const SVPath& path = spotSV->path(iAsset);
                int endIdx = path.end();
                if ((db_iWDG ==iWDGroup) && (Idx == endIdx -1) ){
                    int db_iAgeG = instGMWB->db_ithAgeGroup - 1;

                    db_accountValueAtT[db_pathCount] = accountValue[db_iWDG][db_iAgeG];
                    db_benefitAtT[db_pathCount] = benefit[db_iWDG][db_iAgeG];
                    db_nonLapseAtT[db_pathCount] = groupSize[db_iWDG][db_iAgeG];
                }
            }               
        }
    }catch (exception& e){
        throw ModelException(e, routine);
    }
}

// control variate is done here
void InsuranceAGMWBProdSV::recordExtraOutput(Control* control, Results* results, const IMCPrices& prices) const {

    const InsuranceAPrices& myprices = static_cast<const InsuranceAPrices&>(prices);

    // pv from last date to value date
//    double discountFactor = instGMWB->instSettle->pv(instGMWB->valueDate,  
//                 instGMWB->monitorDates[instGMWB->monitorDates.size()-1],
//                 instGMWB->discount.get(),
//                 0/*instGMWB->asset.get()*/);

    if (control && control->isPricing()) {
        OutputRequest* request =
        control->requestsOutput(OutputRequest::OPTION_PRICE);
        if (request) {
            double optionPrice = myprices.getOptionPrice();
            // need to present value here
//            optionPrice *= discountFactor;
            results->storeRequestResult(request, optionPrice);
        }
        request = control->requestsOutput(OutputRequest::FEE_PRICE);
        if (request) {
            double feePrice = myprices.getFeePrice();
            // need to present value here
//            feePrice *= discountFactor;
            results->storeRequestResult(request, feePrice);
        }

        // ---debug output---

        if ((control && control->isPricing()) && (instGMWB->DEBUG_OUTPUT)){
            int i;
            int nbVar = 12;
            int nbVar1 = 3;
            int numSmth = db_accountValue.size();
            DateTimeArraySP tmp(new DateTimeArray(numSmth));
            CStringArraySP tmp_name(new StringArray( nbVar));
            CStringArraySP tmp_name1(new StringArray( nbVar1));

            if (numSmth > 0){
                CDoubleMatrixSP smthMatrix(new DoubleMatrix(nbVar, numSmth));

                //CStringMatrixSP tmp_path(new StringMatrix(1, nbVar));

                (*tmp_name)[0] = "Path";
                (*tmp_name)[1] = "AV at Beg";
                (*tmp_name)[2] = "AV at End ";
                (*tmp_name)[3] = "benefit at Beg";
                (*tmp_name)[4] = "benefit at End";
                (*tmp_name)[5] = "db_currentWithdrawals";
                (*tmp_name)[6] = "db_lastStepUp";
                (*tmp_name)[7] = "db_negativeAV";
                (*tmp_name)[8] = "db_payoff";
                (*tmp_name)[9] = "groupSize at Beg";
                (*tmp_name)[10] = "groupSize at End";
                (*tmp_name)[11] = "JPM Fees for this period";

                for (i=0; i<numSmth; i++){

                    (*tmp)[i] = instGMWB->monitorDates[i];

                    (*smthMatrix)[0][i] = db_index[i];

                    (*smthMatrix)[1][i] = db_accountValueBeg[i];

                    (*smthMatrix)[2][i] = db_accountValue[i];

                    (*smthMatrix)[3][i] = db_benefitBeg[i];

                    //(*tmp_name)[0][1] = "db_accountValue";
                    (*smthMatrix)[4][i] = db_benefit[i];

                    (*smthMatrix)[5][i] = db_currentWithdrawals[i];

                    //(*tmp_name)[3] = "db_lastStepUp";
                    (*smthMatrix)[6][i] = db_lastStepUp[i];
                    
                    (*smthMatrix)[7][i] = db_negativeAV[i];

                    (*smthMatrix)[8][i] = db_payoff[i];

                    (*smthMatrix)[9][i] = db_groupSizeBeg[i];

                    (*smthMatrix)[10][i] = db_groupSize[i];

                    (*smthMatrix)[11][i] = db_JPMFees[i];
                }

                results->storeGreek(tmp_name, Results::DEBUG_PACKET, OutputNameSP(new OutputName("PATH_INFO_NAME")));
                results->storeGreek(tmp, Results::DEBUG_PACKET, OutputNameSP(new OutputName("SIM_DATE")));
                results->storeGreek(smthMatrix, Results::DEBUG_PACKET, OutputNameSP(new OutputName("PATH_INFO")));
            }

            if (numSmth > 0 && instGMWB->DEBUG_OUTPUT_ESTIMATE){
                CDoubleMatrixSP pathMatrix(new DoubleMatrix(numSmth, instGMWB->db_maxIte));  //db_maxIte is the maximum nb of simulation
                CDoubleMatrixSP tempMatrix(new DoubleMatrix(nbVar1, instGMWB->db_maxIte ));  //db_maxIte is the maximum nb of simulation
    
                for (i=0; i<numSmth; i++){
                    for (int j=0; j < instGMWB->db_maxIte; j++){
                        (*pathMatrix)[i][j] = db_Mindex[i][j];
                    }
                }

                (*tmp_name1)[0] = "AV";
                (*tmp_name1)[1] = "Ben";
                (*tmp_name1)[2] = "GroupSize";

                for (int j=0; j < instGMWB->db_maxIte; j++){
                    (*tempMatrix)[0][j] = db_accountValueAtT[j];
                    (*tempMatrix)[1][j] = db_benefitAtT[j];
                    (*tempMatrix)[2][j] = db_nonLapseAtT[j];
                }

                results->storeGreek(pathMatrix, Results::DEBUG_PACKET, OutputNameSP(new OutputName("PATH_INFO1")));
                results->storeGreek(tempMatrix, Results::DEBUG_PACKET, OutputNameSP(new OutputName("PATH_INFO2")));
            }
        }
    }
}

/////////////////////////////////////////////////////

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* InsuranceAGMWB::createProduct(const MonteCarlo* model) const {
    static const string routine = "InsuranceAGMWB::createProduct";
    try {
        if(model->stateVarUsed()) {
            // we need to create a SimSeries object which says which assets need
            // which dates to be simulated
            SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                        one */
            simSeries->addDates(monitorDates);

            // State variables
            return new InsuranceAGMWBProdSV(this, simSeries);
        }else{
            throw ModelException(routine, "only useStateVars = True is supported for GMWB products.");
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

CClassConstSP const InsuranceAGMWB::TYPE = CClass::registerClassLoadMethod(
    "InsuranceAGMWB", typeid(InsuranceAGMWB), InsuranceAGMWB::load);

// force linker to include this file (avoid having header file) */
extern bool InsuranceAGMWBLoad()
{
    return true && InsuranceAGMWB::TYPE;
}

DRLIB_END_NAMESPACE
