
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAnnuityGMB.cpp
//
//   Author      : Romain Garagnon
//
//   Description   Insurance annuity payoff
//                 1. Guaranteed Minimum Death Benefit (GMDB)
//                 2. Guaranteed Minimum Withdrawal Benefit (GMWB)
//                 3. Guaranteed Minimum Accumulation Benefit (GMAB)
//
//
//   $Log: InsuranceAnnuityGMB.cpp,v $
//   Revision 1.9  2005/05/04 18:59:08  rgaragnon
//   added dynamic suumptions (related to ITM level) for
//   withdrawals and dynamix utilization (GMWB benefit only)
//   lapse rate (all benefits)
//
//   Revision 1.8  2005/02/07 18:25:20  rgaragnon
//   re-designed the product implementation to enable combined GMAB + GMDB benefit
//
//   Revision 1.7  2005/01/28 14:05:20  rgaragnon
//   added step ups for GMWB benefit
//   corrected bug for pv of fees for the GMWB benefit
//
//   Revision 1.6  2005/01/24 17:12:36  rgaragnon
//   factorized update for account value and added step-ups variables
//
//   Revision 1.5  2005/01/10 19:19:38  rgaragnon
//   added some features (GMDB and GMAB) following new clients requests
//
//   Revision 1.4  2004/12/14 22:49:36  rgaragnon
//   modified some comments
//
//   Revision 1.3  2004/12/14 17:06:19  rgaragnon
//   modified some comments
//
//   Revision 1.2  2004/12/14 15:56:35  rgaragnon
//   added check against struck and protected currency treatments
//
//   Revision 1.1  2004/12/13 21:35:59  rgaragnon
//   Guaranteed Minimum Death / Accumulation / Withdrawal Benefits
//   Initial Version
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/Random.hpp"

DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////
class InsuranceAnnuityGMB: public Generic1Factor, 
                           virtual public LastSensDate,
                           virtual public IMCIntoProduct {
protected:
    //////// fields ////////
    
    // historical spot samples
    SampleListSP    spotSamples; // sample dates and samples (when the sample date is in the past)

    // mortality and lapse rate definitions
    DoubleArray     mortalityRate; // mortality rate as function of the client age
    DoubleArray     lapseRate; // contract lapse rate as function of the time from start date
    bool            dieOption; // if yes, the mortality rates table is taken into account
                           // if no, the mortality rate is 0%
    bool            lapseOption; // if yes, the lapse rate table is taken into account
                             // if no, the mortality rate is 0%
    bool            isDynamicLapseRate; // if yes, the lapse rate is dynamic (i.e.determined with respect to ITM)
    int             lapseRateMultiplierType; //default: 1 (old version, power ft); 2: AV/Benefits if Benefit/AV >=1
    double          withdrawMultiplier; //additional multiplier depends on the current status of withdrawal.
    double          noWithdrawMultiplier;


    // age groups with percentages and age limit for benefit accumulation   
    IntArray        ageAtStart; // age of the people in the age group at start date
    DoubleArray     groupSizeAtStart; // percentage of people in the corresponding group

    double          initialAmount; //initial amount invested in the insurance policy
    
    
    
    bool            isStepUp; // is there a step up feature
    bool            isDynamicStepUp; // if yes, step ups are dynamic
    int             initialWaitPeriod; // initial waiting period
    int             stepUpPeriod; // step up period
    double          stepUpITMLevel; // ITM level for step up
    double          feeRate; // fee rate as percentage of the account value

    // Guaranteed Minimum Death Benefit feature
    bool            isGMDB; // is there a Guaranteed Minimum Death Benefit feature
    double          GMDBWithdrawRate; // GMDB withdrawal rate
    int             GMDBAgeCap; // age at which the GMDB benefit gets capped
    int             GMDBAgeLimit; // age limit for GMDB benefit
    bool            isAnnualRatchet; // is there an annual ratchet for the GMWD
    bool            isPercentageRollUp; // is there a percentage roll up feature
    double          percentageRollUp; // percentage of annual roll-up

    // Guaranteed Minimum Accumulation Benefit feature
    bool            isGMAB; // is there a Guaranteed Minimum Accumulation Benefit feature
    double          GMABWithdrawRate; // GMAB withdrawal rate
    int             GMABAgeLimit; // age limit for GMAB benefit
    DoubleArray     GMABPercentages; // percentages of initial amount guaranteed
    IntArray        GMABMaturities; // maturities of the GMAB benefit

    // Guaranteed Minimum Withdrawal Benefit feature
    bool            isGMWB; // is there a Guaranteed Minimum Withdrawal Benefit feature
    double          maxGMWBPercentage; // maximal percentage of initial amount for annual withdrawal
    bool            isDynamicUtilization; // if yes, the withdrawal rate in GMWB benefit is dynamic
                                  // (i.e.determined with respect to ITM)
    
    // Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature 
    DoubleArray     WDGroupPercentage;  //withdrawal group percentage.
    DoubleArray     WDITMLevels;  //withdrawal ITM Levels.
    CDoubleMatrixSP startWDYears; // diff groups start to withdrawal at diff. policy year and ITM levels         

    int             deterministicWDType; //withdraw starts at 
                                            //1: 0; 
                                            //2: different withdraw group may start WD at 0 or some time in the future; 
                                            //3: different withdraw group may start WD at 0 or some time in the future depending on ITM

    double adjustmentAV;  // a buffer to ajust the AV

    // Debug output
    bool DEBUG_OUTPUT_ESTIMATE ; //output info to do the estimation between 10 years and 30 years
    int db_maxIte; // hard caoded, for max iteration that the DEBUG_OUTPUT_ESTIMATE can handle $unregistered

    bool DEBUG_OUTPUT ;
    int db_ithPath ;
    int db_ithAgeGroup ;
    int db_ithWDGroup;  //1 if deterministicWDType =1

public:
    static CClassConstSP const TYPE;
    friend class InsuranceAnnuityGMBProd;

    virtual void Validate() {
        static const string routine("InsuranceAnnuityGMB::Validate");
        // NB don't call parent's validate - issues with fwdStart
        if (fwdStarting)
            throw ModelException(routine, "Fwd starting flag on "
                      "Generic1Factor is not used");
        if (oneContract)
            throw ModelException(routine, "oneContract flag on "
                      "Generic1Factor is not used");

        // checks that the asset is not currency struck
        if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
            throw ModelException(routine, "ccy struck type not supported.");

        // checks that the asset is not currency protected
        if (ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED)
            throw ModelException(routine, "ccy protected type not supported.");

        // checks that there is at least one age group
        if (ageAtStart.size() <= 0)
            throw ModelException(routine, "ageAtStart array size must be > 0.");

        // checks that the number of ages and the number of age percentages are consistent
        if (ageAtStart.size() != groupSizeAtStart.size())
            throw ModelException(routine, "ageAtStart and groupSizeAtStart arrays must be the same size.");

        // computes the sum of age percentages
        double SumPercentageGroup = 0.0;
        int iGroup;
        for (iGroup = 0 ; iGroup < groupSizeAtStart.size() ; iGroup++)
        {
            SumPercentageGroup += groupSizeAtStart[iGroup];
        }

        // checks that the sum of age percentages is 100%
        if (!Maths::equals(1.0, SumPercentageGroup))
            throw ModelException(routine, "the sum of the age percentages has to be 1.0.");

        // checks there is at least one past sample
        if (spotSamples->getDates().size() == 0)
            throw ModelException(routine, "there must be at least one sample date.");

        // checks the first sample date is on or before value date
//to review, do we really need this?
//        if (spotSamples->getFirstDate() > valueDate)
//            throw ModelException(routine, "there must be one sample date in the past (or = valueDate).");

        // computes the maximum age among the groups
        int maxAge = 0;
        int maxAgeIndex = 0;
        for (iGroup = 0 ; iGroup < groupSizeAtStart.size() ; iGroup++)
        {
            if (ageAtStart[iGroup] > maxAge)
            {
                maxAge = max(maxAge, ageAtStart[iGroup]);
                maxAgeIndex = iGroup;
            }
        }

        // checks that the mortality rates are provided up to the last simulation date
        if (dieOption && (mortalityRate.size() < maxAge + spotSamples->getDates().size()))
            throw ModelException(routine, "the mortality rates for the last simulation dates are missing.");

        // checks that the lapse rates are provided up to the last simulation date
        if (lapseOption && (lapseRate.size() < spotSamples->getDates().size()))
            throw ModelException(routine, "the lapse rates for the last simulation dates are missing.");

        // checks that GMDB and GMWB flags are not both set to Yes 
        if (!isGMDB && !isGMWB && !isGMAB)
            throw ModelException(routine, "isGMDB or isGMWB or isGMAB has to be set to Yes.");

        // checks that GMDB and GMWB flags are not both set to Yes 
        if (isGMWB && (isGMDB || isGMAB))
            throw ModelException(routine, "GMWB benefit can't be combined with another benefit.");

        // if the benefits GMWB and GMAB are combined
        // we have to make some checks on the consistency of the inputs
        if (isGMDB && isGMAB)
        {
            if (GMDBAgeLimit != GMABAgeLimit)
                throw ModelException(routine, "if the benefits GMDB and GMAB are combined "
                              "GMDB and GMAB age limits have to be the same.");

            if (!Maths::equals(GMDBWithdrawRate, GMABWithdrawRate))
                throw ModelException(routine, "if the benefits GMDB and GMAB are combined "
                              "GMDB and GMAB withdrawal rates have to be the same.");
        }

        // checks that the maximum age is strictly lower than GMDB age limit at start date
        if (isGMDB && (maxAge >= GMDBAgeLimit))
        {
            throw ModelException(routine, "GMDB age limit is greater or equal to the age of group " + 
                                                          Format::toString(maxAgeIndex) + " at start date.");
        }

        // checks that GMDB age limit is greater or equal to GMDB age cap
        if (isGMDB && (GMDBAgeCap >= GMDBAgeLimit))
        {
            throw ModelException(routine, "GMDB age cap is greater or equal to GMDB age limit");
        }

        // checks that the maximum age + GMAB maturity is strictly lower than GMAB age limit at start date
        if (isGMAB && (maxAge + GMABMaturities[GMABMaturities.size() - 1] > GMABAgeLimit))
        {
            throw ModelException(routine, "age of group " + Format::toString(maxAgeIndex) +
                                                          " at maturity is strictly greater than GMAB age limit.");

        }

        if ((lapseRateMultiplierType != 1) && (lapseRateMultiplierType != 2)){

            throw ModelException(routine, "Wrong lapseRateMultiplierType! should be 1 or 2 ");
        }

        if ((deterministicWDType != 1) && (deterministicWDType != 2) && (deterministicWDType != 3)){

            throw ModelException(routine, "Wrong deterministicWDType! should be 1 or 2 or 3! ");
        }
        
        // checks that the inputs are consistent with step up feature
        if (isStepUp)
        {
            if (!isGMWB)
                throw ModelException(routine, "step ups only allowed for GMWB benefit.");

            if (!isDynamicStepUp && stepUpITMLevel < 0.0)
                throw ModelException(routine, "ITM level for step up has to be >= 0.0.");

            if (initialWaitPeriod < 0)
                throw ModelException(routine, "if step ups allowed, the initial waiting period has to be >= 0.");

            if (stepUpPeriod <= 0)
                throw ModelException(routine, "if step ups allowed, the step up period has to be >= 1.");
        }


        int nbWDGroupPercentage = WDGroupPercentage.size();
        int nbWDITMLevels = WDITMLevels.size();
                               
        if (nbWDGroupPercentage == 1 && WDGroupPercentage[0] != 1.0){
            throw ModelException(routine, "withdrawal group percentage should be sum to 100%.");
        }

        if (deterministicWDType == 2){ //withdraw only depends on Policy Year
            if (nbWDITMLevels > 1)
            {
                throw ModelException(routine, "If withdrawal depends only on policy year not on ITM level, then WDITMLevels should be have only 1 col.");
            }
        }

        if (isGMWB && (deterministicWDType == 2 || deterministicWDType == 3)) //withdraw depends on Policy Year and/or ITM
        {

            if (adjustmentAV==0.0){
                throw ModelException(routine, "adjustment Buffer for AV can't be 0.0 ");
            }

            //check the sum of WDGroupPercentage is 100%
            double sum = 0.0;
            int i;
            for (i = 0; i < nbWDGroupPercentage; i++){
                sum += WDGroupPercentage[i];
            }

            if(sum != 1.0 )
            {
                throw ModelException(routine, "withdrawal group percentage should be sum to 100%.");
            }

            for (i = 0; i < nbWDITMLevels - 1; i++){
                if (WDITMLevels[i] > WDITMLevels[i+1] ){
                    throw ModelException(routine, "ITM levels for withdrawal should be increasing ");
                }                
            }
            
            if (startWDYears.get())
            {
                //check size
                int nbCols = startWDYears->numCols();
                if (nbCols != nbWDITMLevels){
                    throw ModelException(routine,
                                         Format::toString("Nb of columns (%d) in startWDYears should be equal to nb of WDITMLevels (%d)",
                                                          nbCols,
                                                          nbWDITMLevels));
                }

                int nbRows = startWDYears->numRows();
                if (nbRows != nbWDGroupPercentage)
                {
                    throw ModelException(routine,
                                         Format::toString("Nb of rows (%d) in startWDYears should be equal to nb of WDGroupPercentage (%d)",
                                                          nbRows,
                                                          nbWDGroupPercentage));
                }

                //decrease on the dimension of nbWDITMLevels
                int j;
                for (i = 0; i < nbWDGroupPercentage; i++){
                    for (j = 0; j < nbWDITMLevels - 1; j++){
                        if ( (*startWDYears)[j][i] < (*startWDYears)[j+1][i]  ){
                            throw ModelException(routine, "start withdraw year should be an decreasing function of ITM levels!");
                        }                
                    }
                }


                //increase on the dimension of nbWDGroupPercentage
                for (j = 0; j < nbWDITMLevels ; j++){
                    for (i = 0; i < nbWDGroupPercentage-1; i++){
                        if ( (*startWDYears)[j][i] > (*startWDYears)[j][i+1]  ){
                            throw ModelException(routine, "start withdraw year should be an increasing function!");
                        }                
                    }
                }

                startWDYears->checkNonNegative();
            }

            // if no weights provided, 
            if (!startWDYears)
            {
                int nbWDGroupPercentage = WDGroupPercentage.size();
                int nbWDITMLevels   = WDITMLevels.size();
                startWDYears = CDoubleMatrixSP(new CDoubleMatrix(nbWDGroupPercentage, nbWDITMLevels));
                int iWDGroupPercentage = 0;
                for (; iWDGroupPercentage < nbWDGroupPercentage; ++iWDGroupPercentage){
                    int iWDITMLevels = 0;
                    for (; iWDITMLevels < nbWDITMLevels; ++iWDITMLevels){
                        (*startWDYears)[iWDGroupPercentage][iWDITMLevels] = 0.0;
                    }
                }
            }
        }else{
            //set the value to 100%
            WDGroupPercentage.resize(1,1.0 );

            if ((DEBUG_OUTPUT) && (db_ithWDGroup > 1)){
                throw ModelException(routine, "only 1 withdraw group supported!");
            }
        }

        if ((DEBUG_OUTPUT) && ((db_ithWDGroup < 1) || db_ithAgeGroup < 1)){
            throw ModelException(routine, "output index for groups >= 1!");
        }


        AssetUtil::assetCrossValidate(asset.get(),
                                      false, //fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);
    }

    void validatePop2Object(){
        static const string routine("InsuranceAnnuityGMB::validatePop2Object");
        }
   
    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** when to stop tweaking - default implementation assumes product can
        be priced with a MC */
    DateTime endDate(const Sensitivity*) const {
        DateTime end = valueDate.rollDate(365 * (mortalityRate.size() - ageAtStart[0]));
        return end;
    }


private:
    InsuranceAnnuityGMB(): Generic1Factor(TYPE),
                            isDynamicLapseRate(false),
                            lapseRateMultiplierType(1),                            
                            withdrawMultiplier(1.0),
                            noWithdrawMultiplier(1.0),
                            isStepUp(false),
                            isDynamicStepUp(false),
                            initialWaitPeriod(-1),
                            stepUpPeriod(-1),
                            stepUpITMLevel(-1.0),
                            isDynamicUtilization(false),
                            deterministicWDType(1),
                            adjustmentAV(1.0)
    {
        WDGroupPercentage.resize(1,1.0);  //withdrawal group percentage.
        WDITMLevels.resize(1, 10000000);  //withdrawal ITM Levels.
        DEBUG_OUTPUT = false;
        db_ithPath = 1 ;
        db_ithAgeGroup = 1;
        db_ithWDGroup = 1;  //1 if deterministicWDType =1

        DEBUG_OUTPUT_ESTIMATE = true;
        db_maxIte = 100001;
    }


    // for reflection
    InsuranceAnnuityGMB(const InsuranceAnnuityGMB& rhs); // not implemented
    InsuranceAnnuityGMB& operator=(const InsuranceAnnuityGMB& rhs); // not implemented

    static IObject* defaultInsuranceAnnuityGMB() {
        return new InsuranceAnnuityGMB();
    }

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta) {
        // use valueDate before it changes
        spotSamples->roll(theta->getUtil(valueDate), 0 /* iAsset */,
        asset.get()); // then roll our past values
        Generic1Factor::sensShift(theta); // and then call parent's method
        return true; // continue to tweak components which implement Theta
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){

        REGISTER(InsuranceAnnuityGMB, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LastSensDate);   
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultInsuranceAnnuityGMB);
        FIELD(spotSamples, "sample dates and samples (when the sample date is in the past)");
        FIELD(mortalityRate, "mortality rate as function of the client age");
        //lapse
        FIELD(lapseOption, "if yes, the lapse rate table is taken into account");
        FIELD(lapseRate, "contract lapse rate as function of the time from start date");
        FIELD(isDynamicLapseRate, "if yes, the lapse rate is dynamic (i.e.determined with respect to ITM)");
        FIELD_MAKE_OPTIONAL(isDynamicLapseRate);
        FIELD(lapseRateMultiplierType, "1 power ft; 2: AV/Benefits, both are active only if Benefit/AV >=1");
        FIELD_MAKE_OPTIONAL(lapseRateMultiplierType);
        FIELD(withdrawMultiplier, "additional multiplier for lapse rate: depends on the current status of withdrawal.");
        FIELD_MAKE_OPTIONAL(withdrawMultiplier);
        FIELD(noWithdrawMultiplier, "additional multiplier for lapse rate: depends on the current status of withdrawal.");
        FIELD_MAKE_OPTIONAL(noWithdrawMultiplier);

        FIELD(dieOption, "if yes, the mortality rates table is taken into account");
        FIELD(ageAtStart, "age of the people in the age group at start date");
        FIELD(groupSizeAtStart, "the percentage of people in the corresponding age group");
        FIELD(initialAmount, "initial amount invested in the insurance policy");
        FIELD(isStepUp, "is there a step up feature");
        FIELD_MAKE_OPTIONAL(isStepUp);
        FIELD(isDynamicStepUp, "if yes, step ups are dynamic");
        FIELD_MAKE_OPTIONAL(isDynamicStepUp);
        FIELD(initialWaitPeriod, "initial waiting period");
        FIELD_MAKE_OPTIONAL(initialWaitPeriod);
        FIELD(stepUpPeriod, "step up period");
        FIELD_MAKE_OPTIONAL(stepUpPeriod);
        FIELD(stepUpITMLevel, "ITM level for step up");
        FIELD_MAKE_OPTIONAL(stepUpITMLevel);
        FIELD(feeRate, "fee rate as percentage of the account value");
        // Guaranteed Minimum Death Benefit feature
        FIELD(isGMDB, "is there a Guaranteed Minimum Death Benefit feature");
        FIELD(GMDBWithdrawRate, "GMDB withdrawal rate");
        FIELD(GMDBAgeCap, "age at which the GMDB benefit gets capped");
        FIELD(GMDBAgeLimit, "age limit for GMDB benefit");
        FIELD(isAnnualRatchet, "is there an annual ratchet for the GMDB");
        FIELD(isPercentageRollUp, "is there a percentage roll up feature");
        FIELD(percentageRollUp, "percentage of annual roll-up");
        // Guaranteed Minimum Accumulation Benefit feature
        FIELD(isGMAB, "is there a Guaranteed Minimum Accumulation Benefit feature");
        FIELD(GMABWithdrawRate, "GMDB withdrawal rate");
        FIELD(GMABAgeLimit, "age limit for GMAB benefit");
        FIELD(GMABPercentages, "percentages of initial amount garanteed");
        FIELD(GMABMaturities, "maturities of the GMAB benefit");

        // Guaranteed Minimum Withdrawal Benefit feature
        FIELD(isGMWB, "is there a Guaranteed Minimum Withdrawal Benefit feature");
        FIELD(maxGMWBPercentage, "maximal percentage of initial amount for annual withdrawal");
        FIELD(isDynamicUtilization, "if yes, the withdrawal rate in GMWB benefit is dynamic "
                                   "(i.e.determined with respect to ITM)");
        FIELD_MAKE_OPTIONAL(isDynamicUtilization);

        // Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature 
        FIELD(deterministicWDType, "withdraw starts at 1: 0; 2 or 3: different withdraw group may start WD at 0 or some time in the futur or depend on ITM");
        FIELD_MAKE_OPTIONAL(deterministicWDType);

        FIELD(WDGroupPercentage, "the percentage of people in the corresponding withdrawal group");
        FIELD_MAKE_OPTIONAL(WDGroupPercentage);

        FIELD(WDITMLevels, "ITM levels for different withdrawal group");
        FIELD_MAKE_OPTIONAL(WDITMLevels);

        FIELD(startWDYears, "diff groups start to withdrawal at diff. policy year and ITM levels ");
        FIELD_MAKE_OPTIONAL(startWDYears);

        FIELD(adjustmentAV, "buffer to adjust the account value");
        FIELD_MAKE_OPTIONAL(adjustmentAV);

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


        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for insurance annuity */
//////// product class //////////
class InsuranceAnnuityGMBProd : public IMCProduct, virtual public IMCProductLN {
private:
    const InsuranceAnnuityGMB* inst; // reference to original instrument
    DateTimeArray              simulationDates; // simulation dates
    //DoubleArray                histSpots; // all historical spots
    DateTimeArray              histDates; // all historical dates
    DoubleArray                growthFactor; // growth factor to last date
    double                     pathValue; // store the path value
    double                     feeValue; // store the present value of the fees

    RandUniformDefaultSP       uniRand; // generator of uniform random variables

    vector<DoubleArray>         accountValue; // current account value;

    vector<bool>                isWaitPeriod; // are we in the initial waiting period ?
    vector<bool>                isWaitPeriodSoFar; // are we in the initial waiting period ?

    vector<IntArray >           lastStepUp; // nber of years since the last step-up for the different groups
    vector<IntArray>            lastStepUpSoFar; // nber of years since the last step-up for the different groups


    DoubleArray                notionalGMADB; // garanteed minimum death and accumulation benefit notionals for the different age groups
    DoubleArray                benefitGMADB; // garanteed minimum death benefit for the different age groups

    vector<DoubleArray>         benefitGMWB; // garanteed minimum withdrawal benefit for the different age groups
    vector<DoubleArray>         lastMaxGMWBWithdrawal; // last maximum GMWB withdrawal

    //DoubleArray                 groupSize; // current percentage of people in the corresponding age group
    vector<DoubleArray>         groupSize; // current percentage of people in the corresponding age group

    vector<DoubleArray>         groupSizeSoFar; // current percentage of people in the corresponding age group

    /* for preservation of the past */
    double                     pathValueSoFar; // store the path value
    double                     feeValueSoFar; // store the present value of the fees

    vector<DoubleArray>        accountValueSoFar; // current account value;

    DoubleArray                notionalGMADBSoFar; // garanteed minimum death and accunulation benefit notionals for the different age groups
    DoubleArray                benefitGMADBSoFar; // garanteed minimum death benefit for the different age groups

    vector<DoubleArray>         benefitGMWBSoFar; // garanteed minimum withdrawal benefit for the different age groups
    vector<DoubleArray>         lastMaxGMWBWithdrawalSoFar; // last maximum GMWB withdrawal

    vector<BoolArray>           isWithdraw; //always true except for GMWB when deterministicWDType =2 or 3 or isDynamicUtilization = true
    vector<BoolArray>           isWithdrawSoFar;


    //for debug output

    int             db_pathCount;

    DoubleArray     db_index; //index value 

    DoubleArray     db_accountValueBeg;     
    DoubleArray     db_benefitBeg;     
    DoubleArray     db_groupSizeBeg;

    DoubleArray     db_accountValue;   
    DoubleArray     db_benefit;     
    DoubleArray     db_currentWithdrawals ;

    //DoubleArray     db_lastStepUp ; 
    DoubleArray     db_groupSize;
    DoubleArray     db_negativeAV;  //account value when option needs to pay, ie AV <0, only for GMWB
    DoubleArray     db_payoff;

    //debug output info for extimation between 10 and 30 years
    vector<DoubleArray>  db_Mindex; //index value  //[Idx, itera] 
    DoubleArray     db_benefitAtT;     //benefits amount at maturity,   [iter] 
    DoubleArray     db_accountValueAtT;   // acount value at maturity,  [iter]
    DoubleArray     db_nonLapseAtT;   //non lapse at maturity, [iter]


    class InsuranceAnnuityGMBPrices;
    typedef refCountPtr<InsuranceAnnuityGMBPrices> InsuranceAnnuityGMBPricesSP;

    class InsuranceAnnuityGMBPrices: public IMCPricesSimple {
    public:

        enum PricesIndices {
            OPTION_PRICE = 0,
            FEE_PRICE,
            NB_PRICES
        };

        /** adds supplied price to this set of IMCPrices */
        void add(double price, int index) {
            simplePrices[index]->add(price);
        }

        /** Returns the averaged result, together with the standard error */
        virtual void getResult(double* result, double* resultStdErr) const {
            simplePrices[OPTION_PRICE]->getResult(result, resultStdErr);
        }

        double getOptionPrice() const {
            double optionPrice, dummy;
            simplePrices[OPTION_PRICE]->getResult(&optionPrice, &dummy);
            return optionPrice;
        }

        double getFeePrice() const {
            double feePrice, dummy;
            simplePrices[FEE_PRICE]->getResult(&feePrice, &dummy);
            return feePrice;
        }

        /** Returns true if the path, identified by pathIdx, should be
        repriced. If false is returned then there will be no add()
        method called for this path and the IMCPrices object must
        take any appropriate action */
        virtual bool repriceForGreek(int pathIdx) {
            return true;
        }

        /** Support product-level caching across tweaks */
        virtual int storagePerPath(IMCProduct* product) const {
            return 0;
        }

        virtual void configureCache(const IntArray& changedAssets) {
        }

        /** Returns a deep copy of this object */
        IMCPrices* clone() const {
            InsuranceAnnuityGMBPricesSP copy(new InsuranceAnnuityGMBPrices());
            copy->simplePrices.resize(simplePrices.size());
            for (unsigned int iPrice = 0 ; iPrice < copy->simplePrices.size() ; iPrice++) {
                IMCPricesSP temp(this->simplePrices[iPrice]->clone());
                copy->simplePrices[iPrice] = DYNAMIC_POINTER_CAST<IMCPricesSimple>(temp);
            }
            return copy.get();
        }

        InsuranceAnnuityGMBPrices(int NbIter, int NbSubSamples):
        simplePrices(NB_PRICES){
            for (unsigned int iPrice = 0 ; iPrice < simplePrices.size() ; iPrice++){
                simplePrices[iPrice] = IMCPricesSimpleSP(new MCPricesSimple(NbIter, NbSubSamples));
            }
        }

    private:
        InsuranceAnnuityGMBPrices() {
        }

        /** adds supplied price to this set of IMCPrices */
        virtual void add(double price) {
            throw ModelException("IMCPrices::add(double)",
                     "internal error");
        }

        /** Returns the last price stored. Undefined behaviour if no
        prices yet stored */
        virtual double lastPrice() const {
            throw ModelException("IMCPrices::lastPrice",
                     "internal error");
        }

        /** On pricing run returns MAX(x, 0.0). It should be used for applying
        the 'final point of optionality' within the product. This allows
        QuickGreeks type of IMCPrices not to apply the max when doing first
        order greeks (this may sound strange but see the docs for why) */
        virtual double maxWithZero(double x) const {
            throw ModelException("IMCPrices::maxWithZero",
                     "internal error");
        }

        /** Reset this object so that it can be used for the same operation
        again. Normally, a new IMCPrices object is created for each
        pricing run. However, for quick x gamma, it is important to use
        the same one for each pair of assets */
        virtual void reset() {
            throw ModelException("IMCPrices::reset",
                     "internal error");
        }

        /** Ease cloning */
        virtual IMCPrices* emptyConstructor() const {
            throw ModelException("IMCPrices::emptyConstructor",
                     "internal error");
        }

        vector<IMCPricesSimpleSP> simplePrices;

    };

    /** calc payoff for each path in order to price diff group withdrawal inside payoff*/
    void calcPayoff(const IPathGenerator*  pathGen, const int iWDGroup);


    void calcDeterministicWDType3(const int& timeFromStart,                                                     
                                                    DoubleArraySP withdrawals,
                                                    const int& iGroup,
                                                    const int& iWDGroup 
                                                   );


public:
    
    virtual IMCPrices* createOrigPrices(int nbIter,
                             int nbSubSamples,
                             int mode) {
        return new InsuranceAnnuityGMBPrices(nbIter, nbSubSamples);
    }

    /** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen) const {
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    InsuranceAnnuityGMBProd(const InsuranceAnnuityGMB* inst,
                            IRefLevelSP refLevel,
                            const SimSeriesSP& simSeries) :
        IMCProduct(inst->asset.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  refLevel,
                  simSeries,
                  inst->spotSamples,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
                  inst(inst) {
            
        simulationDates = simSeries->getAllDates();

        histDates = inst->spotSamples->getDates();
        //histSpots = inst->spotSamples->getValues();

        // compute growth factors from payment dates to last date
        DateTime lastDate = histDates[histDates.size()-1];

        growthFactor.resize(histDates.size());
        
        int iDate;
        for (iDate = 0 ; iDate < histDates.size() ; iDate++)
        {
            DateTime settlementDate = inst->instSettle->settles(histDates[iDate], 0);
            if (inst->valueDate <= settlementDate)
                growthFactor[iDate] = 1.0 / inst->instSettle->pv(settlementDate,  
                                     lastDate,
                                     inst->discount.get(),
                                     inst->asset.get());
            else
                growthFactor[iDate] = 9999.9; // this number should never be used
        }

        int nberGroup = inst->groupSizeAtStart.size();        
        int nbWDGroupPercentage = 1;
        int iWD ;

        if ( (inst->isGMWB) && 
            ( (inst->deterministicWDType == 2) //withdraw depends on policy year
            || (inst->deterministicWDType == 3) ) ) //withdraw depends on policy year and ITM levels
        {
            nbWDGroupPercentage = inst->WDGroupPercentage.size();
            //may start withdraw later

            isWithdraw.resize(nbWDGroupPercentage);
            isWithdrawSoFar.resize(nbWDGroupPercentage);
            for (iWD = 0 ; iWD < nbWDGroupPercentage ; iWD++){
                isWithdraw[iWD].resize(nberGroup,false);
                isWithdrawSoFar[iWD].resize(nberGroup,false);
            }
        }else{

            //isWithdraw.resize(nbWDGroupPercentage, true); //start withdraw from the beginning
            isWithdraw.resize(nbWDGroupPercentage);
            isWithdrawSoFar.resize(nbWDGroupPercentage);
            for (iWD = 0 ; iWD < nbWDGroupPercentage ; iWD++){
                isWithdraw[iWD].resize(nberGroup,true);
                isWithdrawSoFar[iWD].resize(nberGroup,true);
            }
        }

        // initialize path value
        pathValue = 0.0;
        pathValueSoFar = 0.0;

        // initialize present value of the fees
        feeValue = 0.0;
        feeValueSoFar = 0.0;

        // create and initialize uniform random generator
        uniRand = RandUniformDefaultSP(new RandUniformDefault());

        int iGroup;

        // resize and initialize account value
        accountValue.resize(nbWDGroupPercentage);
        accountValueSoFar.resize(nbWDGroupPercentage);
        for (iWD = 0 ; iWD < nbWDGroupPercentage ; iWD++){
            accountValueSoFar[iWD].resize(nberGroup, 1.0);
            accountValue[iWD].resize(nberGroup,1.0);
        }
       
        // initialize isWaitPeriod flag
        isWaitPeriod .resize(nbWDGroupPercentage, true);
        isWaitPeriodSoFar.resize(nbWDGroupPercentage, true);

        // resize and initialize last step-up flag
        lastStepUp.resize(nbWDGroupPercentage);
        lastStepUpSoFar.resize(nbWDGroupPercentage);
        for (iWD = 0 ; iWD < nbWDGroupPercentage ; iWD++){
            lastStepUpSoFar[iWD].resize(nberGroup,-1);
            lastStepUp[iWD].resize(nberGroup,-1);
//            lastStepUpSoFar[iWD].resize(nberGroup,0); //to review
//            lastStepUp[iWD].resize(nberGroup,0);
        }

        // resize and initialize notionals and benefits for GMDB and GMAB
        notionalGMADB.resize(nberGroup);
        notionalGMADBSoFar.resize(nberGroup);
        benefitGMADB.resize(nberGroup);
        benefitGMADBSoFar.resize(nberGroup);

        for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
        {
            notionalGMADB[iGroup] = (inst->isGMDB || inst->isGMAB) ? 1.0 : 0.0;
            notionalGMADBSoFar[iGroup] = (inst->isGMDB || inst->isGMAB) ? 1.0 : 0.0;
            benefitGMADB[iGroup] = (inst->isGMDB || inst->isGMAB) ? 1.0 : 0.0;
            benefitGMADBSoFar[iGroup] = (inst->isGMDB || inst->isGMAB) ? 1.0 : 0.0;
        }

        // resize and initialize benefits and maximum withdrawal for GMWB
        // initialize and resize age group sizes
        benefitGMWB.resize(nbWDGroupPercentage);
        lastMaxGMWBWithdrawal.resize(nbWDGroupPercentage);
        benefitGMWBSoFar.resize(nbWDGroupPercentage);
        lastMaxGMWBWithdrawalSoFar.resize(nbWDGroupPercentage);
        groupSizeSoFar.resize(nbWDGroupPercentage);
        groupSize.resize(nbWDGroupPercentage);

        for (iWD = 0 ; iWD < nbWDGroupPercentage ; iWD++){
            groupSize[iWD].resize(nberGroup);
            groupSizeSoFar[iWD].resize(nberGroup);

            benefitGMWBSoFar[iWD].resize(nberGroup);
            benefitGMWB[iWD].resize(nberGroup);
            lastMaxGMWBWithdrawal[iWD].resize(nberGroup);
            lastMaxGMWBWithdrawalSoFar[iWD].resize(nberGroup);           
        }
        
        for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
        {
            for (iWD = 0 ; iWD < nbWDGroupPercentage ; iWD++){
                groupSize[iWD][iGroup] = inst->groupSizeAtStart[iGroup] 
                                        * inst->WDGroupPercentage[iWD]; //WDGroupPercentage[0] = 1 except deterministicWDType =2 or 3
                benefitGMWB[iWD][iGroup] = inst->isGMWB ? 1.0: 0.0;
                lastMaxGMWBWithdrawal[iWD][iGroup] = 
                    inst->isGMWB ? inst->maxGMWBPercentage : 0.0; // account value = 1.0 initially

                groupSizeSoFar[iWD][iGroup] = groupSize[iWD][iGroup];
                benefitGMWBSoFar[iWD][iGroup]= benefitGMWB[iWD][iGroup];
                lastMaxGMWBWithdrawalSoFar[iWD][iGroup]= lastMaxGMWBWithdrawal[iWD][iGroup] ;
            }
        }        


        //for Debug output
        if (inst->DEBUG_OUTPUT){            
            db_pathCount = -1; //0;
            int size = histDates.size();

            db_index.resize(size);

            db_accountValueBeg.resize(size); 
            db_benefitBeg.resize(size);
            db_groupSizeBeg.resize(size);

            db_accountValue.resize(size); 
            db_benefit.resize(size);
            db_currentWithdrawals.resize(size);
            //db_lastStepUp.resize(size);
            db_groupSize.resize(size);
            db_payoff.resize(size);
            db_negativeAV.resize(size);

            if (inst->DEBUG_OUTPUT_ESTIMATE){            

                db_Mindex.resize(histDates.size());            
                for(int i = 0; i <histDates.size(); i++ ){
                    db_Mindex[i].resize(inst->db_maxIte, 0.0);
                }

                db_benefitAtT.resize(inst->db_maxIte, 0.0);
                db_accountValueAtT.resize(inst->db_maxIte, 0.0);
                db_nonLapseAtT.resize(inst->db_maxIte, 0.0);
            }
        }
    }

    // compute fees (supposed to be paid at the begining of the period)
    void computeFees(double* fees, const int& timeFromStart, const int& iWDGroup);

    // update step up
    void updateStepUp(const int& timeFromStart, const int& iWDGroup);

    // compute coupons taken by clients
    void computeWithdrawals(const int& timeFromStart,
                            DoubleArraySP withdrawals,
                            const int& iWDGroup);

    // update account value
    void updateAccountValue(DoubleArraySP withdrawals, bool alreadySettled, const int& iWDGroup);

    // update benefits Guaranteed Minimum Death Benefit (GMDB)
    // update benefits Guaranteed Minimum Death Benefit (GMAB) // not used
    void updateBenefitGMADB(const int& timeFromStart, const double& performance);

    // update notionals for Guaranteed Minimum Death Benefit (GMDB)
    // update notionals for Guaranteed Minimum Accumulation Benefit (GMAB)
    void updateNotionalGMADB(DoubleArraySP withdrawals);

    // update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)
    void updateBenefitGMWB(DoubleArraySP withdrawals, const int& iWDGroup);

    // update number of clients by taking into account the lapse rate
    void updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup);

    // update number of clients by taking into account the mortality rate
    void updateDeadClients(const int& timeFromStart, const int& iWDGroup);

    void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices);

        /** print extra output **/
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices&) const;
    
        /** vol interp for LN */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;

};
// end of class InsuranceAnnuityGMBProd

// compute fees (supposed to be paid at the begining of the period)
void InsuranceAnnuityGMBProd::computeFees(double* fees, const int& timeFromStart, const int& iWDGroup) {

    int iGroup;
    for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size(); iGroup++)
    {
        *fees +=  groupSize[iWDGroup][iGroup] * accountValue[iWDGroup][iGroup] * inst->feeRate *
            growthFactor[timeFromStart]; // the fees are forward valued to last date
    }
}

// decide if the client steps up or not and update benefit
void InsuranceAnnuityGMBProd::updateStepUp(const int& timeFromStart, const int& iWDGroup){
    int iGroup;
    int nberGroup = inst->groupSizeAtStart.size();
    
    if (isWaitPeriod[iWDGroup] && (timeFromStart >= inst->initialWaitPeriod)) 
    {
        isWaitPeriod[iWDGroup] = false;

        for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
        {
            // authorize step-up from the end of the initial waiting period
            lastStepUp[iWDGroup][iGroup] = inst->stepUpPeriod;
        }
    }

    if (inst->isGMWB && inst->isStepUp) // step-ups only supported by GMWB benefit
    {
        for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
        {
            if (!isWaitPeriod[iWDGroup] && (lastStepUp[iWDGroup][iGroup] >= inst->stepUpPeriod))   /// to review

            {
                if (!Maths::equals(0.0, benefitGMWB[iWDGroup][iGroup]))
                {                                                        
                    // if step ups are static, compare to stepUpITMLevel
                    if (!inst->isDynamicStepUp)
                    {
                        if (accountValue[iWDGroup][iGroup] >= inst->stepUpITMLevel * benefitGMWB[iWDGroup][iGroup])
                        {
                            benefitGMWB[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup];
                            lastMaxGMWBWithdrawal[iWDGroup][iGroup] =
                                max(lastMaxGMWBWithdrawal[iWDGroup][iGroup], inst->maxGMWBPercentage * accountValue[iWDGroup][iGroup]);
                            // set lastStepUp[iGroup] to -1 since lastStepUp[iGroup]
                            // is going to be incremented at the end of the loop
                            //lastStepUp[iWDGroup][iGroup] = 0;  //to review
                            lastStepUp[iWDGroup][iGroup] = -1;

                        }
                    }
                    // if step ups are dynamic, use formula for step up utilization
                    else
                    {
                        // compute the step up ratio
                        double stepUpRatio = 0.0; // ratio of clients who steps up
                        if (!Maths::isZero(benefitGMWB[iWDGroup][iGroup]))
                        {
                            stepUpRatio = 2.0 * ((accountValue[iWDGroup][iGroup] / benefitGMWB[iWDGroup][iGroup]) - 1.0);
                            stepUpRatio = min(max(stepUpRatio, 0.0), 1.0);
                        }

                        // determine if the group is going to step up or not
                        double uniform; // uniform variable
                        uniRand->draw(1, &uniform);
                        bool stepUpUtilization = (uniform < stepUpRatio); // flag indicating if the group steps up or not

                        if (stepUpUtilization)
                        {
                            benefitGMWB[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup];
                            lastMaxGMWBWithdrawal[iWDGroup][iGroup] =
                                max(lastMaxGMWBWithdrawal[iWDGroup][iGroup], inst->maxGMWBPercentage * accountValue[iWDGroup][iGroup]);
                            // set lastStepUp[iGroup] to -1 since lastStepUp[iGroup]
                            // is going to be incremented at the end of the loop
                            lastStepUp[iWDGroup][iGroup] = -1;                                }
                    }
                }
            }
        }
    }

    for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
    {
        ++lastStepUp[iWDGroup][iGroup];
    }
}


// update number of clients by taking into account the lapse rate
void InsuranceAnnuityGMBProd::updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup) {

    double lapseRate = 0.0; // base lapse rate (the same for all clients)
    double dynamicLapseRate = 0.0; // same as base lapse rate if no dynamic assumptions
    double lapseRateMultiplier = 0.0; // lapse rate multiplier if dynamic lapse rate
    double floorLapseRate = 0.0; //inst->isDynamicLapseRate ? 0.02 : 0.0; // floor for lapse rate if dynamic lapse rate

    if (inst->lapseOption && (timeFromStart != 0)) // timeFromStart = 0 at start date
    {
        lapseRate = inst->lapseRate[timeFromStart-1]; // base lapse rate depending on the time from start date only
        double ITMLevel; // in the money level

        int iGroup;
        for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size(); iGroup++)
        {
            // ITM level for GMAB GMWB and GMDB benefits if dynamic assumptions for lapse rate
            if (inst->isDynamicLapseRate)
            {
                if (inst->lapseRateMultiplierType == 1){ //old version
                    
                    floorLapseRate = 0.02; // floor for lapse rate if dynamic lapse rate

                    // account value = 0.0, ITM level is set to 1.0
                    if (!Maths::isPositive(accountValue[iWDGroup][iGroup]))
                    {
                        ITMLevel = 1.0;
                    }
                    // account value != 0.0 and GMAB or GMWB benefit
                    else if ((inst->isGMDB || inst->isGMAB) && Maths::isPositive(accountValue[iWDGroup][iGroup]))
                    {
                        ITMLevel = (notionalGMADB[iGroup] / accountValue[iWDGroup][iGroup]) - 1.0;
                    }
                    // account value != 0.0 and GMWB benefit
                    else
                    {
                        ITMLevel = (benefitGMWB[iWDGroup][iGroup] / accountValue[iWDGroup][iGroup]) - 1.0;
                    }

                    // the lapse rate at a floor set to 2% (when lapse rate is dynamic only)
                    lapseRateMultiplier = pow(max((1.0 - 2.0 * max(0.0, ITMLevel)), 0.000001), 1.5);

                }else if((inst->lapseRateMultiplierType == 2) && (inst->isGMWB)){ //only available for GMWB

                    //floorLapseRate = 0.0; //reset floor to zero for this contract



                    double ITMAv_ben = 0.0;

                    if (Maths::isPositive(benefitGMWB[iWDGroup][iGroup])){
                        ITMAv_ben = accountValue[iWDGroup][iGroup] / benefitGMWB[iWDGroup][iGroup] /inst->adjustmentAV;
                    }
                    lapseRateMultiplier = Maths::min(1.0, ITMAv_ben);

                    //applying additional multiplier depends on the current status of withdrawal.
                    if ((*withdrawals)[iGroup] > 0.0){ //current withdraw status
                        lapseRateMultiplier *= inst->withdrawMultiplier; //only active if ITM > 0
                    }else{
                        lapseRateMultiplier *= inst->noWithdrawMultiplier; //only active if ITM > 0
                    }
                }else{

                }
            }
            // ITM level set to 0.0 for GMAB GMWB and GMDB benefits if no dynamic assumptions for lapse rate
            else
            {
                lapseRateMultiplier = 1.0; //just use the basic lapse rate
            }

            dynamicLapseRate = max(lapseRateMultiplier * lapseRate, floorLapseRate);
            groupSize[iWDGroup][iGroup] *= (1.0 - dynamicLapseRate);
        }
    }
}

// update number of clients by taking into account the mortality rate
void InsuranceAnnuityGMBProd::updateDeadClients(const int& timeFromStart, const int& iWDGroup) {

    double mortalityRate = 0.0; // mortality rate

    if (inst->dieOption && (timeFromStart != 0))
    {
        int iGroup;
        for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size() ; iGroup++)
        {   
            // current age of the group
            int currentAge = inst->ageAtStart[iGroup] + timeFromStart;
            mortalityRate = inst->mortalityRate[currentAge-1]; // mortality rate depending on the current age of the group

            groupSize[iWDGroup][iGroup] *= (1.0 - mortalityRate);
        }
    }
}

// update benefits for Guaranteed Minimum Death Benefit (GMDB)
// update benefits Guaranteed Minimum Death Benefit (GMAB) // not used
void InsuranceAnnuityGMBProd::updateBenefitGMADB(const int& timeFromStart, const double& performance) {

    if (!inst->isAnnualRatchet && !inst->isPercentageRollUp)
    // do nothing
    return;

    int iGroup;
    for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size(); iGroup++)
    {
        // current age of the group
        int currentAge = inst->ageAtStart[iGroup] + timeFromStart;

        // no increase after the benefir age limit
        if (currentAge <= inst->GMDBAgeCap)
        {
            if (inst->isAnnualRatchet)
                // take the max with current performance
                benefitGMADB[iGroup] = max(benefitGMADB[iGroup], performance);

            if (inst->isPercentageRollUp)
            {
                // take the max with the percentage roll-up
                double rollup = ::pow((1.0 + inst->percentageRollUp), timeFromStart);
                benefitGMADB[iGroup] = max(benefitGMADB[iGroup], rollup);
            }
        }
    }
}

// update notionals for Guaranteed Minimum Death Benefit (GMDB)
// update notionals for Guaranteed Minimum Accumulation Benefit (GMAB)
void InsuranceAnnuityGMBProd::updateNotionalGMADB(DoubleArraySP withdrawals) {

    int iGroup;
    for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size();  iGroup++)
    {
        notionalGMADB[iGroup] *= (1.0 -  (*withdrawals)[iGroup]);
    }
}

// update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)
void InsuranceAnnuityGMBProd::updateBenefitGMWB(DoubleArraySP withdrawals, const int& iWDGroup) {

    int iGroup;
    for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size();  iGroup++)
    {
        benefitGMWB[iWDGroup][iGroup] -= (*withdrawals)[iGroup];

        if (Maths::isNegative(benefitGMWB[iWDGroup][iGroup]))
            // should never reach this line
            benefitGMWB[iWDGroup][iGroup] = 0.0;
    }
}

// compute coupons taken by clients
void InsuranceAnnuityGMBProd::computeWithdrawals(const int& timeFromStart,                                                 
                                                DoubleArraySP withdrawals,
                                                const int& iWDGroup
                                                 ) {
        
    int iGroup;
    for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size();  iGroup++)
    {       
        if (inst->isGMDB || inst->isGMAB)
        {
            if (inst->isGMDB)
                // for GMDB or GMDB + GMAB, withdrawals expressed as percentage of the current account value
                (*withdrawals)[iGroup] = (timeFromStart == 0) ? 0.0 : inst->GMDBWithdrawRate;
            else
                // for GMAB, withdrawals expressed as percentage of the current account value
                (*withdrawals)[iGroup] = (timeFromStart == 0) ? 0.0 : inst->GMABWithdrawRate;
        }
        else
        { //for GMWB, different group withdrawal.

            //isWithdraw is flag indicating if the withdrawal takes place or not
            // always true when no dynamic assumption for utilization
            if (inst->isDynamicUtilization) //stochastic withdraw
            {
                // if the current GMWB benefit is not zero
                // determine if withdrawal happens using a uniform variable
                if (!Maths::isZero(accountValue[iWDGroup][iGroup]))
                {
                    // compute the in the money level
                    double ITMLevel = (benefitGMWB[iWDGroup][iGroup] / accountValue[iWDGroup][iGroup]) - 1.0;
                    // compute the dynamic utilization factor
                    double utilizationFactor = 0.4 * pow(max((1.0 + 2.0 * ITMLevel), 0.0), 1.5);
                    utilizationFactor = min(max(utilizationFactor, 0.2), 1.0);

                    double uniform; // uniform variable
                    uniRand->draw(1, &uniform);

                    isWithdraw[iWDGroup][iGroup] = (uniform < utilizationFactor);
                }
                // if the current GMWB benefit is not zero, the withdrawal always happens
                else
                {
                    isWithdraw[iWDGroup][iGroup] = true;
                }
            }else{//deterministic withdraw, but may depends on ITM 
                double timeToWithdrawal;
                               
                switch (inst->deterministicWDType) {
                    case 1:{
                        //withdraw from beginning
                    }        
                    break;
                    case 2: { //timeToWithdrawal didn't depend on ITMLevel, iWDITMLevels ==0
                        if (!isWithdraw[iWDGroup][iGroup]){//haven't start to withdraw
                            timeToWithdrawal = (*inst->startWDYears)[0][iWDGroup];  
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
                ((timeFromStart == 0) ? 0.0 : max(min(lastMaxGMWBWithdrawal[iWDGroup][iGroup], benefitGMWB[iWDGroup][iGroup]), 0.0)) ;
        }
    }
}

void InsuranceAnnuityGMBProd::calcDeterministicWDType3(const int& timeFromStart, 
                                                DoubleArraySP withdrawals,
                                                const int& iGroup,
                                                const int& iWDGroup 
                                               ){

    double ITMLevel;

    if ( (Maths::isNegative(accountValue[iWDGroup][iGroup])) || (Maths::isZero(accountValue[iWDGroup][iGroup])))
    {
        ITMLevel = 100000000000000.0;
    }else{
        ITMLevel = inst->adjustmentAV * (benefitGMWB[iWDGroup][iGroup] / accountValue[iWDGroup][iGroup]);
    }      

    int iFirst;
    int nbWDITMLevel = inst->WDITMLevels.size();

    //looking for iFirst
    iFirst = nbWDITMLevel -1;
    while ((iFirst >= 0) && (timeFromStart >= (*inst->startWDYears)[iFirst][iWDGroup])){
        iFirst -= 1;
    }
    iFirst = (iFirst == nbWDITMLevel -1)? iFirst: (iFirst + 1);

    double timeStartWD = (*inst->startWDYears)[iFirst][iWDGroup];
    if(timeFromStart >= timeStartWD){ //check ITM cond only if timeFromStart >= scheduled WD year
        //check the ITM levels
        if (iFirst ==0){
            isWithdraw[iWDGroup][iGroup] = true;                          
        }else{
            isWithdraw[iWDGroup][iGroup] = (ITMLevel >= inst->WDITMLevels[iFirst-1]);
        }
    }


    //looking for iLast
//        int iLast = nbWDITMLevel -1;
//        for (int i = nbWDITMLevel -1; i>=0; i--){
//            if (( (*inst->startWDYears)[i][iWDGroup] <= timeFromStart) 
//                && ((*inst->startWDYears)[i][iWDGroup] > (*inst->startWDYears)[iLast][iWDGroup] ) ){
//                iLast = i;
//            }
//        }
//                            
//        double timeStartWD = (*inst->startWDYears)[iFirst][iWDGroup];
//        if(timeFromStart >= timeStartWD){ //check ITM cond only if timeFromStart >= scheduled WD year
//            //check the ITM levels
//            if (iFirst ==0){
//                isWithdraw[iWDGroup][iGroup] = (ITMLevel < inst->WDITMLevels[iLast]);                          
//            }else{
//                isWithdraw[iWDGroup][iGroup] = ((ITMLevel < inst->WDITMLevels[iLast]) 
//                                && ((ITMLevel >= inst->WDITMLevels[iFirst-1])));
//            }
//        }

}

// update account value
void InsuranceAnnuityGMBProd::updateAccountValue(DoubleArraySP withdrawals, 
                                                 bool alreadySettled,
                                                 const int& iWDGroup) {
        
    int iGroup;
    for (iGroup = 0 ; iGroup < inst->groupSizeAtStart.size();   iGroup++)
    {
        if (inst->isGMDB || inst->isGMAB)
            // for GMDB and GMAB, withdrawals expressed as percentage of the current account value
            accountValue[iWDGroup][iGroup] *= (1.0 - (*withdrawals)[iGroup]);

        else
        {
            // for GMWB, withdrawals expressed as a percentage of the initial amount
            accountValue[iWDGroup][iGroup] -= (*withdrawals)[iGroup];
            if (alreadySettled)
            {
                // if the account value is inferior to 0 in the past, we set it at 0
                if (Maths::isNegative(accountValue[iWDGroup][iGroup]))
                    accountValue[iWDGroup][iGroup] = 0.0;

                // if the GMWB benefit is inferior or equal to 0.0, we set the account value at 0.0
                if (!Maths::isPositive(benefitGMWB[iWDGroup][iGroup]))
                    accountValue[iWDGroup][iGroup] = 0.0;
            }
        }
    }
}

void InsuranceAnnuityGMBProd::payoff(const IPathGenerator*  pathGen, IMCPrices& prices) {

    static const string routine("InsuranceAnnuityGMBProd::payoff");
    try {
        InsuranceAnnuityGMBPrices& myprices = static_cast<InsuranceAnnuityGMBPrices&>(prices);

        const double* path = pathGen->Path(0,0); // access path

        int iGroup;
        int iWDGroup ;

        // set a flag isDoingPast
        bool doingPast = pathGen->doingPast();

        //int nberGroup = groupSize.size();
        int nberGroup = inst->groupSizeAtStart.size();        
        int nbWDGroupPercentage = 1; //1 if is not GMWB and deterministicWDType 2 or 3

        if ((inst->isGMWB) && 
            ((inst->deterministicWDType == 2) || 
            (inst->deterministicWDType == 3) ) )
        {
            nbWDGroupPercentage = inst->WDGroupPercentage.size(); 
        }

        // preservation of the past
        if (!doingPast)
        {
            for (iWDGroup = 0; iWDGroup < nbWDGroupPercentage; iWDGroup++){

                isWaitPeriod[iWDGroup] = isWaitPeriodSoFar[iWDGroup]; // isWaitingPeriod flag

                for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
                {
                    pathValue = 0.0; //pathValueSoFar; // path value
                    feeValue = 0.0; //feeValueSoFar; // present value of the fees
                    accountValue[iWDGroup][iGroup] = accountValueSoFar[iWDGroup][iGroup]; // current account value
                    lastStepUp[iWDGroup][iGroup] = lastStepUpSoFar[iWDGroup][iGroup]; // nber of years since the last step-up for the different groups
                    notionalGMADB[iGroup] = notionalGMADBSoFar[iGroup]; // garanteed minimum death and accumulation benefit notional for the different age groups
                    benefitGMADB[iGroup] = benefitGMADBSoFar[iGroup]; // garanteed minimum death benefit for the different age groups
                    benefitGMWB[iWDGroup][iGroup] = benefitGMWBSoFar[iWDGroup][iGroup]; // garanteed minimum withdrawal benefit for the different age groups                                
                    lastMaxGMWBWithdrawal[iWDGroup][iGroup] = lastMaxGMWBWithdrawalSoFar[iWDGroup][iGroup]; // last maximum GMWB withdrawal
                    groupSize[iWDGroup][iGroup] = groupSizeSoFar[iWDGroup][iGroup];
                    isWithdraw[iWDGroup][iGroup]=isWithdrawSoFar[iWDGroup][iGroup];
                }
            }
        }

        //GMWB, loop on different withdrawal groups

        //debug output, count the nb of simulation        
        if (inst->DEBUG_OUTPUT){            
            db_pathCount += 1;
        }

        for (iWDGroup = 0; iWDGroup < nbWDGroupPercentage; iWDGroup++){
            calcPayoff(pathGen, iWDGroup);
        }

        // preservation of the past 
        if (doingPast)
        {
            for (iWDGroup = 0; iWDGroup < nbWDGroupPercentage; iWDGroup++){

                isWaitPeriodSoFar[iWDGroup] = isWaitPeriod[iWDGroup]; // isWaitingPeriod flag

                for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
                {
                    pathValueSoFar = pathValue; // path value
                    feeValueSoFar = feeValue; // present value of the fees
                    accountValueSoFar[iWDGroup][iGroup] = accountValue[iWDGroup][iGroup]; // current account value;
                    notionalGMADBSoFar[iGroup] = notionalGMADB[iGroup]; // garanteed minimum death and accumulation benefit notional for the different age groups
                    benefitGMADBSoFar[iGroup] = benefitGMADB[iGroup]; // garanteed minimum death benefit for the different age groups
                    benefitGMWBSoFar[iWDGroup][iGroup] = benefitGMWB[iWDGroup][iGroup]; // garanteed minimum withdrawal benefit for the different age groups
                    lastMaxGMWBWithdrawalSoFar[iWDGroup][iGroup] = lastMaxGMWBWithdrawal[0][iGroup]; // last maximum GMWB withdrawal
                    groupSizeSoFar[iWDGroup][iGroup] = groupSize[iWDGroup][iGroup];
                    isWithdrawSoFar[iWDGroup][iGroup]= isWithdraw[iWDGroup][iGroup];
                    lastStepUpSoFar[iWDGroup][iGroup] =lastStepUp[iWDGroup][iGroup]; // nber of years since the last step-up for the different groups
                }
            }
        }


        if (!doingPast)
        {
            myprices.add((inst->initialAmount * pathValue), InsuranceAnnuityGMBPrices::OPTION_PRICE); // option price
            myprices.add((inst->initialAmount * feeValue), InsuranceAnnuityGMBPrices::FEE_PRICE); // fee price
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


void InsuranceAnnuityGMBProd::calcPayoff(const IPathGenerator*  pathGen, const int iWDGroup){

    static const string routine("InsuranceAnnuityGMBProd::calcPayoff");
    try {

        const double* path = pathGen->Path(0,0); // access path

        int beginIdx = pathGen->begin(0);
        int endIdx = pathGen->end(0);

        int iGroup;
        int nberGroup = inst->groupSizeAtStart.size();        

        DateTime settlementDate; // settlement date for the current start date anniversary
        DateTime valueDate = inst->valueDate; // value date

        double indexValue; // index value (equal to 1$ at start date)
        DoubleArraySP currentWithdrawals = DoubleArraySP(new DoubleArray(nberGroup)); // current coupon taken by clients
        //bool isCurrentWithdrawal = false; //additional lapse rate Multiplier depends on the current status of withdrawal

        //for debug output
        int db_iWDG = (inst->DEBUG_OUTPUT )? inst->db_ithWDGroup - 1 : 0;

        int Idx;

        for (Idx = beginIdx ; Idx < endIdx ; Idx++)
        {
            //for debug output
            if (inst->DEBUG_OUTPUT ){ // not doing past

                bool doingPast = pathGen->doingPast();
              
                if (doingPast || ( (db_pathCount == inst->db_ithPath) && (db_iWDG ==iWDGroup) ) ){
                    int db_iAgeG = inst->db_ithAgeGroup - 1;
                    db_accountValueBeg[Idx] = accountValue[db_iWDG][db_iAgeG];
                    db_benefitBeg[Idx] = benefitGMWB[db_iWDG][db_iAgeG];
                    db_groupSizeBeg[Idx] = groupSize[db_iWDG][db_iAgeG];
                    db_index[Idx] = path[Idx]; //output path for tests purpose
                }

                if (inst->DEBUG_OUTPUT_ESTIMATE ){ // not doing past
                    if(doingPast || db_pathCount < inst->db_maxIte){
                        db_Mindex[Idx][db_pathCount] = path[Idx]; //output path for tests purpose
                    }else{
                        throw ModelException(routine, "nb of iteration is bigger than maximum iteration allowed for DEBUG_OUTPUT_ESTIMATE.");
                    }
                }
            }

            // update index value
            indexValue = 1.0 * (path[Idx] / path[0]);

            // update settlement date
            settlementDate = inst->instSettle->settles(histDates[Idx], 0);
            bool alreadySettled = valueDate.isGreater(settlementDate);

            // update account value (before coupon is paid)
            for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
            {
                if (Idx != 0) // start date
                    accountValue[iWDGroup][iGroup] *= (path[Idx] / path[Idx-1]);
            }

            // compute coupons taken by clients
            // calculated as percentage of the current account value for GMDB
            // calculated as percentage of the initial amount for GMWB
            // calculated as percentage of the current account value for GMAB
            computeWithdrawals(Idx, currentWithdrawals, iWDGroup);

            //update step up for GMWB
            updateStepUp(Idx, iWDGroup);

            if (inst->isGMDB || inst->isGMAB)
            {
                if (inst->isGMDB)
                {
                    // update benefits for Guaranteed Minimum Death Benefit (GMDB)
                    updateBenefitGMADB(Idx, indexValue);
                }

                // update notionals for Guaranteed Minimum Death Benefit (GMDB)
                // update notionals for Guaranteed Minimum Death Benefit (GMAB)
                updateNotionalGMADB(currentWithdrawals);
                // update account value
                updateAccountValue(currentWithdrawals, alreadySettled, iWDGroup);

                if (alreadySettled)
                {
                    for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
                    {
                        // no benefit after GMDB age limit for GMDB and GMDB + GMAB
                        if (inst->isGMDB)
                        {
                            int groupAge = inst->ageAtStart[iGroup] + Idx;
                            if (groupAge >= inst->GMDBAgeLimit)
                            {
                                notionalGMADB[iGroup] = 0.0; // no accumulation benefit after age limit
                                accountValue[iWDGroup][iGroup] = 0.0;
                                groupSize[iWDGroup][iGroup] = 0.0;
                            }
                        }
                        // no benefit after last GMAB maturity for GMAB
                        else
                        {
                            if (Idx == inst->GMABMaturities[inst->GMABMaturities.size() - 1])
                            {
                                notionalGMADB[iGroup] = 0.0; // no accumulation benefit after maturity
                                accountValue[iWDGroup][iGroup] = 0.0;
                                groupSize[iWDGroup][iGroup] = 0.0;
                            }
                        }
                    }
                }
            }
            else
            {
                // update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)
                updateBenefitGMWB(currentWithdrawals, iWDGroup);
                // update account value
                updateAccountValue(currentWithdrawals, alreadySettled, iWDGroup);
            }

            // compute fees (supposed to be paid at the begining of the period)
            if (!alreadySettled && (Idx !=0))
                computeFees(&feeValue, Idx, iWDGroup);

            // remove the clients which have lapsed during the year
            // additional lapse rate Multiplier may depend on the current status of withdrawal
            updateLapseClients(Idx, currentWithdrawals, iWDGroup); 

            if (!alreadySettled && (Idx != 0))
            {
                // add value of Guaranteed Minimum Death Benefit (GMDB) option
                for (iGroup = 0 ; iGroup < nberGroup ; iGroup++)
                {
                    // Guaranteed Minimum Death Benefit (GMDB) or / and
                    // Guaranteed Minimum Accumulation Benefit (GMAB)
                    if (inst->isGMDB || inst->isGMAB)
                    {
                        int groupAge = inst->ageAtStart[iGroup] + Idx;
                        double mortalityRate = inst->dieOption ? inst->mortalityRate[groupAge-1] : 0.0;

                        // Guaranteed Minimum Accumulatiom Benefit (GMAB)
                        if (inst->isGMAB)
                        {
                            int iMat;
                            for (iMat = 0 ; iMat < inst->GMABMaturities.size() ; iMat++)
                            {
                                if (Idx == inst->GMABMaturities[iMat])
                                {
                                    double benefit = notionalGMADB[iGroup] * inst->GMABPercentages[iMat];

                                    pathValue += growthFactor[Idx] *  // the fees are forward valued to last date
                                        groupSize[iWDGroup][iGroup] * 
                                        max(benefit - accountValue[iWDGroup][iGroup], 0.0);

                                    if (Maths::isPositive(benefit - accountValue[iWDGroup][iGroup]))
                                        accountValue[iWDGroup][iGroup] = benefit;
                                }
                            }

                            // no benefit after last GMAB maturity for GMAB
                            if (!inst->isGMDB && (Idx == inst->GMABMaturities[inst->GMABMaturities.size() - 1]))
                            {
                                notionalGMADB[iGroup] = 0.0; // no accumulation benefit after last GMAB maturity
                                accountValue[iWDGroup][iGroup] = 0.0;
                                groupSize[iWDGroup][iGroup] = 0.0;
                            }
                        }

                        // Guaranteed Minimum Death Benefit (GMDB)
                        if (inst->isGMDB)
                        {
                            pathValue += growthFactor[Idx] *  // the fees are forward valued to last date
                                groupSize[iWDGroup][iGroup] * mortalityRate *
                                max(notionalGMADB[iGroup] * benefitGMADB[iGroup] - accountValue[iWDGroup][iGroup], 0.0);

                            // no benefit after GMDB age limit for GMDB and GMDB + GMAB
                            if (groupAge >= inst->GMDBAgeLimit)
                            {
                                notionalGMADB[iGroup] = 0.0; // no accumulation benefit after age limit
                                accountValue[iWDGroup][iGroup] = 0.0;
                                groupSize[iWDGroup][iGroup] = 0.0;
                            }
                        }
                    }
                    else
                    {
                    if (Maths::isNegative(accountValue[iWDGroup][iGroup]))
                    {
                        pathValue +=  growthFactor[Idx] *  // the fees are forward valued to last date
                            groupSize[iWDGroup][iGroup] * (-accountValue[iWDGroup][iGroup]);
                        
                        // for output
                        if (inst->DEBUG_OUTPUT){
//                            if( (db_pathCount == inst->db_ithPath)){
                            if (( db_pathCount == inst->db_ithPath) && (db_iWDG ==iWDGroup) ){
                                db_negativeAV[Idx] = accountValue[iWDGroup][iGroup];
                                db_payoff[Idx] = groupSize[iWDGroup][iGroup] * (-accountValue[iWDGroup][iGroup]);
                            }
                        }

                        accountValue[iWDGroup][iGroup] = 0.0;
                    }
                    // if the GMWB benefit is inferior to 0.0, we set the account value at 0.0
                    if (!Maths::isPositive(benefitGMWB[iWDGroup][iGroup]))
                        accountValue[iWDGroup][iGroup] = 0.0;
                    }
                }
            }

            // remove the clients which are dead during the year
            // (we assume that the clients die at the end of the year)
            updateDeadClients(Idx, iWDGroup);

            //for debug output
            if (inst->DEBUG_OUTPUT ){ // not doing past
//                if( db_pathCount == inst->db_ithPath){
                if ( pathGen->doingPast() || (( db_pathCount == inst->db_ithPath) && (db_iWDG ==iWDGroup)) ){
                    int db_iAgeG = inst->db_ithAgeGroup - 1;
                    db_accountValue[Idx] = accountValue[db_iWDG][db_iAgeG];
                    db_benefit[Idx] = benefitGMWB[db_iWDG][db_iAgeG];
                    //db_isWithdraw[Idx] = isWithdraw[db_ithWDGroup][db_ithAgeGroup]
                    db_currentWithdrawals[Idx] = (*currentWithdrawals)[db_iAgeG]; 
                    //db_lastStepUp[Idx]  = lastStepUp[db_iWDG][db_iAgeG];
                    db_groupSize[Idx] = groupSize[db_iWDG][db_iAgeG];
                }

                if(inst->DEBUG_OUTPUT_ESTIMATE){
                    if ((db_iWDG ==iWDGroup) && (Idx == endIdx -1) ){
                        int db_iAgeG = inst->db_ithAgeGroup - 1;

                        db_accountValueAtT[db_pathCount] = accountValue[db_iWDG][db_iAgeG];
                        db_benefitAtT[db_pathCount] = benefitGMWB[db_iWDG][db_iAgeG];
                        db_nonLapseAtT[db_pathCount] = groupSize[db_iWDG][db_iAgeG];
                    }
                }               
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}


// for the LogNormal path generator
CVolRequestLNArray InsuranceAnnuityGMBProd::getVolInterp(const IMCPathGenerator* pathGen,
                                                            int   iAsset) const {
    DateTime imntStartDate = inst->fwdStarting ? inst->startDate : inst->valueDate;
    CVolRequestLNArray reqarr(1); // one interp level per asset here

    reqarr[0] = CVolRequestLNSP(
    new LinearStrikeVolRequest(inst->asset->getSpot(), //ATM
                               imntStartDate, 
                               inst->endDate(0),
                               inst->fwdStarting));
    return reqarr;
}

// control variate is done here
void InsuranceAnnuityGMBProd::recordExtraOutput(Control* control, Results* results, const IMCPrices& prices) const {

    const InsuranceAnnuityGMBPrices& myprices = static_cast<const InsuranceAnnuityGMBPrices&>(prices);

    // pv from last date to value date
    double discountFactor = inst->instSettle->pv(inst->valueDate,  
                 histDates[histDates.size()-1],
                 inst->discount.get(),
                 inst->asset.get());

    if (control && control->isPricing()) {
        OutputRequest* request =
        control->requestsOutput(OutputRequest::OPTION_PRICE);
        if (request) {
            double optionPrice = myprices.getOptionPrice();
            // need to present value here
            optionPrice *= discountFactor;
            results->storeRequestResult(request, optionPrice);
        }
        request = control->requestsOutput(OutputRequest::FEE_PRICE);
        if (request) {
            double feePrice = myprices.getFeePrice();
            // need to present value here
            feePrice *= discountFactor;
            results->storeRequestResult(request, feePrice);
        }

        // ---debug output---

        if ((control && control->isPricing()) && (inst->DEBUG_OUTPUT)){
            int i;
            int nbVar = 10;
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
                (*tmp_name)[6] = "db_negativeAV";
                (*tmp_name)[7] = "db_payoff";
                (*tmp_name)[8] = "groupSize at Beg";
                (*tmp_name)[9] = "groupSize at End";

                for (i=0; i<numSmth; i++){

                    (*tmp)[i] = histDates[i];

                    (*smthMatrix)[0][i] = db_index[i];

                    (*smthMatrix)[1][i] = db_accountValueBeg[i];

                    (*smthMatrix)[2][i] = db_accountValue[i];

                    (*smthMatrix)[3][i] = db_benefitBeg[i];

                    //(*tmp_name)[0][1] = "db_accountValue";
                    (*smthMatrix)[4][i] = db_benefit[i];

                    (*smthMatrix)[5][i] = db_currentWithdrawals[i];

                    //(*tmp_name)[3] = "db_lastStepUp";
                    //(*smthMatrix)[3][i] = db_lastStepUp[i];
                    
                    (*smthMatrix)[6][i] = db_negativeAV[i];

                    (*smthMatrix)[7][i] = db_payoff[i];

                    (*smthMatrix)[8][i] = db_groupSizeBeg[i];

                    (*smthMatrix)[9][i] = db_groupSize[i];

                }

                results->storeGreek(tmp_name, Results::DEBUG_PACKET, OutputNameSP(new OutputName("PATH_INFO_NAME")));
                results->storeGreek(tmp, Results::DEBUG_PACKET, OutputNameSP(new OutputName("SIM_DATE")));
                results->storeGreek(smthMatrix, Results::DEBUG_PACKET, OutputNameSP(new OutputName("PATH_INFO")));
            }

            if (numSmth > 0 && inst->DEBUG_OUTPUT_ESTIMATE){
                CDoubleMatrixSP pathMatrix(new DoubleMatrix(numSmth, inst->db_maxIte));  //db_maxIte is the maximum nb of simulation
                CDoubleMatrixSP tempMatrix(new DoubleMatrix(nbVar1, inst->db_maxIte ));  //db_maxIte is the maximum nb of simulation
    
                for (i=0; i<numSmth; i++){
                    for (int j=0; j < inst->db_maxIte; j++){
                        (*pathMatrix)[i][j] = db_Mindex[i][j];
                    }
                }

                (*tmp_name1)[0] = "AV";
                (*tmp_name1)[1] = "Ben";
                (*tmp_name1)[2] = "GroupSize";

                for (int j=0; j < inst->db_maxIte; j++){
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


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* InsuranceAnnuityGMB::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(spotSamples->getAllDates());

    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(spotSamples->getDates()[0]));
    return new InsuranceAnnuityGMBProd(this, refLevel, simSeries);
}

CClassConstSP const InsuranceAnnuityGMB::TYPE = CClass::registerClassLoadMethod(
    "InsuranceAnnuityGMB", typeid(InsuranceAnnuityGMB), InsuranceAnnuityGMB::load);

// force linker to include this file (avoid having header file) */
extern bool InsuranceAnnuityGMBLoad()
{
    return true && InsuranceAnnuityGMB::TYPE;
}

DRLIB_END_NAMESPACE
