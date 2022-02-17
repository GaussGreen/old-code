
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAGMWB.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 common features for GMWB
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/InsuranceA.hpp"

DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////
class InsuranceAGMWB: public InsuranceA{

public:
    static CClassConstSP const TYPE;
    friend class InsuranceAGMWBProd;
    friend class InsuranceAGMWBProdSV;

    virtual void Validate(); 
   
    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

protected:
    /** constructor*/
    InsuranceAGMWB(CClassConstSP const type=TYPE): 
                                InsuranceA(type), 
                                adjustmentAV(1.0),
                                isStepUp(false),
                                initialWaitPeriod(-1),
                                stepUpPeriod(-1){
        DEBUG_OUTPUT_ESTIMATE = true;
        db_maxIte = 100001;

        wDITMLevels.resize(1, -1000000.0);  //withdrawal ITM Levels.

        DEBUG_OUTPUT = false;
        db_ithPath = 1 ;
        db_ithAgeGroup = 1;
        db_ithWDGroup = 1;  //1 if deterministicWDType =1
        db_ithAsset = 0;    //0 for the mono underlying or the first asset in Basket
    }

    // for reflection
    InsuranceAGMWB(const InsuranceAGMWB& rhs); // not implemented
    InsuranceAGMWB& operator=(const InsuranceAGMWB& rhs); // not implemented

    static IObject* defaultInsuranceAGMWB() {
        return new InsuranceAGMWB();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

protected:
    //////// fields ////////
       
    // Guaranteed Minimum Withdrawal Benefit feature
    double          maxGMWBPercentage;            // maximal percentage of initial 
                                                  //amount for annual withdrawal
    double          adjustmentAV;                 // a buffer to ajust the AV

    /** step up */
    bool            isStepUp;                     // is there a step up feature
//    bool            isDynamicStepUp;              // if yes, step ups are dynamic
    int             initialWaitPeriod;            // initial waiting period
    int             stepUpPeriod;                 // step up period
    //double          stepUpITMLevel;               // ITM level for step up

    /** withdraw */
    //bool            isDynamicUtilization;         // if yes, the withdrawal rate in GMWB benefit is dynamic
                                                  // (i.e.determined with respect to ITM)

    //int             deterministicWDType;          //withdraw starts at 
                                                  //1: 0; 
                                                  //2: different withdraw group may start WD at 0 or some time in the future; 
                                                  //3: different withdraw group may start WD at 0 or some time in the future depending on ITM

    /** Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature */
    DoubleArray     wDITMLevels;                  // withdrawal ITM Levels.
    //CDoubleMatrixSP startWDYears;                 // diff groups start to withdrawal 
                                                  // at diff. policy year and ITM levels         


    /** Debug output */
    bool            DEBUG_OUTPUT_ESTIMATE ;        //output info to do the estimation between 10 years and 30 years
    int             db_maxIte;                     // hard caoded, for max iteration that the DEBUG_OUTPUT_ESTIMATE can handle $unregistered

    bool            DEBUG_OUTPUT ;
    int             db_ithPath ;
    int             db_ithAgeGroup ;
    int             db_ithWDGroup;                  //1 if deterministicWDType =1
    int             db_ithAsset;                    // output iAsset path
};

//////////////// product class /////////////////////////////////

/* MC product class for insurance annuity */
class InsuranceAGMWBProdSV :   public InsuranceAProdSV{

public:    
    /** equivalent to InstIntoMCProduct. */
    InsuranceAGMWBProdSV(const InsuranceAGMWB* inst,
                            const SimSeriesSP& simSeries) ;

    /** compute fees (supposed to be paid at the begining of the period) */
    /** return this period fees */
    virtual double computeFees(double* fees, const int& timeFromStart, const int& iWDGroup);

    /** update benefits for Guaranteed Minimum Withdrawal Benefit (GMWB)*/
    virtual void updatebenefit(DoubleArraySP withdrawals, const int& iWDGroup);

    /** update number of clients by taking into account the mortality rate */
    //virtual void updateDeadClients(const int& timeFromStart, const int& iWDGroup);

    /** update account value */
    virtual void updateAccountValue(DoubleArraySP withdrawals, bool alreadySettled, const int& iWDGroup);


    virtual void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices);

    /** print extra output **/
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices&) const;


    /** update step up */
    virtual void updateStepUp(const int& timeFromStart, const int& iWDGroup){}

    /** compute coupons taken by clients */
    virtual void computeWithdrawals(const int& timeFromStart,
                            DoubleArraySP withdrawals,
                            const int& iWDGroup){}


    /** update number of clients by taking into account the lapse rate */
    virtual void updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup){}

   
    /** calc payoff for each path in order to price diff group withdrawal inside payoff*/
    virtual  void calcPayoff(const int iWDGroup,const int db_pathCount){};

    //for output purpose
    void Debug_outputBegin(const bool doingPast,  
                            const int iAsset,   //0 if mono        
                            const int iWDGroup, 
                            const int db_pathCount, 
                            const int Idx);

    void Debug_outputMid( 
                            const int iWDGroup, 
                            const int iGroup, 
                            const int db_pathCount, 
                            const int Idx);

    void Debug_outputEnd(const bool doingPast, 
                        const DoubleArraySP withdrawals,
                        const double JPMFees_Idx,
                        const int iAsset,   //0 if mono
                        const int iWDGroup, 
                        const int db_pathCount, 
                        const int Idx);
protected:

    const InsuranceAGMWB* instGMWB;                     // reference to original instrument

    vector<bool>                isWaitPeriod;              // are we in the initial waiting period ?
    vector<IntArray >           lastStepUp;                // nber of years since the last step-up for the different groups
    //vector<DoubleArray>         benefit;               // garanteed minimum withdrawal benefit for the different age groups
    vector<DoubleArray>         lastMaxGMWBWithdrawal;     // last maximum GMWB withdrawal
    vector<BoolArray>           isWithdraw;                //always true except for GMWB when deterministicWDType =2 or 3 or isDynamicUtilization = true

    /** for preservation of the past */
    vector<bool>                isWaitPeriodSoFar;         // are we in the initial waiting period ?
    vector<IntArray>            lastStepUpSoFar;           // nber of years since the last step-up for the different groups
    //vector<DoubleArray>         benefitSoFar;          // garanteed minimum withdrawal benefit for the different age groups
    vector<DoubleArray>         lastMaxGMWBWithdrawalSoFar;// last maximum GMWB withdrawal
    vector<BoolArray>           isWithdrawSoFar;

    /**for debug output*/

    int                         db_pathCount;
    DoubleArray                 db_index; //index value 
    DoubleArray                 db_accountValueBeg;     
    DoubleArray                 db_benefitBeg;     
    DoubleArray                 db_groupSizeBeg;
    DoubleArray                 db_accountValue;   
    DoubleArray                 db_benefit;     
    DoubleArray                 db_currentWithdrawals ;

    DoubleArray                 db_lastStepUp ; 
    DoubleArray                 db_groupSize;
    DoubleArray                 db_negativeAV;               //account value when option needs to pay, ie AV <0, only for GMWB
    DoubleArray                 db_payoff;
    DoubleArray                 db_JPMFees;

    /**debug output info for extimation between 10 and 30 years*/
    vector<DoubleArray>         db_Mindex;                   //index value  //[Idx, itera] 
    DoubleArray                 db_benefitAtT;               //benefits amount at maturity,   [iter] 
    DoubleArray                 db_accountValueAtT;          // acount value at maturity,  [iter]
    DoubleArray                 db_nonLapseAtT;              //non lapse at maturity, [iter]

};
// end of class InsuranceAGMWBProd

DRLIB_END_NAMESPACE
