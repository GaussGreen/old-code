
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAGMWB1.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 common features for GMWB
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/InsuranceAGMWB.hpp"

DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////
class InsuranceAGMWB1: public InsuranceAGMWB{

public:
    static CClassConstSP const TYPE;
    friend class InsuranceAGMWB1Prod;
    friend class InsuranceAGMWB1ProdSV;

    virtual void Validate(); 
   
    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

private:
    /** constructor*/
    InsuranceAGMWB1(CClassConstSP const type=TYPE): 
                                InsuranceAGMWB(type), 
                                isDynamicStepUp(false),
                                stepUpITMLevel(-1.0),
                                isDynamicUtilization(false),
                                deterministicWDType(1){}

    // for reflection
    InsuranceAGMWB1(const InsuranceAGMWB1& rhs); // not implemented
    InsuranceAGMWB1& operator=(const InsuranceAGMWB1& rhs); // not implemented

    static IObject* defaultInsuranceAGMWB1() {
        return new InsuranceAGMWB1();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

protected:
    //////// fields ////////
       
    // Guaranteed Minimum Withdrawal Benefit feature
//    /** step up */
    bool            isDynamicStepUp;              // if yes, step ups are dynamic
    double          stepUpITMLevel;               // ITM level for step up

    /** withdraw */
    bool            isDynamicUtilization;         // if yes, the withdrawal rate in GMWB benefit is dynamic
                                                  // (i.e.determined with respect to ITM)

    int             deterministicWDType;          //withdraw starts at 
                                                  //1: 0; 
                                                  //2: different withdraw group may start WD at 0 or some time in the future; 
                                                  //3: different withdraw group may start WD at 0 or some time in the future depending on ITM

    /** Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature */
//    DoubleArray     wDITMLevels;                  // withdrawal ITM Levels.
    CDoubleMatrixSP startWDYears;                 // diff groups start to withdrawal 
                                                  // at diff. policy year and ITM levels         
};

//////////////// product class /////////////////////////////////
///State variable version//////////////////////////////////

/* MC product class for insurance annuity */
class InsuranceAGMWB1ProdSV :   public InsuranceAGMWBProdSV{

public:    
    /** equivalent to InstIntoMCProduct. */
    InsuranceAGMWB1ProdSV(const InsuranceAGMWB1* inst,
                            const SimSeriesSP& simSeries) ;


    /** update step up */
    virtual void updateStepUp(const int& timeFromStart, const int& iWDGroup);

    /** compute coupons taken by clients */
    virtual void computeWithdrawals(const int& timeFromStart,
                            DoubleArraySP withdrawals,
                            const int& iWDGroup);

    /** update number of clients by taking into account the mortality rate */
    virtual void updateDeadClients(const int& timeFromStart, const int& iWDGroup);

    /** update number of clients by taking into account the lapse rate */
    virtual void updateLapseClients(const int& timeFromStart, DoubleArraySP withdrawals, const int& iWDGroup);

    

    /** calc payoff for each path in order to price diff group withdrawal inside payoff*/
    virtual  void calcPayoff(const int iWDGroup, const int db_pathCount);

private:
    const InsuranceAGMWB1* instGMWB1;                     // reference to original instrument

    void calcDeterministicWDType3(const int& timeFromStart,                                                     
                                    DoubleArraySP withdrawals,
                                    const int& iGroup,
                                    const int& iWDGroup 
                                    );

};
// end of class InsuranceAGMWB1Prod

DRLIB_END_NAMESPACE
