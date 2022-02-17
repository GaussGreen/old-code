
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InsuranceAGMWB2.cpp
//
//   Author      : Xiaolan Zhang
//
//   Description   Insurance annuity payoff
//                 common features for GMDB/GMAB/GMWB
//   Date        : Sept 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/InsuranceAGMWB.hpp"

DRLIB_BEGIN_NAMESPACE

//////// instrument class ////////
class InsuranceAGMWB2: public InsuranceAGMWB{

public:
    static CClassConstSP const TYPE;
    friend class InsuranceAGMWB2Prod;
    friend class InsuranceAGMWB2ProdSV;


    virtual void Validate(); 
   
    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

private:
    /** constructor*/
    InsuranceAGMWB2(CClassConstSP const type=TYPE): 
                                InsuranceAGMWB(type){

        stepUpITMLevels.resize(1,10000000.0);
        wDTime.resize(1, 10000000);  //withdrawal time.
        wDPercent.resize(1, 0);  //withdrawal %.
    }

    // for reflection
    InsuranceAGMWB2(const InsuranceAGMWB2& rhs); // not implemented
    InsuranceAGMWB2& operator=(const InsuranceAGMWB2& rhs); // not implemented

    static IObject* defaultInsuranceAGMWB2() {
        return new InsuranceAGMWB2();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

protected:
    //////// fields ////////
       
    // Guaranteed Minimum Withdrawal Benefit feature

    /** step up */
    DoubleArray     stepUpITMLevels;              // ITM level for step up= AV/Ben

    /** withdraw */
    /** Additional withdrawal features for Guaranteed Minimum Withdrawal Benefit feature */
    DoubleArray     wDTime;                       // withdrawal time.
    DoubleArray     wDPercent;                    // withdrawal time.

    ExpirySP        maturity;			          // maturity of the traded product
    double          longVol;                      // use flat longvol if idx > maturity
};

//////////////// product class /////////////////////////////////
//state variable

/* MC product class for insurance annuity */
class InsuranceAGMWB2ProdSV :   public InsuranceAGMWBProdSV{

public:    
    /** equivalent to InstIntoMCProduct. */
    InsuranceAGMWB2ProdSV(const InsuranceAGMWB2* inst,
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
    const InsuranceAGMWB2* instGMWB2;                     // reference to original instrument

    void decideWithdrawal(const int& timeFromStart,                                                     
                                    DoubleArraySP withdrawals,
                                    const int& iGroup,
                                    const int& iWDGroup 
                                    );

};
// end of class InsuranceAGMWB2Prod


DRLIB_END_NAMESPACE
