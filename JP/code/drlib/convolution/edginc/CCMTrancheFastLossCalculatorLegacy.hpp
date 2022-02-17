//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CCMTRANCHEFASTLOSSCALCULATORLEGACY_HPP
#define QR_CCMTRANCHEFASTLOSSCALCULATORLEGACY_HPP

#include "edginc/TrancheLossCalculatorLegacy.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(CreditTrancheLossConfig);
FORWARD_DECLARE(CounterPartyCredit);

//// Need to scrap this class
class CONVOLUTION_DLL CCMTrancheFastLossCalculatorLegacy : 
    public ITrancheLossCalculatorLegacy
{
public:
    CCMTrancheFastLossCalculatorLegacy () {}    
    virtual ~CCMTrancheFastLossCalculatorLegacy() {}

    /** Builds calculator skipping all variant data ie data that changes
        across time */
    CCMTrancheFastLossCalculatorLegacy(
        bool                           doCreditMetrics,
        bool                           doRFL,
        const bool                     recoverNotional,
        CreditTrancheLossConfigConstSP tranche,
        CounterPartyCreditConstSP      counterParty);

    /** creates a copy */
    CCMTrancheFastLossCalculatorLegacy(
        CCMTrancheFastLossCalculatorLegacy* copy);

    /** This setup is the same as the setup in CDO.cpp but
        takes a ConvolutionProduct (the goal is to scrap the one in
        CDO.cpp and break the dependency between it and CDO). The betaOverride
        must be either empty (=> no override) or the size of the other arrays */
    void setupCM(
        const DoubleArray&        survivalProb,      /* name survival proba */
        double                    counterPartyProb, /* name survival proba */
        const DoubleArray&        betaOverride); // optional per name override

    /** Same as setupCM but for CCM rather than Credit Metrics. The betaOverride
        must be either empty (=> no override) or the size of the other arrays */
    void setupCCM(
        const DoubleArray&        survivalProb,      /* name survival proba */
        const DoubleArray&        floorSurvivalProb, /* from senior curve */
        double                    counterPartyProb, /* name survival proba */
        const DoubleArray&        betaOverride); // optional per name override

    virtual void loss(
        double        K1,    /* (I) lower strike      */
        double        K2,    /* (I) upper strike      */
        double       &L,     /* (O) tranche loss amt  */
        double       &Lcond) /* (O) tranche loss amt cond on cpty surviving */
        const;

    CCMConvolution::NameParamArray& getBasketInfo();
        
private:
    double                         pastLoss;
    CCMConvolution::NameParamArray basketInfo;
    CCMConvolution::NameParamSP    cpty;
};

typedef refCountPtr<
    CCMTrancheFastLossCalculatorLegacy> CCMTrancheFastLossCalculatorLegacySP;


DRLIB_END_NAMESPACE
#endif
