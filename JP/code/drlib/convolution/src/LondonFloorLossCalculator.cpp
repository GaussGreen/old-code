//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson (from CreditMetricsLondonFloor.cpp/ccm2 lib)
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LondonFloorLossCalculator.hpp"
#include "edginc/ICDSParSpreads.hpp"

DRLIB_BEGIN_NAMESPACE

LondonFloorLossCalculator::~LondonFloorLossCalculator(){}

#define TINY 1e-10

/**
 * partition the tranche [K1, K2] into 3 tranches:
 * an equity, a mezz and a senior part
 * if any tranche is of size 0, its
 */
static void ccmSplitStrike(
    double   K1,   /* (I) lower strike    */
    double   K2,   /* (I) upper strike    */
    double   Keq,  /* (I) equity strike   */
    double   Kse,  /* (I) senior strike   */
    vector<double>& Teq,  /* (O) tranche equity strikes (size 0 or 2) */
    vector<double>& Tme,  /* (O) tranche mezz   strikes (size 0 or 2) */
    vector<double>& Tse)  /* (O) tranche senior strikes (size 0 or 2) */
{
    static char routine[] = "ccmSplitStrike";

    Teq.clear();
    Tme.clear();
    Tse.clear();

    if (K2 < K1) {
        throw ModelException(routine, "strike K2 must be greater than K1");
    }
    if (Kse < Keq) {
        throw ModelException(routine, "strike Kse must be greater than Keq");
    }
    /* set Teq */
    if (K1 < Keq) {
        Teq.resize(2);
        Teq[0] = K1;
        Teq[1] = K2 > Keq ? Keq : K2;
    }

    /* set Tse */
    if (K2 > Kse) {
        Tse.resize(2);
        Tse[0] = K1 < Kse ? Kse : K1;
        Tse[1] = K2;
    }

    /* set Tme */
    if (K2 > Keq && K1 < Kse ) {
        Tme.resize(2);
        Tme[0] = K1 < Keq ? Keq : K1;
        Tme[1] = K2 > Kse ? Kse : K2;
    }

    ASSERT(Keq <= Kse);
    ASSERT(Teq.empty() || (Teq[0] <= Teq[1] && Teq[1] <= Keq) );
    ASSERT(Tme.empty() || (Tme[0] <= Tme[1] && Keq <= Tme[0] && Tme[1] <= Kse));
    ASSERT(Tse.empty() || (Tse[0] <= Tse[1] && Kse <= Tse[0]));
}

/** Constructor - takes in the loss calculator which does the unadjusted  
    loss. */
LondonFloorLossCalculator::LondonFloorLossCalculator(
        ITrancheLossCalculatorConstSP originalLossCalculator, /* (I) */
        const DateTimeArray&          timeline,               /* (I) */
        const DateTime&               today,                  /* (I) */
        bool                          computeCondCurve,       /* (I) */
        ICDSParSpreadsConstSP         londonFloor,            /* (I) */
        double                        portfolioNotional,      /* (I) */
        double                        equityStrike,           /* (I) */
        double                        seniorStrike):          /* (I) */
    originalLossCalculator(originalLossCalculator), 
    computeCondCurve(computeCondCurve),
    londonFloorProb(timeline.size()), totalNotional(portfolioNotional),
    equityStrike(equityStrike), seniorStrike(seniorStrike)
{
    try {
        // compute probabilities
        DefaultRatesSP curve(londonFloor->defaultRates());
        for (int i = 0; i < timeline.size(); ++i) {
            londonFloorProb[i] = curve->calcDefaultPV(today, timeline[i]);
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "LondonFloorLossCalculator::"
                             "LondonFloorLossCalculator");
    }
}

/** Calculate the expected loss for specified time point and strikes
    using 'London Floor' methodology */
void LondonFloorLossCalculator::loss(
    int     timePoint,        // (I) do the calculation for this timepoint
    double  k1,               /* (I) lower strike      */
    double  k2,               /* (I) upper strike      */
    double& loss,             /* (O) tranche loss amt  */
    double& lossCond) const{  /* (O) tranche loss amt cond on cpty
                                     surviving */
    IKeySP key(lossKey(timePoint));
    key->loss(k1, k2, loss, lossCond);
}

// implementation of IKey when using a CCMTrancheCalculatorLegacy
class LondonFloorLossCalculator::Key: public virtual IKey{
    IKeySP                           originalKey;
    int                              timePoint;
    LondonFloorLossCalculatorConstSP londonFloorCalculator;
public:
    Key(IKeySP                           originalKey, 
        int                              timePoint,
        LondonFloorLossCalculatorConstSP londonFloorCalculator):
        originalKey(originalKey), timePoint(timePoint),
        londonFloorCalculator(londonFloorCalculator){}
    
    /** Calculates the appropriate rate/factor between the two dates */
    virtual void loss(
        double  k1,                   /* (I) lower strike      */
        double  k2,                   /* (I) upper strike      */
        double& loss,                 /* (O) tranche loss amt  */
        double& lossCond) const{      /* (O) tranche loss amt cond on cpty
                                         surviving */
        londonFloorCalculator->loss(originalKey, timePoint,
                                    k1, k2, loss, lossCond);
    }
};

/** Returns a key used to optimise repeated calculations of losses
    at the same time point. */
ITrancheLossCalculator::IKey* LondonFloorLossCalculator::lossKey(
    int timePoint) const // (I) do the calculation for this timepoint
{
    IKeySP originalKey(originalLossCalculator->lossKey(timePoint));
    return new Key(originalKey, timePoint,
                   LondonFloorLossCalculatorConstSP::attachToRef(this));
}

/** Same as loss() method above but takes an IKey rather than a 
    ITrancheLossCalculator */
void LondonFloorLossCalculator::loss(
    IKeySP  originalKey,          // (I) access to original calculator
    int     timePoint,            // (I) for accessing londonFloorProb
    double  k1,                   /* (I) lower strike      */
    double  k2,                   /* (I) upper strike      */
    double& loss,                 /* (O) tranche loss amt  */
    double& lossCond) const{      /* (O) tranche loss amt cond on cpty
                                         surviving */
    static char routine[] = "LondonFloorLossCalculator::loss";
    try{
        // first calculate unadjusted loss
        originalKey->loss(k1, k2, loss, lossCond);

        double seniorProb = londonFloorProb[timePoint]; // for ease
        /* optimization */
        if (seniorProb > 1.-1e-14){
            return;
        }

        /* input checks */
        if (loss < 0.0) {
            throw ModelException (routine, "Unadjusted loss < 0.");
        }
        if (computeCondCurve && lossCond < 0.0) {
            throw ModelException(routine, "Unadjusted loss conditional on "
                                 "counterparty < 0.");
        }
        /* identify what calculations need be done */
        vector<double> equityTranche;
        vector<double> mezzanineTranche;
        vector<double> seniorTranche;
        ccmSplitStrike(k1, k2, equityStrike, seniorStrike, 
                       equityTranche, mezzanineTranche, seniorTranche);

        // initialise our adjustment
        double adj = 0.0;
        double adjCond = 0.0;
        /* mezzanine part */
        if (!mezzanineTranche.empty()) {
            double LTme[2];
            double Lmezz[2];
            
            /* calculate loss for full mezz area */
            if (fabs(equityStrike-k1) < TINY && fabs(seniorStrike-k2) < TINY) {
                Lmezz[0] = loss;
                Lmezz[1] = lossCond;
            } else {
                originalKey->loss(equityStrike, seniorStrike, 
                                  Lmezz[0], Lmezz[1]);
            }
            
            /* calculate loss for mezz part of the tranche */
            if (fabs(mezzanineTranche[0]-k1)<TINY && 
                fabs(mezzanineTranche[1]-k2)<TINY) {
                LTme[0] = loss;
                LTme[1] = lossCond;
            } else {
                originalKey->loss(mezzanineTranche[0], mezzanineTranche[1], 
                                  LTme[0], LTme[1]);
            }
            
            /* calculate mezz adjustment */
            adj = - (totalNotional-seniorStrike) * 
                (LTme[0]/Lmezz[0]) * (1.-seniorProb);
            adjCond = - (totalNotional-seniorStrike) * 
                (LTme[1]/Lmezz[1]) * (1.-seniorProb);
        }
        
        /* senior part */
        if (!seniorTranche.empty()) {
            adj     += (seniorTranche[1]-seniorTranche[0]) * (1.-seniorProb);
            adjCond += (seniorTranche[1]-seniorTranche[0]) * (1.-seniorProb);
        }

        /* for correlator v1 compatibility */
        if (adj < -loss) {
            adj = -loss;
        }
        if (adjCond < -lossCond) {
            adjCond = -lossCond;
        }
        // then return results
        loss += adj;
        lossCond += adjCond;
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

DRLIB_END_NAMESPACE
